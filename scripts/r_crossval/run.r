#!/usr/bin/env Rscript
# DrugOS R-side cross-validation verdict.
# Independent re-computation in R (stats/utils only) of AUC, terminal t1/2 and
# the tail lambda_z from a DrugOS simulated plasma curve, plus PASS/WARN/FAIL
# against the vendored published bands.
#
# Usage:
#   Rscript run.r --sim <sim.tsv> --bands <bands.tsv> --out <report.md>
# sim.tsv:    time_h \t plasma_total_mg_l  (header included, no # comments)
# bands.tsv:  quantity \t lo \t hi  where quantity is a published_pk.json key,
#             e.g. {cl_plasma_l_h, t_half_h} (other rows are ignored)
suppressPackageStartupMessages({
  library(stats)   # nls, lm      (base R)
  library(utils)   # read.table  (base R)
})

args <- commandArgs(trailingOnly = TRUE)
argv <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1L]
}
sim_f   <- argv("--sim")
bands_f <- argv("--bands")
out_f   <- argv("--out", paste0("report_", gsub("\\W", "", Sys.time()), ".md"))
if (is.na(sim_f) || is.na(bands_f)) {
  stop("--sim and --bands are required")
}

dat <- read.delim(sim_f, header = TRUE, check.names = FALSE)
t    <- dat$time_h
c    <- dat$plasma_total_mg_l
if (!is.finite(t[1]) || t[1] < 0 || max(t) <= 0 || length(c) < 4) {
  stop("sim.tsv must contain a numeric time_h and plasma_total_mg_l column")
}

# --- 1) AUC_t by trapezoid on the full curve (relaxed, independent of python) ---
auc_t <- sum(0.5 * (c[-1] + c[-length(c)]) * diff(t))
auc_t <- if (is.finite(auc_t) && auc_t >= 0) auc_t else NA_real_

# --- 2) terminal phase on the tail: last 25% of times above 5% Cmax ---
cmax  <- max(c, na.rm = TRUE)
tail_idx <- which(t > 0.25 * max(t) & c > 0.05 * cmax)
if (length(tail_idx) >= 3 && c[tail_idx[length(tail_idx)]] > 0) {
  l <- lm(log(c[tail_idx]) ~ t[tail_idx])
  k  <- -coef(l)[2]
  t_half <- log(2) / k
  k_se <- summary(l)$coefficients[2, 2]
} else {
  tail_idx <- integer()
  k <- NA_real_; t_half <- NA_real_; k_se <- NA_real_
}

# --- 3) independent nls mono-exponential alternative on the tail ---
nls_fit <- NULL
if (length(tail_idx) >= 3) {
  tryCatch({
    tl <- t[tail_idx]; cl <- c[tail_idx]
    nls_fit <- nls(cl ~ A * exp(-kk * tl),
                   start = list(A = cl[1], kk = max(k, 0.01)),
                   control = nls.control(maxiter = 100, warnOnly = TRUE))
  }, error = function(e) NULL)
}

# --- 4) bands verdict ---
bands <- read.delim(bands_f, header = TRUE, check.names = FALSE)
v_band <- function(quantity, value) {
  row <- bands[as.character(bands$quantity) == quantity, , drop = FALSE]
  if (nrow(row) == 0 || !is.finite(value)) return(list(status = "N/A"))
  lo <- min(row$lo); hi <- max(row$hi)
  status <- if (value < lo) "FAIL (below band)" else if (value > hi) "FAIL (above band)" else "PASS"
  list(status = status, lo = lo, hi = hi)
}

cl_v    <- v_band("cl_plasma_l_h", NA_real_)  # cl needs dose/AUC_inf; not computed pointwise here
th_v    <- if (is.finite(t_half)) v_band("t_half_h", t_half) else list(status = "N/A")
overall <- "PASS"
if (grepl("FAIL", th_v$status)) overall <- "FAIL"
if (all(!is.finite(c(t_half, k)))) overall <- "WARN (no estimable terminal phase)"

md <- c(
  "# R cross-validation report",
  "",
  sprintf("- Simulated curve: `%s`", sim_f),
  sprintf("- Points: %d, time span %.2f h, Cmax %.4f mg/L", length(c), max(t), cmax),
  sprintf("- AUC0-t (R trapezoid): %s", ifelse(is.finite(auc_t), sprintf("%.3f mg.h/L", auc_t), "N/A")),
  sprintf("- Terminal t1/2 (R log-linear): %s", ifelse(is.finite(t_half), sprintf("%.3f h", t_half), "N/A")),
  sprintf("- lambda_z: %s /h (SE %s)", ifelse(is.finite(k), sprintf("%.5f", k), "N/A"),
          ifelse(is.finite(k_se), sprintf("%.5f", k_se), "N/A")),
)
if (!is.null(nls_fit)) {
  cf <- coef(nls_fit)
  md <- c(md, sprintf("- nls mono-exponential: A=%.4f mg/L, k=%.5f /h, t1/2=%.3f h",
                      cf[1], cf[2], log(2) / cf[2]))
}
md <- c(md,
  sprintf("- Verdict within published band (t1/2): %s (band [%.3f, %.3f] h)",
          th_v$status, th_v$lo, th_v$hi),
  "",
  sprintf("**OVERALL: %s**", overall),
  "",
  "*Generated independently in R; input bands sourced from data/benchmarks/published_pk.json.*"
)
writeLines(md, out_f)
cat("Wrote", out_f, "overall:", overall, "\n")