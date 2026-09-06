# DrugOS literature pharmacokinetic reference model (base R).
#
# Implements the classic structural equations cited across doc/02 and doc/05:
#
#   - Bateman / macro-rate bi-exponential disposition with absorption
#     (Wagner 1976, "Linear pharmacokinetic equations"; Gibaldi & Perrier 1982,
#     "Pharmacokinetics", 2nd ed.):
#
#         C(t) = A*exp(-alpha*t) + B*exp(-beta*t)
#
#   - Trapezoidal AUC with tail correction to infinity via the terminal slope
#     lambda_z (Rowland & Tozer 2010, "Clinical Pharmacokinetics", ch. 4):
#
#         AUC_inf = AUC_t + C_last / lambda_z       CL = F*D / AUC_inf
#
# This R source is the *same* estimator DrugOS re-implements in numpy; the
# bridge checks both implementations agree, so the Python report value is
# independently reproduced by the literature code (doc/06, R bridge).

verify_pk <- function(t, c, dose) {
  stopifnot(is.numeric(t), is.numeric(c), is.numeric(dose),
            length(dose) == 1L, dose > 0,
            length(t) == length(c), length(t) >= 4L)

  keep <- is.finite(t) & is.finite(c) & t >= 0 & c >= 0
  t <- t[keep]
  c <- c[keep]
  n <- length(t)
  cmax <- max(c)

  # --- trapezoid AUC to last sample (identical rule as the numpy estimator) ---
  auc_t <- 0.0
  if (n >= 2L) auc_t <- sum(0.5 * (c[-1L] + c[-n]) * diff(t))

  # --- terminal phase on the last 25% of the time span -----------------------
  sel <- which(t >= 0.75 * max(t) & c > 0.05 * cmax)
  tail_ok <- length(sel) >= 3L && c[sel[length(sel)]] > 0

  lam <- NaN
  auc_inf <- auc_t
  if (tail_ok) {
    fit_l <- lm(log(c[sel]) ~ t[sel])
    lam <- -coef(fit_l)[2L]
    auc_inf <- auc_t + c[sel[length(sel)]] / lam
  }

  cl <- dose / auc_inf

  # --- terminal half-life: from lambda_z + classic method-of-residuals peel -----
  # (Greenblatt & Koch-Weser 1975; Rowland & Tozer 2010, ch.4): subtract the
  # terminal exponential and re-regress the residual early phase -> alpha.
  t_half_beta <- NaN
  two_comp <- 0.0
  if (tail_ok && is.finite(lam) && lam > 0) {
    t_half_beta <- log(2) / lam
    b0 <- coef(fit_l)[1L]
    b_term <- exp(b0) * exp(-lam * t)
    residual <- c - b_term
    early <- which(t <= quantile(t, 0.5) & residual > 0 & is.finite(residual))
    if (lam > 0 && length(early) >= 3L) {
      peel_l <- lm(log(residual[early]) ~ t[early])
      al <- -coef(peel_l)[2L]
      if (is.finite(al) && al > lam) two_comp <- 1.0
    }
  }

  c(n_obs = unname(n), auc_inf_mg_hl = unname(auc_inf), cl_l_h = unname(cl),
    t_half_beta_h = unname(t_half_beta), two_comp = unname(two_comp))
}