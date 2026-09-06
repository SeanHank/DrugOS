"""Licensing-safe literature references for the benchmark corpus."""

CITATIONS: dict[str, str] = {
    "midazolam": (
        "IV PK summary: Greenblatt D.J. et al. (Clin Pharmacokinet 1984); "
        "Goodman & Gilman's Pharmacological Basis of Therapeutics, 13th ed.; "
        "Versed (midazolam) Injection USPI PK table."
    ),
    "acetaminophen": (
        "Goodman & Gilman 13th ed.; Tylenol (acetaminophen) USPI: CL ~0.24 "
        "L/h/kg, Vss ~0.7-0.9 L/kg, t1/2 ~1.5-3 h, Fa ~0.85-0.95."
    ),
    "warfarin": (
        "Goodman & Gilman 13th ed.; Coumadin (warfarin) USPI: high plasma "
        "protein binding (fup ~1%), small Vd (8-14 L), long t1/2."
    ),
    "ciprofloxacin": (
        "Cipro (ciprofloxacin) USPI; Bergan T. (J Antimicrob Chemother "
        "1986): CL_oral ~21.6 L/h, ~60-70% renal, Fa ~0.7, t1/2 ~3-6 h."
    ),
    "dofetilide": (
        "Tikosyn (dofetilide) USPI: fup 0.36, CL ~16.3 L/h (~80% renal), "
        "Vss 3.3-4.7 L/kg, t1/2 ~7-9 h, Fa ~0.9; Smith D.A. et al. "
        "(Br J Clin Pharmacol 1992)."
    ),
    "partition": (
        "Rodgers T., Rowland M. (J Pharm Sci 2006;95:1238-57) tissue:plasma "
        "partition model; Ye M., Nagar S., Korzekwa K. (Biopharm Drug "
        "Dispos 2016;37:123-141) perfusion-limited whole-body PBPK graph."
    ),
    "occupancy": (
        "Daryaee F., Tonge P.J. (Annu Rev Pharmacol Toxicol 2019;59:507-529) "
        "pharmacometric target-turnover occupancy model; see doc/05 2.4."
    ),
    "pathway": (
        "Huang C.-Y., Ferrell J.E. Jr. (Proc Natl Acad Sci USA 1996;93:10078-83) "
        "ultrasensitivity and signal amplification in the MAPK cascade; "
        "DrugOS 3-tier MAPK with basal+signal arms, doc/05 3.2-3.5."
    ),
}

__all__ = ["CITATIONS"]
