"""Coverage for the single-nephron tubular transport cascade (G3, 100 % branch).

doc/12 L17: the segmental micropuncture-style model turns a filtered load into
urine through the proximal tubule, loop of Henle and distal/collecting duct,
with saturable Na+ and xenobiotic transport.  These unit tests pin every
branch of the cascade; the analytic anchors (mass conservation, FE = 1 / > 1 /
< 1, transporter saturability) are asserted again at suite level in
``case_nephron_tubular_transport``.
"""

from __future__ import annotations

import pytest

from drugos.organ.nephron import (
    NephronHandling,
    NephronParams,
    _pt_xenobiotic_reabsorption,
    _secretory_flux,
    nephron_handling,
)

_CFREE_UM = 1.0e-3  # 1 uM as mmol/L


def test_filtered_na_is_c_times_gfr() -> None:
    res = nephron_handling(125.0, plasma_na_mm=140.0)
    # 140 mM * 125 ml/min -> 17.5e6 nmol/min (checked in the module doc).
    assert res.filtered_na_nmol_min == pytest.approx(17.5e6, rel=1e-9)
    assert res.excreted_na_nmol_min == pytest.approx(0.009 * 17.5e6, rel=1e-6)
    assert res.fe_na == pytest.approx(0.009, rel=1e-6)
    assert res.urine_flow_ml_min == pytest.approx(125.0 * 0.009, rel=1e-6)


def test_na_mass_balance_holds() -> None:
    res = nephron_handling(90.0, plasma_na_mm=135.0)
    assert (
        abs(res.filtered_na_nmol_min - res.reabsorbed_na_nmol_min - res.excreted_na_nmol_min) < 1e-6
    )


def test_pure_filtration_fex_one_cl_equals_gfr() -> None:
    res = nephron_handling(
        125.0,
        c_free_mm=_CFREE_UM,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=0.0),
    )
    assert res.fe_x == pytest.approx(1.0, rel=1e-9)
    assert res.cl_r_ml_min == pytest.approx(125.0, rel=1e-9)


def test_net_reabsorption_fex_below_one() -> None:
    res = nephron_handling(
        125.0,
        c_free_mm=_CFREE_UM,
        params=NephronParams(
            pas_reab_fraction=0.5,
            sec_vmax_nmol_min_per_nephron=0.0,
        ),
    )
    assert res.fe_x == pytest.approx(0.5, rel=1e-9)
    assert res.cl_r_ml_min == pytest.approx(62.5, rel=1e-9)


def test_net_secretion_fex_above_one() -> None:
    # A saturable secretory flux (OAT/OCT style) above the filtered load yields
    # fractional excretion > 1 and a renal clearance above the GFR.
    res = nephron_handling(
        125.0,
        c_free_mm=_CFREE_UM,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=4.0e-3),
    )
    assert res.fe_x > 1.0
    assert res.cl_r_ml_min > 125.0
    assert res.secreted_x_nmol_min > 0.0


def test_secretion_saturates_with_km() -> None:
    # At high free concentration the MM flux approaches vmax: doubling a large
    # vmax no longer changes the excreted amount, and the flux tracks the
    # analytic MM value exactly.
    low = _secretory_flux(0.05, 2.0e-3, 0.05)
    assert low == pytest.approx(1.0e-3, rel=1e-9)
    assert _secretory_flux(1e9, 2.0e-3, 0.05) == pytest.approx(2.0e-3, rel=1e-9)
    assert _secretory_flux(0.0, 2.0e-3, 0.05) == pytest.approx(0.0)


def test_full_reabsorption_zero_excretion() -> None:
    res = nephron_handling(
        125.0,
        c_free_mm=_CFREE_UM,
        params=NephronParams(pas_reab_fraction=1.0, sec_vmax_nmol_min_per_nephron=0.0),
    )
    assert res.excreted_x_nmol_min == pytest.approx(0.0)
    assert res.fe_x == pytest.approx(0.0)
    assert res.cl_r_ml_min == pytest.approx(0.0)


def test_no_drug_returns_none_fex() -> None:
    res = nephron_handling(125.0)
    assert res.fe_x is None
    assert res.cl_r_ml_min == 0.0
    assert res.excreted_x_nmol_min == 0.0
    assert isinstance(res, NephronHandling)
    d = res.to_dict()
    assert d["fe_x"] is None and d["cl_r_ml_min"] == 0.0


def test_plasma_na_override() -> None:
    res = nephron_handling(125.0, plasma_na_mm=100.0)
    assert res.filtered_na_nmol_min == pytest.approx(100.0 * 125.0 * 1000.0)


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError, match="gfr"):
        nephron_handling(0.0)
    with pytest.raises(ValueError, match="n_nephrons"):
        nephron_handling(125.0, n_nephrons=0.0)
    with pytest.raises(ValueError, match="plasma_na_mm"):
        nephron_handling(125.0, plasma_na_mm=0.0)
    with pytest.raises(ValueError, match="c_free_mm"):
        nephron_handling(125.0, c_free_mm=-1.0)
    with pytest.raises(ValueError, match="reabsorption fractions"):
        nephron_handling(125.0, params=NephronParams(pt_na_reab_fraction=1.5))
    with pytest.raises(ValueError, match="pas_reab_fraction"):
        _pt_xenobiotic_reabsorption(10.0, 1.5)
    with pytest.raises(ValueError, match="vmax"):
        _secretory_flux(1.0, -1.0, 0.05)
    with pytest.raises(ValueError, match="km"):
        _secretory_flux(1.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="c_free_mm"):
        _secretory_flux(-1.0, 1.0, 0.05)


def test_passive_reabsorption_half() -> None:
    assert _pt_xenobiotic_reabsorption(100.0, 0.0) == pytest.approx(0.0)
    assert _pt_xenobiotic_reabsorption(100.0, 0.25) == pytest.approx(25.0)


def test_fex_saturates_with_concentration() -> None:
    # Marginal digestion proof that the cascade's secretion is genuinely
    # Michaelis-Menten: at the same vmax, a ten-fold higher free concentration
    # no longer keeps fractional excretion constant — the saturable flux rolls
    # off while the filtered load scales linearly, so FE and the clearance
    # both fall (transport saturation, a real drug-drug-interaction probe).
    low = nephron_handling(
        100.0,
        c_free_mm=0.001,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=2.0e-3),
    )
    high = nephron_handling(
        100.0,
        c_free_mm=0.010,
        params=NephronParams(sec_vmax_nmol_min_per_nephron=2.0e-3),
    )
    assert high.fe_x is not None and low.fe_x is not None
    assert 0.0 < high.fe_x < low.fe_x
    assert 0.0 < high.cl_r_ml_min < low.cl_r_ml_min
