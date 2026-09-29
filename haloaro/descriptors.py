"""Molecular and site-specific descriptors, and predicted product distributions.

Inputs are the scraped tables: `tab[(level, name)]` (scalar rows from results/raw/*.csv)
and `atoms[(level, name)]` (geometry, Hirshfeld/CM5 charges and spins, Wiberg matrix from
results/raw/*_atoms.json). Single-point job names are <species>__<role> (see stages.py).

Molecular (one row per level x parent)
    Koopmans:  mu = (HOMO + LUMO)/2, eta = LUMO - HOMO, omega = mu^2 / (2 eta)
    Delta-SCF (sp stage, electronic energies):
        VEA = E(N@N) - E(A@N)      VIE = E(C@N) - E(N@N)      AEA = E(N@N) - E(A@A)
        mu = -(VIE + VEA)/2, eta = VIE - VEA, omega = mu^2 / (2 eta)
    Nelsen four-point inner-sphere reorganisation energy (ArX / ArX.- couple):
        lambda_N = E(N@A) - E(N@N),  lambda_A = E(A@N) - E(A@A),  lambda_i = lambda_N + lambda_A
    N = neutral, A = radical anion, C = cation; X@Y = state X at the optimised geometry of Y.

Site (one row per level x symmetry-unique C-X bond)
    structure: halogens ortho / meta / para to the site
    parent:    r(C-X), CM5 and Hirshfeld charges on C and X, Wiberg C-X bond index
    vertical anion at the parent geometry: condensed Fukui f+ = q_H(N) - q_H(N+1) on C, X
               and C+X; Hirshfeld spin density on C and X
    relaxed radical anion: r(C-X), elongation vs parent, out-of-plane angle of X, spin density
               and CM5 charge on C and X, Wiberg C-X (alpha + beta)
    energetics: dG_1e, dG_frag, dG_2e (carbanion, HDH), BDE, Saveant intrinsic barrier,
               stepwise TS barrier, pKa of the ArH C-H formed at that site
    rank_*:    rank within the parent (1 = most reactive by that descriptor)

Product distribution (one row per level x parent x site)
    fraction of mono-dehalogenation at each site, from Boltzmann-weighted rates
    (degeneracy x exp(-dG‡/RT)) for the stepwise TS, the concerted Saveant barrier and the
    combined (stepwise + concerted) rate at each potential in config.DET_POTENTIALS_V.
"""
from __future__ import annotations

import math

import config as C
from . import thermo
from .geomtools import distance, oop_angle

K = C.HARTREE_TO_KCAL
EV = C.HARTREE_TO_EV
RT = 0.0019872036 * C.TEMPERATURE          # kcal/mol
LN10 = math.log(10)


def _f(x):
    if x in (None, ""):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _r(x, nd=3):
    return None if x is None else round(x, nd)


def _usable(row):
    return bool(row) and row.get("status") in ("ok", "warn")


def _E(tab, level, name):
    row = tab.get((level, name))
    return _f(row.get("E_scf")) if _usable(row) else None


def _atoms(atoms, level, name):
    a = atoms.get((level, name))
    return a if a and a.get("status") in ("ok", "warn") else None


def _diff(a, b, scale=1.0):
    return None if a is None or b is None else (a - b) * scale


# --------------------------------------------------------------------------- #
# Molecular descriptors
# --------------------------------------------------------------------------- #
def molecular_table(tab, species, rx_rows):
    rx_by = {}
    for r in rx_rows:
        rx_by.setdefault((r["level"], r["parent"]), []).append(r)
    out, lam_i = [], {}
    for level in C.LEVELS:
        for name, sp in species.items():
            if sp.kind != "parent":
                continue
            row = tab.get((level, name), {})
            ra = f"{name}_RA"
            ra_row = tab.get((level, ra), {})
            H, L = _f(row.get("homo")), _f(row.get("lumo"))
            d = {"level": level, "parent": name, "halogen": sp.halogen or "", "n_X": sp.n_hal,
                 "status": row.get("status", "missing")}
            if H is not None and L is not None:
                mu, eta = (H + L) / 2 * EV, (L - H) * EV
                d.update({"HOMO_eV": _r(H * EV), "LUMO_eV": _r(L * EV), "gap_eV": _r(eta),
                          "mu_Koopmans_eV": _r(mu), "eta_Koopmans_eV": _r(eta),
                          "omega_Koopmans_eV": _r(mu * mu / (2 * eta)) if eta else None})
            d["dipole_D"] = _f(row.get("dipole_D"))
            d["polar_iso_A3"] = _f(row.get("polar_iso_A3"))
            d["radius_A"] = _f(row.get("radius_A"))

            eNN = _E(tab, level, f"{name}__pop") or _E(tab, level, name)
            eAN = _E(tab, level, f"{name}__vA")
            eCN = _E(tab, level, f"{name}__vC")
            ra_ok = _usable(ra_row) and "RA_dissociated" not in (ra_row.get("flags") or "")
            eAA = _E(tab, level, ra) if ra_ok else None
            eNA = _E(tab, level, f"{ra}__vN") if ra_ok else None
            vea, vie, aea = _diff(eNN, eAN, EV), _diff(eCN, eNN, EV), _diff(eNN, eAA, EV)
            d.update({"VEA_eV": _r(vea), "AEA_elec_eV": _r(aea), "VIE_eV": _r(vie)})
            if vea is not None and vie is not None:
                mu, eta = -(vie + vea) / 2, vie - vea
                d.update({"mu_dSCF_eV": _r(mu), "eta_dSCF_eV": _r(eta),
                          "omega_dSCF_eV": _r(mu * mu / (2 * eta)) if eta else None})
            lN, lA = _diff(eNA, eNN, K), _diff(eAN, eAA, K)
            li = None if lN is None or lA is None else lN + lA
            d.update({"lambda_N_kcal": _r(lN, 2), "lambda_A_kcal": _r(lA, 2),
                      "lambda_i_kcal": _r(li, 2), "lambda_i_eV": _r(li / 23.0605 if li is not None else None)})
            # lambda_i is an electron-transfer reorganisation energy only for an intact pi
            # radical anion; for a bent sigma-type RA it mostly measures C-X cleavage
            d["RA_state"] = ra_row.get("RA_state")
            lam_i[(level, name)] = li if ra_row.get("RA_state") == "pi" else None
            d["RA_bound"] = None if not ra_row else ra_ok
            d["RA_S2"] = _f(ra_row.get("S2"))

            sites = rx_by.get((level, name), [])
            if sites:
                d["dG_ET_kcal"] = sites[0].get("dG_ET_kcal")
                d["E_ET_V"] = sites[0].get("E_ET_V")
                for col, out_col in (("dG_1e_kcal", "min_dG_1e_kcal"),
                                     ("dG_2e_carbanion_kcal", "min_dG_2e_carbanion_kcal"),
                                     ("dG_2e_HDH_kcal", "min_dG_2e_HDH_kcal")):
                    vals = [s[col] for s in sites if s.get(col) is not None]
                    d[out_col] = min(vals) if vals else None
                d["n_unique_sites"] = len(sites)
            out.append(d)
    return out, lam_i


# --------------------------------------------------------------------------- #
# Site descriptors
# --------------------------------------------------------------------------- #
def site_table(tab, atoms, species, reactions, rx_rows, ts_rows, det_rows):
    key = lambda r: (r["level"], r["parent"], r["site"])  # noqa: E731
    rx = {key(r): r for r in rx_rows}
    ts = {key(r): r for r in ts_rows}
    det = {key(r): r for r in det_rows}
    out = []
    for level in C.LEVELS:
        solv = C.LEVELS[level]["solv"] is not None
        for r in reactions:
            P = species[r.parent]
            ci, xi = P.atom_index(r.site)
            o, m, p = P.neighbours(r.site, r.halogen)
            ring = P.ring_atoms(r.site)
            d = {"level": level, "halogen": r.halogen, "parent": r.parent, "n_X": P.n_hal,
                 "site": r.site, "degeneracy": r.degeneracy, "aryl_radical": r.aryl_radical,
                 "ArH_product": r.hydro_product, "n_ortho_X": o, "n_meta_X": m, "n_para_X": p,
                 "ortho_to_ether": (str(r.site).rstrip("p") in ("2", "6")) if P.skeleton == "dpe" else None}

            ga = _atoms(atoms, level, r.parent)
            if ga:
                d["r_CX_parent_A"] = _r(distance(ga["geometry"], ci, xi), 4)
            pop = _atoms(atoms, level, f"{r.parent}__pop")
            hN = pop.get("hirshfeld") if pop else None
            if hN:
                d.update({"qCM5_C": hN[ci - 1]["q_CM5"], "qCM5_X": hN[xi - 1]["q_CM5"],
                          "qH_C": hN[ci - 1]["q_H"], "qH_X": hN[xi - 1]["q_H"]})
            if pop and pop.get("wiberg"):
                d["wiberg_CX_parent"] = _r(pop["wiberg"][ci - 1][xi - 1], 4)
            va = _atoms(atoms, level, f"{r.parent}__vA")
            hA = va.get("hirshfeld") if va else None
            if hN and hA:
                fC = hN[ci - 1]["q_H"] - hA[ci - 1]["q_H"]
                fX = hN[xi - 1]["q_H"] - hA[xi - 1]["q_H"]
                d.update({"fplus_C": _r(fC, 4), "fplus_X": _r(fX, 4), "fplus_CX": _r(fC + fX, 4)})
            if hA:
                d.update({"spin_vA_C": hA[ci - 1]["s_H"], "spin_vA_X": hA[xi - 1]["s_H"]})

            ra_row = tab.get((level, r.radical_anion), {})
            bound = _usable(ra_row) and "RA_dissociated" not in (ra_row.get("flags") or "")
            d["RA_bound"] = bound if ra_row else None
            gr = _atoms(atoms, level, r.radical_anion)
            if gr:
                rr = distance(gr["geometry"], ci, xi)
                d["r_CX_RA_A"] = _r(rr, 4)
                if "r_CX_parent_A" in d:
                    d["dr_CX_RA_A"] = _r(rr - d["r_CX_parent_A"], 4)
                d["oop_X_RA_deg"] = _r(oop_angle(gr["geometry"], ci, xi, ring), 2)
            rp = _atoms(atoms, level, f"{r.radical_anion}__pop")
            hR = rp.get("hirshfeld") if rp else None
            if hR:
                d.update({"spin_RA_C": hR[ci - 1]["s_H"], "spin_RA_X": hR[xi - 1]["s_H"],
                          "qCM5_C_RA": hR[ci - 1]["q_CM5"], "qCM5_X_RA": hR[xi - 1]["q_CM5"]})
            if rp and rp.get("wiberg"):
                d["wiberg_CX_RA"] = _r(rp["wiberg"][ci - 1][xi - 1], 4)

            k = (level, r.parent, r.site)
            x = rx.get(k, {})
            for col in ("dG_1e_kcal", "E_1e_V", "dG_frag_kcal", "dG_2e_carbanion_kcal",
                        "dG_2e_HDH_kcal", "dG_rad_red_kcal"):
                d[col] = x.get(col)
            y = det.get(k, {})
            d["BDE_kcal"] = y.get("BDE_kcal (dH, this level)")
            d["dG0_act_concerted_kcal"] = y.get("dG0_act_kcal (intrinsic)")
            t = ts.get((f"{level}", r.parent, r.site), {})
            if t.get("ts_status") in ("ok", "warn"):
                d["dG_act_stepwise_frag_kcal"] = t.get("dG_act_kcal (TS - ArX.-)")
            if solv:
                gA = thermo._val(tab, level, r.aryl_anion)
                gH = thermo._val(tab, level, r.hydro_product)
                if gA is not None and gH is not None:
                    d["pKa_ArH_site"] = _r((gA + thermo._proton(level) - gH) / (RT * LN10), 2)
            out.append(d)

    # ranks within each parent (1 = most reactive / most favourable)
    rank_spec = [("dG_1e_kcal", False), ("dG_2e_carbanion_kcal", False), ("dG_2e_HDH_kcal", False),
                 ("BDE_kcal", False), ("dG0_act_concerted_kcal", False),
                 ("dG_act_stepwise_frag_kcal", False), ("fplus_CX", True), ("spin_RA_X", True),
                 ("dr_CX_RA_A", True), ("wiberg_CX_parent", False)]
    groups = {}
    for d in out:
        groups.setdefault((d["level"], d["parent"]), []).append(d)
    for g in groups.values():
        for col, high_first in rank_spec:
            vals = [(float(d[col]), d) for d in g if _f(d.get(col)) is not None]
            vals.sort(key=lambda t: -t[0] if high_first else t[0])
            for i, (_, d) in enumerate(vals, 1):
                d[f"rank_{col}"] = i
    return out


# --------------------------------------------------------------------------- #
# Product distribution
# --------------------------------------------------------------------------- #
def _fractions(items):
    """items: list of (key, degeneracy, barrier kcal or None) -> {key: fraction}."""
    ok = [(k, g, b) for k, g, b in items if b is not None]
    if not ok:
        return {}
    bmin = min(b for _, _, b in ok)
    w = {k: g * math.exp(-(b - bmin) / RT) for k, g, b in ok}
    tot = sum(w.values())
    return {k: v / tot for k, v in w.items()}


def product_distribution(site_rows, det_rows, path_rows):
    key = lambda r: (r["level"], r["parent"], r["site"])  # noqa: E731
    det = {key(r): r for r in det_rows}
    path = {key(r): r for r in path_rows}
    groups = {}
    for d in site_rows:
        groups.setdefault((d["level"], d["parent"]), []).append(d)
    out = []
    for (level, parent), g in groups.items():
        if len(g) < 1:
            continue
        cols = {}
        cols["frac_stepwise_TS"] = _fractions(
            [(d["site"], d["degeneracy"], _f(d.get("dG_act_stepwise_frag_kcal"))) for d in g])
        for E in C.DET_POTENTIALS_V:
            tag = f"@ {E:+.2f} V"
            cols[f"frac_concerted {tag}"] = _fractions(
                [(d["site"], d["degeneracy"], _f(det.get(key(d), {}).get(f"dG_act_kcal {tag}"))) for d in g])
            # combined: rate_i = exp(-step_i/RT) + exp(-conc_i/RT)
            items = []
            for d in g:
                p = path.get(key(d), {})
                bs = [_f(p.get(f"dG_act_stepwise_kcal {tag}")), _f(p.get(f"dG_act_concerted_kcal {tag}"))]
                bs = [b for b in bs if b is not None]
                if bs:
                    b0 = min(bs)
                    eff = b0 - RT * math.log(sum(math.exp(-(b - b0) / RT) for b in bs))
                    items.append((d["site"], d["degeneracy"], eff))
            cols[f"frac_combined {tag}"] = _fractions(items)
        for d in g:
            row = {"level": level, "parent": parent, "site": d["site"], "degeneracy": d["degeneracy"],
                   "aryl_radical": d["aryl_radical"], "ArH_product": d["ArH_product"]}
            for c, fr in cols.items():
                row[c] = _r(fr.get(d["site"]), 4) if fr else None
            out.append(row)
    return out
