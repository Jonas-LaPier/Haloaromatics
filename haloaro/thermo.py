"""Reaction free energies, reduction potentials and TS barriers from scraped data.

Conventions (all per mole, T = 298.15 K)
    G(e-)   = config.G_ELECTRON_KCAL (gas phase, Fermi-Dirac)
    SMD levels: every solute gets +config.STD_STATE_CORR_KCAL (1 atm -> 1 M);
                G(H+, aq, 1 M) = G_gas + dG_solv + std-state corr
    Gas levels: G(H+, g, 1 atm) = config.G_PROTON_GAS_KCAL; no potentials reported
    E (V vs SHE) = -dG / (n F) - E_abs(SHE)

Reactions (one row per symmetry-unique C-X bond, degeneracy listed):
    ET      ArX + e-         -> ArX.-
    1e      ArX + e-         -> Ar.  + X-        (dissociative 1e reduction)
    frag    ArX.-            -> Ar.  + X-
    rad     Ar. + e-         -> Ar-
    2e_an   ArX + 2e-        -> Ar-  + X-        (to carbanion)
    2e_HDH  ArX + H+ + 2e-   -> ArH  + X-        (hydrodehalogenation)
    TS      dG‡ = G(TS) - G(ArX.-);  dG‡(from ArX + e-) = G(TS) - G(ArX) - G(e-)
"""
from __future__ import annotations

import config as C

K = C.HARTREE_TO_KCAL


def _solvated(level):
    return C.LEVELS[level]["solv"] is not None


def _G(tab, level, name):
    """G in kcal/mol incl. std-state correction, or None."""
    row = tab.get((level, name))
    if not row or row.get("G") in (None, "") or row.get("status") not in ("ok", "warn"):
        return None
    g = float(row["G"]) * K
    if _solvated(level):
        g += C.STD_STATE_CORR_KCAL
    return g


def _proton(level):
    if _solvated(level):
        return C.G_PROTON_GAS_KCAL + C.DG_SOLV_PROTON_KCAL + C.STD_STATE_CORR_KCAL
    return C.G_PROTON_GAS_KCAL


def _dg(products, reactants):
    if any(v is None for v in products + reactants):
        return None
    return sum(products) - sum(reactants)


def _E(dg, n, level):
    if dg is None or not _solvated(level):
        return None
    return -dg / (n * C.FARADAY_KCAL) - C.E_ABS_SHE_V


def _r(x, nd=2):
    return None if x is None else round(x, nd)


def reaction_table(tab, reactions, species):
    """tab: {(level, name): row dict from scrape}. Returns list of row dicts."""
    ge = C.G_ELECTRON_KCAL
    out = []
    for level in C.LEVELS:
        gH = _proton(level)
        for rx in reactions:
            g = {k: _G(tab, level, v) for k, v in dict(
                ArX=rx.parent, RA=rx.radical_anion, Ar_rad=rx.aryl_radical,
                Ar_an=rx.aryl_anion, X=f"{rx.halogen}_anion", ArH=rx.hydro_product).items()}
            prow = tab.get((level, rx.parent), {})
            lumo = prow.get("lumo_eV")
            e = {}
            e["ET"] = _dg([g["RA"]], [g["ArX"], ge])
            e["1e"] = _dg([g["Ar_rad"], g["X"]], [g["ArX"], ge])
            e["frag"] = _dg([g["Ar_rad"], g["X"]], [g["RA"]])
            e["rad"] = _dg([g["Ar_an"]], [g["Ar_rad"], ge])
            e["2e_an"] = _dg([g["Ar_an"], g["X"]], [g["ArX"], 2 * ge])
            e["2e_HDH"] = _dg([g["ArH"], g["X"]], [g["ArX"], gH, 2 * ge])
            n = {"ET": 1, "1e": 1, "rad": 1, "2e_an": 2, "2e_HDH": 2}
            missing = [species_name for key, species_name in dict(
                ArX=rx.parent, RA=rx.radical_anion, Ar_rad=rx.aryl_radical, Ar_an=rx.aryl_anion,
                X=f"{rx.halogen}_anion", ArH=rx.hydro_product).items() if g[key] is None]
            ra_flags = tab.get((level, rx.radical_anion), {}).get("flags", "")
            row = {
                "level": level, "halogen": rx.halogen, "parent": rx.parent,
                "n_X": species[rx.parent].n_hal, "site": rx.site, "degeneracy": rx.degeneracy,
                "aryl_radical": rx.aryl_radical, "aryl_anion": rx.aryl_anion,
                "ArH_product": rx.hydro_product,
                "LUMO_parent_eV": lumo,
                "dG_ET_kcal": _r(e["ET"]), "E_ET_V": _r(_E(e["ET"], 1, level), 3),
                "dG_1e_kcal": _r(e["1e"]), "E_1e_V": _r(_E(e["1e"], 1, level), 3),
                "dG_frag_kcal": _r(e["frag"]),
                "dG_rad_red_kcal": _r(e["rad"]), "E_rad_red_V": _r(_E(e["rad"], 1, level), 3),
                "dG_2e_carbanion_kcal": _r(e["2e_an"]), "E_2e_carbanion_V": _r(_E(e["2e_an"], 2, level), 3),
                "dG_2e_HDH_kcal": _r(e["2e_HDH"]), "E_2e_HDH_V": _r(_E(e["2e_HDH"], 2, level), 3),
                "RA_flags": ra_flags,
                "missing_species": ";".join(missing),
            }
            out.append(row)
    return out


def ts_table(tab, reactions):
    ge = C.G_ELECTRON_KCAL
    out = []
    for level in C.TS_LEVELS:
        tsl = f"ts_{level}"
        for rx in reactions:
            ts = tab.get((tsl, rx.ts), {})
            gts = _G({(level, rx.ts): ts} if ts else {}, level, rx.ts)
            gra = _G(tab, level, rx.radical_anion)
            garx = _G(tab, level, rx.parent)
            gprod = [_G(tab, level, rx.aryl_radical), _G(tab, level, f"{rx.halogen}_anion")]
            ets = ts.get("E_scf")
            era = tab.get((level, rx.radical_anion), {}).get("E_scf")
            dE = (float(ets) - float(era)) * K if ets not in (None, "") and era not in (None, "") else None
            out.append({
                "level": level, "halogen": rx.halogen, "parent": rx.parent, "site": rx.site,
                "degeneracy": rx.degeneracy, "ts": rx.ts,
                "dG_act_kcal (TS - ArX.-)": _r(_dg([gts], [gra])),
                "dE_act_kcal (TS - ArX.-)": _r(dE),
                "dG_act_from_ArX_kcal (TS - ArX - e-)": _r(_dg([gts], [garx, ge])),
                "dG_frag_kcal": _r(_dg(gprod, [gra])),
                "imag_freq_cm-1": ts.get("lowest_freq"),
                "ts_status": ts.get("status", "missing"),
                "ts_flags": ts.get("flags", ""),
                "RA_flags": tab.get((level, rx.radical_anion), {}).get("flags", ""),
            })
    return out
