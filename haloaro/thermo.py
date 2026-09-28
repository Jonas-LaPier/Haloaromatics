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
            (only if config.RUN_RA_TS)

Saveant concerted dissociative electron transfer (det_table), from the neutral parent:
    ArX + e-(E) -> Ar. + X-          E°_DET = E_1e above
    D       = H(Ar.) + H(X.) - H(ArX)            (or G for BDFE; X. includes spin-orbit)
    lambda0 = N_A e^2/(8 pi eps0 a) (1/eps_op - 1/eps_s)    a = vdW sphere radius of ArX
    dG0‡    = (D + lambda0) / 4                  intrinsic barrier (at E = E°_DET)
    dG°(E)  = F (E - E°_DET)
    dG‡(E)  = dG0‡ (1 + dG°(E) / (4 dG0‡))^2      alpha(E) = 0.5 (1 + dG°(E) / (4 dG0‡))
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


def _val(tab, level, name, key="G"):
    """G or H (kcal/mol) with std-state (G, solvated) and spin-orbit (X. atoms) corrections."""
    row = tab.get((level, name))
    if not row or row.get(key) in (None, "") or row.get("status") not in ("ok", "warn"):
        return None
    v = float(row[key]) * K
    if key == "G" and _solvated(level):
        v += C.STD_STATE_CORR_KCAL
    if row.get("kind") == "halogen_atom":
        v -= C.SPIN_ORBIT_KCAL.get(row.get("halogen"), 0.0)
    return v


def gas_counterpart(level):
    L = C.LEVELS[level]
    for k, v in C.LEVELS.items():
        if v["method"] == L["method"] and v["basis"] == L["basis"] and v["solv"] is None:
            return k
    return None


def lambda0_kcal(radius_A):
    """Marcus-Hush one-sphere outer-sphere reorganisation energy (kcal/mol)."""
    if C.LAMBDA0_KCAL is not None:
        return C.LAMBDA0_KCAL
    if radius_A in (None, ""):
        return None
    a = float(radius_A) + C.RADIUS_PROBE_A
    coul = 332.0637  # N_A e^2 / (4 pi eps0) in kcal A / mol
    lam = coul / (2 * a) * (1 / C.EPS_OPTICAL - 1 / C.EPS_STATIC)
    return 2 * lam if C.LAMBDA0_MODEL == "homogeneous" else lam


def det_table(tab, reactions):
    ge = C.G_ELECTRON_KCAL
    out = []
    for level in C.LEVELS:
        gas = gas_counterpart(level)
        for rx in reactions:
            X, Xr = f"{rx.halogen}_anion", f"{rx.halogen}_rad"

            def bond(lv, key):
                return _dg([_val(tab, lv, rx.aryl_radical, key), _val(tab, lv, Xr, key)],
                           [_val(tab, lv, rx.parent, key)])

            D_H, D_G = bond(level, "H"), bond(level, "G")
            D_H_gas = bond(gas, "H") if gas else None
            D_G_gas = bond(gas, "G") if gas else None
            row = {"level": level, "halogen": rx.halogen, "parent": rx.parent, "site": rx.site,
                   "degeneracy": rx.degeneracy, "aryl_radical": rx.aryl_radical,
                   "BDE_kcal (dH, this level)": _r(D_H), "BDFE_kcal (dG, this level)": _r(D_G),
                   "BDE_gas_kcal": _r(D_H_gas), "BDFE_gas_kcal": _r(D_G_gas)}
            dg1 = _dg([_val(tab, level, rx.aryl_radical), _val(tab, level, X)],
                      [_val(tab, level, rx.parent), ge])
            E0 = _E(dg1, 1, level)
            row["dG_DET_kcal (ArX + e- -> Ar. + X-)"] = _r(dg1)
            row["E0_DET_V"] = _r(E0, 3)
            radius = tab.get((level, rx.parent), {}).get("radius_A")
            lam = lambda0_kcal(radius) if _solvated(level) else None
            D = D_H if C.SAVEANT_D == "H" else D_G
            row["radius_A"] = radius
            row["lambda0_kcal"] = _r(lam)
            flags = []
            if _solvated(level) and D is not None and lam is not None:
                g0 = (D + lam) / 4
                row["dG0_act_kcal (intrinsic)"] = _r(g0)
                for E in C.DET_POTENTIALS_V:
                    if E0 is None:
                        continue
                    dgE = C.FARADAY_KCAL * (E - E0)
                    f = 1 + dgE / (4 * g0)
                    if f < 0:  # beyond the activationless limit of the quadratic law
                        flags.append(f"E={E}: past activationless limit")
                        f = 0.0
                    row[f"dG_act_kcal @ {E:+.2f} V"] = _r(g0 * f * f)
                    row[f"alpha @ {E:+.2f} V"] = _r(0.5 * f, 3)
            else:
                row["dG0_act_kcal (intrinsic)"] = None
                if not _solvated(level):
                    flags.append("gas level: bond energies only (no solvent reorganisation)")
            miss = [n for n in (rx.parent, rx.aryl_radical, Xr, X) if _val(tab, level, n) is None]
            row["missing_species"] = ";".join(miss)
            row["flags"] = "; ".join(flags)
            out.append(row)
    return out


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
            # a dissociated "radical anion" is an Ar...X- complex, not ArX.-: no ET/frag values
            if "RA_dissociated" in (tab.get((level, rx.radical_anion), {}).get("flags") or ""):
                g["RA"] = None
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


def pathway_table(rx_rows, ts_rows, det_rows, lam_i=None):
    """Stepwise (via ArX.-) vs concerted (Saveant) reduction barriers at each potential.

    Stepwise, at electrode/donor potential E (SMD TS levels only):
        dG_ET(E)    = F (E - E°_ET)                       ArX + e- -> ArX.-
        lambda_ET   = lambda0 + lambda_i                  outer sphere (Marcus-Hush, as in
                                                          det_table) + inner sphere (Nelsen
                                                          4-point, sp stage)
        dG‡_ET(E)   = (lambda_ET/4)(1 + dG_ET(E)/lambda_ET)^2   Marcus
        dG‡_frag    = G(TS) - G(ArX.-)                    from the ts_* stage
        dG‡_step(E) = max( dG‡_ET(E),  max(dG_ET(E), 0) + dG‡_frag )
    i.e. the rate-limiting step is either the electron transfer or the C-X cleavage of the
    radical anion (pre-equilibrium with ArX if dG_ET(E) > 0).
    Concerted:  dG‡_conc(E) from det_table.
    The pathway with the lower barrier at E is reported as 'favored'.
    """
    key = lambda r: (r["level"], r["parent"], r["site"])  # noqa: E731
    rx = {key(r): r for r in rx_rows}
    ts = {key(r): r for r in ts_rows}
    out = []
    for d in det_rows:
        level = d["level"]
        if level not in C.TS_LEVELS or not _solvated(level):
            continue
        k = key(d)
        r, t = rx.get(k, {}), ts.get(k, {})
        lam0, e_et = d.get("lambda0_kcal"), r.get("E_ET_V")
        li_level = gas_counterpart(level) if C.LAMBDA_I_FROM_GAS else level
        li = (lam_i or {}).get((li_level, d["parent"]))
        lam = None if lam0 is None else lam0 + (li or 0.0)
        frag = t.get("dG_act_kcal (TS - ArX.-)")
        ra_unbound = "RA_dissociated" in (r.get("RA_flags") or "")
        ts_ok = t.get("ts_status") in ("ok", "warn")
        row = {"level": level, "halogen": d["halogen"], "parent": d["parent"], "site": d["site"],
               "degeneracy": d["degeneracy"], "E0_ET_V (ArX/ArX.-)": e_et,
               "E0_DET_V (ArX/Ar.+X-)": d.get("E0_DET_V"), "lambda0_kcal": lam0,
               "lambda_i_kcal": _r(li), "lambda_ET_total_kcal": _r(lam),
               "dG0_act_ET_kcal (lambda_ET/4)": _r(lam / 4) if lam is not None else None,
               "dG_act_frag_kcal (TS - ArX.-)": frag if ts_ok else None,
               "dG0_act_concerted_kcal": d.get("dG0_act_kcal (intrinsic)")}
        notes = []
        if li is None and not ra_unbound:
            notes.append("lambda_i not used (sp stage missing, or RA is not a pi radical anion): "
                         "ET barrier uses lambda0 only")
        if ra_unbound:
            notes.append("radical anion unbound at this level: concerted only")
        elif not ts_ok:
            notes.append(f"stepwise TS unavailable ({t.get('ts_status', 'missing')})")
        for E in C.DET_POTENTIALS_V:
            conc = d.get(f"dG_act_kcal @ {E:+.2f} V")
            step = None
            if not ra_unbound and ts_ok and None not in (lam, e_et, frag):
                dget = C.FARADAY_KCAL * (E - e_et)
                x = max(0.0, 1 + dget / lam)
                step = max(lam / 4 * x * x, max(dget, 0.0) + frag)
            row[f"dG_act_stepwise_kcal @ {E:+.2f} V"] = _r(step)
            row[f"dG_act_concerted_kcal @ {E:+.2f} V"] = conc
            if step is None and conc is None:
                fav = None
            elif step is None:
                fav = "concerted"
            elif conc is None:
                fav = "stepwise"
            else:
                fav = "stepwise" if step < conc else "concerted"
            row[f"favored @ {E:+.2f} V"] = fav
        row["notes"] = "; ".join(notes)
        out.append(row)
    return out
