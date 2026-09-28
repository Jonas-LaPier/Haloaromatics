"""Correlate measured electrochemical rate constants with computed descriptors.

For every level of theory, each numeric descriptor x of the parent compounds is regressed as

    ln(kobs / h^-1) = slope * x + intercept

Two fits are reported per descriptor:
    fit_set = <CORR_FIT_GROUP> (default: the bromobenzenes) - the training set; the fitted
              line is then used to predict the remaining compounds (e.g., BDE-47, BDE-99)
    fit_set = all           - every compound with data

Site-resolved quantities are reduced to one value per compound in two ways:
    min_* / max_*  the most favourable site
    eff_*          statistically weighted effective barrier over all sites,
                   dG_eff = -RT ln sum_i g_i exp(-dG_i/RT)   (g_i = degeneracy)
For barriers, transition-state theory predicts slope = -1/RT (= -1.69 per kcal/mol at
298 K); `slope_x_RT` reports slope * (-RT) so a value of 1 means ideal TST behaviour.

Potential-dependent barriers are evaluated at config.EXP_POTENTIAL_V (the potential at which
the rate constants were measured).
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import config as C

RT = 0.0019872036 * C.TEMPERATURE


def _f(x):
    if x in (None, ""):
        return None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------- #
# statistics (no numpy/scipy dependency)
# --------------------------------------------------------------------------- #
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > FPMIN else FPMIN)
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > FPMIN else FPMIN)
        c = 1 + aa / c; c = c if abs(c) > FPMIN else FPMIN
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > FPMIN else FPMIN)
        c = 1 + aa / c; c = c if abs(c) > FPMIN else FPMIN
        de = d * c
        h *= de
        if abs(de - 1) < EPS:
            break
    return h


def _betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def t_pvalue(t, df):
    """Two-sided p-value of Student's t."""
    if df <= 0:
        return None
    return _betai(df / 2, 0.5, df / (df + t * t))


def linfit(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0:
        return None
    slope = sxy / sxx
    icpt = my - slope * mx
    r = sxy / math.sqrt(sxx * syy)
    df = n - 2
    t = r * math.sqrt(df / max(1e-300, 1 - r * r)) if abs(r) < 1 else float("inf")
    res = [y - (slope * x + icpt) for x, y in zip(xs, ys)]
    rmse = math.sqrt(sum(e * e for e in res) / n)
    # leave-one-out
    loo = []
    for i in range(n):
        xi = xs[:i] + xs[i + 1:]
        yi = ys[:i] + ys[i + 1:]
        mxi, myi = sum(xi) / (n - 1), sum(yi) / (n - 1)
        sxxi = sum((x - mxi) ** 2 for x in xi)
        if sxxi == 0:
            return None
        b = sum((x - mxi) * (y - myi) for x, y in zip(xi, yi)) / sxxi
        loo.append(ys[i] - (b * xs[i] + myi - b * mxi))
    press = sum(e * e for e in loo)
    return {"n": n, "slope": slope, "intercept": icpt, "r": r, "R2": r * r,
            "p": t_pvalue(t, df) if math.isfinite(t) else 0.0, "RMSE": rmse,
            "LOO_RMSE": math.sqrt(press / n), "Q2": 1 - press / syy}


def mlr2(x1, x2, ys):
    """y = b0 + b1 x1 + b2 x2 by normal equations."""
    n = len(ys)
    if n < 4:
        return None
    X = [[1.0, a, b] for a, b in zip(x1, x2)]
    A = [[sum(X[k][i] * X[k][j] for k in range(n)) for j in range(3)] for i in range(3)]
    v = [sum(X[k][i] * ys[k] for k in range(n)) for i in range(3)]
    M = [row[:] + [vv] for row, vv in zip(A, v)]
    for i in range(3):
        piv = max(range(i, 3), key=lambda r: abs(M[r][i]))
        if abs(M[piv][i]) < 1e-12:
            return None
        M[i], M[piv] = M[piv], M[i]
        for r in range(3):
            if r != i:
                fct = M[r][i] / M[i][i]
                M[r] = [a - fct * b for a, b in zip(M[r], M[i])]
    b = [M[i][3] / M[i][i] for i in range(3)]
    pred = [b[0] + b[1] * a + b[2] * c for a, c in zip(x1, x2)]
    my = sum(ys) / n
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - p) ** 2 for y, p in zip(ys, pred))
    R2 = 1 - ss_res / ss_tot if ss_tot else None
    adj = 1 - (1 - R2) * (n - 1) / (n - 3) if R2 is not None and n > 3 else None
    return {"n": n, "b0": b[0], "b1": b[1], "b2": b[2], "R2": R2, "adj_R2": adj,
            "RMSE": math.sqrt(ss_res / n)}


# --------------------------------------------------------------------------- #
# descriptor table
# --------------------------------------------------------------------------- #
def _eff(values):
    """[(degeneracy, barrier)] -> -RT ln sum g exp(-b/RT)."""
    vals = [(g, b) for g, b in values if b is not None]
    if not vals:
        return None
    b0 = min(b for _, b in vals)
    return b0 - RT * math.log(sum(g * math.exp(-(b - b0) / RT) for g, b in vals))


MOL_COLS = ["n_X", "HOMO_eV", "LUMO_eV", "gap_eV", "mu_Koopmans_eV", "eta_Koopmans_eV",
            "omega_Koopmans_eV", "VEA_eV", "AEA_elec_eV", "VIE_eV", "mu_dSCF_eV", "eta_dSCF_eV",
            "omega_dSCF_eV", "lambda_i_kcal", "dipole_D", "polar_iso_A3", "radius_A",
            "dG_ET_kcal", "E_ET_V"]


def descriptor_table(mol_rows, rx_rows, det_rows, path_rows, ts_rows, site_rows):
    """{(level, parent): {descriptor: value}}."""
    E = C.EXP_POTENTIAL_V
    tag = f"@ {E:+.2f} V"
    D = {}

    def put(level, parent, k, v):
        if v is not None:
            D.setdefault((level, parent), {})[k] = v

    for r in mol_rows:
        for c in MOL_COLS:
            put(r["level"], r["parent"], c, _f(r.get(c)))

    def grouped(rows):
        g = {}
        for r in rows:
            g.setdefault((r["level"], r["parent"]), []).append(r)
        return g

    for (lv, par), rs in grouped(rx_rows).items():
        for c in ("dG_ET_any_kcal",):
            put(lv, par, c, _f(rs[0].get(c)))
        for c in ("dG_1e_kcal", "dG_frag_kcal", "dG_frag_any_kcal", "dG_2e_carbanion_kcal",
                  "dG_2e_HDH_kcal", "dG_rad_red_kcal"):
            vals = [_f(r.get(c)) for r in rs if _f(r.get(c)) is not None]
            put(lv, par, f"min_{c}", min(vals) if vals else None)
            put(lv, par, f"eff_{c}", _eff([(r["degeneracy"], _f(r.get(c))) for r in rs]))
    for (lv, par), rs in grouped(det_rows).items():
        for c, out in (("BDE_kcal (dH, this level)", "min_BDE_kcal"),
                       ("dG0_act_kcal (intrinsic)", "min_dG0_act_concerted_kcal")):
            vals = [_f(r.get(c)) for r in rs if _f(r.get(c)) is not None]
            put(lv, par, out, min(vals) if vals else None)
        put(lv, par, f"eff_dG_act_concerted_kcal {tag}",
            _eff([(r["degeneracy"], _f(r.get(f"dG_act_kcal {tag}"))) for r in rs]))
    for (lv, par), rs in grouped(path_rows).items():
        put(lv, par, f"eff_dG_act_stepwise_kcal {tag}",
            _eff([(r["degeneracy"], _f(r.get(f"dG_act_stepwise_kcal {tag}"))) for r in rs]))
        comb = []
        for r in rs:
            bs = [b for b in (_f(r.get(f"dG_act_stepwise_kcal {tag}")),
                              _f(r.get(f"dG_act_concerted_kcal {tag}"))) if b is not None]
            if bs:
                comb.append((r["degeneracy"], _eff([(1, b) for b in bs])))
        put(lv, par, f"eff_dG_act_combined_kcal {tag}", _eff(comb))
    for (lv, par), rs in grouped(ts_rows).items():
        put(lv, par, "eff_dG_act_frag_TS_kcal",
            _eff([(r["degeneracy"], _f(r.get("dG_act_kcal (TS - ArX.-)"))) for r in rs]))
    for (lv, par), rs in grouped(site_rows).items():
        for c, fn in (("fplus_CX", max), ("spin_vA_X", max), ("spin_RA_X", max),
                      ("dr_CX_RA_A", max), ("wiberg_CX_parent", min), ("qCM5_X", max),
                      ("pKa_ArH_site", min)):
            vals = [_f(r.get(c)) for r in rs if _f(r.get(c)) is not None]
            put(lv, par, f"{'max' if fn is max else 'min'}_{c}", fn(vals) if vals else None)
    return D


def load_experimental(path=None):
    p = Path(path or C.EXPERIMENTAL_KOBS)
    if not p.exists():
        return []
    rows = []
    with p.open(newline="") as f:
        for r in csv.DictReader(f):
            k = _f(r["kobs_per_h"])
            if k and k > 0:
                r["ln_k"] = math.log(k)
                rows.append(r)
    return rows


def run(mol_rows, rx_rows, det_rows, path_rows, ts_rows, site_rows):
    """Returns (correlation rows, data rows, pair rows, predictions)."""
    exp = load_experimental()
    if not exp:
        return [], [], [], []
    D = descriptor_table(mol_rows, rx_rows, det_rows, path_rows, ts_rows, site_rows)
    levels = list(C.LEVELS)
    data_rows = []
    for lv in levels:
        for e in exp:
            d = D.get((lv, e["species"]), {})
            data_rows.append({"level": lv, "compound": e["compound"], "species": e["species"],
                              "group": e["group"], "kobs_per_h": e["kobs_per_h"],
                              "ln_kobs": round(e["ln_k"], 4), **{k: round(v, 5) for k, v in d.items()}})
    descs = sorted({k for d in D.values() for k in d})
    out, preds = [], []
    for lv in levels:
        for desc in descs:
            for fit_set in (C.CORR_FIT_GROUP, "all"):
                pts = [(D.get((lv, e["species"]), {}).get(desc), e["ln_k"], e) for e in exp
                       if fit_set == "all" or e["group"] == fit_set]
                pts = [p for p in pts if p[0] is not None]
                fit = linfit([p[0] for p in pts], [p[1] for p in pts])
                if not fit:
                    continue
                row = {"level": lv, "descriptor": desc, "fit_set": fit_set,
                       **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in fit.items()},
                       "slope_x_RT": round(-fit["slope"] * RT, 3) if "act" in desc else None,
                       "compounds": ";".join(p[2]["species"] for p in pts)}
                if fit_set != "all":
                    others = [e for e in exp if e["group"] != fit_set]
                    errs = []
                    for e in others:
                        x = D.get((lv, e["species"]), {}).get(desc)
                        if x is None:
                            continue
                        pred = fit["slope"] * x + fit["intercept"]
                        errs.append(pred - e["ln_k"])
                        preds.append({"level": lv, "descriptor": desc, "compound": e["compound"],
                                      "ln_kobs": round(e["ln_k"], 3), "ln_k_pred": round(pred, 3),
                                      "k_pred_per_h": round(math.exp(pred), 5),
                                      "residual_ln": round(pred - e["ln_k"], 3)})
                    row["pred_RMSE_other_groups"] = round(math.sqrt(sum(x * x for x in errs) / len(errs)), 3) if errs else None
                out.append(row)
    out.sort(key=lambda r: (r["level"], r["fit_set"] != C.CORR_FIT_GROUP, -(r["R2"] or 0)))

    pair_rows = []
    for lv in levels:
        for a, b in C.CORR_PAIRS:
            for fit_set in (C.CORR_FIT_GROUP, "all"):
                pts = [(D.get((lv, e["species"]), {}).get(a), D.get((lv, e["species"]), {}).get(b), e["ln_k"])
                       for e in exp if fit_set == "all" or e["group"] == fit_set]
                pts = [p for p in pts if p[0] is not None and p[1] is not None]
                fit = mlr2([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts])
                if fit:
                    pair_rows.append({"level": lv, "x1": a, "x2": b, "fit_set": fit_set,
                                      **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in fit.items()}})
    return out, data_rows, pair_rows, preds


def plots(corr_rows, data_rows, outdir, top=6):
    """Scatter plots of the best descriptors per level (needs matplotlib; skipped otherwise)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = []
    for lv in C.LEVELS:
        best = [r for r in corr_rows if r["level"] == lv and r["fit_set"] == C.CORR_FIT_GROUP][:top]
        if not best:
            continue
        rows = [d for d in data_rows if d["level"] == lv]
        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        for ax, r in zip(axes.flat, best):
            for grp, mk in ((C.CORR_FIT_GROUP, "o"), (None, "s")):
                pts = [(d.get(r["descriptor"]), d["ln_kobs"], d["compound"]) for d in rows
                       if (d["group"] == grp if grp else d["group"] != C.CORR_FIT_GROUP)
                       and d.get(r["descriptor"]) is not None]
                if pts:
                    ax.scatter([p[0] for p in pts], [p[1] for p in pts], marker=mk,
                               facecolor="k" if grp else "none", edgecolor="k")
                    for x, y, lab in pts:
                        ax.annotate(lab.replace("Tribromobenzene", "TBB").replace("Dibromobenzene", "DBB")
                                    .replace("Tetrabromobenzene", "TeBB").replace("Hexabromobenzene", "HBB")
                                    .replace("Bromobenzene", "BB"), (x, y), fontsize=7,
                                    xytext=(3, 3), textcoords="offset points")
            xs = [d[r["descriptor"]] for d in rows if d.get(r["descriptor"]) is not None]
            if xs:
                x0, x1 = min(xs), max(xs)
                ax.plot([x0, x1], [r["slope"] * x0 + r["intercept"], r["slope"] * x1 + r["intercept"]], "k--", lw=1)
            ax.set_xlabel(r["descriptor"], fontsize=8)
            ax.set_ylabel("ln(k_obs / h$^{-1}$)", fontsize=8)
            ax.set_title(f"R$^2$ = {r['R2']:.2f}, n = {r['n']}", fontsize=9)
        for ax in list(axes.flat)[len(best):]:
            ax.axis("off")
        fig.suptitle(f"{lv}: top descriptors (filled: fit set, open: predicted)", fontsize=11)
        fig.tight_layout()
        f = outdir / f"correlations_{lv}.png"
        fig.savefig(f, dpi=150)
        plt.close(fig)
        files.append(str(f))
    return files
