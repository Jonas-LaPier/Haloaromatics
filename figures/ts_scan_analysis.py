"""Classify the relaxed C-X scans of the radical anions (tsscan stages) and plot their profiles.

Run from the project root:
    .venv/bin/python figures/ts_scan_analysis.py [outdir]        (default results/)

For every scan log, using converged scan points only (energies are electronic, relative to the
first point, kcal/mol):
    barrier          first interior local maximum after which the energy falls by > 1.5 kcal/mol
                     (cleavage into the bent sigma complex Ar.···X-); "then uphill" if the energy
                     rises again afterwards (separation of Ar. and X-)
    smooth barrier   that maximum >= 1.5 kcal/mol and the next step falls by <= CROSS_KCAL
    ~barrierless     maximum < 1.5 kcal/mol
    state crossing   the step after the maximum falls by > CROSS_KCAL: the scan jumps
                     between electronic states (e.g. the extra electron moves to the cleaving bond),
                     so its maximum is a crossing point, not a saddle point
    uphill to end    energy rises to the last point (no cleavage within the scanned range)
    downhill         energy falls from the start (no barrier on the radical-anion surface)
Scans that were never generated (radical anion dissociated at that level) are listed as such.
Writes <outdir>/ts_scan_summary.csv and <outdir>/plots/ts_scans_<level>.png.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, ".")
import config as C  # noqa: E402
from haloaro.geomtools import distance  # noqa: E402
from haloaro.parse import scan_points  # noqa: E402
from haloaro.stages import REACTIONS, SPECIES, STAGES  # noqa: E402

CROSS_KCAL = 5.0
FLAT_KCAL = 1.5
K = C.HARTREE_TO_KCAL


def classify(e):
    """(class, index of the barrier point, tail) from converged scan energies (kcal/mol).

    The barrier is the first interior local maximum after which the energy falls by more than
    FLAT_KCAL: the barrier to C-X cleavage (into the bent sigma complex Ar.···X-). tail says what
    the energy does after the following minimum ("then uphill" = separating Ar. and X-)."""
    n = len(e)
    if n < 3:
        return "too few points", None, ""
    for i in range(1, n - 1):
        if e[i] >= e[i - 1] and e[i] > e[i + 1] and e[i] - min(e[i + 1:]) > FLAT_KCAL:
            j = i + 1 + min(range(n - i - 1), key=lambda k: e[i + 1 + k])
            tail = "then uphill" if e[-1] - e[j] > FLAT_KCAL else ""
            if e[i] < FLAT_KCAL:
                return "~barrierless", i, tail
            cls = "state crossing" if e[i] - e[i + 1] > CROSS_KCAL else "smooth barrier"
            return cls, i, tail
    k = max(range(n), key=e.__getitem__)
    if e[k] < FLAT_KCAL:
        return "~barrierless", k, ""
    if k == 0:
        return "downhill", k, ""
    return "uphill to end", k, ""


def ra_states():
    out = {}
    for lv in C.TS_LEVELS:
        f = Path("results/raw") / f"{lv}.csv"
        if f.exists():
            for r in csv.DictReader(open(f)):
                if r["kind"] == "radical_anion":
                    out[(lv, r["name"])] = (r.get("RA_state"), r.get("RA_longest_CX_A"), r.get("flags"))
    return out


def analyse():
    states = ra_states()
    rows, profiles = [], {}
    for lv in C.TS_LEVELS:
        st = STAGES[f"tsscan_{lv}"]
        for rx in REACTIONS:
            name = rx.ts
            ra = states.get((lv, rx.radical_anion), (None, None, ""))
            base = {"level": lv, "parent": rx.parent, "site": rx.site, "scan": name, "halogen": rx.halogen,
                    "RA_state": ra[0], "RA_longest_CX_A": ra[1]}
            lp = st.log(name)
            if not lp.exists():
                rows.append({**base, "class": "no scan (radical anion dissociated)" if ra[0] == "dissociated"
                             else "no scan"})
                continue
            txt = lp.read_text(errors="ignore")
            pts = [p for p in scan_points(txt.splitlines()) if p["energy"] is not None]
            conv = [p for p in pts if p["converged"]]
            if not conv:
                rows.append({**base, "class": "no converged points"})
                continue
            ci, xi = SPECIES[name].atom_index(rx.site)
            e = [(p["energy"] - conv[0]["energy"]) * K for p in conv]
            r = [distance(p["geometry"], ci, xi) for p in conv]
            cls, k, tail = classify(e)
            end = ("normal" if "Normal termination" in txt else
                   "error" if "Error termination" in txt else "incomplete")
            drops = [e[i] - e[i + 1] for i in range(len(e) - 1)]
            j = max(range(len(drops)), key=drops.__getitem__) if drops else None
            rows.append({**base, "class": cls, "after_barrier": tail, "termination": end, "points": len(pts), "unconverged": len(pts) - len(conv),
                         "r0_A": round(r[0], 3),
                         "E_barrier_kcal": round(e[k], 2) if k is not None else None,
                         "r_at_barrier_A": round(r[k], 3) if k is not None else None,
                         "E_min_after_kcal": round(min(e[k:]), 2) if k is not None else None,
                         "largest_drop_kcal": round(drops[j], 2) if drops else None,
                         "r_drop_A": round(r[j], 3) if drops else None,
                         "E_last_kcal": round(e[-1], 2), "r_last_A": round(r[-1], 3)})
            profiles[(lv, name)] = (r, e, cls)
    return rows, profiles


def plot(profiles, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from haloaro import plotting as P
    except ImportError:
        return []
    P.use_style()
    files = []
    for lv in C.TS_LEVELS:
        items = sorted((k, v) for k, v in profiles.items() if k[0] == lv)
        if not items:
            continue
        ncol = 6
        nrow = math.ceil(len(items) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(12, 1.9 * nrow), squeeze=False)
        for ax, ((_, name), (r, e, cls)) in zip(axes.flat, items):
            ax.plot(r, e, "o-", ms=3, lw=0.8, color="black", mfc="#a6a6a6", mec="#a6a6a6")
            ax.set_title(f"{name.replace('TS_', '')}\n{cls[0] if isinstance(cls, tuple) else cls}", fontsize=7)
            ax.tick_params(labelsize=6)
        for ax in list(axes.flat)[len(items):]:
            ax.axis("off")
        fig.supxlabel("r(C–X) (Å)", fontsize=9)
        fig.supylabel("ΔE relative to the radical anion (kcal/mol)", fontsize=9)
        fig.tight_layout()
        f = Path(outdir) / "plots" / f"ts_scans_{lv}.png"
        f.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(f, dpi=200, facecolor="white")
        plt.close(fig)
        files.append(str(f))
    return files


def main(outdir="results"):
    rows, profiles = analyse()
    out = Path(outdir) / "ts_scan_summary.csv"
    keys = {}
    for r in rows:
        keys.update(dict.fromkeys(r))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(keys))
        w.writeheader()
        w.writerows(rows)
    figs = plot(profiles, outdir)
    print(f"wrote {out} ({len(rows)} scans) and {len(figs)} figures")
    return rows


if __name__ == "__main__":
    main(*sys.argv[1:])
