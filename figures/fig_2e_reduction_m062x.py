"""ln kobs vs dG(ArX + 2e- -> Ar- + X-) at M06-2X/6-311++G(d), (a) gas and (b) SMD, in Jonas's style.

Run from the project root after compile:  .venv/bin/python figures/fig_2e_reduction_m062x.py figures/fig_2e_reduction_m062x
Label offsets for crowded points are set by hand (`fixed`); check them if the data change.
"""
import csv, sys
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from haloaro import plotting as P
import config as C

desc = "min_dG_2e_carbanion_kcal"
data = list(csv.DictReader(open("results/correlation_data.csv")))
q = {(r["level"], r["fit_set"]): r for r in csv.DictReader(open("results/qsar_summary.csv")) if r["descriptors"] == desc}
font = P.use_style()
yt, yd = P.nice_ticks([float(d["ln_kobs"]) for d in data], nbins=6)
yticks = (yt + [yt[-1] + (yt[1] - yt[0]) * k for k in (1, 2)], yd)   # headroom for the legend
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
panels = []
for i, (ax, lv) in enumerate(zip(axes, ("m062x_gas", "m062x_smd"))):
    pts = P._points(data, lv, desc)
    fit = q[(lv, C.CORR_FIT_GROUP)]
    drawn, line = P._panel(ax, fit, pts, desc, "dG ArX + 2e- -> Ar- + X-", yticks, ms=7)
    if lv == "m062x_smd":   # room on the right for the labels of the crowded low-k points
        P._set_ticks(ax, "x", [-190, -180, -170, -160, -150, -140, -130], 0)
    ax.text(0.03, 0.97, f"({'ab'[i]})", transform=ax.transAxes, ha="left", va="top")
    panels.append((ax, pts, line, drawn))
fig.tight_layout(w_pad=2.0)
for i, (ax, pts, line, drawn) in enumerate(panels):
    letter = [t for t in ax.texts][-1]
    leg = None
    if i == 1:   # top band of (b), above the data and fit line
        from matplotlib.lines import Line2D
        leg = ax.legend([Line2D([], [], linestyle="none", ms=6, mew=0.8, **{k: v for k, v in st.items() if k != "label"})
                         for _, st in drawn], [st["label"] for _, st in drawn], loc="upper right",
                        frameon=False, handletextpad=0.3, fontsize=9)
    fixed = ({"1,2-DiCl": (0, -10), "1,3-DiBr": (7, 0), "1,4-DiBr": (6, -8), "BDE-99": (-7, 0),
              "1,2,4,5-TetraCl": (-10, -30)} if i == 1 else {"PentaCl": (7, -10)})
    P._place_labels(fig, ax, pts, line, avoid=[leg, letter], fixed=fixed)
out = sys.argv[1]
for ext, dpi in (("png", 600), ("pdf", None)):
    fig.savefig(f"{out}.{ext}", dpi=dpi, facecolor="white")
g, s = q[("m062x_gas", "bromobenzene")], q[("m062x_smd", "bromobenzene")]
ga, sa = q[("m062x_gas", "halobenzene")], q[("m062x_smd", "halobenzene")]
print("font", font)
print(f"gas: bromobenzene fit R2={float(g['R2']):.2f} slope={float(g['slope']):.3f} p={float(g['p']):.4f} n={g['n']}, "
      f"pred RMSE (Cl + PBDE)={g['pred_RMSE_other_groups']}; halobenzene R2={float(ga['R2']):.2f} n={ga['n']}")
print(f"smd: bromobenzene fit R2={float(s['R2']):.2f} slope={float(s['slope']):.3f} p={float(s['p']):.4f} n={s['n']}, "
      f"pred RMSE (Cl + PBDE)={s['pred_RMSE_other_groups']}; halobenzene R2={float(sa['R2']):.2f} n={sa['n']}")
