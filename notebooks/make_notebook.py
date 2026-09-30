"""Write notebooks/Haloaromatics_QSAR_workflow.ipynb (run from the project root).

    .venv/bin/python notebooks/make_notebook.py
    .venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/Haloaromatics_QSAR_workflow.ipynb
"""
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
    md("""# Haloaromatics: computational workflow and QSAR analysis

This notebook recreates the analysis of the project in Python, step by step, from the scraped Gaussian results to the QSARs:

1. The pipeline on Sherlock (commands only; the calculations are not run here)
2. Species and reactions
3. Scraped data and quality control
4. Reaction free energies, Savéant barriers and descriptors
5. Experimental rate constants
6. Single-descriptor QSARs
7. Example figure: two-electron reduction free energy
8. Multiple-linear-regression QSARs with training/test splits
9. Rebuilding the workbooks

Rate constants were measured at −2.0 V vs SHE (uncompensated potential, without iR-drop compensation). All methods, constants and thresholds come from `config.py`; the code lives in `haloaro/`. Run it from the project's `.venv` (`.venv/bin/jupyter lab`) after `git pull`, so that `results/raw/` holds the latest scrape."""),
    code("""import os, sys, csv, json
from pathlib import Path
from collections import Counter

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
os.chdir(ROOT); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "figures"))

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import config as C
from haloaro import thermo, descriptors, correlations, mlr, plotting
from haloaro.stages import STAGES, SPECIES, REACTIONS, jobs
pd.set_option("display.max_columns", 30); pd.set_option("display.width", 200)
print("project root:", ROOT)
print("levels:", list(C.LEVELS))"""),
    md("""## 1. Pipeline on Sherlock

Gaussian 16 jobs run on Sherlock through `hx.py`, from `$GROUP_HOME/Haloaromatics` (the checkpoints are too large for `$HOME`):

| stage | what | command |
|---|---|---|
| `am1` | AM1 pre-optimisation | `python3 hx.py generate am1 && python3 hx.py submit am1 --serial --time 0-00:30:00` |
| `optfreq` | Opt + Freq at B3LYP and M06-2X, 6-311++G(d), gas and SMD(water) | `python3 hx.py generate optfreq && python3 hx.py submit optfreq` |
| `sp` | vertical anion/cation/neutral single points, CM5/Hirshfeld, Wiberg | `python3 hx.py generate sp && python3 hx.py submit sp` |
| `tsscan`, `ts` | relaxed C–X scans of the radical anion and TS searches (M06-2X) | `python3 hx.py generate tsscan && python3 hx.py submit tsscan` |
| status / retry | QC and automatic fixes (3-try cap) | `python3 hx.py status all`, `python3 hx.py retry <stage> --only <job>` |
| scrape | logs → `results/raw/<stage>.csv` | `python3 hx.py scrape all` |

`compile` (reactions, barriers, descriptors, correlations, figures, workbook) runs on the Mac with `.venv/bin/python hx.py compile`. The rest of this notebook does what `compile` does, one step at a time."""),
    md("## 2. Species and reactions"),
    code("""kinds = Counter(s.kind for s in SPECIES.values())
print(pd.Series(kinds, name="species").to_string())
print(f"\\n{len(REACTIONS)} symmetry-unique dehalogenation reactions (C–X bonds)")
pd.DataFrame([{"stage": n, "jobs": len(jobs(st))} for n, st in STAGES.items()]).set_index("stage").T"""),
    code("""pd.DataFrame([vars(r) for r in REACTIONS]).head(8)"""),
    md("## 3. Scraped data and quality control\n`results/raw/<stage>.csv` holds one row per job with energies, orbital energies, ⟨S²⟩, imaginary frequencies and the QC status (`ok`, `warn`, `fail`). Radical anions are classified by their longest C–X bond (π, bent σ, or dissociated)."),
    code("""RAW = Path("results/raw")
raw = {st: pd.read_csv(f) for st in STAGES if (f := RAW / f"{st}.csv").exists() and f.stat().st_size}
qc = pd.DataFrame({st: df["status"].value_counts() for st, df in raw.items()}).fillna(0).astype(int).T
qc"""),
    code("""ra = pd.concat([df.assign(level=lv) for lv, df in raw.items() if lv in C.LEVELS], ignore_index=True)
ra = ra[ra["kind"] == "radical_anion"]
pd.crosstab(ra["name"].str.replace("_RA", ""), ra["level"], values=ra["RA_state"], aggfunc="first")"""),
    md("## 4. Reaction free energies, Savéant barriers and descriptors\nThe same functions `hx.py compile` uses: `thermo.reaction_table` (ΔG and E° of each step per C–X bond), `thermo.det_table` (C–X bond energies and Savéant concerted barriers), `descriptors.molecular_table` and `descriptors.site_table`."),
    code("""tab, species_rows = {}, []
for sn, st in STAGES.items():
    f = RAW / f"{sn}.csv"
    if not f.exists() or f.stat().st_size == 0:
        continue
    for r in csv.DictReader(open(f)):
        tab[(sn if st.kind in ("ts", "tsscan") else r["level"], r["name"])] = r
atoms = {}
for sn, st in STAGES.items():
    f = RAW / f"{sn}_atoms.json"
    if f.exists() and st.kind in ("optfreq", "sp"):
        for name, v in json.loads(f.read_text()).items():
            atoms[(st.level, name)] = v

rx = thermo.reaction_table(tab, REACTIONS, SPECIES)
det = thermo.det_table(tab, REACTIONS)
ts = thermo.ts_table(tab, REACTIONS)
mol, lam_i = descriptors.molecular_table(tab, SPECIES, rx)
path = thermo.pathway_table(rx, ts, det, lam_i)
sites = descriptors.site_table(tab, atoms, SPECIES, REACTIONS, rx, ts, det)
rx_df = pd.DataFrame(rx)
cols = ["level", "parent", "site", "RA_state", "dG_ET_any_kcal", "dG_frag_any_kcal", "dG_1e_kcal",
        "dG_2e_carbanion_kcal", "E_2e_carbanion_V"]
rx_df[rx_df["parent"].isin(["BrBz_1", "BrBz_123456", "ClBz_123456", "BDE_245_24"])][cols].head(12)"""),
    code("""mol_df = pd.DataFrame(mol)
mol_df[mol_df["level"] == "m062x_gas"][["parent", "n_X", "LUMO_eV", "VEA_eV", "omega_dSCF_eV", "RA_state",
                                        "dG_ET_kcal", "min_dG_2e_carbanion_kcal"]].head(16)"""),
    md("## 5. Experimental rate constants\n`data/experimental_kobs.csv`: bromobenzenes and PBDEs from LaPier et al., *Environ. Sci. Technol.* 2026, 60, 1346 (Table 1; 1,3-dibromobenzene corrected to 0.0329 h⁻¹); chlorobenzenes unpublished."),
    code("""kobs = pd.read_csv(C.EXPERIMENTAL_KOBS)
kobs["ln_kobs"] = np.log(kobs["kobs_per_h"])
kobs[["compound", "label", "group", "kobs_per_h", "se_per_h", "ln_kobs", "half_life", "source"]]"""),
    md("""## 6. Single-descriptor QSARs
`correlations.run` fits ln k_obs = slope·x + intercept for every descriptor, level and compound set, and reports R², p, RMSE and leave-one-out Q² (predictive ability within the fitted set), plus the RMSE of predicting the compounds not used in the fit (`pred_RMSE_other_groups`). The bromobenzenes (n = 7) are the training set."""),
    code("""corr, cdata, cpairs, cpred, qsar = correlations.run(mol, rx, det, path, ts, sites)
q = pd.DataFrame(qsar)
q = q[(q["fit_set"] == C.CORR_FIT_GROUP) & (q["descriptor_set"] != "MLR") & q["R2"].notna()]
for c in ("R2", "Q2", "RMSE", "pred_RMSE_other_groups", "p", "n"):
    q[c] = pd.to_numeric(q[c])
q.sort_values("Q2", ascending=False)[["level", "model", "n", "R2", "Q2", "RMSE", "pred_RMSE_other_groups", "p"]].head(12)"""),
    code("""grid = q.pivot_table(index="model", columns="level", values="Q2")
grid.sort_values("m062x_gas", ascending=False).round(2)"""),
    md("Figures in the project's plot style (`haloaro/plotting.py`): the six best models per level, fitted on the bromobenzenes (filled gray), with the chlorobenzenes (open) and PBDEs (squares) predicted."),
    code("""plotting.use_style()
figs = plotting.correlation_figures(qsar, cdata, "results/plots")   # writes PNG/PDF files
%matplotlib inline
from IPython.display import Image, display
display(Image("results/plots/correlations_m062x_gas.png", width=800))"""),
    md("## 7. Example: ln k_obs vs ΔG of ArX + 2e⁻ → Ar⁻ + X⁻ (M06-2X, gas and SMD)"),
    code("""fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
panels = []
for ax, lv in zip(axes, ("m062x_gas", "m062x_smd")):
    fit = q[(q["level"] == lv) & (q["descriptors"] == "min_dG_2e_carbanion_kcal")].iloc[0].to_dict()
    pts = plotting._points(cdata, lv, "min_dG_2e_carbanion_kcal")
    drawn, line = plotting._panel(ax, fit, pts, "min_dG_2e_carbanion_kcal", "dG ArX + 2e- -> Ar- + X-",
                                  ([-8, -6, -4, -2, 0, 2, 4], 0), ms=6, fontsize=9)
    ax.set_title(f"{plotting.level_text(lv)}: R² = {fit['R2']:.2f}", fontsize=9)
    panels.append((ax, pts, line))
fig.tight_layout()
for ax, pts, line in panels:          # labels after the layout is final
    plotting._place_labels(fig, ax, pts, line, fontsize=7)
plt.show()"""),
    md("""## 8. Multiple-linear-regression QSARs with training/test splits
`haloaro/mlr.py` fits every combination of 1–3 descriptors that have values for all 14 compounds, and scores each on 200 random stratified splits (10 training compounds, 4 test compounds: 2 bromobenzenes, 1 chlorobenzene, 1 PBDE). The test RMSE is pooled over all splits, so a model is judged only on compounds it was not fitted on. Models whose descriptors correlate with |r| > 0.9 are flagged as collinear. A transfer test fits on the bromobenzenes only and predicts the chlorobenzenes and PBDEs. (This cell takes about a minute.)"""),
    code("""models, info = mlr.explore(cdata)
m = pd.DataFrame(models)
m["model"] = m["descriptors"].apply(lambda d: " + ".join(d))
ok = m[~m["collinear"]]
print(f"{len(m)} models, {len(ok)} non-collinear; mean-model test RMSE {info['m062x_gas']['RMSE_test_mean_model']:.2f}")
ok[["level", "k", "model", "RMSE_test", "RMSE_test_sd", "Q2_F1", "R2_train", "R2", "RMSE_transfer", "max_abs_r"]].head(10).round(2)"""),
    code("""print("best per number of descriptors:")
display(ok.groupby("k").head(1)[["level", "k", "model", "RMSE_test", "RMSE_test_sd", "Q2_F1"]].round(2))
print("descriptor frequency in the top 30 models:")
Counter(d for ds in ok.head(30)["descriptors"] for d in ds).most_common(8)"""),
    code("""best = ok.iloc[0]
pred = mlr.oos_predictions(cdata, best["level"], best["descriptors"])
fig, ax = plt.subplots(figsize=(3.5, 3.3))
pts = [{"x": yo, "y": yp, "group": r["group"], "label": r["label"], "yerr": None} for r, yo, yp, n in pred]
plotting._series(ax, pts, 6)
one = ax.plot([-8, 8], [-8, 8], color="#a6a6a6", lw=0.75, zorder=1)[0]      # 1:1 line
for a in ("x", "y"):
    plotting._set_ticks(ax, a, [-8, -4, 0, 4, 8], 0)
ax.set_xlabel("observed ln k$_\\\\mathrm{obs}$ (h$^\\\\mathrm{−1}$)"); ax.set_ylabel("out-of-sample predicted ln k$_\\\\mathrm{obs}$")
ax.set_title(best["model"].replace(" + ", " +\\n"), fontsize=7)
fig.tight_layout(); plotting._place_labels(fig, ax, pts, one, fontsize=6.5); plt.show()"""),
    md("""Interpretation (details and sources in `figures/QSAR_MLR_exploration.xlsx`, Interpretation tab): the Wiberg C–X bond order recurs in the best mixed-training models because it tracks ln k_obs within each halogen family but is offset between them, so it pairs with a halogen-sensitive term (Savéant intrinsic barrier, C–X bond energy, LUMO or a reduction free energy). Two descriptors are the most these 14 compounds support. Trained on bromobenzenes alone, frontier-orbital and Fukui descriptors transfer best to the chlorobenzenes and PBDEs."""),
    md("## 9. Rebuilding the workbooks\nThe Excel deliverables are built by scripts in `figures/` from the compiled results (uncomment to run):"),
    code("""# !python hx.py compile
# !python figures/qsar_workbook.py figures/QSAR_workbook.xlsx
# !python figures/qsar_mlr_explore.py figures/QSAR_MLR_exploration.xlsx
print(sorted(p.name for p in Path("figures").glob("*.xlsx")))"""),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3"}})
nbf.write(nb, "notebooks/Haloaromatics_QSAR_workflow.ipynb")
print("wrote notebooks/Haloaromatics_QSAR_workflow.ipynb")
