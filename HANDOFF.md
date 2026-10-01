# Session handoff — 2026-09-30

Read this, README.md, config.py and hx.py at the start of a new session.

## Where things stand

- **All Gaussian work is finished.** No jobs are queued on Sherlock (`squeue --me` empty).
- **Git:** Mac, Sherlock (`$GROUP_HOME/Haloaromatics`) and GitHub were all in sync at the end of
  the session. Uncommitted on purpose:
  - Mac: `figures/QSAR_workbook.xlsx`, `figures/QSAR_MLR_exploration.xlsx` show as modified,
    probably saved from Excel by Jonas. Ask before committing or discarding.
  - Sherlock: `calcs/b3lyp_smd/logs/BDE_245_24_RA.log` (failed job, waiting for its retry, see To do 1).
- **Compile on the Mac** (`.venv/bin/python hx.py compile`); Sherlock has no openpyxl/matplotlib
  for python 3.12. Sherlock only scrapes and commits `results/raw`.
- The Mac's `results/` (compiled CSVs, workbook, `results/plots/`) has **not** yet been recompiled
  with the final PBDE data: run `hx.py compile` and `figures/ts_scan_analysis.py` after `git pull`.

## Done this session (most recent first)

1. **TS scans interpreted** — `figures/ts_scan_analysis.py` (classifies all 94 M06-2X scans,
   converged points only: smooth barrier / ~barrierless / state crossing / uphill to end;
   writes `results/ts_scan_summary.csv` and `results/plots/ts_scans_<level>.png`) and the draft
   `figures/TS_scan_interpretation.md`. Main points: no stepwise barrier for the measured
   compounds (bromobenzene and BDE-47 radical anions dissociate; chlorobenzene π radical anions
   cleave over ≤ ~5 kcal/mol); other sites show state crossings, not saddle points; gas-phase
   scans of the most halogenated compounds run uphill (unsolvated halide); BDE-99: computed C2
   cleavage (→ BDE-66) vs observed BDE-47 (C5) is unresolved. **No TS jobs were submitted, by
   decision**; only compounds without rate constants have smooth 4–9 kcal/mol barriers.
2. **Parser / QC fixes** (`haloaro/parse.py`, `haloaro/qc.py`, tested; no change to the 1234
   existing logs): error reason = last known message after the last "Normal termination" (was:
   last 800 lines, which missed PBDE errors); partial terminations without an error line are
   `incomplete` instead of `fail`, so plain `retry` cannot touch running jobs.
3. **Comparison with LaPier et al. 2026 SI (Table S5)** at M06-2X/gas: LUMO and the 2e ΔG
   (eq 3) reproduce for the bromobenzenes (R² 0.915 vs 0.92). Radical-anion descriptors differ for
   bromobenzene, 1,2,4,5-tetra-, hexabromo and PBDEs (radical-anion energies 5–25 kcal/mol lower
   here). Table S5 eq 2 entries for 1,2,4-TriBr (+11.1, same as its eq 5) and 1,2,4,5-TetrBr (−23.7)
   look wrong; with ours, eq 2 R² 0.06 → 0.61. Paper's PBDE values match our C4 (BDE-47) and C2
   (BDE-99) sites, not our most favourable site; with ours the PBDEs are overpredicted by ~3 ln
   units. The paper's "MSE" column is SSE. Not written up in a file yet (results are in the chat
   summary of 2026-09-30).
4. **MLR exploration** — `haloaro/mlr.py`, `figures/qsar_mlr_explore.py` →
   `figures/QSAR_MLR_exploration.xlsx` (Summary, Interpretation with sourced literature,
   Methods, All models, Candidates, parity charts). 2928 models with 1–3 descriptors, 200
   stratified train/test splits (`MLR_*` in config.py). Best: LUMO + Savéant intrinsic barrier +
   Wiberg C–X (M06-2X/SMD), test RMSE 1.13 ± 0.37; best single 1.31; n_X 1.67; mean 2.68. Wiberg
   bond order in 27/30 top models (halogen-family offset); two descriptors are the most the data
   support; for extrapolation from bromobenzenes, LUMO + Fukui f⁺ transfers best.
5. **Jupyter notebook** — `notebooks/Haloaromatics_QSAR_workflow.ipynb` (executed; regenerate
   with `notebooks/make_notebook.py`).
6. **Main QSAR workbook** — `figures/qsar_workbook.py` → `figures/QSAR_workbook.xlsx` (ReadMe
   with every parameter and statistic, Performance ranked by LOO Q², Data, 69 live Excel charts;
   5 × 4 in charts, smaller type, legend below the plot). Chart builder:
   `figures/make_excel_chart.py` (copy of the plot-like-jonas skill's, extended).
7. Earlier: chlorobenzene k_obs added (unpublished); QSAR summary and per-group fits; styled
   matplotlib figures (`haloaro/plotting.py`); project moved to `$GROUP_HOME` after `$HOME` quota
   filled; PBDE AM1 → optfreq → sp → scans run; 2-day limit for PBDE scans.

## To do

1. **B3LYP/SMD BDE-99 radical anion** (`BDE_245_24_RA`): failed at the opt step limit
   (dissociating, C4–Br 2.72 Å). Now that the parser reads `opt_maxcycles`, on Sherlock:
   `python3 hx.py retry b3lyp_smd --only BDE_245_24_RA && python3 hx.py submit b3lyp_smd --only BDE_245_24_RA`
   (restart from its chk, `Opt=(MaxCycles=300,CalcFC)`). Then its sp jobs if it ends bound.
   Needs Jonas's OK (it was put on hold).
2. **Recompile on the Mac** after `git pull`: `hx.py compile`, `figures/ts_scan_analysis.py`,
   `figures/qsar_workbook.py`, `figures/qsar_mlr_explore.py`; check figures against the
   plot-like-jonas checklist; commit results.
3. **Decisions for Jonas:**
   - PBDE site convention for descriptors: most favourable site (current), degeneracy-weighted,
     or experimentally observed cleavage site (BDE-99 → BDE-47 is C5).
   - Whether to locate TSs for the few smooth-barrier compounds without rate constants.
   - Whether to add the LaPier-2026 comparison as a tab or document.
4. **Literature:** add Costentin, Robert & Savéant, J. Am. Chem. Soc. 2004, 126, 16051 (bending
   of the cleaving bond in aryl-halide radical anions) to Zotero, then extend
   `TS_scan_interpretation.md`. Ask Jonas to check Table S5 eq 2 entries (1,2,4-TriBr, 1,2,4,5-TetrBr).
5. **Open science question:** BDE-99 regioselectivity (computed C2 vs observed C5). Ideas: PBDE
   conformer search (CREST/xTB, recommended in README), other radical-anion states, BDE-66 standard.
6. **Task 3 from the start of the project** (Slurm "advance" job chaining stages) was never built;
   the pipeline was babysat interactively instead. Probably not needed now.
7. **Housekeeping:** 11 orphaned text files in the Zotero index (`build_report.json`
   `orphan_text_files`), deletion only with Jonas's OK. Index updates on this Mac must run with
   `.venv/bin/python` (has pypdfium2/pypdf; the system python3 has none).

## Working rules that held all session

- Sherlock: Jonas opens a ControlMaster session; connect with
  `ssh -o ControlPath=~/.ssh/cm-%r@%h:%p -o BatchMode=yes jlapier@login.sherlock.stanford.edu`.
  Never enter passwords/Duo. Work in `$GROUP_HOME/Haloaromatics` with `module load python/3.12.1`.
  Read `/etc/agents/AGENTS.md` rules (no polling loops, aggregate short jobs, etc.).
- Commits: Claude commits locally (Mac and Sherlock); Jonas pushes. Never commit logs of running
  jobs. Use explicit `git add <files>` (never `commit -a`; Jonas may have edited workbooks).
- Ask before: >~50 jobs, deleting files, changing methods/thresholds, pushing.
- Scientific conventions: kobs at −2.0 V vs SHE (uncompensated potential, without iR-drop
  compensation) at first mention; 1,3-DBB kobs 0.0329 h⁻¹; EXP_POTENTIAL_V = −2.0;
  chlorobenzene kobs source "unpublished"; don't change computational methods to match the paper.
- Literature claims only from passages actually read in the Zotero index (`use-zotero-index`).
