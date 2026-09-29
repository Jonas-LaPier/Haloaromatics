# Haloaromatics

Electronic-structure workflow (Gaussian 16 on Sherlock) for the reductive
dehalogenation of chloro- and bromobenzenes. Everything runs through one
driver script, `hx.py`. Methods, solvent, resources and constants are set in `config.py`.

## Chemical scope

The code enumerates every symmetry-unique isomer and every symmetry-unique C–X bond
(`python3 hx.py species` writes `results/species_list.csv` and `results/reaction_list.csv`):

| species | charge / mult | count |
|---|---|---|
| parents ArX: benzene + 12 chloro + 12 bromo + BDE-47 + BDE-99 | 0 / 1 | 27 |
| PBDE hydrodebromination products (BDE-17, -28, -48, -49, -66, -74) | 0 / 1 | 6 |
| radical anions ArX•⁻ (stepwise electron transfer) | −1 / 2 | 26 |
| aryl radicals Ar• (one per unique C–X bond; phenyl shared by both series) | 0 / 2 | 46 |
| aryl carbanions Ar⁻ | −1 / 1 | 46 |
| Cl⁻, Br⁻ | −1 / 1 | 2 |
| Cl•, Br• atoms (for C–X bond energies) | 0 / 2 | 2 |
| C–X cleavage transition states of ArX•⁻ | −1 / 2 | 47 |

That gives 47 unique dehalogenation reactions: 20 per halogen for the benzenes, plus 2
for BDE-47 and 5 for BDE-99. The reaction table also lists each reaction's degeneracy,
meaning how many equivalent C–X bonds it covers.

**PBDEs** (`PBDES` in `config.py`) use a diphenyl-ether skeleton.
- **Rings and sites:** rings A and B are each numbered from the ether carbon. Sites are
  `2`–`6` in ring A and `2p`–`6p` in ring B.
- **Names:** canonical names put the more substituted ring first (`BDE_24_24` is BDE-47,
  `BDE_245_24` is BDE-99). An `r` marks the radical carbon (`BDE_24_2r4_rad`).
- **Starting geometry:** a twisted C₂-like conformation (C–O–C 120°, both rings rotated
  50°). No conformer search is done; a CREST/xTB search is recommended before
  interpreting small energy differences.
- **Resources:** PBDE jobs use `RESOURCES_DPE` (16 CPUs, 48 GB, 12 h for opt/freq, scan and TS),
  and `submit` sends them as a separate array.

**Names.** `ClBz_124` is 1,2,4-trichlorobenzene and `ClBz_124_RA` is its radical anion.
`ClPh_24_rad` and `ClPh_24_anion` are the 2,4-dichlorophenyl radical and carbanion,
with the radical/anion carbon numbered C1. `TS_ClBz_124_x4` is the TS for breaking the
C4–Cl bond in `ClBz_124_RA`. Atom order in every input is C1…C6 followed
by the substituents in site order. For PBDEs the order is ring A C1–C6, ring B C1–C6, O,
then the substituents of ring A and then ring B.

## Stages

| stage | what | geometry from |
|---|---|---|
| `am1` | Opt AM1 (every molecular species) | idealised planar ring |
| `b3lyp_gas`, `b3lyp_smd`, `m062x_gas`, `m062x_smd` (group `optfreq`) | Opt + Freq, 6-311++G(d), gas and SMD(water) | `Geom=Check` from AM1 chk |
| `tsscan_m062x_gas`, `tsscan_m062x_smd` (group `tsscan`) | relaxed C–X scan of ArX•⁻ (16 × 0.08 Å) | DFT radical-anion minimum, X tilted 15° out of plane |
| `ts_m062x_gas`, `ts_m062x_smd` (group `ts`) | Opt=(TS,CalcFC,NoEigenTest) + Freq | highest point of the scan |
| `sp_<level>` for each opt/freq level (group `sp`) | single points with `Pop=(Hirshfeld,NBORead)` and Wiberg bond indices; 123 jobs per level (below) | `Geom=Check` from the opt/freq chk |

Single-point jobs, run for each parent P:

| job | state | geometry | gives |
|---|---|---|---|
| `P__pop` | neutral | neutral | CM5/Hirshfeld charges, Wiberg C–X |
| `P__vA` | anion (−1, doublet) | neutral | vertical EA, Fukui f⁺, spin density |
| `P__vC` | cation (+1, doublet) | neutral | vertical IE |
| `P_RA__pop` | radical anion | radical anion | spin density, charges, Wiberg C–X |
| `P_RA__vN` | neutral | radical anion | four-point λᵢ |

Benzene has only the first three. Radical anions that dissociated have no `P_RA__*` jobs.

Barriers are computed for two pathways. `compile` compares them in `pathway_comparison.csv`.
## Pathway 1, stepwise: TS for C–X cleavage in the radical anion

ArX + e⁻ → ArX•⁻ → [Ar···X]‡•⁻ → Ar• + X⁻. The TS is located at M06-2X/6-311++G(d),
gas and SMD, starting from the maximum of a relaxed C–X scan. The halogen is tilted out
of the plane at the start of the scan. This breaks the planar mirror symmetry so the π*
and σ* states can mix; otherwise the planar path crosses between the two states at a
cusp rather than a smooth saddle point. If a radical anion has no bound minimum at a
level (flag `RA_dissociated`), no scan is generated there, because only the concerted
pathway exists at that level. To change the TS level, edit `TS_LEVELS` in `config.py`.

## Pathway 2, concerted: Savéant dissociative electron transfer

Barriers are computed for the concerted pathway starting from the **neutral parent**,
ArX + e⁻ → [ArX]‡ → Ar• + X⁻, using Savéant's model (J. Am. Chem. Soc. 1987, 109, 6788;
Acc. Chem. Res. 1993, 26, 455). No extra Gaussian jobs are needed beyond `optfreq`, which
now includes the Cl• and Br• atoms.

- **D**, the C–X bond energy of the neutral parent: D = H(Ar•) + H(X•) − H(ArX). The X•
  energies include the ²P₃/₂ spin–orbit correction (Cl 0.84, Br 3.51 kcal/mol). Set
  `SAVEANT_D = "G"` in `config.py` to use the bond dissociation free energy (BDFE) instead.
- **λ₀**, the solvent reorganization energy, from the Marcus–Hush one-sphere model for an
  electrode reaction: λ₀ = e²/(8πε₀a)·(1/ε_op − 1/ε_s). Here a is the radius of a sphere
  with the parent's van der Waals volume (Bondi radii, computed from the optimized
  geometry). Set `LAMBDA0_MODEL = "homogeneous"` for a molecular electron donor, or
  `LAMBDA0_KCAL` to impose a value.
- **Intrinsic barrier** ΔG₀‡ = (D + λ₀)/4, the barrier at E = E°_DET.
- **Barrier at potential E:** ΔG‡(E) = ΔG₀‡(1 + ΔG°(E)/4ΔG₀‡)², with ΔG°(E) = F(E − E°_DET).
  E°_DET is the one-electron potential `E_1e`. It is evaluated at each potential in
  `DET_POTENTIALS_V`, along with the transfer coefficient α.
- The full barriers are reported for the SMD levels only. The gas levels give bond energies only.

## Stepwise vs concerted comparison (SMD TS levels)

At each potential E in `DET_POTENTIALS_V`:

- **ET step:** ΔG_ET(E) = F(E − E°_ET). Its barrier ΔG‡_ET(E) = (λ/4)(1 + ΔG_ET/λ)², from
  Marcus theory, with λ = λ₀ + λᵢ. λ₀ is the outer-sphere term from the Savéant section.
  λᵢ is the inner-sphere term from Nelsen's four-point method (`sp` stage), taken from the
  gas-phase level of the same functional by default (`LAMBDA_I_FROM_GAS`).
- **Stepwise barrier:** ΔG‡_step(E) = max(ΔG‡_ET(E), max(ΔG_ET(E), 0) + ΔG‡_frag). The
  rate-limiting step is either the electron transfer or the C–X cleavage of a
  pre-equilibrated radical anion.
- **Concerted barrier:** ΔG‡_conc(E), from Savéant's model.
- **Favored pathway:** whichever of the two barriers is lower.

Set `RUN_RA_TS = False` in `config.py` to skip the TS stages and report only the concerted analysis.

## Workflow on Sherlock

```bash
cd $GROUP_HOME/Haloaromatics  # project root; all paths are relative to here
module load python/3.12.1     # the system python3 (3.6) is too old for this code
python3 hx.py generate am1
python3 hx.py submit am1 --serial --time 0-00:30:00   # AM1 jobs take seconds: run them in one job
python3 hx.py submit m062x_gas  # one job array; add --dry-run to only print the sbatch line
python3 hx.py status am1      # QC summary: ok / warn / fail / incomplete / missing
python3 hx.py retry am1       # rebuild failed inputs with automatic fixes, then submit again

python3 hx.py generate optfreq      # only species whose AM1 job passed QC
python3 hx.py submit optfreq
python3 hx.py status optfreq

python3 hx.py generate tsscan && python3 hx.py submit tsscan   # needs the m062x_* radical anions
python3 hx.py generate ts     && python3 hx.py submit ts       # from the scan maxima
python3 hx.py generate sp     && python3 hx.py submit sp       # single points (any time after optfreq)

python3 hx.py scrape all      # logs -> results/raw/<stage>.csv
python3 hx.py compile         # -> results/*.csv and results/Haloaromatics_results.xlsx
```

**Tip:** before submitting all 129 jobs of a new level, test one bromine species first:
`python3 hx.py submit m062x_gas --only BrBz_1`. This confirms that 6-311++G(d) is
defined for Br in your Gaussian build.

Run jobs from `$GROUP_HOME/Haloaromatics` (1 TB), not `$HOME` (15 GB). The checkpoint
files (about 70 MB each for the PBDEs, over 1,000 in total) filled the `$HOME` quota and
killed every running job. `~/Haloaromatics` is kept as a copy of the logs and results
only, without checkpoints, so it cannot run jobs that read a `.chk` (`Geom=Check`,
retries).

Details:

- `submit` writes `calcs/<stage>/manifest_<time>.txt`, which lists every input that has
  no log yet, and submits `slurm/g16_array.sbatch` as one array (throttle `%50`). This
  means resubmitting never reruns finished jobs. It prints the Slurm job ID.
- `submit --serial --time <limit>` runs the whole manifest one input after another in a
  single job instead of an array. Use it for short jobs such as AM1 (seconds each):
  Sherlock asks that the work inside a job last at least 10 minutes. The time limit must
  be given explicitly and covers the whole batch.
- Gaussian scratch goes to `$L_SCRATCH` (node-local). Set `FORMCHK=1` before `sbatch`
  to also write `.fchk` files.
- `retry` moves the failed log and input to `calcs/<stage>/logs/failed/<name>.tryN.*`
  and writes a new input:
  - SCF failure → `SCF=(XQC,VShift=400,MaxCycle=512) NoSymm`. If it failed partway through
    an optimization, the job restarts from the geometry it reached, but not its wavefunction.
  - optimization ran out of cycles or time → restart from its own chk with
    `Geom=Check Guess=Read`, `Opt=(MaxCycles=300,CalcFC)`, `NoSymm`
  - internal-coordinate errors → `Opt=Cartesian`
  - an imaginary frequency at a minimum → displace along that mode (0.25 Å for modes
    above 60 cm⁻¹, otherwise 0.15 Å plus `Tight`) and reoptimize with `Opt=CalcFC NoSymm`
  - `retry` refuses a job that has already failed `MAX_RETRIES` times (default 3); use
    `--force` to override.
- Radical anions are optimized without symmetry, starting from a slightly perturbed AM1
  geometry. From symmetric starting rings, the symmetry-constrained optimization converges
  to planar saddle points with an imaginary out-of-plane or Jahn–Teller mode.
- `--only NAME ...` limits any stage command to specific jobs. `--force` regenerates inputs.

### Quality checks (`haloaro/qc.py`)

- Normal termination count (2 for Opt+Freq)
- Charge, multiplicity and atom count match the species
- A stationary point was found and the free energy is present
- There are no imaginary frequencies at minima and exactly 1 at TSs. Extra imaginary modes
  smaller than `IMAG_TOL_CM` (20 cm⁻¹) are flagged `small_imag` as a warning rather than a
  failure, since they are integration-grid noise on floppy radical anions.
- The TS imaginary mode contains the C–X stretch (flag `ts_mode_not_CX_stretch` otherwise)
- ⟨S²⟩ is within 0.10 of 0.75 for doublets
- **Radical-anion state** (`RA_state` column; set by the longest C–X bond and the
  out-of-plane angle of X):
  - `pi`: C–X ≤ 2.00 Å (Cl) or 2.15 Å (Br)
  - `sigma_bent` (flag `RA_sigma`): X bent ≥ 10° out of the plane, with C–X ≤ 2.8 Å (Cl)
    or 3.0 Å (Br). This is a loose, bent σ-type radical anion.
  - `dissociated` (flag `RA_dissociated`): a planar C–X longer than 2.3 Å (Cl) or 2.5 Å
    (Br), or any C–X beyond the σ limit. It is an Ar•···X⁻ complex, so there is no
    ArX•⁻ intermediate at that level. Its ET and fragmentation ΔG are left blank, no
    stepwise TS is generated, and the pathway comparison reports "concerted only".
  - `elongated`: anything else.

## Output (`results/`)

- `species_energies.csv`: every species × level with E, H, G (Hartree), HOMO/LUMO, ⟨S²⟩, imaginary frequencies, QC status
- `parent_LUMO.csv`: parent LUMO energies (Hartree and eV) at every level
- `reactions.csv`: one row per level × unique C–X bond:
  - `ET` ArX + e⁻ → ArX•⁻
  - `1e` ArX + e⁻ → Ar• + X⁻ (dissociative one-electron reduction)
  - `frag` ArX•⁻ → Ar• + X⁻
  - `rad_red` Ar• + e⁻ → Ar⁻
  - `2e_carbanion` ArX + 2e⁻ → Ar⁻ + X⁻
  - `2e_HDH` ArX + H⁺ + 2e⁻ → ArH + X⁻ (hydrodehalogenation)
  - `RA_2e` ArX•⁻ + e⁻ → Ar⁻ + X⁻ (second electron to the radical anion)

  Each gives ΔG in kcal/mol and E° in V vs SHE (SMD levels only).
- `ts_barriers.csv`: stepwise ΔG‡ and ΔE‡ = TS − ArX•⁻, ΔG‡ relative to ArX + e⁻, imaginary frequency, QC flags
- `pathway_comparison.csv`: E°_ET, E°_DET, λ₀, ΔG‡_frag and ΔG₀‡, plus the stepwise and concerted ΔG‡ and the favored pathway at each potential
- `det_barriers.csv`: C–X bond enthalpy and free energy (at this level and in gas phase), E°_DET, radius a, λ₀, ΔG₀‡, and ΔG‡ and α at each potential in `DET_POTENTIALS_V`
- `molecular_descriptors.csv`: one row per level × parent:
  - HOMO, LUMO and gap
  - Koopmans and ΔSCF values of μ, η and ω (electrophilicity)
  - vertical and adiabatic EA, vertical IE
  - LUMO of the relaxed radical anion (`LUMO_RA_eV`, any state)
  - λ_N, λ_A and λᵢ
  - dipole moment, isotropic polarizability and van der Waals radius
  - ΔG_ET and E°_ET
  - the most favorable site's ΔG for each reaction
- `site_descriptors.csv`: one row per level × unique C–X bond, for regioselectivity:
  - structure: halogens ortho, meta and para to the site
  - parent: r(C–X), CM5 and Hirshfeld charges, Wiberg C–X
  - vertical anion: Fukui f⁺ on C, X and C+X, spin density
  - relaxed radical anion: r(C–X), elongation, out-of-plane angle of X, spin density,
    charges, Wiberg C–X
  - site energetics: ΔG_1e, ΔG_frag, ΔG_2e, BDE, Savéant ΔG₀‡, stepwise ΔG‡_frag, and the
    pKa of the ArH C–H formed at the site (SMD levels, direct scheme)
  - `rank_*` columns: each site's rank within its parent (1 = most reactive by that descriptor)
- `product_distribution.csv`: predicted fraction of mono-dehalogenation at each site, from
  degeneracy × exp(−ΔG‡/RT) using the stepwise TS, the concerted barrier, and the combined
  stepwise + concerted rate at each potential
- `qc_issues.csv`: every job whose status is not `ok`
- `qsar_summary.csv` (first sheet of the workbook): the correlation models to check, for
  every level and fit set, in the order of `QSAR_DESCRIPTORS` in `config.py`: set `thermo`
  (parent and radical-anion LUMO, ΔG of ArX + e⁻ → Ar• + X⁻, ArX + 2e⁻ → Ar⁻ + X⁻,
  ArX + e⁻ → ArX•⁻, ArX•⁻ → Ar• + X⁻ and ArX•⁻ + e⁻ → Ar⁻ + X⁻), set `extended` (EA,
  electrophilicity, E°, λᵢ, BDE, Savéant and stepwise barriers, Fukui f⁺, Wiberg bond
  orders, spin densities), and the two-descriptor models in `CORR_PAIRS`.
- `correlations.csv`: ln(k_obs) against every compound-level descriptor at every level.
  Each fit reports n, slope, intercept, r, R², p, SSE, MSE (= SSE/n), RMSE, leave-one-out
  RMSE and Q². Each descriptor is fitted on every set in `CORR_FIT_SETS`: the bromobenzenes
  (`CORR_FIT_GROUP`, whose fit is then used to predict all other compounds), the
  chlorobenzenes, the halobenzenes (bromo- and chlorobenzenes) and all compounds. Site quantities are reduced to the
  most favorable site (`min_*`/`max_*`) or to a degeneracy-weighted effective barrier
  (`eff_*` = −RT ln Σ gᵢ exp(−ΔG‡ᵢ/RT)). Potential-dependent barriers are evaluated at
  `EXP_POTENTIAL_V` (−2.0 V vs SHE; uncompensated potential, without iR-drop compensation). For barriers, `slope_x_RT` = 1 would be ideal
  transition-state-theory behavior.
- `correlation_data.csv`: ln(k_obs) and every descriptor, one row per compound and level
- `correlation_predictions.csv`: predicted vs observed ln(k) for compounds outside the fit set
- `correlation_pairs.csv`: two-descriptor models (`CORR_PAIRS`: ΔG_ET + ΔG_frag and ΔG_ET + ΔG(ArX•⁻ + e⁻ → Ar⁻ + X⁻))
- `plots/correlations_<level>.png`: the six best single-descriptor fits for each level
  (needs matplotlib)

Experimental rate constants are in `data/experimental_kobs.csv`, from LaPier et al.,
*Environ. Sci. Technol.* 2026, 60, 1346, Table 1. Two corrections were applied:
- The measurements were made at −2.0 V vs SHE; all calculations and correlations use this
  potential.
- The published 1,3-dibromobenzene value (0.79 ± 0.074) is per day; it was divided by 24
  to give 0.0329 h⁻¹, consistent with the 21 h half-life.

Five chlorobenzenes (group `chlorobenzene`: hexa-, penta-, 1,2,4,5-tetra-, 1,2- and
1,4-dichlorobenzene), also measured at −2.0 V vs SHE, were added on 2026-09-29 from
unpublished data (no standard errors yet). They enter the `all`
fits and are predicted by the bromobenzene fit, together with the PBDEs.

The `dG_ET_any` and `dG_frag_any` columns use the lowest radical-anion energy whatever
its structure (π, bent σ or dissociated), as in the published QSAR. The plain `dG_ET` and
`dG_frag` columns leave dissociated radical anions blank.
- `Haloaromatics_results.xlsx`: all of the above as sheets, plus the constants used

**Conventions (edit in `config.py`):**

- G(e⁻) = −0.867 kcal/mol (Fermi–Dirac)
- E_abs(SHE) = 4.281 V
- G(H⁺,aq) = −6.28 − 265.9 + 1.89 kcal/mol
- For SMD levels, +1.89 kcal/mol is added to every solute (1 atm → 1 M)
- E° = −ΔG/(nF) − E_abs(SHE)
- Potentials are not reported for gas-phase levels (only ΔG)

## Layout

```
config.py            methods, solvent, resources, constants
hx.py                command-line driver
haloaro/             species.py (enumeration/naming), stages.py (input writers, retry),
                     parse.py (log parser), qc.py, thermo.py, descriptors.py,
                     correlations.py, geomtools.py
data/experimental_kobs.csv          measured rate constants used by compile
slurm/g16_array.sbatch
calcs/<stage>/{inputs,logs,chks}/   manifests in calcs/<stage>/
results/
```
# Haloaromatics
