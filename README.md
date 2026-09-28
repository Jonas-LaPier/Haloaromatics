# Haloaromatics

Electronic-structure workflow (Gaussian 16 on Sherlock) for the reductive
dehalogenation of chloro- and bromobenzenes. Everything runs through one
driver script, `hx.py`. Methods, solvent, resources and constants are set in `config.py`.

## Chemical scope

The code enumerates every symmetry-unique isomer and every symmetry-unique C–X bond
(`python3 hx.py species` writes `results/species_list.csv` and `results/reaction_list.csv`):

| species | charge / mult | count |
|---|---|---|
| parents ArX: benzene + 12 chloro + 12 bromo | 0 / 1 | 25 |
| radical anions ArX•⁻ (reactant for the C–X cleavage TS) | −1 / 2 | 24 |
| aryl radicals Ar• (one per unique C–X bond; phenyl shared by both series) | 0 / 2 | 39 |
| aryl carbanions Ar⁻ | −1 / 1 | 39 |
| Cl⁻, Br⁻ | −1 / 1 | 2 |
| C–X cleavage transition states | −1 / 2 | 40 |

That gives 40 unique dehalogenation reactions (20 per halogen). The reaction table
also lists each reaction's degeneracy, meaning how many equivalent C–X bonds it covers.

**Names.** `ClBz_124` is 1,2,4-trichlorobenzene and `ClBz_124_RA` is its radical anion.
`ClPh_24_rad` and `ClPh_24_anion` are the 2,4-dichlorophenyl radical and carbanion,
with the radical/anion carbon numbered C1. `TS_ClBz_124_x4` is the TS for breaking the
C4–Cl bond in `ClBz_124_RA`. Atom order in every input is C1…C6 followed by the
substituents in site order.

## Stages

| stage | what | geometry from |
|---|---|---|
| `am1` | Opt AM1 (every molecular species) | idealised planar ring |
| `b3lyp_gas`, `b3lyp_smd`, `m062x_gas`, `m062x_smd` (group `optfreq`) | Opt + Freq, 6-311++G(d), gas and SMD(water) | `Geom=Check` from AM1 chk |
| `tsscan_m062x_gas`, `tsscan_m062x_smd` (group `tsscan`) | relaxed C–X scan of ArX•⁻ (16 × 0.08 Å) | DFT radical-anion minimum, X tilted 15° out of plane |
| `ts_m062x_gas`, `ts_m062x_smd` (group `ts`) | Opt=(TS,CalcFC,NoEigenTest) + Freq | highest point of the scan |

**TS level.** M06-2X/6-311++G(d) was chosen because it is well benchmarked for barrier
heights, the diffuse functions are needed for the anions, and the reactant radical anions
are already computed at that level in the `optfreq` stage. The out-of-plane tilt breaks
planar symmetry so the π* and σ* states can mix. Otherwise the planar scan crosses
between the two states at a cusp and never reaches a smooth saddle point. To change
the TS level, edit `TS_LEVELS` in `config.py`.

## Workflow on Sherlock

```bash
cd ~/Haloaromatics            # project root; all paths are relative to here
python3 hx.py generate am1
python3 hx.py submit am1      # one job array; add --dry-run to only print the sbatch line
python3 hx.py status am1      # QC summary: ok / warn / fail / incomplete / missing
python3 hx.py retry am1       # rebuild failed inputs with automatic fixes, then submit again

python3 hx.py generate optfreq      # only species whose AM1 job passed QC
python3 hx.py submit optfreq
python3 hx.py status optfreq

python3 hx.py generate tsscan && python3 hx.py submit tsscan
python3 hx.py generate ts     && python3 hx.py submit ts

python3 hx.py scrape all      # logs -> results/raw/<stage>.csv
python3 hx.py compile         # -> results/*.csv and results/Haloaromatics_results.xlsx
```

**Tip:** before submitting all 129 jobs of a new level, test one bromine species first:
`python3 hx.py submit m062x_gas --only BrBz_1`. This confirms that 6-311++G(d) is
defined for Br in your Gaussian build.

Details:

- `submit` writes `calcs/<stage>/manifest_<time>.txt`, which lists every input that has
  no log yet, and submits `slurm/g16_array.sbatch` as one array (throttle `%50`). This
  means resubmitting never reruns finished jobs.
- Gaussian scratch goes to `$L_SCRATCH` (node-local). Set `FORMCHK=1` before `sbatch`
  to also write `.fchk` files.
- `retry` moves the failed log and input to `calcs/<stage>/logs/failed/<name>.tryN.*`
  and writes a new input:
  - SCF failure → `SCF=(XQC,MaxCycle=512)`
  - optimisation ran out of cycles or time → restart from its own chk with `Geom=Check Guess=Read`
  - internal-coordinate errors → `Opt=Cartesian`
  - an imaginary frequency at a minimum → displace along that mode and reoptimise with `Opt=(CalcFC,Tight)`
- `--only NAME ...` limits any stage command to specific jobs. `--force` regenerates inputs.

### Quality checks (`haloaro/qc.py`)

- Normal termination count (2 for Opt+Freq)
- Charge, multiplicity and atom count match the species
- A stationary point was found and the free energy is present
- The number of imaginary frequencies is 0 for minima and 1 for TSs
- ⟨S²⟩ is within 0.10 of 0.75 for doublets
- **Radical anion dissociated**: C–X > 2.3 Å (Cl) or 2.5 Å (Br). If a radical anion has
  no bound minimum at a level, cleavage is barrierless there, so no TS scan is generated
  and the reaction row is flagged.
- The TS imaginary mode contains the C–X stretch

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

  Each gives ΔG in kcal/mol and E° in V vs SHE (SMD levels only).
- `ts_barriers.csv`: ΔG‡ and ΔE‡ = TS − ArX•⁻, ΔG‡ relative to ArX + e⁻, imaginary frequency, QC flags
- `qc_issues.csv`: every job whose status is not `ok`
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
                     parse.py (log parser), qc.py, thermo.py, geomtools.py
slurm/g16_array.sbatch
calcs/<stage>/{inputs,logs,chks}/   manifests in calcs/<stage>/
results/
```
