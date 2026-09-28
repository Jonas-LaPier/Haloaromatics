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
| radical anions ArX•⁻ (stepwise electron transfer) | −1 / 2 | 24 |
| aryl radicals Ar• (one per unique C–X bond; phenyl shared by both series) | 0 / 2 | 39 |
| aryl carbanions Ar⁻ | −1 / 1 | 39 |
| Cl⁻, Br⁻ | −1 / 1 | 2 |
| Cl•, Br• atoms (for C–X bond energies) | 0 / 2 | 2 |

That gives 40 unique dehalogenation reactions (20 per halogen). The reaction table
also lists each reaction's degeneracy, meaning how many equivalent C–X bonds it covers.

**Names.** `ClBz_124` is 1,2,4-trichlorobenzene and `ClBz_124_RA` is its radical anion.
`ClPh_24_rad` and `ClPh_24_anion` are the 2,4-dichlorophenyl radical and carbanion,
with the radical/anion carbon numbered C1. Atom order in every input is C1…C6 followed
by the substituents in site order.

## Stages

| stage | what | geometry from |
|---|---|---|
| `am1` | Opt AM1 (every molecular species) | idealised planar ring |
| `b3lyp_gas`, `b3lyp_smd`, `m062x_gas`, `m062x_smd` (group `optfreq`) | Opt + Freq, 6-311++G(d), gas and SMD(water) | `Geom=Check` from AM1 chk |
## Activation barriers: Savéant concerted dissociative electron transfer

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

The stepwise radical-anion TS stages (a C–X scan of ArX•⁻, then Opt=TS) are still in the
code but switched off. Set `RUN_RA_TS = True` in `config.py` to bring back the
`tsscan_*` and `ts_*` stages.

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
- There are no imaginary frequencies at minima (exactly 1 for TSs, if `RUN_RA_TS` is on)
- ⟨S²⟩ is within 0.10 of 0.75 for doublets
- **Radical anion dissociated**: C–X > 2.3 Å (Cl) or 2.5 Å (Br). The radical anion has no
  bound minimum at that level, so reduction there is concerted. Its row in the reaction
  table is flagged.

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
- `det_barriers.csv`: C–X bond enthalpy and free energy (at this level and in gas phase), E°_DET, radius a, λ₀, ΔG₀‡, and ΔG‡ and α at each potential in `DET_POTENTIALS_V`
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
# Haloaromatics
