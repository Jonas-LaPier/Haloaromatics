# Interpretation of the radical-anion C–X scans (tsscan, M06-2X/6-311++G(d))

Draft, 2026-09-30. Data: relaxed C–X scans of every bound radical anion (16 points × 0.08 Å
from the radical-anion minimum, X tilted 15° out of the ring plane), gas phase and SMD(water).
Classification and profiles: `figures/ts_scan_analysis.py` → `results/ts_scan_summary.csv`,
`results/plots/ts_scans_m062x_<level>.png`. Energies are electronic, relative to the radical
anion, from converged scan points only; maxima of relaxed scans are upper-bound estimates, not
transition-state energies. Rate constants were measured at −2.0 V vs SHE (uncompensated
potential, without iR-drop compensation).

Literature statements below are paraphrased from passages read in the cited sources (reference
keys at the end); everything else is our result or our inference, and is marked as such.

## 1. What the scans show

**Most radical anions of the measured compounds have no barrier to C–X cleavage.**

| compounds with k_obs | gas phase | SMD (water) |
|---|---|---|
| bromobenzene, 1,3-, 1,4-dibromo, 1,2,4-, 1,3,5-tribromo, 1,2,4,5-tetrabromo | radical anion dissociates on optimisation (no scan) | dissociates, except 1,4-dibromo: π, cleaves with ≈1.2 kcal/mol |
| hexabromobenzene | π; energy rises to the end of the scan (12 kcal/mol at 3.2 Å) | bent σ; 4.0 kcal/mol barrier at 2.06 Å into an Ar•···Br⁻ complex |
| BDE-47 | dissociates | dissociates |
| BDE-99 | bent σ (C2–Br 2.59 Å); C2 cleavage essentially barrierless; other sites cross states (17–50 kcal/mol) or rise to 70 kcal/mol | bent σ (C2–Br 2.64 Å); C2 flat (≤ 2.3 kcal/mol to 3.9 Å); other sites cross (23–47) or rise (56–68) |
| 1,2-dichloro | dissociates | π; barrierless |
| 1,4-dichloro | π; state crossing at 2.03 Å, 3.4 kcal/mol | π; crossing at 2.01 Å, 5.1 kcal/mol |
| 1,2,4,5-tetrachloro | π; crossing at 2.02 Å, 3.0 kcal/mol | π; ≈1.0 kcal/mol |
| pentachloro | bent σ; 4–5 kcal/mol at C1/C2, then uphill | bent σ; crossings at 12–21 kcal/mol (C1, C2); C3 already elongated |
| hexachloro | π; rises to the end (11.7 kcal/mol) | π; ≈0.1 kcal/mol |

Across all 47 reactions at both levels (94 scans), only a handful of scans show a smooth
barrier of 4–9 kcal/mol (e.g. 1,2,3,4- and 1,2,3,4,5-tetra/pentabromobenzene in SMD,
1,2,3- and 1,2,3,4-chlorobenzenes in the gas phase), and none of them belongs to a compound
with a measured rate constant except hexabromobenzene (SMD, 4.0 kcal/mol).

**Three recurring shapes** (our classification):
1. *Barrierless or ≤ ~5 kcal/mol* at the most favourable site: cleavage follows electron uptake
   almost immediately.
2. *State crossings* at the other sites of polyhalogenated compounds: the energy rises steeply
   and then drops by 10–40 kcal/mol within one 0.08 Å step at r(C–X) ≈ 2.0–2.8 Å. We interpret
   this as the extra electron jumping to the bond being stretched (a diabatic crossing), so the
   maximum is a crossing point, not a saddle point; a TS search started there is not meaningful.
3. *Uphill to the end* in the gas phase for hexabromo-, hexachloro- and pentachlorobenzene:
   separating Ar• and X⁻ without solvation costs energy throughout the scan; in SMD the same
   cleavages are barrierless or small.

## 2. Interpretation

**(a) Electron transfer, not C–X cleavage, should limit the rate.** If the radical anion
dissociates on formation, or cleaves with ≤ ~5 kcal/mol, the cleavage step is fast and cannot
explain rate differences of more than three orders of magnitude between the compounds (our
inference). Andrieux et al. (1986) describe chlorobenzene and bromobenzene as fast-cleaving
aryl halides whose direct and mediated electrochemical reductions are kinetically controlled by
the forward electron transfer (measurements in DMF). This is consistent with the QSAR result
that thermodynamic descriptors of electron uptake and of the overall two-electron reaction
predict ln k_obs, and with the absence of usable stepwise barriers for the measured compounds.

**(b) Bromine vs chlorine.** The computed bromobenzene radical anions mostly dissociate,
whereas several chlorobenzene radical anions are bound π species that cleave over small
barriers or crossings (our result). Aromatic chlorides are described as forming relatively
stable, long-lived radical anions through their π orbitals and as tending to a stepwise
mechanism (Yu et al. 2025), and aryl halides in general can accommodate the incoming electron
in a low-energy π* orbital and usually follow a stepwise mechanism (Mazzucato et al. 2023).
The faster cleavage of C–Br than C–Cl is in line with reports that C–Br and C–I cleavage is
favoured over C–Cl (higher bromide and iodide yields; King & Mitch 2024) and that alkane
dehalogenation rates follow chlorine < bromine < iodine on activated-carbon cathodes (King &
Mitch 2022). The literature comparisons involve other substrates or solvents, as noted.

**(c) Stepwise vs concerted.** For the bromobenzenes and BDE-47 at these levels there is no
radical-anion minimum, which in the terminology of the main workbook means that only the
concerted pathway exists (or a stepwise path with a vanishingly short-lived intermediate). A
weak C–X bond favours the concerted mechanism (Mazzucato et al. 2023). For the chlorobenzenes,
bound π radical anions with small cleavage barriers correspond to a stepwise mechanism. The
LaPier et al. (2026) conclusion that bromobenzene reduction is stepwise rests on the QSAR
(radical-anion formation was the most predictive individual step there); our scans suggest the
radical anion of most bromobenzenes is at best a very shallow intermediate. These two views are
compatible if electron uptake into the π* system is the kinetically relevant step and cleavage
follows without a significant barrier (our inference).

**(d) Solvation.** The gas-phase scans of the most halogenated compounds rise throughout,
whereas SMD makes the same cleavages nearly barrierless. King & Mitch (2022) note that the
halide products must ultimately be solvated, a factor gas-phase reduction free energies do not
include, while gas-phase descriptors correlated best for alkanes sorbed on activated carbon,
possibly reflecting the non-aqueous character of the carbon phase. For cleavage barriers, the
SMD scans are therefore the more relevant of the two (our inference), even though gas-phase
thermodynamic descriptors give the better QSARs.

**(e) BDE-99 regioselectivity — an open question.** The computed BDE-99 radical anion holds the
extra electron in the C2–Br bond, and C2 cleavage is essentially barrierless at both levels,
which would give 2,3',4,4'-tetrabromodiphenyl ether (BDE-66) (our result). LaPier et al.
(2026) observed BDE-47 (loss of the C5 bromine) as a debromination product of BDE-99, together
with an unidentified tetrabrominated congener. Our C5 scan from the C2-localised radical anion
shows a state crossing at about 50 kcal/mol, so it does not explain BDE-47 formation. Possible
reasons (not tested): other radical-anion states or conformers (no conformer search was done
for the PBDEs), cleavage at the electrode surface, or a route not captured by a single-bond
scan. A BDE-66 standard would test whether the unidentified product is the C2-cleavage product.

## 3. What this means for the TS stage

- No TS optimisations are warranted for the measured compounds: their scans are barrierless,
  state crossings, or rise to the end. The few smooth barriers (4–9 kcal/mol) belong to
  compounds without rate constants; their TSs could be located if wanted (10–12 h each).
- Stepwise TS barriers therefore cannot serve as QSAR descriptors for this data set; the
  concerted (Savéant) barriers and the thermodynamic descriptors remain the usable kinetic
  descriptors.
- State crossings call for a different treatment than a single-surface TS search (e.g. a
  crossing-point or two-state model). The bending of the cleaving bond in aryl-halide π radical
  anions is the subject of Costentin, Robert & Savéant, J. Am. Chem. Soc. 2004, 126, 16051
  (cited by Mazzucato et al. 2023); that paper is not in the library, so its content is not used
  here. Adding it would allow a sourced discussion of the 15° tilt and of the crossings.

## References (read in the project's Zotero library)

- **Andrieux1986** Andrieux, C. P.; Savéant, J. M.; Su, K. B. J. Phys. Chem. 1986, 90, 3815. https://doi.org/10.1021/j100407a059
- **Mazzucato2023** Mazzucato, M.; Isse, A. A.; Durante, C. Curr. Opin. Electrochem. 2023. https://doi.org/10.1016/j.coelec.2023.101254
- **Yu2025** Yu et al. Electrochemical reduction for chlorinated hydrocarbons contaminated groundwater remediation: Mechanisms, challenges, and perspectives. Water Res. 2025 (Zotero key 47SIUCWM).
- **KingMitch2022** King, J. F.; Mitch, W. A. Environ. Sci. Technol. 2022, 56, 17965. https://doi.org/10.1021/acs.est.2c05608
- **KingMitch2024** King, J. F.; Mitch, W. A. Crit. Rev. Environ. Sci. Technol. 2024. https://doi.org/10.1080/10643389.2023.2239130
- **LaPier2026** LaPier, J. K. et al. Environ. Sci. Technol. 2026, 60, 1346. https://doi.org/10.1021/acs.est.5c03324
