"""Project configuration: levels of theory, stages, cluster resources, constants.

Edit this file (not the code in haloaro/) to change methods, solvents or
resources. Every stage writes to calcs/<stage>/{inputs,logs,chks}.
"""

# --------------------------------------------------------------------------- #
# Chemistry
# --------------------------------------------------------------------------- #
HALOGENS = ("Cl", "Br")

# Polybrominated diphenyl ethers: label -> (ring A Br locants, ring B Br locants).
# Their hydrodebromination products (e.g., BDE-28, BDE-17 from BDE-47) are added as
# neutral species automatically. Conformers are not searched: each starts from a twisted
# C2-like geometry, so consider a CREST/xTB conformer search before trusting small
# energy differences.
PBDES = {"BDE-47": ("24", "24"), "BDE-99": ("245", "24")}
BASIS = "6-311++G(d)"
SOLVENT = "Water"                      # used for every *_smd level
SMD = f"SCRF=(SMD,Solvent={SOLVENT})"

# Level-of-theory definitions. `solv` is None (gas) or an SCRF keyword string.
LEVELS = {
    "b3lyp_gas": {"method": "B3LYP",  "basis": BASIS, "solv": None},
    "b3lyp_smd": {"method": "B3LYP",  "basis": BASIS, "solv": SMD},
    "m062x_gas": {"method": "M062X",  "basis": BASIS, "solv": None},
    "m062x_smd": {"method": "M062X",  "basis": BASIS, "solv": SMD},
}

# --------------------------------------------------------------------------- #
# Activation barriers: Saveant concerted dissociative electron transfer (DET)
#     ArX + e-  ->  [ArX]‡  ->  Ar. + X-          (one step, from the NEUTRAL parent)
#     dG0‡ = (D + lambda0) / 4                        intrinsic barrier
#     dG‡(E) = dG0‡ * (1 + dG°(E) / (4 dG0‡))^2      dG°(E) = F (E - E°_DET) ... see thermo.py
# D = homolytic C-X bond energy of the neutral parent (ArX -> Ar. + X.), from the
# opt/freq stage (needs the X. atoms). lambda0 = solvent reorganisation (Marcus-Hush
# one-sphere model, electrode reaction). Computed for the SMD levels only.
# J.-M. Saveant, J. Am. Chem. Soc. 1987, 109, 6788; Acc. Chem. Res. 1993, 26, 455.
# --------------------------------------------------------------------------- #
SAVEANT_D = "H"              # "H": bond dissociation enthalpy (Saveant's D); "G": BDFE
SPIN_ORBIT_KCAL = {"Cl": 0.84, "Br": 3.51}   # 2P3/2 ground-state stabilisation of X.
                                             # (1/3 of the 2P splitting), subtracted from G, H of X.
EPS_STATIC = 78.36           # water (must match SOLVENT)
EPS_OPTICAL = 1.777          # n^2 of water
LAMBDA0_MODEL = "electrode"  # "electrode": e^2/(8 pi eps0 a)(1/eps_op - 1/eps_s) (image term neglected)
                             # "homogeneous": 2x electrode value (self-exchange-like, donor of equal size)
LAMBDA0_KCAL = None          # set a number to override the computed lambda0 for every compound
RADIUS_PROBE_A = 0.0         # added to the van der Waals sphere-equivalent radius a
DET_POTENTIALS_V = [-1.0, -1.5, -2.0]        # electrode/donor potentials (V vs SHE) for dG‡(E)

# Experimental rate constants for correlation (compile -> results/correlations.csv)
EXPERIMENTAL_KOBS = "data/experimental_kobs.csv"
EXP_POTENTIAL_V = -2.0          # potential at which the kobs were measured (V vs SHE; uncompensated, no iR-drop compensation)
CORR_FIT_GROUP = "bromobenzene" # fit on this group, then predict the others (e.g. PBDEs)
# Compound sets fitted separately (group names from EXPERIMENTAL_KOBS; "halobenzene" means
# bromobenzene + chlorobenzene, "all" means every compound with data)
CORR_FIT_SETS = ["bromobenzene", "chlorobenzene", "halobenzene", "all"]
# Two-descriptor MLR models: electron transfer to ArX combined with the radical-anion step
CORR_PAIRS = [("dG_ET_any_kcal", "min_dG_frag_any_kcal"),
              ("dG_ET_any_kcal", "min_dG_RA_2e_any_kcal")]

# QSAR summary (results/qsar_summary.csv, sheet QSAR_summary): (descriptor, label, set), fitted
# against ln(kobs) at every level and fit set. Set "thermo" holds the reduction thermochemistry
# and frontier orbitals; set "extended" the other parameters expected to matter for reduction
# kinetics. Site quantities use the most favourable site; *_any uses the lowest radical-anion
# energy whatever its state. Barriers "@ V" are evaluated at EXP_POTENTIAL_V.
_AT = f"@ {EXP_POTENTIAL_V:+.2f} V"      # suffix of potential-dependent descriptor names
QSAR_DESCRIPTORS = [
    ("LUMO_eV",                  "LUMO of ArX",                        "thermo"),
    ("LUMO_RA_eV",               "LUMO of the radical anion ArX.-",              "thermo"),
    ("min_dG_1e_kcal",           "dG ArX + e- -> Ar. + X-",                "thermo"),
    ("min_dG_2e_carbanion_kcal", "dG ArX + 2e- -> Ar- + X-",               "thermo"),
    ("dG_ET_any_kcal",           "dG ArX + e- -> ArX.-",                   "thermo"),
    ("min_dG_frag_any_kcal",     "dG ArX.- -> Ar. + X-",                   "thermo"),
    ("min_dG_RA_2e_any_kcal",    "dG ArX.- + e- -> Ar- + X-",              "thermo"),
    ("n_X",                      "number of halogens",                           "extended"),
    ("VEA_eV",                   "vertical electron affinity",                   "extended"),
    ("AEA_elec_eV",              "adiabatic electron affinity (electronic)",     "extended"),
    ("omega_dSCF_eV",            "electrophilicity index (dSCF)",                "extended"),
    ("E_ET_V",                   "E deg ArX/ArX.- (SMD levels only)",            "extended"),
    ("dG_ET_kcal",               "dG ArX + e- -> ArX.-, bound RA only",           "extended"),
    ("min_dG_frag_kcal",         "dG ArX.- -> Ar. + X-, bound RA only",           "extended"),
    ("min_dG_RA_2e_kcal",        "dG ArX.- + e- -> Ar- + X-, bound RA only",           "extended"),
    ("min_dG_2e_HDH_kcal",       "dG ArX + H+ + 2e- -> ArH + X-",                "extended"),
    ("min_dG_rad_red_kcal",      "dG Ar. + e- -> Ar-",                           "extended"),
    ("lambda_i_kcal",            "inner-sphere reorganisation energy",           "extended"),
    ("min_BDE_kcal",             "weakest C-X bond enthalpy",                    "extended"),
    ("min_dG0_act_concerted_kcal", "Saveant intrinsic barrier (weakest C-X)",    "extended"),
    (f"eff_dG_act_concerted_kcal {_AT}", "Saveant concerted barrier",         "extended"),
    (f"eff_dG_act_stepwise_kcal {_AT}",  "stepwise ET + cleavage barrier",    "extended"),
    (f"eff_dG_act_combined_kcal {_AT}",  "combined stepwise + concerted",     "extended"),
    ("eff_dG_act_frag_TS_kcal",  "C-X cleavage TS barrier of ArX.-",             "extended"),
    ("max_fplus_CX",             "largest Fukui f+ on a C-X bond",               "extended"),
    ("min_wiberg_CX_parent",     "weakest Wiberg C-X bond order",                "extended"),
    ("max_spin_RA_X",            "largest halogen spin density in ArX.-",        "extended"),
    ("max_dr_CX_RA_A",           "largest C-X elongation in ArX.-",              "extended"),
]

# Stepwise radical-anion TS workflow (tsscan_* / ts_* stages): ArX.- -> [Ar...X]‡.- -> Ar. + X-.
# Run alongside the Saveant concerted analysis; compile compares the two pathways.
RUN_RA_TS = True

# Levels used for the radical-anion C-X cleavage TS (only if RUN_RA_TS)
TS_LEVELS = ["m062x_gas", "m062x_smd"] if RUN_RA_TS else []

# Relaxed C-X scan used to locate the TS guess
SCAN_STEPS = 16          # number of steps
SCAN_STEP_SIZE = 0.08    # Angstrom per step (total +1.28 A from the RA minimum)
OOP_ANGLE_DEG = 15.0     # out-of-plane tilt of X applied to the start geometry
                         # (breaks planarity so pi*/sigma* states can mix)

# Single-point stage (vertical anion/cation/neutral + population analysis) on the
# optimised geometries of these levels. At SMD levels the vertical energies use
# equilibrium solvation; lambda_i for the pathway comparison is taken from the gas-phase
# counterpart of each functional when LAMBDA_I_FROM_GAS is True.
SP_LEVELS = list(LEVELS)
LAMBDA_I_FROM_GAS = True

# Extra route keywords common to all DFT jobs
DFT_EXTRA = "Int=UltraFine"
AM1_SCF = "SCF=XQC"                    # route: #p Opt AM1 <AM1_SCF>

# --------------------------------------------------------------------------- #
# Stage -> resources on Sherlock.  mem is total GB per job; %mem in the .gjf is
# set to MEM_FRACTION of that to leave headroom for Gaussian overhead.
# --------------------------------------------------------------------------- #
PARTITION = None          # e.g. "normal" or your group partition; None = cluster default
MEM_FRACTION = 0.80
ARRAY_THROTTLE = 50       # max simultaneously running array tasks (%N)

RESOURCES = {
    "am1":     {"cpus": 1, "mem_gb": 2,  "time": "0-00:20:00"},
    "optfreq": {"cpus": 8, "mem_gb": 16, "time": "0-12:00:00"},
    "tsscan":  {"cpus": 8, "mem_gb": 16, "time": "1-00:00:00"},
    "ts":      {"cpus": 8, "mem_gb": 16, "time": "1-00:00:00"},
    "sp":      {"cpus": 8, "mem_gb": 16, "time": "0-04:00:00"},
}
# Larger resources for diphenyl ethers (23-25 atoms, 4-5 Br)
RESOURCES_DPE = {
    "am1":     {"cpus": 1, "mem_gb": 2,  "time": "0-01:00:00"},
    "optfreq": {"cpus": 16, "mem_gb": 48, "time": "0-12:00:00"},
    "tsscan":  {"cpus": 16, "mem_gb": 48, "time": "0-12:00:00"},
    "ts":      {"cpus": 16, "mem_gb": 48, "time": "0-12:00:00"},
    "sp":      {"cpus": 16, "mem_gb": 32, "time": "0-08:00:00"},
}

# --------------------------------------------------------------------------- #
# Thermochemistry constants
# --------------------------------------------------------------------------- #
HARTREE_TO_KCAL = 627.509474
HARTREE_TO_EV = 27.211386
FARADAY_KCAL = 23.060548          # kcal mol-1 V-1
TEMPERATURE = 298.15

# 1 atm -> 1 M standard-state correction applied to every solute in SMD levels
STD_STATE_CORR_KCAL = 1.894

# Electron: gas-phase free energy, Fermi-Dirac convention (Bartmess 1994)
G_ELECTRON_KCAL = -0.867
# Absolute potential of the SHE consistent with the above conventions
# (Isse & Gennaro, J. Phys. Chem. B 2010, 114, 7894). IUPAC value is 4.44 V.
E_ABS_SHE_V = 4.281

# Proton: G(H+, gas, 1 atm) and aqueous solvation free energy (Tissandier 1998)
G_PROTON_GAS_KCAL = -6.28
DG_SOLV_PROTON_KCAL = -265.9

# QC thresholds
S2_TOL = 0.10                     # |<S^2> - s(s+1)| above this is flagged
CX_DISSOCIATED_A = {"Cl": 2.30, "Br": 2.50}   # planar C-X beyond this => Ar...X- complex
RA_PI_MAX_A = {"Cl": 2.00, "Br": 2.15}        # C-X up to this => intact pi radical anion
RA_SIGMA_MAX_A = {"Cl": 2.80, "Br": 3.00}     # bent C-X up to this => sigma-type radical anion
RA_SIGMA_OOP_DEG = 10.0                       # X out-of-plane angle marking a bent sigma RA
IMAG_TOL_CM = 20.0                            # |imaginary| below this at a minimum => warn, not fail
MAX_RETRIES = 3                               # retry refuses jobs with this many failed tries
