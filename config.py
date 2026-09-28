"""Project configuration: levels of theory, stages, cluster resources, constants.

Edit this file (not the code in haloaro/) to change methods, solvents or
resources. Every stage writes to calcs/<stage>/{inputs,logs,chks}.
"""

# --------------------------------------------------------------------------- #
# Chemistry
# --------------------------------------------------------------------------- #
HALOGENS = ("Cl", "Br")
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

# Levels used for the C-X cleavage transition states (ArX.- -> [Ar...X]‡.- -> Ar. + X-).
# M06-2X is a well benchmarked meta-hybrid for barrier heights; diffuse functions
# are required for the anions. Reactant radical anions come from the same level
# in the opt/freq stage so barriers are internally consistent.
TS_LEVELS = ["m062x_gas", "m062x_smd"]

# Relaxed C-X scan used to locate the TS guess
SCAN_STEPS = 16          # number of steps
SCAN_STEP_SIZE = 0.08    # Angstrom per step (total +1.28 A from the RA minimum)
OOP_ANGLE_DEG = 15.0     # out-of-plane tilt of X applied to the start geometry
                         # (breaks planarity so pi*/sigma* states can mix)

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
CX_DISSOCIATED_A = {"Cl": 2.30, "Br": 2.50}   # C-X distance => RA fell apart
