"""Gaussian 16 log file parser (pure Python, no dependencies)."""
from __future__ import annotations

import re
from pathlib import Path

SYMBOLS = {1: "H", 6: "C", 8: "O", 17: "Cl", 35: "Br"}

RE_SCF = re.compile(r"SCF Done:\s+E\((\S+)\)\s+=\s+(-?\d+\.\d+)")
RE_G = re.compile(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)")
RE_H = re.compile(r"Sum of electronic and thermal Enthalpies=\s+(-?\d+\.\d+)")
RE_ZPE = re.compile(r"Zero-point correction=\s+(-?\d+\.\d+)")
RE_GCORR = re.compile(r"Thermal correction to Gibbs Free Energy=\s+(-?\d+\.\d+)")
RE_S2 = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+),\s+after\s+(\d+\.\d+)")
RE_CHMULT = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
RE_NUM = re.compile(r"-?\d+\.\d+")
RE_ELAPSED = re.compile(r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+(\d+\.?\d*) seconds")

KNOWN_ERRORS = [
    ("Convergence failure -- run terminated", "scf_convergence"),
    ("No lower point found -- run aborted", "scf_convergence"),   # l508 (QC/XQC) failure
    ("SCF has not converged", "scf_convergence"),                 # l701/l801 after SCF failure
    ("Number of steps exceeded", "opt_maxcycles"),
    ("Optimization stopped", "opt_not_converged"),
    ("FormBX had a problem", "internal_coords"),
    ("Error in internal coordinate system", "internal_coords"),
    ("Linear angle in", "internal_coords"),
    ("Bend failed for angle", "internal_coords"),
    ("galloc:  could not allocate memory", "memory"),
    ("Erroneous write", "disk_or_scratch"),
    ("No such file or directory", "missing_file"),
    ("Problem with the distance matrix", "geometry"),
    ("Small interatomic distances", "geometry"),
    ("The combination of multiplicity", "charge_mult"),
    ("Wanted a", "missing_chk_or_route"),
    ("DUE TO TIME LIMIT", "time_limit"),
]


def _last(regex, text, group=1, cast=float):
    m = None
    for m in regex.finditer(text):
        pass
    return cast(m.group(group)) if m else None


def _geom_block(lines, start):
    """Parse an orientation table starting at header line index `start`."""
    atoms = []
    i = start + 5
    while i < len(lines) and not lines[i].strip().startswith("-----"):
        p = lines[i].split()
        z = int(p[1])
        atoms.append((SYMBOLS.get(z, str(z)), float(p[3]), float(p[4]), float(p[5])))
        i += 1
    return atoms


def last_geometry(lines):
    idx = None
    for i, l in enumerate(lines):
        if "Standard orientation:" in l or "Input orientation:" in l:
            idx = i
    return _geom_block(lines, idx) if idx is not None else None


def orbital_energies(lines):
    """HOMO/LUMO (Hartree) from the last population analysis block."""
    last = None
    for i, l in enumerate(lines):
        if l.startswith(" Alpha  occ. eigenvalues"):
            last = i
    if last is None:
        return {}
    start = last
    while start > 0 and lines[start - 1].startswith(" Alpha  occ. eigenvalues"):
        start -= 1
    vals = {"ao": [], "av": [], "bo": [], "bv": []}
    keys = {" Alpha  occ.": "ao", " Alpha virt.": "av", "  Beta  occ.": "bo", "  Beta virt.": "bv"}
    i = start
    while i < len(lines):
        k = next((v for p, v in keys.items() if lines[i].startswith(p)), None)
        if k is None:
            break
        vals[k] += [float(x) for x in RE_NUM.findall(lines[i].split("--", 1)[1])]
        i += 1
    out = {}
    homos = [v[-1] for v in (vals["ao"], vals["bo"]) if v]
    lumos = [v[0] for v in (vals["av"], vals["bv"]) if v]
    if homos:
        out["homo"] = max(homos)
    if lumos:
        out["lumo"] = min(lumos)
    if vals["ao"]:
        out["homo_alpha"] = vals["ao"][-1]
    if vals["av"]:
        out["lumo_alpha"] = vals["av"][0]
    if vals["bo"]:
        out["homo_beta"] = vals["bo"][-1]
    if vals["bv"]:
        out["lumo_beta"] = vals["bv"][0]
    return out


def frequencies(lines):
    """(list of freqs, displacement vectors of mode 1) from the last freq job."""
    start = None
    for i, l in enumerate(lines):
        if "Harmonic frequencies (cm**-1)" in l:
            start = i
    if start is None:
        return [], None
    freqs, mode1 = [], None
    i = start
    while i < len(lines):
        l = lines[i]
        if "Thermochemistry" in l:
            break
        if l.strip().startswith("Frequencies --"):
            freqs += [float(x) for x in RE_NUM.findall(l.split("--", 1)[1])]
            if mode1 is None:
                j = i
                while j < len(lines) and "Atom  AN" not in lines[j]:
                    j += 1
                mode1 = []
                j += 1
                while j < len(lines):
                    p = lines[j].split()
                    if len(p) < 5 or not p[0].isdigit():
                        break
                    mode1.append((float(p[2]), float(p[3]), float(p[4])))
                    j += 1
        i += 1
    return freqs, mode1


def scan_points(lines):
    """Relaxed scan: list of dicts {energy, geometry, converged}."""
    pts, geom, e = [], None, None
    for i, l in enumerate(lines):
        if "Standard orientation:" in l or "Input orientation:" in l:
            geom = _geom_block(lines, i)
        elif "SCF Done:" in l:
            m = RE_SCF.search(l)
            if m:
                e = float(m.group(2))
        elif "Optimization completed" in l or "Optimization stopped" in l:
            pts.append({"energy": e, "geometry": geom,
                        "converged": "completed" in l})
    return pts


def hirshfeld(lines):
    """Last Hirshfeld/CM5 table -> list of dicts {sym, q_H, s_H, q_CM5}."""
    start = None
    for i, l in enumerate(lines):
        if l.startswith(" Hirshfeld charges, spin densities, dipoles, and CM5 charges"):
            start = i
    if start is None:
        return None
    out = []
    for l in lines[start + 2:]:
        p = l.split()
        if len(p) < 8 or not p[0].isdigit():
            break
        out.append({"sym": p[1], "q_H": float(p[2]), "s_H": float(p[3]), "q_CM5": float(p[7])})
    return out


def wiberg(lines, natoms):
    """Wiberg bond index matrices (NAO basis), in the order printed: total density first,
    then (open shell only) alpha and beta spin orbitals."""
    mats, i = [], 0
    starts = [k for k, l in enumerate(lines) if "Wiberg bond index matrix in the NAO basis" in l]
    for k in starts:
        m = [[0.0] * natoms for _ in range(natoms)]
        cols = None
        j = k + 1
        while j < len(lines):
            l = lines[j]
            if "Totals by atom" in l or "Atom-atom overlap" in l or "Wiberg bond index matrix" in l:
                break
            t = l.split()
            if t and t[0] == "Atom" and all(x.isdigit() for x in t[1:]):
                cols = [int(x) - 1 for x in t[1:]]
            elif cols and t and t[0].endswith(".") and t[0][:-1].isdigit():
                r = int(t[0][:-1]) - 1
                vals = [float(x) for x in t[2:2 + len(cols)]]
                for c, v in zip(cols, vals):
                    if r < natoms and c < natoms:
                        m[r][c] = v
            j += 1
        mats.append(m)
    return mats


def parse_log(path):
    path = Path(path)
    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    d = {"file": str(path), "name": path.stem}
    d["n_normal_term"] = text.count("Normal termination of Gaussian")
    d["error_term"] = "Error termination" in text
    if d["error_term"]:
        # Reason of an error termination: the known message that occurs last after the last
        # successful step. l9999 prints a long summary before terminating (over 800 lines for
        # the PBDEs), so the message can sit far above the end; searching only after the last
        # "Normal termination" keeps messages of an earlier, successful step out.
        start = text.rfind("Normal termination of Gaussian")
        seg = text[start:] if start >= 0 else text
        hits = [(seg.rfind(pat), code) for pat, code in KNOWN_ERRORS if pat in seg]
        d["error_type"] = max(hits)[1] if hits else "unknown"
    else:
        # jobs without an error line (running or killed): messages near the end only
        tail = "\n".join(lines[-800:])
        d["error_type"] = next((code for pat, code in KNOWN_ERRORS if pat in tail), None)
    m = None
    for m in RE_SCF.finditer(text):
        pass
    d["scf_method"] = m.group(1) if m else None
    d["E_scf"] = float(m.group(2)) if m else None
    d["G"] = _last(RE_G, text)
    d["H"] = _last(RE_H, text)
    d["ZPE"] = _last(RE_ZPE, text)
    d["G_corr"] = _last(RE_GCORR, text)
    d["S2"] = _last(RE_S2, text, 1)
    d["S2_annihilated"] = _last(RE_S2, text, 2)
    cm = RE_CHMULT.search(text)
    d["charge"], d["mult"] = (int(cm.group(1)), int(cm.group(2))) if cm else (None, None)
    d["stationary"] = "Stationary point found" in text
    d.update(orbital_energies(lines))
    freqs, mode1 = frequencies(lines)
    d["freqs"] = freqs
    d["n_imag"] = sum(1 for f in freqs if f < 0)
    d["lowest_freq"] = min(freqs) if freqs else None
    d["mode1"] = mode1
    d["geometry"] = last_geometry(lines)
    d["hirshfeld"] = hirshfeld(lines)
    nat = len(d["geometry"]) if d["geometry"] else 0
    wb = wiberg(lines, nat) if nat else []
    # NBO prints the total-density Wiberg matrix first; open-shell runs then repeat it for
    # the alpha and beta spin orbitals (spin-resolved values, not halves of the total)
    d["wiberg"] = wb[0] if wb else None
    dip = None
    for k, l in enumerate(lines):
        if l.startswith(" Dipole moment (field-independent basis, Debye)") and k + 1 < len(lines):
            m = re.search(r"Tot=\s*(-?\d+\.\d+)", lines[k + 1])
            if m:
                dip = float(m.group(1))
    d["dipole_D"] = dip
    pol = None
    for l in lines:
        if l.strip().startswith("Exact polarizability:"):
            v = [float(x) for x in RE_NUM.findall(l.split(":", 1)[1])]
            if len(v) >= 6:
                pol = (v[0] + v[2] + v[5]) / 3
    d["polar_iso_bohr3"] = pol
    secs = 0.0
    for m in RE_ELAPSED.finditer(text):
        dd, hh, mm, ss = m.groups()
        secs += int(dd) * 86400 + int(hh) * 3600 + int(mm) * 60 + float(ss)
    d["wall_hours"] = round(secs / 3600, 3) if secs else None
    d["_lines"] = lines
    return d
