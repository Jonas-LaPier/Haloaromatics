"""Quality checks on parsed Gaussian logs.

status:  ok | warn | fail | incomplete | missing
flags:   list of short machine-readable strings explaining the status
"""
from __future__ import annotations

import config as C
from .geomtools import distance, dot, norm, sub, unit, xyz


def expected_terminations(stage_kind, species):
    if stage_kind in ("am1", "tsscan", "sp"):
        return 1
    if species.is_atom:
        return 1
    return 2  # Opt + Freq (Link1)


def evaluate(parsed, species, stage_kind, site=None):
    """Return (status, flags). parsed may be None (no log)."""
    if parsed is None:
        return "missing", []
    flags = []
    fail = warn = False

    if parsed["error_term"]:
        flags.append(f"error:{parsed['error_type']}")
        fail = True
    need = expected_terminations(stage_kind, species)
    if not parsed["error_term"] and parsed["n_normal_term"] < need:
        if parsed["n_normal_term"] == 0:
            return "incomplete", ["no_termination (running, killed or time limit)"]
        flags.append(f"terminations {parsed['n_normal_term']}/{need}")
        fail = True

    if parsed["charge"] is not None and (parsed["charge"], parsed["mult"]) != (species.charge, species.mult):
        flags.append("charge_mult_mismatch")
        fail = True

    if parsed.get("geometry") and len(parsed["geometry"]) != len(species.atoms()):
        flags.append("atom_count_mismatch")
        fail = True

    is_opt = stage_kind in ("am1", "optfreq", "ts") and not species.is_atom
    if is_opt and not parsed["stationary"] and not parsed["error_term"]:
        flags.append("no_stationary_point")
        fail = True

    if stage_kind in ("optfreq", "ts"):
        if parsed["G"] is None:
            flags.append("no_free_energy")
            fail = True
        want = 1 if stage_kind == "ts" else 0
        if parsed["freqs"] and parsed["n_imag"] != want:
            flags.append(f"n_imag={parsed['n_imag']} (want {want}; lowest {parsed['lowest_freq']:.1f})")
            fail = True

    # spin contamination
    # (skipped for AM1: semi-empirical UHF radicals are routinely contaminated and the
    #  AM1 stage only supplies starting geometries)
    if species.mult == 2 and parsed["S2"] is not None and stage_kind != "am1":
        if abs(parsed["S2"] - 0.75) > C.S2_TOL:
            flags.append(f"spin_contam S2={parsed['S2']:.3f}")
            warn = True

    geom = parsed.get("geometry")
    # radical anion fell apart during optimisation?
    if species.kind == "radical_anion" and geom and stage_kind != "am1":
        lim = C.CX_DISSOCIATED_A[species.halogen]
        for s in range(1, 7):
            ci, xi = species.atom_index(s)
            if xi and geom[xi - 1][0] == species.halogen and distance(geom, ci, xi) > lim:
                flags.append(f"RA_dissociated(C{s}-X {distance(geom, ci, xi):.2f} A)")
                warn = True

    # TS: imaginary mode should be the C-X stretch
    if species.kind == "ts" and stage_kind == "ts" and geom and parsed["mode1"] and site:
        ci, xi = species.atom_index(site)
        u = unit(sub(xyz(geom[xi - 1]), xyz(geom[ci - 1])))
        rel = sub(parsed["mode1"][xi - 1], parsed["mode1"][ci - 1])
        big = max(norm(v) for v in parsed["mode1"]) or 1
        if abs(dot(rel, u)) / big < 0.3:
            flags.append("ts_mode_not_CX_stretch")
            warn = True
        flags.append(f"r(C-X)={distance(geom, ci, xi):.3f}")

    status = "fail" if fail else ("warn" if warn else "ok")
    return status, flags


def usable(status):
    return status in ("ok", "warn")
