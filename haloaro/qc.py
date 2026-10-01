"""Quality checks on parsed Gaussian logs.

status:  ok | warn | fail | incomplete | missing
flags:   list of short machine-readable strings explaining the status
"""
from __future__ import annotations

import config as C
from .geomtools import distance, dot, norm, oop_angle, sub, unit, xyz


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
        # an earlier step finished (e.g. Opt of Opt+Freq) and no error line: the next step is
        # still running or was killed, so this is not a failure that retry should act on
        return "incomplete", [f"terminations {parsed['n_normal_term']}/{need} "
                              "(later step running, killed or time limit)"]

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
            imag = sorted(f for f in parsed["freqs"] if f < 0)
            # tiny imaginary modes of floppy radical anions are integration-grid noise
            extra = imag[want:] if parsed["n_imag"] > want else imag
            if parsed["n_imag"] > want and all(abs(f) < C.IMAG_TOL_CM for f in extra):
                flags.append(f"small_imag {', '.join(f'{f:.1f}' for f in extra)} cm-1")
                warn = True
            else:
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
    # radical-anion character from its longest C-X bond
    if species.kind == "radical_anion" and geom and stage_kind != "am1":
        state, site_, r, oop = ra_state(geom, species)
        if state == "dissociated":
            flags.append(f"RA_dissociated(C{site_}-X {r:.2f} A)")
            warn = True
        elif state == "sigma_bent":
            flags.append(f"RA_sigma(C{site_}-X {r:.2f} A, {oop:.0f} deg)")
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


def ra_state(geom, species):
    """Classify a radical-anion geometry from its longest C-X bond.

    pi          r <= RA_PI_MAX_A                    (intact pi radical anion)
    sigma_bent  X bent >= RA_SIGMA_OOP_DEG out of plane and r <= RA_SIGMA_MAX_A
                                                    (loose, bent sigma-type radical anion)
    dissociated r > CX_DISSOCIATED_A (planar) or r > RA_SIGMA_MAX_A  (Ar...X- complex)
    elongated   anything else (planar, pi_max < r <= dissociation threshold)
    Returns (state, site, r, oop_deg).
    """
    X = species.halogen
    best = None
    for s, ci, xi, ring in species.cx_sites(X):
        if geom[xi - 1][0] == X:
            r = distance(geom, ci, xi)
            if best is None or r > best[2]:
                best = (s, ci, r, oop_angle(geom, ci, xi, ring))
    if best is None:
        return "pi", None, None, None
    s, ci, r, oop = best
    if r <= C.RA_PI_MAX_A[X]:
        state = "pi"
    elif oop >= C.RA_SIGMA_OOP_DEG and r <= C.RA_SIGMA_MAX_A[X]:
        state = "sigma_bent"
    elif r > C.CX_DISSOCIATED_A[X] or r > C.RA_SIGMA_MAX_A[X]:
        state = "dissociated"
    else:
        state = "elongated"
    return state, s, r, oop


def usable(status):
    return status in ("ok", "warn")
