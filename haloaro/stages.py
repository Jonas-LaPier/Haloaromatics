"""Stage definitions and Gaussian input (.gjf) generation.

Stages
------
am1                 AM1 pre-optimisation of every molecular species
<level>             Opt+Freq at each level in config.LEVELS (geometry read from AM1 chk)
tsscan_<level>      relaxed C-X scan on the radical anion (TS guess), per config.TS_LEVELS
ts_<level>          Opt=TS + Freq from the scan maximum
sp_<level>          single points on the optimised geometries of <level> (config.SP_LEVELS):
                      <P>__pop     parent at its own geometry, population analysis
                      <P>__vA      anion (-1, doublet) at the neutral geometry  (vertical EA, Fukui f+)
                      <P>__vC      cation (+1, doublet) at the neutral geometry (vertical IE)
                      <P>_RA__pop  radical anion at its own geometry, population analysis
                      <P>_RA__vN   neutral at the radical-anion geometry        (4-point lambda_i)
                    every job: Pop=(Hirshfeld,NBORead) with $NBO BNDIDX $END (CM5 charges,
                    Hirshfeld spin densities, Wiberg bond indices)
"""
from __future__ import annotations

import os
from dataclasses import dataclass, replace
from pathlib import Path

import config as C
from . import qc
from .geomtools import distance, displace_along, fmt_atoms, perturb, tilt_out_of_plane
from .parse import parse_log, scan_points
from .species import build_all

ROOT = Path(__file__).resolve().parent.parent
CALCS = "calcs"   # relative to project root; .gjf paths are relative to ROOT


@dataclass
class Stage:
    name: str
    kind: str          # am1 | optfreq | tsscan | ts | sp
    level: str | None

    @property
    def dir(self):
        return Path(CALCS) / self.name

    def path(self, sub, job, ext):
        return self.dir / sub / f"{job}.{ext}"

    def inp(self, job):
        return self.path("inputs", job, "gjf")

    def log(self, job):
        return self.path("logs", job, "log")

    def chk(self, job):
        return self.path("chks", job, "chk")


def all_stages():
    st = {"am1": Stage("am1", "am1", None)}
    for lvl in C.LEVELS:
        st[lvl] = Stage(lvl, "optfreq", lvl)
    for lvl in C.TS_LEVELS:
        st[f"tsscan_{lvl}"] = Stage(f"tsscan_{lvl}", "tsscan", lvl)
        st[f"ts_{lvl}"] = Stage(f"ts_{lvl}", "ts", lvl)
    for lvl in C.SP_LEVELS:
        st[f"sp_{lvl}"] = Stage(f"sp_{lvl}", "sp", lvl)
    return st


STAGES = all_stages()
GROUPS = {
    "optfreq": list(C.LEVELS),
    "tsscan": [f"tsscan_{l}" for l in C.TS_LEVELS],
    "ts": [f"ts_{l}" for l in C.TS_LEVELS],
    "sp": [f"sp_{l}" for l in C.SP_LEVELS],
}
GROUPS["all"] = list(STAGES)


def resolve(names):
    out = []
    for n in names:
        out += GROUPS.get(n, [n])
    for n in out:
        if n not in STAGES:
            raise SystemExit(f"Unknown stage '{n}'. Known: {', '.join(list(STAGES) + list(GROUPS))}")
    return out


SPECIES, REACTIONS = build_all(C.HALOGENS, include_ts=C.RUN_RA_TS, pbdes=C.PBDES)
RXN_BY_TS = {r.ts: r for r in REACTIONS}

NBO_TAIL = "$NBO BNDIDX $END"
SP_ROLES = {  # role: (source kind, charge, mult)
    "pop": (None, None, None),
    "vA": ("parent", -1, 2),
    "vC": ("parent", +1, 2),
    "vN": ("radical_anion", 0, 1),
}


def _sp_jobs():
    out = {}
    for n, s in SPECIES.items():
        if s.kind == "parent" and s.meta.get("product_only"):
            continue   # hydrodehalogenation products that are not themselves studied
        if s.kind == "parent":
            roles = ["pop", "vA", "vC"]
        elif s.kind == "radical_anion":
            roles = ["pop", "vN"]
        else:
            continue
        for role in roles:
            _, q, m = SP_ROLES[role]
            q = s.charge if q is None else q
            m = s.mult if m is None else m
            out[f"{n}__{role}"] = replace(s, name=f"{n}__{role}", kind="sp", charge=q, mult=m,
                                          meta={**s.meta, "src": n, "src_kind": s.kind, "role": role})
    return out


SP_JOBS = _sp_jobs()


def jobs(stage):
    """List of (job_name, Species, site)."""
    if stage.kind == "am1":
        return [(n, s, None) for n, s in SPECIES.items() if s.kind != "ts" and not s.is_atom]
    if stage.kind == "optfreq":
        return [(n, s, None) for n, s in SPECIES.items() if s.kind != "ts"]
    if stage.kind == "sp":
        return [(n, s, None) for n, s in SP_JOBS.items()]
    return [(r.ts, SPECIES[r.ts], r.site) for r in REACTIONS]


# --------------------------------------------------------------------------- #
# Input writing
# --------------------------------------------------------------------------- #
class Skip(Exception):
    pass


def resources(stage, job=None):
    """Slurm/Gaussian resources for a job; diphenyl ethers use config.RESOURCES_DPE."""
    sp = SPECIES.get(job) or SP_JOBS.get(job) if job else None
    if sp is not None and sp.skeleton == "dpe" and stage.kind in C.RESOURCES_DPE:
        return C.RESOURCES_DPE[stage.kind]
    return C.RESOURCES[stage.kind]


def _res(stage, job=None):
    return resources(stage, job)


def _link0(stage, job, oldchk=None):
    r = _res(stage, job)
    mem = max(1, int(r["mem_gb"] * C.MEM_FRACTION))
    s = f"%nprocshared={r['cpus']}\n%mem={mem}GB\n"
    if oldchk:
        s += f"%oldchk={oldchk}\n"
    s += f"%chk={stage.chk(job)}\n"
    return s


def _kw(name, opts):
    return name if not opts else f"{name}=({','.join(opts)})"


def _route(stage, species, opt_opts=(), extra=(), geom_check=False, guess_read=False, scf=None):
    parts = ["#p"]
    if stage.kind == "am1":
        parts += [_kw("Opt", list(opt_opts)), "AM1", scf or C.AM1_SCF]
    else:
        lvl = C.LEVELS[stage.level]
        if species.is_atom:
            parts.append("Freq")
        elif stage.kind == "optfreq":
            parts += [_kw("Opt", list(opt_opts)), "Freq"]
        elif stage.kind == "tsscan":
            parts.append(_kw("Opt", ["ModRedundant", "MaxCycles=100"] + list(opt_opts)))
        elif stage.kind == "ts":
            parts += [_kw("Opt", ["TS", "CalcFC", "NoEigenTest", "MaxCycles=150"] + list(opt_opts)), "Freq"]
        elif stage.kind == "sp":
            parts.append("SP")
        parts.append(f"{lvl['method']}/{lvl['basis']}")
        if lvl["solv"]:
            parts.append(lvl["solv"])
        if C.DFT_EXTRA:
            parts.append(C.DFT_EXTRA)
        if stage.kind in ("tsscan", "ts"):
            parts.append("NoSymm")
        if stage.kind == "sp":
            parts.append("Pop=(Hirshfeld,NBORead)")
            if species.meta.get("role") != "pop":
                # vertical states of symmetric rings put the extra/missing electron in a
                # degenerate orbital: let the SCF break symmetry, with a robust fallback
                parts += ["NoSymm", scf or "SCF=XQC"]
                scf = None
        if scf:
            parts.append(scf)
    if geom_check:
        parts.append("Geom=Check")
    if guess_read:
        parts.append("Guess=Read")
    parts += list(extra)
    return " ".join(dict.fromkeys(parts))   # drop duplicate keywords, keep order


def _gjf(link0, route, title, species, coords=None, tail=""):
    s = f"{link0}{route}\n\n{title}\n\n{species.charge} {species.mult}\n"
    if coords:
        s += fmt_atoms(coords) + "\n"
    s += "\n"
    if tail:
        s += tail + "\n\n"
    return s


def _load_ok(stage, job, species, site=None):
    lp = ROOT / stage.log(job)
    if not lp.exists():
        raise Skip(f"prerequisite {stage.name}/{job} has no log")
    p = parse_log(lp)
    status, flags = qc.evaluate(p, species, stage.kind, site)
    if not qc.usable(status):
        raise Skip(f"prerequisite {stage.name}/{job} status={status} {flags}")
    return p, flags


def build_input(stage, job, species, site, retry=None, ignore_prereq=False):
    """Return gjf text. `retry` is a dict of fixes produced by retry_fixes()."""
    retry = retry or {}
    title = f"{job} {stage.name}"
    opt_opts = list(retry.get("opt_opts", []))
    scf = retry.get("scf")
    coords = retry.get("coords")
    own_chk = retry.get("from_own_chk", False)
    extra = tuple(retry.get("extra", ()))

    def _rt(sp, oo=(), **kw):
        kw.setdefault("scf", scf)
        return _route(stage, sp, oo, extra=extra, **kw)

    if own_chk:  # restart from this job's own checkpoint (last geometry + wfn)
        link0 = _link0(stage, job)
        route = _rt(species, opt_opts, geom_check=True, guess_read=not retry.get("no_guess"))
        tail = retry.get("tail", NBO_TAIL if stage.kind == "sp" else "")
        return _gjf(link0, route, title, species, tail=tail)

    if stage.kind == "sp":
        src_stage = STAGES[stage.level]
        src = species.meta["src"]
        p, flags = _load_ok(src_stage, src, SPECIES[src])
        if species.meta["src_kind"] == "radical_anion" and any(
                f.startswith("RA_dissociated") for f in flags):
            raise Skip(f"{src} is not bound at {stage.level}: no radical-anion geometry")
        role = species.meta["role"]
        return _gjf(_link0(stage, job, oldchk=src_stage.chk(src)),
                    _rt(species, geom_check=True, guess_read=(role == "pop")),
                    title, species, tail=NBO_TAIL)

    if stage.kind == "am1":
        return _gjf(_link0(stage, job), _rt(species, opt_opts),
                    title, species, coords or species.atoms())

    if stage.kind == "optfreq":
        if species.is_atom:
            return _gjf(_link0(stage, job), _rt(species), title, species, species.atoms())
        if coords:
            return _gjf(_link0(stage, job), _rt(species, opt_opts), title, species, coords)
        am1 = STAGES["am1"]
        if species.kind == "radical_anion":
            # Symmetric starting rings converge to symmetry-constrained saddle points (planar
            # pi radical anions with an imaginary out-of-plane or Jahn-Teller mode). Start
            # from a slightly perturbed AM1 geometry and optimise without symmetry.
            p, _ = _load_ok(am1, job, species)
            return _gjf(_link0(stage, job), _route(stage, species, opt_opts, extra=("NoSymm",) + extra, scf=scf),
                        title, species, perturb(p["geometry"], 0.03, seed=len(job)))
        if not ignore_prereq:
            _load_ok(am1, job, species)
        return _gjf(_link0(stage, job, oldchk=am1.chk(job)),
                    _rt(species, opt_opts, geom_check=True), title, species)

    rxn = RXN_BY_TS[job]
    ci, xi = species.atom_index(site)

    if stage.kind == "tsscan":
        src = STAGES[stage.level]
        ra = SPECIES[rxn.radical_anion]
        p, flags = _load_ok(src, rxn.radical_anion, ra)
        if any(f.startswith("RA_dissociated") for f in flags):
            raise Skip(f"{rxn.radical_anion} is not bound at {stage.level} ({flags}); "
                       "C-X cleavage is barrierless from the RA at this level")
        geom = tilt_out_of_plane(p["geometry"], ci, xi, C.OOP_ANGLE_DEG, species.ring_atoms(site))
        tail = f"B {ci} {xi} S {C.SCAN_STEPS} {C.SCAN_STEP_SIZE:.3f}"
        return _gjf(_link0(stage, job), _rt(species, opt_opts),
                    title, species, geom, tail=tail)

    if stage.kind == "ts":
        if coords:
            return _gjf(_link0(stage, job), _rt(species, opt_opts), title, species, coords)
        scan = STAGES[f"tsscan_{stage.level}"]
        lp = ROOT / scan.log(job)
        if not lp.exists():
            raise Skip(f"no scan log {scan.log(job)}")
        text_ok = qc.evaluate(parse_log(lp), species, "tsscan", site)
        pts = [q for q in scan_points(lp.read_text(errors="ignore").splitlines()) if q["energy"] is not None]
        if len(pts) < 3:
            raise Skip(f"scan has only {len(pts)} points ({text_ok})")
        es = [q["energy"] for q in pts]
        k = max(range(len(es)), key=es.__getitem__)
        if k == 0:
            raise Skip("scan energy decreases monotonically: no barrier on the RA surface")
        if k == len(es) - 1:
            raise Skip("scan maximum at last point: extend SCAN_STEPS in config.py")
        g = pts[k]["geometry"]
        title += f" (scan point {k}/{len(es)-1}, r(C-X)={distance(g, ci, xi):.3f} A)"
        return _gjf(_link0(stage, job), _rt(species, opt_opts), title, species, g)

    raise ValueError(stage.kind)


# --------------------------------------------------------------------------- #
# Retry logic
# --------------------------------------------------------------------------- #
def _scf_failed_at_start(parsed):
    """True if the SCF never converged at the first geometry (no opt step taken)."""
    lines = parsed.get("_lines", [])
    return not any("Step number" in l for l in lines)


def retry_fixes(stage, parsed, status, flags, species=None):
    """Map a failure to input modifications."""
    et = parsed.get("error_type") if parsed else None
    fx = {}
    nosymm = [] if stage.kind in ("tsscan", "ts") else ["NoSymm"]
    if et == "scf_convergence":
        # level shift damps orbital flipping between near-degenerate pi* orbitals;
        # NoSymm stops the occupation being locked to one irreducible representation
        fx["scf"] = "SCF=(XQC,VShift=400,MaxCycle=512)"
        fx["extra"] = nosymm
        if (stage.kind == "am1" and species is not None and species.kind == "radical_anion"
                and _scf_failed_at_start(parsed)):
            # start from the converged neutral parent's AM1 geometry instead of the ideal ring
            pp = ROOT / stage.log(species.meta["parent"])
            if pp.exists():
                par = parse_log(pp)
                if qc.usable(qc.evaluate(par, SPECIES[species.meta["parent"]], "am1")[0]):
                    fx["coords"] = par["geometry"]
        elif stage.kind in ("optfreq", "ts") and not _scf_failed_at_start(parsed) and parsed.get("geometry"):
            # SCF failed mid-optimisation: keep the geometry reached, but not the bad wavefunction
            fx.update({"from_own_chk": True, "no_guess": True, "opt_opts": ["MaxCycles=300"]})
        return fx
    if et in ("opt_maxcycles", "opt_not_converged", "time_limit", None) and parsed and parsed.get("geometry"):
        fx["from_own_chk"] = True
        fx["opt_opts"] = ["MaxCycles=300", "CalcFC"] if stage.kind in ("am1", "optfreq") else []
        if stage.kind == "optfreq":
            fx["extra"] = nosymm
    if et == "internal_coords":
        fx["opt_opts"] = ["Cartesian", "MaxCycles=300"]
        fx["from_own_chk"] = True
    if status == "incomplete" and parsed and parsed.get("geometry"):
        fx["from_own_chk"] = True
    imag = [f for f in flags if f.startswith("n_imag=")]
    if imag and stage.kind == "optfreq" and parsed.get("mode1") and parsed.get("geometry"):
        # push the minimum off the saddle along the imaginary mode, without symmetry so the
        # optimiser cannot re-symmetrise the displaced structure
        big = abs(parsed.get("lowest_freq") or 0) > 60
        g = displace_along(parsed["geometry"], parsed["mode1"], 0.25 if big else 0.15)
        fx = {"coords": perturb(g, 0.01, seed=7),
              "opt_opts": ["CalcFC", "MaxCycles=300"] + ([] if big else ["Tight"]),
              "extra": nosymm}
    if stage.kind == "tsscan" and fx.get("from_own_chk"):
        fx.pop("from_own_chk")  # restart scans from scratch (modredundant is not kept)
    return fx
