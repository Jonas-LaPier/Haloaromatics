"""Species enumeration, canonical naming and starting geometries.

A benzene ring is represented as a 6-tuple of substituents on C1..C6:
    'H', 'Cl', 'Br'  -> substituent atom
    '*'              -> vacant site (aryl radical / aryl carbanion carbon)

Canonical forms are taken over the 12 operations of the D6h ring permutation
group (6 rotations x 2 reflections), so only structurally unique isomers are
generated.

Naming
------
Parents          Bz, ClBz_12, BrBz_135, ...   (lowest halogen locant set)
Radical anions   ClBz_12_RA                   (ArX.-, charge -1, doublet)
Aryl radicals    Ph_rad, ClPh_2_rad, ...      (vacant carbon = C1; lowest locants)
Aryl carbanions  Ph_anion, ClPh_2_anion, ...  (charge -1, singlet)
Halides          Cl_anion, Br_anion
Halogen atoms    Cl_rad, Br_rad                (X., doublet; for C-X bond energies)
Transition state TS_ClBz_12_x1                (C-X bond at parent locant 1 breaking
                                                in the radical anion)

Polybrominated diphenyl ethers (skeleton "dpe")
----------------------------------------------
The two rings are A and B, each numbered from the ether carbon (C1). A substitution pattern
is a 10-tuple: ring A positions 2-6, then ring B positions 2-6. Canonical forms use the 8
symmetry operations (flip of either ring, 2<->6 and 3<->5, and exchange of the rings);
the more substituted ring is written first.
Parents          BDE_24_24 (BDE-47), BDE_245_24 (BDE-99), BDE_24_4 (BDE-28), ...
Aryl radicals    BDE_2r4_24_rad   (token "2r4" = radical carbon at 2, Br at 4 in that ring)
Sites            "2", "4", "5" (ring A) and "2p", "4p" (ring B, primed)
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

HALOGENS = ("Cl", "Br")

# Bond lengths (Angstrom) used only for the crude starting geometries.
BOND = {"H": 1.08, "Cl": 1.74, "Br": 1.90}
CC = 1.39


# --------------------------------------------------------------------------- #
# Ring symmetry helpers
# --------------------------------------------------------------------------- #
def _ring_perms():
    """All 12 index permutations p such that new[i] = old[p[i]]."""
    perms = []
    for r in range(6):
        perms.append(tuple((i + r) % 6 for i in range(6)))
        perms.append(tuple((r - i) % 6 for i in range(6)))
    return perms


RING_PERMS = _ring_perms()


def _apply(ring, p):
    return tuple(ring[p[i]] for i in range(6))


def locants(ring, atom):
    return tuple(i + 1 for i, s in enumerate(ring) if s == atom)


def canonical_parent(ring):
    """Orientation of a parent ring with the lowest halogen locant set."""
    best = None
    for p in RING_PERMS:
        r = _apply(ring, p)
        key = tuple(i + 1 for i, s in enumerate(r) if s != "H")
        if best is None or key < best[0]:
            best = (key, r)
    return best[1]


def canonical_aryl(ring):
    """Orientation with the vacant carbon at C1 and lowest halogen locants."""
    best = None
    for p in RING_PERMS:
        r = _apply(ring, p)
        if r[0] != "*":
            continue
        key = tuple(i + 1 for i, s in enumerate(r) if s not in ("H", "*"))
        if best is None or key < best[0]:
            best = (key, r)
    return best[1]


def _halogen_of(ring):
    hs = {s for s in ring if s in HALOGENS}
    if len(hs) > 1:
        raise ValueError("Mixed halogens are not supported by the naming scheme")
    return hs.pop() if hs else None


def parent_name(ring):
    ring = canonical_parent(ring)
    x = _halogen_of(ring)
    if x is None:
        return "Bz"
    return f"{x}Bz_{''.join(map(str, locants(ring, x)))}"


def aryl_name(ring, kind):
    """kind = 'rad' or 'anion'."""
    ring = canonical_aryl(ring)
    x = _halogen_of(ring)
    if x is None:
        return f"Ph_{kind}"
    return f"{x}Ph_{''.join(map(str, locants(ring, x)))}_{kind}"


# --------------------------------------------------------------------------- #
# Data classes
# --------------------------------------------------------------------------- #
@dataclass
class Species:
    name: str
    kind: str            # parent | radical_anion | aryl_radical | aryl_anion | halide | halogen_atom | ts
    charge: int
    mult: int
    ring: tuple | None   # substituent tuple (None for atoms)
    halogen: str | None  # series (Cl/Br) or None for benzene / phenyl
    n_hal: int
    element: str | None = None  # for single atoms
    meta: dict = field(default_factory=dict)
    skeleton: str = "benzene"   # benzene | dpe (diphenyl ether)

    @property
    def is_atom(self):
        return self.ring is None

    def atoms(self):
        """List of (symbol, x, y, z) starting geometry.

        benzene: C1..C6, then substituents in site order (vacant sites skipped).
        dpe:     ring A C1..C6, ring B C1'..C6', O, then substituents of A (2-6), B (2-6).
        """
        if self.is_atom:
            return [(self.element, 0.0, 0.0, 0.0)]
        if self.skeleton == "dpe":
            return dpe_geometry(self.ring)
        return ring_geometry(self.ring)

    # ---- site bookkeeping (site = int 1..6 for benzene, "2".."6"/"2p".."6p" for dpe)
    def _slots(self):
        """List of (site_label, ring_carbon_index, ring_atom_indices, symbol) in atom order."""
        if self.skeleton == "dpe":
            out = []
            for j, sym in enumerate(self.ring):
                ringB = j >= 5
                pos = j % 5 + 2
                out.append((f"{pos}p" if ringB else str(pos), (6 if ringB else 0) + pos,
                            list(range(7, 13)) if ringB else list(range(1, 7)), sym))
            return out, 13
        return [(i + 1, i + 1, list(range(1, 7)), sym) for i, sym in enumerate(self.ring)], 6

    def atom_index(self, site):
        """1-based Gaussian atom indices (ring C, substituent) for a site."""
        slots, n = self._slots()
        for label, c, _, sym in slots:
            if sym != "*":
                n += 1
            if str(label) == str(site):
                return c, (None if sym == "*" else n)
        raise ValueError(site)

    def ring_atoms(self, site):
        slots, _ = self._slots()
        return next(r for label, _, r, _ in slots if str(label) == str(site))

    def cx_sites(self, halogen=None):
        """Sites that carry a halogen: list of (site, c_idx, x_idx, ring_atom_indices)."""
        X = halogen or self.halogen
        out = []
        for label, c, r, sym in self._slots()[0]:
            if sym == X:
                out.append((label, *self.atom_index(label), r))
        return out

    def neighbours(self, site, halogen=None):
        """Halogens (ortho, meta, para) to a site within the same ring."""
        X = halogen or self.halogen
        if self.skeleton == "dpe":
            ringB = str(site).endswith("p")
            pos = int(str(site).rstrip("p"))
            sub = self.ring[5:] if ringB else self.ring[:5]
            at = lambda q: q in range(2, 7) and sub[q - 2] == X  # noqa: E731
            o = sum(at(pos + d) for d in (1, -1))
            m = sum(at(pos + d) for d in (2, -2))
            p = int(at(pos + 3) or at(pos - 3))
            return o, m, p
        i = int(site) - 1
        ring = self.ring
        return (sum(ring[(i + d) % 6] == X for d in (1, -1)),
                sum(ring[(i + d) % 6] == X for d in (2, -2)), int(ring[(i + 3) % 6] == X))


@dataclass
class Dehalogenation:
    """ArX (+e-) -> ArX.- -> Ar. + X-  (one unique symmetry class of C-X bond)."""
    parent: str
    radical_anion: str
    site: int | str      # lowest parent locant of this class (str for dpe sites)
    degeneracy: int      # number of equivalent C-X bonds
    halogen: str
    aryl_radical: str
    aryl_anion: str
    hydro_product: str   # ArH formed by 2e- + H+ hydrodehalogenation
    ts: str


def ring_geometry(ring):
    atoms = []
    # Regular hexagon in the xy plane; C1 at angle 90 deg, numbering clockwise
    for i in range(6):
        a = math.radians(90 - 60 * i)
        atoms.append(("C", CC * math.cos(a), CC * math.sin(a), 0.0))
    for i, s in enumerate(ring):
        if s == "*":
            continue
        a = math.radians(90 - 60 * i)
        r = CC + BOND[s]
        atoms.append((s, r * math.cos(a), r * math.sin(a), 0.0))
    return atoms


# --------------------------------------------------------------------------- #
# Diphenyl ethers
# --------------------------------------------------------------------------- #
CO = 1.39


def _dpe_ops():
    flip = lambda r: tuple(r[4 - k] for k in range(5))  # noqa: E731  positions 2..6 -> 6..2
    ops = []
    for swap in (False, True):
        for fa in (False, True):
            for fb in (False, True):
                def op(t, swap=swap, fa=fa, fb=fb):
                    a, b = t[:5], t[5:]
                    a = flip(a) if fa else a
                    b = flip(b) if fb else b
                    return (b + a) if swap else (a + b)
                ops.append(op)
    return ops


DPE_OPS = _dpe_ops()


def _ring_token(r, X):
    tok = ""
    for k, sym in enumerate(r):
        if sym == "*":
            tok += f"{k + 2}r"
        elif sym == X:
            tok += str(k + 2)
    return tok or "0"


def _dpe_key(t, X):
    a, b = t[:5], t[5:]
    na = sum(c != "H" for c in a)
    nb = sum(c != "H" for c in b)
    return (-na, _ring_token(a, X), -nb, _ring_token(b, X))


def canonical_dpe(t, X="Br"):
    return min((op(t) for op in DPE_OPS), key=lambda u: _dpe_key(u, X))


def dpe_name(t, X="Br", kind=None):
    c = canonical_dpe(t, X)
    base = f"BDE_{_ring_token(c[:5], X)}_{_ring_token(c[5:], X)}"
    return base if kind is None else f"{base}_{kind}"


def dpe_from_locants(ringA, ringB, X="Br"):
    t = ["H"] * 10
    for q in ringA:
        t[int(q) - 2] = X
    for q in ringB:
        t[5 + int(q) - 2] = X
    return tuple(t)


def dpe_geometry(t, twist_deg=50.0):
    """Twisted diphenyl ether starting geometry (C-O 1.39 A, C-O-C 120 deg)."""
    import math as m

    def add(a, b): return tuple(x + y for x, y in zip(a, b))
    def sc(a, k): return tuple(x * k for x in a)
    def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    def unit(a):
        n = m.sqrt(sum(x * x for x in a)); return sc(a, 1 / n)

    O = (0.0, 0.0, 0.0)
    rings, subs = [], []
    for side, sub in ((-1, t[:5]), (1, t[5:])):
        u = unit((side * m.sin(m.radians(60)), -m.cos(m.radians(60)), 0.0))
        c1 = sc(u, CO)
        center = add(c1, sc(u, CC))
        v0 = (0.0, 0.0, 1.0)
        phi = m.radians(twist_deg)
        v = add(sc(v0, m.cos(phi)), sc(cross(u, v0), m.sin(phi)))
        ring = []
        for k in range(6):
            th = m.radians(60 * k)
            ring.append(add(center, add(sc(u, -CC * m.cos(th)), sc(v, CC * m.sin(th)))))
        rings.append(ring)
        for k, sym in enumerate(sub):
            if sym == "*":
                continue
            pc = ring[k + 1]
            d = unit(tuple(p - q for p, q in zip(pc, center)))
            subs.append((sym, *add(pc, sc(d, BOND[sym]))))
    atoms = [("C", *p) for p in rings[0]] + [("C", *p) for p in rings[1]] + [("O", *O)]
    return atoms + subs


# known congener numbers (IUPAC/BZ numbering) for labelling
BDE_CONGENERS = {"BDE_24_24": "BDE-47", "BDE_245_24": "BDE-99", "BDE_24_4": "BDE-28",
                 "BDE_24_2": "BDE-17", "BDE_24_34": "BDE-66", "BDE_24_25": "BDE-49",
                 "BDE_245_2": "BDE-48", "BDE_245_4": "BDE-74"}


def add_pbdes(add, rxns, pbdes, include_ts=True, X="Br"):
    """Add PBDE parents (dict label -> (ringA locants, ringB locants)) and their reactions."""
    for label, (la, lb) in pbdes.items():
        t = canonical_dpe(dpe_from_locants(la, lb, X), X)
        pname = dpe_name(t, X)
        n = t.count(X)
        meta = {"congener": BDE_CONGENERS.get(pname, label)}
        add(Species(pname, "parent", 0, 1, t, X, n, meta=meta, skeleton="dpe"))
        add(Species(f"{pname}_RA", "radical_anion", -1, 2, t, X, n,
                    meta={"parent": pname, **meta}, skeleton="dpe"))
        par = Species(pname, "parent", 0, 1, t, X, n, skeleton="dpe")
        classes = {}
        for site, _, _, _ in par.cx_sites():
            j = (5 if str(site).endswith("p") else 0) + int(str(site).rstrip("p")) - 2
            vac = tuple("*" if k == j else t[k] for k in range(10))
            classes.setdefault(dpe_name(vac, X, "rad"), []).append((site, j))
        for rad_name, members in classes.items():
            site, j = members[0]
            vac = canonical_dpe(tuple("*" if k == j else t[k] for k in range(10)), X)
            add(Species(rad_name, "aryl_radical", 0, 2, vac, X, n - 1, skeleton="dpe"))
            an_name = dpe_name(vac, X, "anion")
            add(Species(an_name, "aryl_anion", -1, 1, vac, X, n - 1, skeleton="dpe"))
            h = canonical_dpe(tuple("H" if k == j else t[k] for k in range(10)), X)
            hname = dpe_name(h, X)
            add(Species(hname, "parent", 0, 1, h, X, n - 1, skeleton="dpe",
                        meta={"product_only": True, "congener": BDE_CONGENERS.get(hname, "")}))
            ts_name = f"TS_{pname}_x{site}"
            if include_ts:
                add(Species(ts_name, "ts", -1, 2, t, X, n, skeleton="dpe",
                            meta={"parent": pname, "site": site}))
            rxns.append(Dehalogenation(
                parent=pname, radical_anion=f"{pname}_RA", site=site, degeneracy=len(members),
                halogen=X, aryl_radical=rad_name, aryl_anion=an_name, hydro_product=hname,
                ts=ts_name))


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #
def unique_parents(halogen):
    """All structurally unique C6H(6-n)X(n), n = 1..6, in canonical form."""
    seen = {}
    for n in range(1, 7):
        for sites in itertools.combinations(range(6), n):
            ring = tuple(halogen if i in sites else "H" for i in range(6))
            c = canonical_parent(ring)
            seen.setdefault(parent_name(c), c)
    return seen


def build_all(halogens=HALOGENS, include_ts=True, pbdes=None):
    """Return (species dict name->Species, list[Dehalogenation])."""
    sp: dict[str, Species] = {}
    rxns: list[Dehalogenation] = []

    def add(s: Species):
        sp.setdefault(s.name, s)

    add(Species("Bz", "parent", 0, 1, ("H",) * 6, None, 0))

    for x in halogens:
        add(Species(f"{x}_anion", "halide", -1, 1, None, x, 1, element=x))
        add(Species(f"{x}_rad", "halogen_atom", 0, 2, None, x, 1, element=x))
        for pname, ring in unique_parents(x).items():
            n = ring.count(x)
            add(Species(pname, "parent", 0, 1, ring, x, n))
            add(Species(f"{pname}_RA", "radical_anion", -1, 2, ring, x, n,
                        meta={"parent": pname}))

            # group halogen sites of this parent into symmetry classes by product
            classes: dict[str, list[int]] = {}
            for i, s in enumerate(ring):
                if s != x:
                    continue
                vac = tuple("*" if j == i else ring[j] for j in range(6))
                classes.setdefault(aryl_name(vac, "rad"), []).append(i + 1)

            for rad_name, sites in classes.items():
                site = min(sites)
                vac = tuple("*" if j == site - 1 else ring[j] for j in range(6))
                aring = canonical_aryl(vac)
                hal = x if n > 1 else None
                an_name = aryl_name(vac, "anion")
                add(Species(rad_name, "aryl_radical", 0, 2, aring, hal, n - 1))
                add(Species(an_name, "aryl_anion", -1, 1, aring, hal, n - 1))
                hring = tuple("H" if j == site - 1 else ring[j] for j in range(6))
                ts_name = f"TS_{pname}_x{site}"
                if include_ts:
                    add(Species(ts_name, "ts", -1, 2, ring, x, n,
                                meta={"parent": pname, "site": site}))
                rxns.append(Dehalogenation(
                    parent=pname, radical_anion=f"{pname}_RA", site=site,
                    degeneracy=len(sites), halogen=x, aryl_radical=rad_name,
                    aryl_anion=an_name, hydro_product=parent_name(hring), ts=ts_name))
    if pbdes:
        add_pbdes(add, rxns, pbdes, include_ts)
    return sp, rxns


if __name__ == "__main__":
    sp, rx = build_all()
    from collections import Counter
    print(Counter(s.kind for s in sp.values()))
    print(len(rx), "unique dehalogenation reactions")
    for r in rx:
        print(r)
