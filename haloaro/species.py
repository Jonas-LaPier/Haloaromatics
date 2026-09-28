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
Transition state TS_ClBz_12_x1                (C-X bond at parent locant 1 breaking
                                                in the radical anion)
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
    kind: str            # parent | radical_anion | aryl_radical | aryl_anion | halide | ts
    charge: int
    mult: int
    ring: tuple | None   # substituent tuple (None for atoms)
    halogen: str | None  # series (Cl/Br) or None for benzene / phenyl
    n_hal: int
    element: str | None = None  # for single atoms
    meta: dict = field(default_factory=dict)

    @property
    def is_atom(self):
        return self.ring is None

    def atoms(self):
        """List of (symbol, x, y, z) starting geometry.

        Atom order: C1..C6, then substituents in site order (vacant sites skipped).
        """
        if self.is_atom:
            return [(self.element, 0.0, 0.0, 0.0)]
        return ring_geometry(self.ring)

    def atom_index(self, site):
        """1-based Gaussian atom indices (C, substituent) for ring site 1..6."""
        c_idx = site
        n = 6
        for i, s in enumerate(self.ring):
            if s != "*":
                n += 1
            if i == site - 1:
                if s == "*":
                    return c_idx, None
                return c_idx, n
        raise ValueError(site)


@dataclass
class Dehalogenation:
    """ArX (+e-) -> ArX.- -> Ar. + X-  (one unique symmetry class of C-X bond)."""
    parent: str
    radical_anion: str
    site: int            # lowest parent locant of this class
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


def build_all(halogens=HALOGENS):
    """Return (species dict name->Species, list[Dehalogenation])."""
    sp: dict[str, Species] = {}
    rxns: list[Dehalogenation] = []

    def add(s: Species):
        sp.setdefault(s.name, s)

    add(Species("Bz", "parent", 0, 1, ("H",) * 6, None, 0))

    for x in halogens:
        add(Species(f"{x}_anion", "halide", -1, 1, None, x, 1, element=x))
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
                add(Species(ts_name, "ts", -1, 2, ring, x, n,
                            meta={"parent": pname, "site": site}))
                rxns.append(Dehalogenation(
                    parent=pname, radical_anion=f"{pname}_RA", site=site,
                    degeneracy=len(sites), halogen=x, aryl_radical=rad_name,
                    aryl_anion=an_name, hydro_product=parent_name(hring), ts=ts_name))
    return sp, rxns


if __name__ == "__main__":
    sp, rx = build_all()
    from collections import Counter
    print(Counter(s.kind for s in sp.values()))
    print(len(rx), "unique dehalogenation reactions")
    for r in rx:
        print(r)
