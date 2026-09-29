"""Small vector helpers (no numpy dependency so it runs anywhere on Sherlock)."""
import math


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    return scale(a, 1.0 / n) if n else a


def xyz(atom):
    return tuple(atom[1:4])


def distance(atoms, i, j):
    """1-based indices."""
    return norm(sub(xyz(atoms[i - 1]), xyz(atoms[j - 1])))


def _ring_pts(atoms, ring_idx):
    return [xyz(atoms[i - 1]) for i in (ring_idx or range(1, 7))]


def tilt_out_of_plane(atoms, c_idx, x_idx, angle_deg, ring_idx=None):
    """Rotate substituent x about carbon c out of the ring plane by angle_deg.
    ring_idx: 1-based indices of the six ring atoms (default atoms 1-6)."""
    ring = _ring_pts(atoms, ring_idx)
    n = unit(cross(sub(ring[2], ring[0]), sub(ring[4], ring[0])))
    c, x = xyz(atoms[c_idx - 1]), xyz(atoms[x_idx - 1])
    v = sub(x, c)
    d = norm(v)
    u = unit(sub(v, scale(n, dot(v, n))))   # in-plane direction
    t = math.radians(angle_deg)
    new = add(c, scale(add(scale(u, math.cos(t)), scale(n, math.sin(t))), d))
    out = list(atoms)
    out[x_idx - 1] = (atoms[x_idx - 1][0],) + new
    return out


def oop_angle(geom, ci, xi, ring_idx=None):
    """Angle (deg) between the C-X bond (1-based atom indices) and the ring plane."""
    ring = _ring_pts(geom, ring_idx)
    n = unit(cross(sub(ring[2], ring[0]), sub(ring[4], ring[0])))
    v = unit(sub(xyz(geom[xi - 1]), xyz(geom[ci - 1])))
    return math.degrees(math.asin(min(1.0, abs(dot(v, n)))))


def perturb(atoms, amp=0.03, seed=0):
    """Deterministic small displacement of every atom (breaks spatial symmetry)."""
    import random
    rnd = random.Random(seed)
    return [(a[0],) + tuple(c + rnd.uniform(-amp, amp) for c in a[1:4]) for a in atoms]


def displace_along(atoms, mode, amp=0.15):
    """Displace geometry along a normal mode (largest atom moves by amp A)."""
    m = max(norm(v) for v in mode) or 1.0
    return [(a[0],) + add(xyz(a), scale(v, amp / m)) for a, v in zip(atoms, mode)]


BONDI = {"H": 1.20, "C": 1.70, "O": 1.52, "Cl": 1.75, "Br": 1.85}


def vdw_volume(atoms, spacing=0.2):
    """Van der Waals volume (A^3) of the union of Bondi spheres, by grid counting."""
    pts = [(xyz(a), BONDI[a[0]]) for a in atoms]
    lo = [min(p[0][k] - p[1] for p in pts) for k in range(3)]
    hi = [max(p[0][k] + p[1] for p in pts) for k in range(3)]
    n = [int((hi[k] - lo[k]) / spacing) + 1 for k in range(3)]
    count = 0
    for i in range(n[0]):
        x = lo[0] + (i + 0.5) * spacing
        near_x = [(c, r * r) for c, r in pts if abs(c[0] - x) < r]
        if not near_x:
            continue
        for j in range(n[1]):
            y = lo[1] + (j + 0.5) * spacing
            near = [(c, r2) for c, r2 in near_x if (c[0] - x) ** 2 + (c[1] - y) ** 2 < r2]
            if not near:
                continue
            for k in range(n[2]):
                z = lo[2] + (k + 0.5) * spacing
                for c, r2 in near:
                    if (c[0] - x) ** 2 + (c[1] - y) ** 2 + (c[2] - z) ** 2 < r2:
                        count += 1
                        break
    return count * spacing ** 3


def sphere_radius(atoms):
    """Radius (A) of the sphere with the same volume as the vdW envelope."""
    return (3 * vdw_volume(atoms) / (4 * math.pi)) ** (1 / 3)


def fmt_atoms(atoms):
    return "\n".join(f"{s:<2} {x:14.8f} {y:14.8f} {z:14.8f}" for s, x, y, z in atoms)
