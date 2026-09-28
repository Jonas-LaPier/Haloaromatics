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


def tilt_out_of_plane(atoms, c_idx, x_idx, angle_deg):
    """Rotate substituent x about carbon c out of the ring plane by angle_deg."""
    ring = [xyz(a) for a in atoms[:6]]
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


def displace_along(atoms, mode, amp=0.15):
    """Displace geometry along a normal mode (largest atom moves by amp A)."""
    m = max(norm(v) for v in mode) or 1.0
    return [(a[0],) + add(xyz(a), scale(v, amp / m)) for a, v in zip(atoms, mode)]


def fmt_atoms(atoms):
    return "\n".join(f"{s:<2} {x:14.8f} {y:14.8f} {z:14.8f}" for s, x, y, z in atoms)
