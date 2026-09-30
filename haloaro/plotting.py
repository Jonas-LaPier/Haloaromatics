"""Correlation figures in Jonas's plot style (plot-like-jonas; style file jonas.mplstyle).

White background, Cambria, full black box, outward ticks, no gridlines, fixed decimals per
axis, one y range for every panel, dotted black least-squares line over the fitted compounds
only, unframed legend in the corner that covers the least data (single figures only: in the 3 x 2
grid the caption defines the markers), compound shorthand labels
(data/experimental_kobs.csv `label`) with gray leader lines when pushed away from their
point, no chart titles (captions go in captions.md).

Series roles:
    bromobenzene   published data, fitted       filled gray circle, no outline
    chlorobenzene  new data, predicted          open circle, black outline
    pbde           published data, predicted    filled gray square, no outline
Error bars: +/- 1 standard error of kobs propagated to ln kobs (SE/kobs), where available.

Written by compile (needs matplotlib):
    plots/correlations_<level>.png/.pdf   six best models (fit on CORR_FIT_GROUP), 3 x 2 panels
    plots/qsar_<level>_<rank>_<descriptor>.png/.pdf   best PLOT_TOP_SINGLE models, labelled
    plots/captions.md                     figure captions
"""
from __future__ import annotations

import logging
import math
import re
import textwrap
from pathlib import Path

import config as C

STYLE = Path(__file__).with_name("jonas.mplstyle")
GRAY = "#a6a6a6"
ROLES = {
    "bromobenzene":  dict(marker="o", mfc=GRAY, mec=GRAY, label="published bromobenzene data"),
    "chlorobenzene": dict(marker="o", mfc="white", mec="black", label="new chlorobenzene data"),
    "pbde":          dict(marker="s", mfc=GRAY, mec=GRAY, label="published PBDE data"),
}
MARKS = {"bromobenzene": "filled gray circles are published bromobenzene data",
         "chlorobenzene": "open circles are new chlorobenzene data",
         "pbde": "gray squares are published PBDE data"}
FIRST_E = "−2.0 V vs SHE (uncompensated potential, without iR-drop compensation)"
E_TEXT = "−2.0 V vs SHE"
Y_LABEL = "ln k$_\\mathrm{obs}$ (h$^\\mathrm{−1}$)"


# --------------------------------------------------------------------------- #
# style
# --------------------------------------------------------------------------- #
def use_style():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    names = {f.name for f in font_manager.fontManager.ttflist}
    if "Cambria" not in names:
        for d in C.PLOT_FONT_DIRS:
            for f in sorted(Path(d).expanduser().glob("Cambria*")):
                try:
                    font_manager.fontManager.addfont(str(f))
                except Exception:
                    pass
        names = {f.name for f in font_manager.fontManager.ttflist}
    font = next((n for n in ("Cambria", "Caladea") if n in names), "DejaVu Serif")
    plt.style.use(str(STYLE))
    mpl.rcParams["font.serif"] = [font, "DejaVu Serif"]
    mpl.rcParams["mathtext.fontset"] = "custom"
    mpl.rcParams["mathtext.rm"] = font
    mpl.rcParams["mathtext.it"] = font + ":italic"
    mpl.rcParams["mathtext.bf"] = font + ":bold"
    mpl.rcParams["pdf.fonttype"] = 42          # embed TrueType so text stays editable
    return font


def nice_ticks(values, nbins=5, pad_frac=0.04):
    """Round, evenly spaced ticks framing the values with a little room, and the number of
    decimals that matches the step (one count per axis)."""
    from matplotlib.ticker import MaxNLocator
    lo, hi = min(values), max(values)
    pad = pad_frac * (hi - lo or abs(hi) or 1)
    ticks = list(MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo - pad, hi + pad))
    step = ticks[1] - ticks[0]
    ticks = [t for t in ticks if lo - pad - step < t < hi + pad + step]
    while ticks[0] > lo - pad:
        ticks.insert(0, ticks[0] - step)
    while ticks[-1] < hi + pad:
        ticks.append(ticks[-1] + step)
    while len(ticks) > 2 and ticks[1] <= lo - pad:
        ticks.pop(0)
    while len(ticks) > 2 and ticks[-2] >= hi + pad:
        ticks.pop()
    dec = max(0, -math.floor(math.log10(step) + 1e-9))
    if round(step * 10 ** dec, 6) % 1:                  # e.g. 2.5, 0.25
        dec += 1
    return ticks, dec


def _set_ticks(ax, which, ticks, dec):
    from matplotlib.ticker import FixedLocator, FormatStrFormatter
    axis = ax.xaxis if which == "x" else ax.yaxis
    (ax.set_xlim if which == "x" else ax.set_ylim)(ticks[0], ticks[-1])
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FormatStrFormatter(f"%.{dec}f"))


# --------------------------------------------------------------------------- #
# text
# --------------------------------------------------------------------------- #
# ASCII chemistry in config labels -> (mathtext for axes, plain Unicode for captions)
_CHEM = [("ArX.-", "ArX$^{•−}$", "ArX•⁻"), ("Ar.", "Ar$^{•}$", "Ar•"), ("Ar-", "Ar$^{−}$", "Ar⁻"),
         ("X-", "X$^{−}$", "X⁻"), ("e-", "e$^{−}$", "e⁻"), ("H+", "H$^{+}$", "H⁺"), ("->", "→", "→"),
         ("E deg", "E°", "E°"), ("dG0", "ΔG$_\\mathrm{0}$", "ΔG₀"), ("dG ", "ΔG ", "ΔG "),
         ("C-X", "C–X", "C–X"), ("Saveant", "Savéant", "Savéant"), ("dSCF", "ΔSCF", "ΔSCF"),
         (" RA ", " radical anion ", " radical anion ")]


def chem(text, plain=False):
    """ArX.- -> ArX•⁻, e- -> e⁻, dG -> ΔG ... (mathtext, or Unicode with plain=True)."""
    text = f" {text} "
    for a, math_, uni in _CHEM:
        text = text.replace(a, uni if plain else math_)
    return text.strip()


def unit(desc):
    d = desc.split(" @")[0]
    for suf, u in (("_A3", "Å$^{3}$"), ("_eV", "eV"), ("_kcal", "kcal/mol"), ("_V", "V vs SHE"),
                   ("_A", "Å"), ("_D", "D")):
        if d.endswith(suf):
            return u
    return ""


def axis_label(desc, label, width=38):
    at = re.search(r"@ ([+-]?\d+\.\d+) V", desc)
    text = label + (f" at {float(at.group(1)):.1f} V" if at else "")
    u = unit(desc)
    text = text + (f" ({u})" if u else "")
    return "\n".join(chem(line) for line in (textwrap.wrap(text, width) or [text]))


def level_text(level):
    lv = C.LEVELS[level]
    m = {"M062X": "M06-2X"}.get(lv["method"], lv["method"])
    return f"{m}/{lv['basis']}, " + ("gas phase" if not lv["solv"] else f"SMD ({C.SOLVENT.lower()})")


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #
def _points(data_rows, level, desc):
    pts = []
    for d in data_rows:
        if d["level"] != level or d.get(desc) in (None, ""):
            continue
        k, se = float(d["kobs_per_h"]), d.get("se_per_h")
        pts.append({"x": float(d[desc]), "y": float(d["ln_kobs"]), "group": d["group"],
                    "label": d.get("label") or d["compound"],
                    "yerr": float(se) / k if se not in (None, "") else None})
    return pts


def _series(ax, pts, ms):
    """Markers by role; returns [(group, style)] of the series drawn."""
    drawn = []
    for grp, st in ROLES.items():
        sel = [p for p in pts if p["group"] == grp]
        if not sel:
            continue
        kw = {k: v for k, v in st.items() if k != "label"}
        xs, ys, es = [p["x"] for p in sel], [p["y"] for p in sel], [p["yerr"] or 0 for p in sel]
        if any(es):
            ax.errorbar(xs, ys, yerr=es, fmt="none", ecolor="black", elinewidth=0.75,
                        capthick=0.75, capsize=2.5, zorder=2.5)
        ax.plot(xs, ys, linestyle="none", ms=ms, mew=0.8, zorder=3, **kw)
        drawn.append((grp, st))
    return drawn


def _fit_line(ax, fit, pts):
    xs = [p["x"] for p in pts if p["group"] == C.CORR_FIT_GROUP]
    if len(xs) < 2 or fit.get("slope") in (None, ""):
        return None
    m, b = float(fit["slope"]), float(fit["intercept"])
    x0, x1 = min(xs), max(xs)
    return ax.plot([x0, x1], [m * x0 + b, m * x1 + b], linestyle=(0, (1, 1.2)),
                   color="black", linewidth=1.5, zorder=2)[0]


def _line_pts(ax, line):
    if line is None:
        return []
    (xa, xb), (ya, yb) = line.get_xdata(), line.get_ydata()
    return [ax.transData.transform((xa + (xb - xa) * t / 150, ya + (yb - ya) * t / 150)) for t in range(151)]


def _legend(fig, ax, drawn, pts, line, fontsize, avoid=(), place=True):
    """Unframed, marker-only legend in the position covering the fewest markers (markers just
    outside it count too, so their labels keep room), the least fit line and none of `avoid`.
    Returns (score, legend)."""
    from matplotlib.lines import Line2D
    proxies = [Line2D([], [], linestyle="none", ms=6, mew=0.8,
                      **{k: v for k, v in st.items() if k != "label"}) for _, st in drawn]
    labels = [st["label"] for _, st in drawn]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    marks = [ax.transData.transform((p["x"], p["y"])) for p in pts]
    lpts = _line_pts(ax, line)
    boxes = [a.get_window_extent(r) for a in avoid]
    room = 1.2 * fontsize * fig.dpi / 72
    best = None
    for loc in ("upper right", "upper left", "lower left", "lower right", "center left",
                "center right", "lower center", "upper center"):
        leg = ax.legend(proxies, labels, loc=loc, frameon=False, handletextpad=0.3, fontsize=fontsize)
        bb = leg.get_window_extent(r)
        score = (10 * sum(bb.padded(3).contains(*m) for m in marks)
                 + 2 * sum(bb.padded(room).contains(*m) for m in marks)
                 + sum(bb.contains(*q) for q in lpts) + 50 * sum(_overlap(bb, o) for o in boxes))
        leg.remove()
        if best is None or score < best[0]:
            best = (score, loc)
    if not place:
        return best[0], None
    return best[0], ax.legend(proxies, labels, loc=best[1], frameon=False, handletextpad=0.3,
                              fontsize=fontsize)


def _overlap(a, b):
    return not (a.x1 < b.x0 or b.x1 < a.x0 or a.y1 < b.y0 or b.y1 < a.y0)


def _place_labels(fig, ax, pts, line, avoid=(), fontsize=9, ms=7, fixed=None, directions=None):
    """Greedy label placement, run after the layout is final: try positions around each
    point and keep the first that clears every marker, label, the fit line, the legend and
    the frame. Positions far from the point get a gray leader line. `fixed` maps a label to
    an (dx, dy) offset in points, placed first, for hand-tuned crowded spots. `directions`
    restricts the search to these (dx, dy) offsets (e.g. Excel's right/left/above/below).
    Returns {label: (dx, dy)} of the positions used."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(r)
    half = 0.5 * ms * fig.dpi / 72 + 1
    marks = [Bbox.from_bounds(*(ax.transData.transform((p["x"], p["y"])) - half), 2 * half, 2 * half)
             for p in pts]
    busy = list(marks) + [a.get_window_extent(r) for a in avoid if a is not None]
    lpts = _line_pts(ax, line)
    # candidate offsets (points): rings of increasing distance, 12 directions each
    cands = [(0.0, 6.0, 0.0), (0.0, -6.0, 0.0)]
    for d in (() if directions else (7, 12, 18, 26, 36)):
        for ang in (0, 180, 90, 270, 30, 150, 210, 330, 60, 120, 240, 300):
            dx, dy = d * math.cos(math.radians(ang)), 0.8 * d * math.sin(math.radians(ang))
            cands.append((d, round(dx, 1), round(dy, 1)))
    if directions:
        cands = [(math.hypot(dx, dy), dx, dy) for dx, dy in directions]
    placed = {}
    # crowded points first, so they get the positions closest to them
    crowd = [sum(abs(a.x0 - b.x0) < 40 and abs(a.y0 - b.y0) < 20 for b in marks) for a in marks]
    fixed = fixed or {}
    for i in sorted(range(len(pts)), key=lambda i: (pts[i]["label"] not in fixed, -crowd[i])):
        p = pts[i]
        best = None
        mine = [(math.hypot(*fixed[p["label"]]),) + tuple(fixed[p["label"]])] if p["label"] in fixed else cands
        for d, dx, dy in mine:
            ha = "left" if dx > 1 else ("right" if dx < -1 else "center")
            va = "center" if abs(dx) > 1 else ("bottom" if dy > 0 else "top")
            t = ax.annotate(p["label"], (p["x"], p["y"]), xytext=(dx, dy), textcoords="offset points",
                            fontsize=fontsize, ha=ha, va=va, zorder=4)
            bb = t.get_window_extent(r).expanded(1.03, 1.08)
            t.remove()
            inside = (bb.x0 >= ax_box.x0 + 2 and bb.x1 <= ax_box.x1 - 2 and
                      bb.y0 >= ax_box.y0 + 2 and bb.y1 <= ax_box.y1 - 2)
            hits = (sum(_overlap(bb, o) for k, o in enumerate(busy) if k != i)
                    + sum(bb.contains(*q) for q in lpts))
            # a label must sit closer to its own point than to any other point
            cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
            gap = lambda m: max(m.x0 - bb.x1, bb.x0 - m.x1, 0) ** 2 + max(m.y0 - bb.y1, bb.y0 - m.y1, 0) ** 2  # noqa: E731
            own = gap(marks[i])
            ambiguous = sum(gap(m) < own for k, m in enumerate(marks) if k != i)
            score = 1e5 * (not inside) + 1e3 * hits + 300 * ambiguous + d
            if best is None or score < best[0]:
                best = (score, d, dx, dy, ha, va, bb)
            if score < 8:
                break
        _, d, dx, dy, ha, va, bb = best
        ax.annotate(p["label"], (p["x"], p["y"]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=fontsize, ha=ha, va=va, zorder=4,
                    arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.75, shrinkA=0,
                                    shrinkB=half * 72 / fig.dpi + 1) if d > 12 else None)
        busy.append(bb)
        placed[p["label"]] = (dx, dy)
    return placed


def _panel(ax, fit, pts, desc, label, yticks, ms, fontsize=None):
    """Markers, fit line and axes; the legend and labels come after the layout is final."""
    drawn = _series(ax, pts, ms)
    line = _fit_line(ax, fit, pts)
    _set_ticks(ax, "x", *nice_ticks([p["x"] for p in pts]))
    _set_ticks(ax, "y", *yticks)
    ax.set_xlabel(axis_label(desc, label), fontsize=fontsize)
    ax.set_ylabel(Y_LABEL, fontsize=fontsize)
    return drawn, line


def _save(fig, stem):
    out = []
    for ext, dpi in (("png", 600), ("pdf", None)):
        fig.savefig(f"{stem}.{ext}", dpi=dpi, facecolor="white")
        out.append(f"{stem}.{ext}")
    return out


def _caption(num_text, level, models, first, labelled, present):
    names = {"bromobenzene": "bromobenzenes", "chlorobenzene": "chlorobenzenes", "pbde": "PBDEs"}
    fitted = names.get(C.CORR_FIT_GROUP, C.CORR_FIT_GROUP)
    legend = [MARKS[g] for g in ROLES if g in present]
    legend = (", ".join(legend[:-1]) + " and " + legend[-1]) if len(legend) > 1 else legend[0]
    fits = "; ".join((f"({chr(97 + i)}) {chem(m['model'], plain=True)}: " if len(models) > 1 else "")
                     + f"R² = {float(m['R2']):.2f}, n = {m['n']}" for i, m in enumerate(models))
    se_note = (" (no standard errors are available for the chlorobenzenes)"
               if "chlorobenzene" in present else "")
    predicted = " The other compounds are predicted by that fit." if present - {C.CORR_FIT_GROUP} else ""
    return (f"{num_text} Natural logarithm of the observed pseudo-first-order rate constants "
            f"(k_obs, h⁻¹) for electrochemical dehalogenation at {FIRST_E if first else E_TEXT} "
            f"[add electrode, cell, electrolyte and pH], plotted against descriptors computed at "
            f"{level_text(level)}. {legend[0].upper() + legend[1:]}. Error bars are ±1 standard "
            f"error of k_obs propagated to ln k_obs{se_note}. The dotted line is the least-squares "
            f"fit to the {fitted} ({fits}).{predicted}"
            + (" Compounds are labelled by shorthand." if labelled else ""))


def correlation_figures(qsar_rows, data_rows, outdir, top_grid=6, top_single=None):
    """Grid figure per level (best `top_grid` models fitted on CORR_FIT_GROUP) and labelled
    single figures for the best `top_single` (config.PLOT_TOP_SINGLE). Returns file list."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    top_single = C.PLOT_TOP_SINGLE if top_single is None else top_single
    font = use_style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for old in list(outdir.glob("correlations_*.*")) + list(outdir.glob("qsar_*.*")):
        old.unlink()
    # one y range for every panel, with a free band at the top for the legend
    yt, ydec = nice_ticks([float(d["ln_kobs"]) for d in data_rows], nbins=6)
    yticks = (yt + [yt[-1] + (yt[1] - yt[0])], ydec)

    def present(lv, models):
        return {d["group"] for r in models for d in data_rows
                if d["level"] == lv and d.get(r["descriptors"]) not in (None, "")}

    files, caps, fig_no, first = [], [f"# Correlation figures (font: {font})\n"], 1, True
    for lv in C.LEVELS:
        best = [r for r in qsar_rows if r["level"] == lv and r["fit_set"] == C.CORR_FIT_GROUP
                and r.get("R2") not in (None, "") and r.get("slope") not in (None, "")
                and int(r["n"] or 0) >= 5]
        best.sort(key=lambda r: -float(r["R2"]))
        if not best:
            continue
        grid = best[:top_grid]
        nrow = math.ceil(len(grid) / 3)
        fig, axes = plt.subplots(nrow, 3, figsize=(7.0, 2.45 * nrow), squeeze=False)
        panels = []
        for i, (ax, r) in enumerate(zip(axes.flat, grid)):
            pts = _points(data_rows, lv, r["descriptors"])
            drawn, line = _panel(ax, r, pts, r["descriptors"], r["model"], yticks, ms=5, fontsize=9)
            ax.tick_params(labelsize=9)
            ax.text(0.97, 0.97, f"({chr(97 + i)})", transform=ax.transAxes, ha="right", va="top", fontsize=9)
            panels.append((ax, pts, line))
        for ax in list(axes.flat)[len(grid):]:
            ax.axis("off")
        fig.tight_layout(w_pad=1.0, h_pad=1.2)
        # no legend in the grid: a legend is wider than a 2.3 in panel; the caption defines the markers
        stem = outdir / f"correlations_{lv}"
        files += _save(fig, stem)
        plt.close(fig)
        caps.append(f"## {stem.name}.png\n\n"
                    + _caption(f"Figure {fig_no}.", lv, grid, first, False, present(lv, grid)) + "\n")
        fig_no, first = fig_no + 1, False
        for k, r in enumerate(best[:top_single], 1):
            fig, ax = plt.subplots(figsize=(3.5, 3.0))
            pts = _points(data_rows, lv, r["descriptors"])
            drawn, line = _panel(ax, r, pts, r["descriptors"], r["model"], yticks, ms=7)
            fig.tight_layout()
            leg = None
            if len(drawn) >= 2:
                if _legend(fig, ax, drawn, pts, line, 9, place=False)[0] >= 10:
                    # the legend would cover a point: one more tick of headroom for this figure
                    t, dec = yticks
                    _set_ticks(ax, "y", t + [t[-1] + t[1] - t[0]], dec)
                leg = _legend(fig, ax, drawn, pts, line, 9)[1]
            _place_labels(fig, ax, pts, line, avoid=[leg])
            slug = re.sub(r"[^A-Za-z0-9]+", "_", r["descriptors"]).strip("_")
            stem = outdir / f"qsar_{lv}_{k}_{slug}"
            files += _save(fig, stem)
            plt.close(fig)
            caps.append(f"## {stem.name}.png\n\n"
                        + _caption(f"Figure {fig_no}.", lv, [r], False, True, present(lv, [r])) + "\n")
            fig_no += 1
    (outdir / "captions.md").write_text("\n".join(caps))
    files.append(str(outdir / "captions.md"))
    return files
