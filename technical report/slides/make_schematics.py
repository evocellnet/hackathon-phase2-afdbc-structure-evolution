#!/usr/bin/env python
"""Regenerate the explanatory schematics used in the slide deck.

Usage:  python make_schematics.py            (all figures)
        python make_schematics.py s_3di      (selected figures)

Output: technical report/slides/assets/<name>.png  (1600x1000 px, 200 dpi)
All panels are illustrative cartoons, not data.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Arc, Circle

OUT = Path(__file__).resolve().parent / "assets"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- palette
DEEP = "#0B3C49"
TEAL = "#1B998B"
LIGHT = "#E8F1F2"
AMBER = "#F4A259"
CORAL = "#E76F51"
GREY = "#8D99AE"
INK = "#1F2933"
MID = "#C9DDE0"  # outline for light-teal shapes

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15,
    "axes.edgecolor": GREY,
    "axes.labelcolor": INK,
    "xtick.color": GREY,
    "ytick.color": GREY,
    "axes.linewidth": 1.2,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})

W, H, DPI = 8.0, 5.0, 200  # 1600 x 1000 px
REGISTRY = {}


def figure(name):
    def deco(fn):
        REGISTRY[name] = fn
        return fn
    return deco


def canvas(xmax=16, ymax=10):
    """Full-bleed axes with data units = 1/10 of figure height."""
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=DPI, transparent=False)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def box(ax, x, y, w, h, text="", fc=LIGHT, ec=None, tc=INK, fs=15, bold=False,
        r=0.25, lw=1.5, ls="-", z=2, **kw):
    """Rounded box with (x, y) = lower-left corner; text centred."""
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec or fc, lw=lw, ls=ls, zorder=z)
    ax.add_patch(p)
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc,
                fontsize=fs, fontweight="bold" if bold else "normal", zorder=z + 1,
                linespacing=1.2, **kw)
    return p


def arrow(ax, p0, p1, color=GREY, lw=2.2, rad=0.0, ms=18, z=1, ls="-", style="-|>"):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                        connectionstyle=f"arc3,rad={rad}", zorder=z, linestyle=ls,
                        shrinkA=0, shrinkB=0)
    ax.add_patch(a)
    return a


def clean_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def star(ax, x, y, s=520, color=AMBER, z=6):
    ax.scatter([x], [y], marker="*", s=s, color=color, edgecolor="white", lw=1.2, zorder=z)


def lock_icon(ax, x, y, s=1.0, color=AMBER):
    """Padlock centred at (x, y)."""
    ax.add_patch(Arc((x, y + 0.18 * s), 0.55 * s, 0.75 * s, theta1=0, theta2=180,
                     lw=4.5 * s, color=color, zorder=6))
    ax.plot([x - 0.275 * s, x - 0.275 * s], [y + 0.18 * s, y], color=color, lw=4.5 * s, zorder=6)
    ax.plot([x + 0.275 * s, x + 0.275 * s], [y + 0.18 * s, y], color=color, lw=4.5 * s, zorder=6)
    ax.add_patch(FancyBboxPatch((x - 0.42 * s, y - 0.5 * s), 0.84 * s, 0.62 * s,
                                boxstyle="round,pad=0,rounding_size=0.08", fc=color, ec=color, zorder=7))
    ax.add_patch(Circle((x, y - 0.17 * s), 0.08 * s, color="white", zorder=8))


def chain(seed=0, n=200, x0=0.0, x1=1.0, amp=1.0):
    """Smooth wiggly 2D backbone."""
    t = np.linspace(0, 1, n)
    rng = np.random.default_rng(seed)
    y = np.zeros_like(t)
    for k, a in zip((2, 3.3, 5.1, 7.4), (0.45, 0.3, 0.18, 0.1)):
        y += a * np.sin(2 * np.pi * k * t + rng.uniform(0, 2 * np.pi))
    x = x0 + (x1 - x0) * t + 0.025 * (x1 - x0) * np.sin(2 * np.pi * 4.2 * t)
    return x, amp * y


# =================================================================== 1
@figure("s_orthology")
def s_orthology():
    fig, ax = canvas()
    R = np.array([1.2, 6.0])
    V = np.array([5.6, 3.7])
    tips = {"fly": np.array([10.6, 8.6]), "mouse": np.array([10.6, 5.1]),
            "human": np.array([10.6, 1.9])}
    tube = [(R, tips["fly"]), (R, V), (V, tips["mouse"]), (V, tips["human"])]
    for lw, c in ((64, MID), (58, LIGHT)):
        for a, b in tube:
            ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=lw, solid_capstyle="round", zorder=1)
    for sp, (lx, ly) in {"fly": (8.0, 9.35), "mouse": (8.4, 6.25), "human": (7.6, 1.5)}.items():
        ax.text(lx, ly, sp, ha="center", va="center", fontsize=17, color=GREY,
                fontweight="bold", zorder=3)

    def seg(a, b, color, lw=3.2, z=4):
        ax.plot([a[0], b[0]], [a[1], b[1]], color="white", lw=lw + 4, zorder=z - 0.5,
                solid_capstyle="round")
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, zorder=z, solid_capstyle="round")

    D = R + 0.45 * (V - R)
    off = np.array([0, 0.3])
    leaves = {
        "fly": tips["fly"], "mouseA": tips["mouse"] + off, "mouseB": tips["mouse"] - off,
        "humanA": tips["human"] + off, "humanB": tips["human"] - off,
    }
    seg(R, tips["fly"], DEEP)
    seg(R, D, DEEP)
    for tag, col, o in (("B", DEEP, -off), ("A", TEAL, off)):
        Vo = V + o
        seg(D, Vo, col)
        seg(Vo, leaves["human" + tag], col)
        seg(Vo, leaves["mouse" + tag], col)
    ax.scatter([R[0]], [R[1]], s=150, color=DEEP, zorder=6)
    ax.text(0.25, 5.0, "ancestral\ngene", ha="left", va="top", fontsize=15, color=DEEP,
            linespacing=1.1)
    star(ax, *D, s=700)
    ax.plot([D[0], 3.0], [D[1] - 0.25, 2.75], color=AMBER, lw=1.5, zorder=3)
    ax.text(3.0, 2.6, "duplication", ha="center", va="top", fontsize=15, color=INK)

    for k, p in leaves.items():
        col = TEAL if k.endswith("A") else DEEP
        ax.scatter([p[0]], [p[1]], s=110, color=col, zorder=7, edgecolor="white", lw=1.2)

    # HOG box
    box(ax, 10.15, 1.15, 0.9, 8.0, fc="none", ec=AMBER, lw=2.5, ls="--", r=0.3, z=8)
    ax.text(11.35, 8.85, "HOG", ha="left", va="center", fontsize=18, color="#C77B30", fontweight="bold")
    ax.text(11.35, 8.3, "all descendants of\none ancestral gene", ha="left", va="top", fontsize=14,
            color=INK, linespacing=1.15)

    # ortholog / paralog arcs
    arrow(ax, leaves["mouseA"] + [0.2, 0], leaves["humanA"] + [0.2, 0], color=TEAL, lw=2.5,
          rad=-0.55, style="<|-|>", ms=14, z=9)
    ax.text(12.1, 4.15, "orthologs", ha="left", va="center", fontsize=16, color=TEAL, fontweight="bold")
    ax.text(12.1, 3.35, "split by\nspeciation", ha="left", va="center", fontsize=14, color=INK)
    arrow(ax, leaves["humanA"] + [0.2, 0], leaves["humanB"] + [0.2, 0], color=DEEP, lw=2.5,
          rad=-1.1, style="<|-|>", ms=12, z=9)
    ax.text(11.75, 1.65, "paralogs", ha="left", va="center", fontsize=16, color=DEEP, fontweight="bold")
    ax.text(11.75, 1.1, "split by duplication", ha="left", va="center", fontsize=14, color=INK)
    save(fig, "s_orthology")


# =================================================================== 2
@figure("s_3di")
def s_3di():
    fig, ax = canvas()
    aa = list("MKVLAGEF")
    tdi = list("DPVLVSCQ")
    t = np.linspace(0, 1, 300)
    xs = 0.9 + 6.0 * t
    ys = 5.0 + 1.1 * np.sin(2 * np.pi * 1.6 * t + 0.4)
    ax.plot(xs, ys, color=GREY, lw=6, solid_capstyle="round", zorder=1)
    idx = np.linspace(8, 291, 8).astype(int)
    px, py = xs[idx], ys[idx]
    # one 3D contact
    ax.plot([px[2], px[6]], [py[2], py[6]], ls=(0, (2, 2)), color=TEAL, lw=2, zorder=0.5)
    ax.text(3.6, 6.35, "3D neighbour", ha="center", va="center", fontsize=14, color=TEAL, zorder=2,
            bbox=dict(fc="white", ec="none", pad=1))
    for i in range(8):
        ax.add_patch(Circle((px[i], py[i]), 0.2, fc=DEEP, ec="white", lw=1.5, zorder=3))
        box(ax, px[i] - 0.3, 8.0, 0.6, 0.65, aa[i], fc=LIGHT, tc=DEEP, fs=16, bold=True, r=0.12)
        box(ax, px[i] - 0.3, 1.6, 0.6, 0.65, tdi[i], fc=TEAL, tc="white", fs=16, bold=True, r=0.12)
        ax.plot([px[i], px[i]], [py[i] + 0.25, 7.95], color=MID, lw=1.2, zorder=0)
        ax.plot([px[i], px[i]], [py[i] - 0.25, 2.3], color=MID, lw=1.2, zorder=0)
    ax.text(0.5, 9.25, "sequence", ha="left", fontsize=15, color=DEEP, fontweight="bold")
    ax.text(0.5, 0.9, "structure (3Di)", ha="left", fontsize=15, color=TEAL, fontweight="bold")

    ax.plot([8.1, 8.1], [0.8, 9.2], color=MID, lw=1.5)

    # right: two aligned proteins
    p1a, p2a = "MKVLAGEFRT", "MRVIAGDYRS"
    p1d, p2d = "DPVLVSCQAV", "DPVLVSCHAV"
    cw, x0 = 0.6, 9.65

    def row(y, s, ref, fc, tc):
        for j, ch in enumerate(s):
            diff = ch != ref[j]
            box(ax, x0 + j * cw, y, cw - 0.07, 0.62, ch, fc=AMBER if diff else fc,
                tc=INK if diff else tc, fs=15, bold=True, r=0.08,
                family="DejaVu Sans Mono")

    for y, lab in ((7.7, "P1"), (6.95, "P2"), (3.9, "P1"), (3.15, "P2")):
        ax.text(x0 - 0.2, y + 0.31, lab, ha="right", va="center", fontsize=14, color=GREY)
    row(7.7, p1a, p1a, LIGHT, DEEP)
    row(6.95, p2a, p1a, LIGHT, DEEP)
    row(3.9, p1d, p1d, TEAL, "white")
    row(3.15, p2d, p1d, TEAL, "white")
    ax.text(x0 + 5 * cw, 6.35, "20 amino acids", ha="center", va="top", fontsize=15, color=DEEP)
    ax.text(x0 + 5 * cw, 2.55, "20 structural letters (3Di)", ha="center", va="top", fontsize=15,
            color=TEAL)
    box(ax, x0 - 0.3, 0.75, 0.42, 0.42, fc=AMBER, r=0.06)
    ax.text(x0 + 0.3, 0.96, "= differs between P1 and P2", ha="left", va="center", fontsize=14, color=INK)
    save(fig, "s_3di")


# =================================================================== 3
@figure("s_saturation")
def s_saturation():
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0.11, 0.15, 0.85, 0.8])
    d = np.linspace(0, 3, 300)
    p = 0.95 * (1 - np.exp(-d / 0.95))
    ax.axvspan(2, 3, color=CORAL, alpha=0.12, lw=0)
    ax.plot(d, d, ls="--", color=GREY, lw=2.2)
    ax.plot(d, p, color=TEAL, lw=3.5)
    ax.axhline(0.95, color=GREY, lw=1.2, ls=":")
    ax.text(0.05, 0.99, "plateau", ha="left", va="bottom", fontsize=14, color=GREY)
    ax.text(0.83, 1.12, "true value", ha="right", va="center", fontsize=15, color=GREY, rotation=0)
    ax.text(2.95, 0.83, "observed\n(p-distance)", ha="right", va="top", fontsize=15, color=TEAL,
            fontweight="bold")
    x = 1.0
    arrow(ax, (x, 0.95 * (1 - np.exp(-x / 0.95)) + 0.02), (x, x - 0.03), color=AMBER, lw=3, ms=20, z=5)
    ax.text(1.2, 1.1, "ML distance corrects\nfor multiple hits", ha="left", va="center",
            fontsize=15, color=INK)
    ax.text(2.5, 0.25, "unreliable\nbeyond ~2", ha="center", va="center", fontsize=15, color=CORAL,
            fontweight="bold")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1.2)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=14, length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("true substitutions per site", fontsize=16)
    ax.set_ylabel("fraction of sites differing", fontsize=16)
    save(fig, "s_saturation")


# =================================================================== 4
@figure("s_structcompare")
def s_structcompare():
    fig, ax = canvas()
    x, y = chain(seed=3, n=400, x0=0.8, x1=8.6, amp=1.5)
    y = y + 5.0
    t = np.linspace(0, 1, 400)
    rng = np.random.default_rng(5)
    dev = 0.12 * np.sin(2 * np.pi * 3 * t) + 0.9 * np.clip(t - 0.82, 0, None) / 0.18
    x2 = x + 0.12 * np.cos(2 * np.pi * 2 * t) + 0.3 * np.clip(t - 0.82, 0, None) / 0.18
    y2 = y + dev
    ax.plot(x, y, color=TEAL, lw=5, solid_capstyle="round", zorder=3)
    ax.plot(x2, y2, color=AMBER, lw=5, solid_capstyle="round", zorder=2, alpha=0.95)
    idx = np.arange(20, 400, 38)
    ax.scatter(x[idx], y[idx], s=50, color=DEEP, zorder=4)
    ax.scatter(x2[idx], y2[idx], s=50, color="#B86E2C", zorder=4)
    # intra-chain distance lines (same pairs in both chains)
    for i, j in ((idx[1], idx[4]), (idx[3], idx[6]), (idx[5], idx[8])):
        ax.plot([x[i], x[j]], [y[i], y[j]], color=DEEP, lw=1.6, ls=(0, (3, 2)), zorder=1)
        ax.plot([x2[i], x2[j]], [y2[i], y2[j]], color="#B86E2C", lw=1.6, ls=(0, (3, 2)), zorder=1)
    tail = t > 0.84
    ax.plot((x[tail] + x2[tail]) / 2, (y[tail] + y2[tail]) / 2, color=CORAL, lw=34, alpha=0.15,
            solid_capstyle="round", zorder=0)
    ax.text(8.4, 3.75, "low pLDDT", ha="center", va="top", fontsize=14,
            color=CORAL, fontweight="bold")
    ax.plot([1.0, 1.6], [1.5, 1.5], color=TEAL, lw=5)
    ax.text(1.8, 1.5, "protein A", va="center", fontsize=14, color=INK)
    ax.plot([4.2, 4.8], [1.5, 1.5], color=AMBER, lw=5)
    ax.text(5.0, 1.5, "protein B", va="center", fontsize=14, color=INK)
    ax.plot([1.0, 1.6], [0.8, 0.8], color=DEEP, lw=1.6, ls=(0, (3, 2)))
    ax.text(1.8, 0.8, "local distances", va="center", fontsize=14, color=INK)

    items = [("TM-score", "global fold match (0–1)", TEAL),
             ("LDDT", "fraction of local\ndistances preserved (0–1)", DEEP),
             ("pLDDT", "AlphaFold confidence\nper residue", CORAL)]
    for k, (h, sub, c) in enumerate(items):
        yy = 8.6 - k * 2.9
        ax.add_patch(Rectangle((9.7, yy - 1.65), 0.14, 1.75, color=c, zorder=2))
        ax.text(10.1, yy, h, fontsize=19, fontweight="bold", color=c, va="top")
        ax.text(10.1, yy - 0.65, sub, fontsize=15, color=INK, va="top", linespacing=1.2)
    save(fig, "s_structcompare")


# =================================================================== 5
@figure("s_engines")
def s_engines():
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0.07, 0.15, 0.33, 0.72])
    rng = np.random.default_rng(1)
    xs = rng.uniform(0.05, 2.0, 140)
    eta, a = 0.08, 0.32
    ys = eta + 0.3 * xs ** 0.75 + rng.normal(0, 0.04, xs.size)
    ax.scatter(xs, ys, s=14, color=GREY, alpha=0.45, lw=0)
    xx = np.linspace(0, 2.1, 200)
    ax.plot(xx, eta + 0.2 * xx, color=TEAL, lw=3)
    ax.plot(xx, eta + 0.4 * xx ** 0.55, color=AMBER, lw=3)
    ax.text(2.15, eta + 0.2 * 2.1 - 0.02, "c = 1\nproportional", fontsize=14, color=TEAL, va="center",
            ha="left", fontweight="bold")
    ax.text(2.15, eta + 0.4 * 2.1 ** 0.55 + 0.05, "c < 1\nsaturating", fontsize=14, color="#C77B30",
            va="center", ha="left", fontweight="bold")
    ax.text(0.05, 0.72, "y = η + a·x$^c$", fontsize=18, color=INK)
    ax.set_xlim(0, 2.1)
    ax.set_ylim(0, 0.8)
    clean_axes(ax)
    ax.set_xlabel("sequence distance (x)", fontsize=15)
    ax.set_ylabel("structure distance (y)", fontsize=15)
    fig.text(0.03, 0.93, "(a) pairwise engine", fontsize=17, fontweight="bold", color=DEEP)

    # (b) branch engine
    fig.text(0.6, 0.93, "(b) branch engine", fontsize=17, fontweight="bold", color=DEEP)
    bx = fig.add_axes([0.6, 0.0, 0.4, 0.9])
    bx.set_xlim(0, 6.4)
    bx.set_ylim(0, 9)
    bx.axis("off")
    # tree: root -> n1 (A,B), n2 (C,D); horizontal branches as paired bars
    ys_leaf = {"A": 7.6, "B": 5.9, "C": 3.9, "D": 2.2}
    n1y, n2y = (7.6 + 5.9) / 2, (3.9 + 2.2) / 2
    rooty = (n1y + n2y) / 2
    edges = [  # (x0, y, length AA, ratio)
        (0.2, n1y, 1.4, 0.30), (0.2, n2y, 1.0, 0.22),
        (1.6, 7.6, 2.4, 0.18), (1.6, 5.9, 1.8, 0.25), (1.2, 3.9, 2.8, 0.20), (1.2, 2.2, 2.0, 0.28),
    ]
    bx.plot([0.2, 0.2], [n2y, n1y], color=GREY, lw=2)
    bx.plot([1.6, 1.6], [5.9, 7.6], color=GREY, lw=2)
    bx.plot([1.2, 1.2], [2.2, 3.9], color=GREY, lw=2)
    bx.plot([0.0, 0.2], [rooty, rooty], color=GREY, lw=2)
    for x0, y, L, r in edges:
        bx.add_patch(FancyBboxPatch((x0, y + 0.04), L, 0.26, boxstyle="round,pad=0,rounding_size=0.06",
                                    fc=TEAL, ec="none"))
        bx.add_patch(FancyBboxPatch((x0, y - 0.3), L * r * 2.0, 0.26,
                                    boxstyle="round,pad=0,rounding_size=0.06", fc=AMBER, ec="none"))
    for k, y in ys_leaf.items():
        end = max(x0 + L for x0, yy, L, r in edges if yy == y)
        bx.text(end + 0.15, y, k, va="center", fontsize=15, color=INK)
    bx.add_patch(Rectangle((3.3, 8.25), 0.4, 0.26, color=TEAL))
    bx.text(3.8, 8.38, "AA", va="center", fontsize=14)
    bx.add_patch(Rectangle((4.6, 8.25), 0.4, 0.26, color=AMBER))
    bx.text(5.1, 8.38, "3Di", va="center", fontsize=14)
    bx.text(0.0, 1.0, "k = ΣAA / Σ3Di", va="center", fontsize=18, fontweight="bold", color=DEEP)
    save(fig, "s_engines")


# =================================================================== 6
@figure("s_clock")
def s_clock():
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = fig.add_axes([0.09, 0.13, 0.62, 0.82])
    t = np.linspace(0, 1, 200)
    ax.plot(t, 1.0 * t, color=DEEP, lw=3.5)
    ax.plot(t, 0.3 * t, color=TEAL, lw=3.5)
    sat = 0.135 * (1 - np.exp(-t / 0.45))
    ax.plot(t, sat, color=CORAL, lw=3, ls="--")
    for tt in (0.4, 0.8):
        ax.plot([tt, tt], [0.3 * tt, tt], color=GREY, lw=1, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    clean_axes(ax)
    ax.set_xlabel("time since divergence", fontsize=16)
    ax.set_ylabel("divergence", fontsize=16)
    fig.text(0.73, 0.9, "sequence (AA)", fontsize=16, color=DEEP, fontweight="bold", va="center")
    fig.text(0.73, 0.37, "structure (3Di)", fontsize=16, color=TEAL, fontweight="bold", va="center")
    fig.text(0.73, 0.42, "constant ratio =\nrescaled clock", fontsize=14, color=INK, va="bottom")
    fig.text(0.73, 0.22, "lineage-specific /\nsaturating\n(violation)", fontsize=14, color=CORAL,
             va="center", fontweight="bold")
    save(fig, "s_clock")


# =================================================================== 7
@figure("s_prereg")
def s_prereg():
    fig, ax = canvas()
    ws = [2.3, 2.6, 2.75, 3.25, 2.9]
    gap = 0.5
    xs = [0.1]
    for w in ws[:-1]:
        xs.append(xs[-1] + w + gap)
    bh = 2.1
    yt, yb = 6.3, 1.7
    my = lambda y: y + bh / 2
    box(ax, xs[0], 4.15, ws[0], bh, "7,000\nfamilies", fc=DEEP, tc="white", fs=17, bold=True)
    box(ax, xs[1], yt, ws[1], bh, "1,400\ntune\n(20%)", fc=LIGHT, ec=TEAL, tc=DEEP, fs=15, bold=True)
    box(ax, xs[1], yb, ws[1], bh, "5,600\nconfirm\n(80%)", fc=LIGHT, ec=TEAL, tc=DEEP, fs=15, bold=True)
    box(ax, xs[2], yt, ws[2], bh, "tune\nparameters", fc=TEAL, tc="white", fs=15, bold=True)
    box(ax, xs[3], yt, ws[3], bh, "", fc="white", ec=AMBER, lw=2.5)
    lock_icon(ax, xs[3] + ws[3] / 2, yt + 1.5, s=0.7)
    ax.text(xs[3] + ws[3] / 2, yt + 0.55, "freeze plan\n(md5 in ledger)", ha="center", va="center",
            fontsize=14, color=INK, linespacing=1.1)
    box(ax, xs[3], yb, ws[3], bh, "run once\non confirm", fc=TEAL, tc="white", fs=15, bold=True)
    box(ax, xs[4], yb, ws[4], bh, "independent\nre-derivation\n42 / 42 ✓", fc=DEEP, tc="white",
        fs=14, bold=False)
    e = 0.06
    arrow(ax, (xs[0] + ws[0] + e, 5.4), (xs[1] - e, my(yt)))
    arrow(ax, (xs[0] + ws[0] + e, 4.9), (xs[1] - e, my(yb)))
    arrow(ax, (xs[1] + ws[1] + e, my(yt)), (xs[2] - e, my(yt)))
    arrow(ax, (xs[2] + ws[2] + e, my(yt)), (xs[3] - e, my(yt)))
    arrow(ax, (xs[3] + ws[3] / 2, yt - e), (xs[3] + ws[3] / 2, yb + bh + e))
    arrow(ax, (xs[1] + ws[1] + e, my(yb)), (xs[3] - e, my(yb)))
    arrow(ax, (xs[3] + ws[3] + e, my(yb)), (xs[4] - e, my(yb)))
    ax.text((xs[1] + ws[1] + xs[3]) / 2, my(yb) - 0.3, "untouched", ha="center", va="top",
            fontsize=14, color=GREY)
    save(fig, "s_prereg")


# =================================================================== 8
@figure("s_selection_bias")
def s_selection_bias():
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    rng = np.random.default_rng(7)
    x = rng.uniform(0.02, 2.2, 260)
    y = 0.05 + 0.3 * x + rng.normal(0, 0.06, x.size)
    thr = 0.42
    axes = [fig.add_axes([0.08, 0.17, 0.4, 0.66]), fig.add_axes([0.57, 0.17, 0.4, 0.66])]
    a, b = axes
    keep = (y < thr) & (x < 1.6 + rng.uniform(0, 0.6, x.size))
    a.scatter(x[~keep], y[~keep], s=14, facecolor="none", edgecolor=GREY, alpha=0.5, lw=0.8)
    a.scatter(x[keep], y[keep], s=16, color=TEAL, alpha=0.8, lw=0)
    a.axhline(thr, color=CORAL, lw=2.2, ls="--")
    a.text(0.03, thr + 0.03, "cluster filter\n(TM / coverage)", color=CORAL, fontsize=14,
           fontweight="bold", va="bottom", zorder=5,
           bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    xx = np.linspace(0, 2.2, 100)
    a.plot(xx, 0.05 + 0.34 * (1 - np.exp(-xx / 0.9)), color=CORAL, lw=3)
    a.text(2.18, 0.3, "apparent\nplateau", color=CORAL, fontsize=14, ha="right", va="top")
    a.text(1.0, 0.68, "removed", color=GREY, fontsize=14, ha="center")
    b.scatter(x, y, s=16, color=TEAL, alpha=0.8, lw=0)
    b.plot(xx, 0.05 + 0.3 * xx, color=DEEP, lw=3)
    b.text(0.05, 0.8, "full range kept", color=DEEP, fontsize=14, ha="left", va="top", fontweight="bold")
    for ax_ in axes:
        ax_.set_xlim(0, 2.25)
        ax_.set_ylim(-0.02, 0.85)
        clean_axes(ax_)
        ax_.set_xlabel("sequence distance", fontsize=15)
    a.set_ylabel("structure distance", fontsize=15)
    fig.text(0.08, 0.89, "in-cluster pairs", fontsize=16, fontweight="bold", color=INK)
    fig.text(0.57, 0.89, "ortholog pairs", fontsize=16, fontweight="bold", color=INK)
    save(fig, "s_selection_bias")


# =================================================================== 9
@figure("s_ldo")
def s_ldo():
    fig, ax = canvas()
    R = (0.7, 6.1)
    lw = 3.5

    def L(p, q, c=DEEP):
        ax.plot([p[0], p[0], q[0]], [p[1], q[1], q[1]], color=c, lw=lw, solid_capstyle="round",
                solid_joinstyle="round")

    n1 = (2.0, 5.2)
    D = (3.4, 4.2)
    L(R, (6.6, 9.0))
    ax.plot([R[0], R[0]], [6.1, 5.2], color=DEEP, lw=lw)
    ax.plot([R[0] - 0.4, R[0]], [6.1, 6.1], color=DEEP, lw=lw)
    L((R[0], 5.2), n1)
    L(n1, (6.6, 7.4))
    ax.plot([n1[0], n1[0]], [5.2, 4.2], color=DEEP, lw=lw)
    L(n1, D)
    L(D, (8.6, 5.3), TEAL)   # P1 long
    L(D, (5.4, 3.1), TEAL)   # P2 short
    ax.plot([D[0], D[0]], [4.2, 5.3], color=TEAL, lw=lw)
    ax.plot([D[0], D[0]], [4.2, 3.1], color=TEAL, lw=lw)
    star(ax, D[0], D[1], s=750)
    ax.text(3.1, 3.75, "duplication", ha="right", va="top", fontsize=14, color=INK)
    for (x, y, lab) in ((6.6, 9.0, "O1"), (6.6, 7.4, "O2")):
        ax.scatter([x], [y], s=120, color=DEEP, zorder=5)
        ax.text(x + 0.3, y, lab, va="center", fontsize=17, fontweight="bold", color=DEEP)
    ax.text(6.6 + 1.1, 8.2, "outgroup\northologs", va="center", fontsize=14, color=GREY)
    ax.scatter([8.6], [5.3], s=120, color=TEAL, zorder=5)
    ax.scatter([5.4], [3.1], s=120, color=TEAL, zorder=5)
    ax.text(8.9, 5.3, "P1 = MDO", va="center", fontsize=17, fontweight="bold", color=TEAL)
    ax.text(5.7, 3.1, "P2 = LDO", va="center", fontsize=17, fontweight="bold", color=TEAL)
    ax.text(8.9, 4.7, "longer sequence branch", va="center", fontsize=13, color=INK)
    ax.text(5.7, 2.5, "shorter branch", va="center", fontsize=13, color=INK)
    ax.text(11.7, 6.4, "?", fontsize=46, color=AMBER, fontweight="bold", ha="center", va="center")
    ax.text(12.3, 6.4, "extra structural\nchange?", fontsize=15, color=INK, va="center")

    # alignment split A / B
    y0 = 0.55
    box(ax, 0.6, y0, 3.2, 0.7, "half A", fc=LIGHT, ec=TEAL, tc=DEEP, fs=14, bold=True, r=0.1)
    box(ax, 3.9, y0, 3.2, 0.7, "half B", fc=AMBER, tc=INK, fs=14, bold=True, r=0.1)
    ax.text(0.6, y0 + 1.0, "alignment columns", fontsize=13, color=GREY, va="bottom")
    ax.text(7.5, 1.75, "A + O1: choose which copy is MDO", fontsize=14, color=INK, va="center")
    ax.text(7.5, 1.1, "B + O2: measure the difference", fontsize=14, color=INK, va="center")
    ax.text(7.5, 0.45, "(avoids regression to the mean)", fontsize=14, color=GREY, va="center")
    save(fig, "s_ldo")


# =================================================================== 10
@figure("s_cluster_vs_hog")
def s_cluster_vs_hog():
    fig, ax = canvas()
    lens = [4.6, 4.3, 4.8, 2.2, 2.4, 2.0]
    cols = [TEAL] * 3 + [AMBER] * 3
    order = [0, 3, 1, 4, 2, 5]  # interleave in the HOG
    box(ax, 0.4, 1.9, 5.8, 6.9, fc="white", ec=DEEP, lw=2.2, r=0.35)
    ax.text(3.3, 8.2, "ortholog family (HOG)", ha="center", fontsize=15, fontweight="bold",
            color=DEEP)
    for k, i in enumerate(order):
        box(ax, 0.9, 7.25 - k * 0.95, lens[i], 0.5, fc=cols[i], r=0.12)
    # clusters
    box(ax, 10.2, 5.9, 5.4, 3.1, fc="white", ec=TEAL, lw=2.2, r=0.35)
    box(ax, 10.2, 1.9, 5.4, 3.1, fc="white", ec="#C77B30", lw=2.2, r=0.35)
    ax.text(12.9, 8.45, "AFDB cluster 1", ha="center", fontsize=15, fontweight="bold", color=TEAL)
    ax.text(12.9, 4.45, "AFDB cluster 2", ha="center", fontsize=15, fontweight="bold", color="#C77B30")
    for k, i in enumerate([0, 1, 2]):
        box(ax, 10.5, 7.65 - k * 0.7, lens[i], 0.42, fc=TEAL, r=0.1)
    for k, i in enumerate([3, 4, 5]):
        box(ax, 10.5, 3.65 - k * 0.7, lens[i], 0.42, fc=AMBER, r=0.1)
    arrow(ax, (6.4, 6.2), (10.0, 7.3), color=TEAL, lw=3, ms=22)
    arrow(ax, (6.4, 4.4), (10.0, 3.3), color="#C77B30", lw=3, ms=22)
    ax.text(8.0, 7.15, "long", ha="center", va="center", fontsize=15, color=TEAL, fontweight="bold",
            rotation=17)
    ax.text(8.4, 3.45, "short", ha="center", va="center", fontsize=15, color="#C77B30",
            fontweight="bold", rotation=-17)
    # coverage illustration
    box(ax, 6.75, 5.25, 2.9, 0.32, fc=TEAL, r=0.08)
    box(ax, 6.75, 4.85, 1.4, 0.32, fc=AMBER, r=0.08)
    ax.text(8.3, 5.0, "≈50%", fontsize=14, color=CORAL, va="center", fontweight="bold")
    ax.text(3.3, 1.0, "Foldseek clustering needs ≥ 90% coverage of both proteins",
            ha="left", va="center", fontsize=15, color=INK)
    ax.text(3.05, 1.0, "✗", ha="right", va="center", fontsize=20, color=CORAL, fontweight="bold")
    save(fig, "s_cluster_vs_hog")


# =================================================================== 11
@figure("s_af_bias")
def s_af_bias():
    fig, ax = canvas()

    def pair(cx, cy, d, c1, c2, seed, w=3.6):
        x, y = chain(seed=seed, n=300, x0=cx - w / 2, x1=cx + w / 2, amp=1.15)
        t = np.linspace(0, 1, 300)
        x2 = x + d * 0.6 * np.sin(2 * np.pi * 2.5 * t + 1)
        y2 = y + d * np.sin(2 * np.pi * 1.7 * t + 0.3)
        ax.plot(x, y + cy, color=c1, lw=4.5, solid_capstyle="round", zorder=3)
        ax.plot(x2, y2 + cy, color=c2, lw=4.5, solid_capstyle="round", zorder=2)

    pair(2.6, 6.0, 0.45, TEAL, AMBER, seed=11)
    pair(8.0, 6.0, 0.16, TEAL, AMBER, seed=11)
    ax.text(2.6, 8.7, "experimental (PDB)", ha="center", fontsize=16, fontweight="bold", color=DEEP)
    ax.text(8.0, 8.7, "AlphaFold models", ha="center", fontsize=16, fontweight="bold", color=DEEP)
    ax.text(2.6, 3.65, "human vs mouse ortholog", ha="center", fontsize=13, color=GREY)
    ax.text(8.0, 3.65, "same two orthologs", ha="center", fontsize=13, color=GREY)
    arrow(ax, (4.3, 2.75), (6.3, 2.75), color=AMBER, lw=3.5, ms=24)
    ax.text(5.3, 2.15, "AF makes orthologs look\n30–40% more similar", ha="center", va="top",
            fontsize=15, color=INK, fontweight="bold")
    ax.plot([10.9, 10.9], [0.6, 9.3], color=MID, lw=1.5)
    # single protein AF vs crystal
    x, y = chain(seed=4, n=300, x0=11.7, x1=15.4, amp=1.15)
    t = np.linspace(0, 1, 300)
    ax.plot(x, y + 6.0, color=DEEP, lw=4.5, solid_capstyle="round", zorder=3)
    ax.plot(x + 0.1 * np.sin(2 * np.pi * 5 * t), y + 6.0 + 0.2 * np.sin(2 * np.pi * 3.1 * t + 2),
            color=CORAL, lw=4.5, solid_capstyle="round", zorder=2)
    ax.text(13.55, 8.7, "one protein", ha="center", fontsize=16, fontweight="bold", color=DEEP)
    ax.plot([11.5, 12.1], [3.6, 3.6], color=DEEP, lw=4.5)
    ax.text(12.3, 3.6, "crystal", va="center", fontsize=14)
    ax.plot([11.5, 12.1], [3.0, 3.0], color=CORAL, lw=4.5)
    ax.text(12.3, 3.0, "AlphaFold", va="center", fontsize=14)
    ax.text(13.55, 1.9, "per-model noise", ha="center", va="center", fontsize=16, color=CORAL,
            fontweight="bold")
    save(fig, "s_af_bias")


# =================================================================== 12
@figure("s_agents")
def s_agents():
    fig, ax = canvas()
    groups = [
        ("Scouts", ["literature", "data", "Ancestra", "environment"]),
        ("Builders", ["panel", "distance library", "pair pipeline", "statistics engines", "HOGs",
                      "floors", "annotations"]),
        ("Reviewers", ["design rounds 5→3→1", "results review", "report review", "retrospective"]),
        ("Verifiers", ["V11 re-derivation", "post-hoc refutation"]),
        ("Analysts", ["case studies", "cluster audit", "AF-vs-PDB"]),
        ("Writers", ["report"]),
    ]
    colw, gap, x0 = 4.9, 0.45, 0.3
    ls = 0.47
    hh = 0.62
    pi_w = 7.0
    box(ax, 8 - pi_w / 2, 8.95, pi_w, 0.9, "PI agent (Claude Opus 5.5)", fc=DEEP, tc="white", fs=16,
        bold=True)
    tops = [8.25, 8.25, 8.25, 3.35, 3.35, 3.35]
    bus1 = 8.6
    ax.plot([8, 8], [8.95, bus1], color=GREY, lw=2, zorder=0)
    centers = [x0 + i % 3 * (colw + gap) + colw / 2 for i in range(6)]
    ax.plot([centers[0], centers[2]], [bus1, bus1], color=GREY, lw=2, zorder=0)
    gx = [x0 + colw + gap / 2, x0 + 2 * colw + 1.5 * gap]
    bus2 = 3.7
    for g in gx:
        ax.plot([g, g], [bus1, bus2], color=GREY, lw=2, zorder=0)
    ax.plot([centers[3], centers[5]], [bus2, bus2], color=GREY, lw=2, zorder=0)
    ax.plot([gx[0], gx[0]], [bus1, bus1], color=GREY)
    ax.plot([centers[0], gx[0]], [bus1, bus1], color=GREY, lw=2, zorder=0)
    for i, (name, members) in enumerate(groups):
        cx = centers[i]
        top = tops[i]
        h = hh + len(members) * ls + 0.3
        ax.plot([cx, cx], [top, bus1 if i < 3 else bus2], color=GREY, lw=2, zorder=0)
        box(ax, cx - colw / 2, top - h, colw, h, fc=LIGHT, ec=MID, r=0.25, z=1)
        box(ax, cx - colw / 2, top - hh, colw, hh, name, fc=TEAL, tc="white", fs=16, bold=True,
            r=0.25, z=2)
        for k, m in enumerate(members):
            ax.text(cx, top - hh - 0.32 - (k + 0.5) * ls, m, ha="center", va="center",
                    fontsize=14, color=DEEP, zorder=3)
    save(fig, "s_agents")


# =================================================================== 13
@figure("s_workflow_decisions")
def s_workflow_decisions():
    fig, ax = canvas()
    steps = ["Scout", "Design", "Critique\n(3 rounds)", "Pilot on\nreal data", "Tune\n(20%)",
             "Freeze", "Confirm\nonce (80%)", "Re-derive", "Try to\nrefute", "Report"]
    notes = {3: "D22\nTM-align\nreplaced", 4: "D31/D32\nestimator\nbiases fixed",
             7: "missing\nsensitivities\nwired", 8: "per-model\nnoise found"}
    bw, bh, gap, x0 = 2.6, 1.45, 0.55, 0.25
    rows = [7.6, 3.2]
    pos = []
    for i, s in enumerate(steps):
        r, c = divmod(i, 5)
        x, y = x0 + c * (bw + gap), rows[r]
        pos.append((x, y))
        fc = DEEP if s in ("Freeze", "Report") else (TEAL if i >= 5 else LIGHT)
        tc = "white" if fc != LIGHT else DEEP
        box(ax, x, y, bw, bh, "" if s == "Freeze" else s, fc=fc, ec=fc if fc != LIGHT else MID,
            tc=tc, fs=14, bold=True)
        if s == "Freeze":
            lock_icon(ax, x + 0.6, y + bh / 2 + 0.05, s=0.5)
            ax.text(x + bw / 2 + 0.35, y + bh / 2, s, ha="center", va="center", fontsize=15,
                    fontweight="bold", color="white", zorder=4)
        if c < 4:
            arrow(ax, (x + bw + 0.04, y + bh / 2), (x + bw + gap - 0.04, y + bh / 2), ms=16)
    # row wrap connector
    xe, ye = pos[4]
    xs_, ys_ = pos[5]
    xr = xe + bw + 0.2
    ax.plot([xe + bw + 0.04, xr, xr], [ye + bh / 2, ye + bh / 2, 5.55], color=GREY, lw=2.2)
    ax.plot([xr, xs_ + bw / 2], [5.55, 5.55], color=GREY, lw=2.2)
    arrow(ax, (xs_ + bw / 2, 5.55), (xs_ + bw / 2, ys_ + bh + 0.04), ms=16)
    for i, txt in notes.items():
        x, y = pos[i]
        nb = 0.55 + 0.47 * txt.count("\n")
        ny = y - 0.25 - nb
        if i < 5:
            ny = y - 0.25 - nb
        box(ax, x + 0.05, ny, bw - 0.1, nb, txt, fc="#FDEBDB", ec=AMBER, tc=INK, fs=14, r=0.15,
            lw=1.8)
    save(fig, "s_workflow_decisions")


if __name__ == "__main__":
    names = sys.argv[1:] or list(REGISTRY)
    for n in names:
        REGISTRY[n]()
