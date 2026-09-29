"""Draw the figures of docs/ (docs/figs/*.svg) from a fresh small simulation.

  python scripts/make_figures.py            # needs bin/fsm2d (make) and matplotlib

futures.svg / score.svg  one bird of a small 2D flock (N = 120, tau = 4, 150 steps of bin/fsm2d): its futures, its
                         view, and the score of one of its moves (docs/NOTES.md sections 1-2)
dag.svg                  the tree of futures and the DAG that merges it (docs/OPTIMIZATIONS.md section 2.2)
sphere.svg               two neighbours painted on the 3D latitude-longitude grid (docs/NOTES.md section 4.3)
empowerment.svg          the decisions of the paper's score, the Gaussian bound and survival counting against
                         empowerment (docs/NOTES.md section 3.6; the numbers of scripts/empowerment_check.py)
"""
import math
import os
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
from fsmtools import load_run  # noqa: E402
from fsmtools.metrics import visual_intervals  # noqa: E402
from fsmtools.reference2d import futures, masks  # noqa: E402
from fsmtools.reference3d import DEFAULTS, Grid, paint  # noqa: E402

OUT = os.path.join(ROOT, "docs", "figs")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"  # categorical slots 1-3
INK, MUTED, FAINT, GRID = "#1f2328", "#59636e", "#b8bec6", "#e8eaed"
FOCAL = "#e34948"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                     "axes.edgecolor": FAINT, "xtick.color": MUTED, "ytick.color": MUTED, "svg.fonttype": "none",
                     "figure.facecolor": "white", "axes.facecolor": "white"})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("wrote docs/figs/" + name)


def flock():
    """A small flocked state from the CPU simulator, and the bird nearest its centre."""
    d = tempfile.mkdtemp()
    subprocess.run([os.path.join(ROOT, "bin", "fsm2d"), "--N", "120", "--tau", "4", "--box", "110", "--steps", "150", "--seed", "3",
                    "--log-every", "0", "--out", os.path.join(d, "r")], check=True, capture_output=True)
    r = load_run(os.path.join(d, "r"))
    frame = np.array(r["traj"][-1])
    c = frame[:, :2].mean(0)
    i = int(np.argmin(((frame[:, :2] - c) ** 2).sum(1)))
    return r["params"], frame, i


def fig_futures(p, frame, i, tau=4, move=1):
    """The bird's frame, heading up: other birds, its futures (those of `move` highlighted) and its current view."""
    th = frame[i, 2]
    rot = lambda x, y: (np.sin(th) * x - np.cos(th) * y, np.cos(th) * x + np.sin(th) * y)  # noqa: E731  heading up, left to the left
    levels = futures(frame, i, p, tau)
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    rel = frame[:, :2] - frame[i, :2]
    ox, oy = rot(rel[:, 0], rel[:, 1])
    for j in range(len(frame)):
        if j != i:
            ax.add_patch(plt.Circle((ox[j], oy[j]), 1.3 * p["R"], color="#3d444d", lw=0, zorder=5))
    for m0, h, px, py, _ in levels:
        fx, fy = -py, px  # futures are in the bird's frame (x ahead, y to the left): ahead is up
        ax.scatter(fx[m0 != move], fy[m0 != move], s=4, color="#c5cbd3", lw=0, zorder=2)
        ax.scatter(fx[m0 == move], fy[m0 == move], s=6, color=BLUE, lw=0, zorder=3)
    ax.add_patch(plt.Circle((0, 0), 1.6 * p["R"], color=FOCAL, lw=0, zorder=4))
    rr = 58
    ax.add_patch(plt.Circle((0, 0), rr, fill=False, color=GRID, lw=6, zorder=1))
    for a, b in visual_intervals(frame, i, p["R"], p["half_angle"]):
        t = np.linspace(a, b, max(2, int((b - a) * 60)))  # relative to heading; heading is up
        ax.plot(-rr * np.sin(t), rr * np.cos(t), color=INK, lw=6, solid_capstyle="butt", zorder=2)
    ax.annotate("", (0, 9), (0, 2.5), arrowprops=dict(arrowstyle="->", color=FOCAL, lw=1.5), zorder=6)
    ax.text(0, -rr - 7, "what it sees now (blocked directions in black)", ha="center", color=MUTED)
    ax.text(rr + 6, 30, "futures that start\nby turning left", color=BLUE, ha="left", va="center")
    ax.text(rr + 6, 10, "other futures", color=MUTED, ha="left", va="center")
    ax.text(rr + 6, -10, "other birds", color="#3d444d", ha="left", va="center")
    ax.text(rr + 6, -30, "the bird, heading up", color=FOCAL, ha="left", va="center")
    ax.set_xlim(-rr - 8, rr + 48)
    ax.set_ylim(-rr - 12, rr + 8)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "futures.svg")


def fig_score(p, frame, i, tau=4, move=1, B=3072):
    X, first = masks(futures(frame, i, p, tau), p, B)
    Xb = X[first == move]
    pb = Xb.mean(0)
    deg = (np.arange(B) + 0.5) * 360 / B
    deg = np.where(deg > 180, deg - 360, deg)
    o = np.argsort(deg)
    deg, pb = deg[o], pb[o]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True, gridspec_kw={"hspace": 0.18})
    a1.plot(deg, pb, color=BLUE, lw=1.2)
    a1.set_ylabel("$p_b(\\theta)$")
    a1.set_ylim(-0.03, 1.03)
    a2.fill_between(deg, pb * (1 - pb), color=BLUE, alpha=0.35, lw=0)
    a2.plot(deg, pb * (1 - pb), color=BLUE, lw=1.0)
    a2.set_ylabel("$p_b(1-p_b)$")
    a2.set_ylim(0, 0.27)
    a2.set_xlabel("direction θ relative to the heading (degrees; 0 = ahead, positive = left)")
    W = len(Xb) / math.pi * (pb * (1 - pb)).sum() * 2 * math.pi / B
    a2.text(0.01, 0.9, f"$n_b$ = {len(Xb)} futures,  $W_b$ = {W:.0f}", transform=a2.transAxes, color=INK, va="top")
    for a in (a1, a2):
        a.grid(axis="y", color=GRID, lw=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    a2.set_xlim(-180, 180)
    a2.set_xticks(range(-180, 181, 45))
    save(fig, "score.svg")


def fig_dag(tau=3, v0=10.0, dv=2.0, dth=math.radians(15)):
    moves = [(-1, 0), (1, 0), (0, -1), (0, 0), (0, 1)]
    tree = [[(0, 0.0, 0.0)]]  # (heading index, x, y) per path
    for k in range(tau):
        tree.append([(h + dh, x + (v0 + ds * dv) * math.cos((h + dh) * dth), y + (v0 + ds * dv) * math.sin((h + dh) * dth))
                     for h, x, y in tree[-1] for dh, ds in moves])
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True, gridspec_kw={"wspace": 0.05})
    colors = [BLUE, ORANGE, AQUA]
    for k in range(1, tau + 1):
        pts = np.array([(x, y) for _, x, y in tree[k]])
        keys, counts = np.unique(np.round(np.array([(h, x, y) for h, x, y in tree[k]]), 6), axis=0, return_counts=True)
        axs[0].scatter(pts[:, 1], pts[:, 0], s=10, color=colors[k - 1], lw=0, alpha=0.6, label=f"after {k} move{'s' if k > 1 else ''}")
        axs[1].scatter(keys[:, 2], keys[:, 1], s=10 * counts, color=colors[k - 1], lw=0, alpha=0.6)
        axs[1].text(0.98, 0.04 + 0.08 * (k - 1), f"level {k}: {len(tree[k])} paths → {len(keys)} nodes", transform=axs[1].transAxes,
                    ha="right", color=INK)
    axs[0].set_title(f"tree: {sum(len(t) for t in tree[1:])} nodes", color=INK, fontsize=10)
    axs[1].set_title("DAG: one node per distinct future", color=INK, fontsize=10)
    axs[0].legend(frameon=False, loc="lower left", labelcolor=INK)
    for a in axs:
        a.set_aspect("equal")
        a.set_xlabel("sideways (left +)")
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    axs[0].set_ylabel("ahead")
    save(fig, "dag.svg")


def fig_sphere():
    g = Grid(96, 48)
    p = dict(DEFAULTS)
    Q = np.array([[3.0, 1.2, 1.4], [0.6, 0.3, 4.2]])  # two neighbours: ahead and up-left, nearly overhead
    _, M = next(paint(np.zeros((1, 3)), Q, np.zeros(1, int), p, g))
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.imshow(M[0][::-1], cmap=matplotlib.colors.ListedColormap(["white", BLUE]), extent=(0, 360, -90, 90), aspect="auto", interpolation="nearest")
    for x in np.linspace(0, 360, g.nlon + 1):
        ax.axvline(x, color=GRID, lw=0.3)
    for y in np.linspace(-90, 90, g.nlat + 1):
        ax.axhline(y, color=GRID, lw=0.3)
    ax.set_xlabel("longitude in the yaw frame (degrees; 0 = ahead)")
    ax.set_ylabel("latitude (degrees)")
    ax.set_xticks(range(0, 361, 60))
    ax.set_yticks(range(-90, 91, 30))
    save(fig, "sphere.svg")


def fig_empowerment():
    sig = [0.05, 0.1, 0.2, 0.3]
    data = {"N = 500": {"$W_b$ (the paper)": [81, 89, 98, 93], "Gaussian bound": [87, 92, 99, 92], "$\\log n_b$ (survival)": [100, 84, 62, 62]},
            "N = 2500": {"$W_b$ (the paper)": [88, 98, 94, 82], "Gaussian bound": [88, 100, 95, 82], "$\\log n_b$ (survival)": [96, 80, 57, 49]}}
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=True, gridspec_kw={"wspace": 0.08})
    for a, (title, series) in zip(axs, data.items()):
        for (name, v), c, mk in zip(series.items(), (BLUE, ORANGE, AQUA), ("o", "s", "^")):
            a.plot(sig, v, color=c, lw=2, marker=mk, ms=6, label=name)
        a.set_title(title, color=INK, fontsize=10)
        a.set_xticks(sig)
        a.set_xlim(0.03, 0.42)
        a.set_xlabel("perception noise σ")
        a.set_ylim(40, 102)
        a.grid(axis="y", color=GRID, lw=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    axs[0].set_ylabel("decisions shared with empowerment (%)")
    axs[0].legend(frameon=False, loc="lower left", labelcolor=INK)
    for (name, v), c, dy in zip(data["N = 2500"].items(), (BLUE, ORANGE, AQUA), (-7, 7, 0)):  # direct labels
        axs[1].annotate(name.split(" (")[0], (sig[-1], v[-1]), (6, dy), textcoords="offset points", color=INK, fontsize=9, va="center")
    save(fig, "empowerment.svg")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    p, frame, i = flock()
    fig_futures(p, frame, i)
    fig_score(p, frame, i)
    fig_dag()
    fig_sphere()
    fig_empowerment()
