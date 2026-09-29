"""Flock measures for 2D and 3D runs. A frame is an array [N, 4] (x, y, theta, v) or [N, 6] (x, y, z, theta, phi, v).

Groups ("fragments") are the connected components of the graph linking agents closer than `link`. The default link
is 3 x the frame's median nearest-neighbour distance: a fixed distance splits large flocks, which are more dilute
than small ones, into many pieces even when they are visibly one flock.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, cKDTree


def dims_of_frame(frame):
    return 3 if frame.shape[-1] == 6 else 2


def positions(frame):
    return frame[..., : dims_of_frame(frame)]


def headings(frame):
    """Unit heading vectors [..., N, dims]."""
    if dims_of_frame(frame) == 2:
        th = frame[..., 2]
        return np.stack([np.cos(th), np.sin(th)], -1)
    th, ph = frame[..., 3], frame[..., 4]
    return np.stack([np.cos(ph) * np.cos(th), np.cos(ph) * np.sin(th), np.sin(ph)], -1)


def order_parameter(frames):
    """Polar order |mean heading| per frame (1: all aligned, ~0: disordered)."""
    return np.linalg.norm(headings(frames).mean(-2), axis=-1)


def nn_distance(frame):
    """Each agent's nearest-neighbour distance."""
    x = positions(frame)
    return cKDTree(x).query(x, k=2)[0][:, 1]


def group_labels(frame, link=None):
    """Group label of every agent (see the module docstring for the default link)."""
    x = positions(frame)
    link = 3 * np.median(nn_distance(frame)) if link is None else link
    n = len(x)
    pairs = cKDTree(x).query_pairs(link, output_type="ndarray")
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(A, directed=False)[1]


def largest_group(frame, link=None):
    """Share of the agents in the largest group."""
    return np.bincount(group_labels(frame, link)).max() / len(frame)


def radius_of_gyration(frame):
    x = positions(frame)
    return float(np.sqrt(((x - x.mean(0)) ** 2).sum(1).mean()))


def density(frame):
    """N / area (2D) or N / volume (3D) of the convex hull."""
    return len(frame) / ConvexHull(positions(frame)).volume


def visual_intervals(frame, i, R=1.0, half_angle="atan"):
    """2D: the merged angular intervals (relative to agent i's heading, in [0, 2 pi)) occupied by the other agents."""
    pos, heading = frame[:, :2], frame[i, 2]
    r = np.delete(pos, i, axis=0) - pos[i]
    d = np.hypot(r[:, 0], r[:, 1])
    ang = np.mod(np.arctan2(r[:, 1], r[:, 0]) - heading, 2 * np.pi)
    half = np.arctan(R / d) if half_angle == "atan" else np.arcsin(np.minimum(1, R / d))
    ivs = []
    for a, b in zip(ang - half, ang + half):
        if a < 0:
            ivs += [(a + 2 * np.pi, 2 * np.pi), (0, b)]
        elif b > 2 * np.pi:
            ivs += [(a, 2 * np.pi), (0, b - 2 * np.pi)]
        else:
            ivs.append((a, b))
    ivs.sort()
    merged = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return np.array(merged)


def opacity(frame, R=1.0, half_angle="atan", agents=None):
    """2D: mean fraction of an agent's visual field that is occupied."""
    agents = range(len(frame)) if agents is None else agents
    vals = []
    for i in agents:
        iv = visual_intervals(frame, i, R, half_angle)
        vals.append((iv[:, 1] - iv[:, 0]).sum() / (2 * np.pi) if len(iv) else 0.0)
    return float(np.mean(vals))


def collisions(traj, R=1.0):
    """Collisions per step: pairs whose straight-line motion during the step brings them within 2R, and that did
    not already overlap at its start."""
    d = dims_of_frame(traj[0])
    counts = np.zeros(len(traj) - 1, dtype=int)
    for t in range(len(traj) - 1):
        a, b = traj[t, :, :d], traj[t + 1, :, :d]
        reach = np.linalg.norm(b - a, axis=1).max() * 2 + 2 * R
        pairs = cKDTree(a).query_pairs(reach, output_type="ndarray")
        if len(pairs) == 0:
            continue
        i, j = pairs[:, 0], pairs[:, 1]
        r0 = a[i] - a[j]
        dr = (b[i] - b[j]) - r0
        dd = (dr ** 2).sum(1)
        s = np.clip(-(r0 * dr).sum(1) / np.where(dd > 0, dd, 1), 0, 1)
        dmin2 = ((r0 + s[:, None] * dr) ** 2).sum(1)
        counts[t] = np.sum((dmin2 < 4 * R * R) & ((r0 ** 2).sum(1) >= 4 * R * R))
    return counts


def summary(run, burn_in=200, every=10):
    """Time-averaged measures of a run after `burn_in` steps (every `every`-th frame)."""
    traj = np.asarray(run["traj"])
    p = run["params"]
    sel = traj[burn_in::every]
    if len(sel) == 0:
        raise ValueError(f"run has {len(traj)} frames, fewer than burn_in = {burn_in}")
    big = np.array([largest_group(f) for f in sel])
    out = {
        "frames": len(traj),
        "order": float(order_parameter(sel).mean()),
        "largest_group": float(big.mean()),
        "largest_group_min": float(big.min()),
        "groups": float(np.mean([group_labels(f).max() + 1 for f in sel])),
        "nn_distance": float(np.mean([np.median(nn_distance(f)) for f in sel])),
        "radius_of_gyration": float(np.mean([radius_of_gyration(f) for f in sel])),
        "collisions_per_agent_step": float(collisions(traj[burn_in:], p["R"]).sum() / (traj.shape[1] * max(1, len(traj) - burn_in - 1))),
    }
    if dims_of_frame(traj[0]) == 2:
        few = sel[:: max(1, len(sel) // 5)]
        out["opacity"] = float(np.mean([opacity(f, p["R"], p.get("half_angle", "atan"), range(0, len(f), max(1, len(f) // 50))) for f in few]))
    return out
