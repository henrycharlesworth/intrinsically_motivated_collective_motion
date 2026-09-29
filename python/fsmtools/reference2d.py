"""Brute-force reference for the 2D scores of single agents: every move sequence (5^tau paths) is enumerated
explicitly (no DAG), pruned by collisions as in the simulators, and each surviving future's visual field is painted
on B angular bins in that future's own frame, as the GPU does. Small tau only (it is exponential).

The futures of first move b (all levels pooled, one per surviving path) give the paper's score
    W_b = (1 / (pi n_b)) sum_bins S_b (n_b - S_b) dtheta,   S_b(bin) = number of futures that see something there.
"""
import numpy as np

MOVES = [(-1, 0), (1, 0), (0, -1), (0, 0), (0, 1)]


def futures(frame, i, p, tau, moves=MOVES):
    """Surviving futures of agent i: list over levels k = 1..tau of (first move, heading index, x, y) arrays,
    positions in agent i's frame (x along its heading), others ballistic at v0."""
    R, dt, v0, dv = p["R"], p["dt"], p["v0"], p["dv"]
    dth = np.radians(p["dtheta_deg"])
    x, y, th = frame[:, 0], frame[:, 1], frame[:, 2]
    others = np.arange(len(frame)) != i
    c, s = np.cos(th[i]), np.sin(th[i])
    rel = lambda X, Y: (c * (X - x[i]) + s * (Y - y[i]), -s * (X - x[i]) + c * (Y - y[i]))  # noqa: E731
    Q = [rel(x[others] + k * v0 * dt * np.cos(th[others]), y[others] + k * v0 * dt * np.sin(th[others]))
         for k in range(tau + 1)]
    D2 = 4 * R * R
    # alive paths: first move, heading index, position
    m0 = np.zeros(1, int); h = np.zeros(1, int); px = np.zeros(1); py = np.zeros(1)
    levels = []
    for k in range(1, tau + 1):
        nm0, nh, nx, ny, ox, oy = [], [], [], [], [], []
        for b, (dh, ds) in enumerate(moves):
            hh = h + dh
            v = v0 + ds * dv
            nm0.append(m0 if k > 1 else np.full(len(h), b)); nh.append(hh)
            nx.append(px + v * dt * np.cos(hh * dth)); ny.append(py + v * dt * np.sin(hh * dth))
            ox.append(px); oy.append(py)
        m0, h, px, py = (np.concatenate(a) for a in (nm0, nh, nx, ny))
        p0x, p0y = np.concatenate(ox), np.concatenate(oy)
        # collisions with every other agent during step k-1 -> k (relative motion, as in the simulators)
        r1x = px[:, None] - Q[k][0][None]; r1y = py[:, None] - Q[k][1][None]
        r0x = p0x[:, None] - Q[k - 1][0][None]; r0y = p0y[:, None] - Q[k - 1][1][None]
        end = r1x ** 2 + r1y ** 2 < D2
        dx, dy = r1x - r0x, r1y - r0y
        dd = dx ** 2 + dy ** 2
        sc = np.where(dd > 0, -(r0x * dx + r0y * dy) / np.where(dd > 0, dd, 1), -1)
        mid = (sc > 0) & (sc < 1) & ((r0x + sc * dx) ** 2 + (r0y + sc * dy) ** 2 < D2)
        start = r0x ** 2 + r0y ** 2 < D2
        hit = (end | (mid & ~start)).any(axis=1)
        m0, h, px, py = m0[~hit], h[~hit], px[~hit], py[~hit]
        levels.append((m0, h, px, py, Q[k]))
    return levels


def masks(levels, p, B):
    """Visual masks [n, B] (bin q set when its centre lies in an occupied interval), in each future's own frame,
    with the first move of each future."""
    dth = np.radians(p["dtheta_deg"])
    bpd = dth / (2 * np.pi / B)
    assert abs(bpd - round(bpd)) < 1e-9, "dtheta must be a whole number of bins"
    out, first = [], []
    for m0, h, px, py, (qx, qy) in levels:
        dx = qx[None] - px[:, None]; dy = qy[None] - py[:, None]
        d = np.hypot(dx, dy)
        half = np.arctan(p["R"] / d) if p["half_angle"] == "atan" else np.arcsin(np.minimum(1, p["R"] / d))
        ang = np.arctan2(dy, dx)
        # first/last bin whose centre is inside [ang - half, ang + half] (unwrapped; bins taken mod B)
        k0 = np.ceil((ang - half) * B / (2 * np.pi) - 0.5).astype(int)
        k1 = np.floor((ang + half) * B / (2 * np.pi) - 0.5).astype(int)
        diff = np.zeros((len(px), 2 * B + 1), np.int32)
        rows = np.broadcast_to(np.arange(len(px))[:, None], k0.shape)
        ok = k1 >= k0
        np.add.at(diff, (rows[ok], k0[ok] + B), 1)
        np.add.at(diff, (rows[ok], k1[ok] + B + 1), -1)
        cov = np.cumsum(diff[:, :-1], axis=1)
        mask = (cov[:, :B] + cov[:, B:]) > 0  # fold the unwrapped range [-B, B) onto [0, B)
        sh = np.rint(h * bpd).astype(int) % B
        idx = (np.arange(B)[None] + sh[:, None]) % B  # body bin q = root-frame bin q + sh
        out.append(np.take_along_axis(mask, idx, axis=1))
        first.append(m0)
    return np.concatenate(out), np.concatenate(first)


def var_scores(X, first, nm=len(MOVES)):
    """The paper's score of every first move from the masks X [n, B] of the futures and their first moves
    (-1 when no future survives)."""
    B = X.shape[1]
    out = []
    for b in range(nm):
        Xb = X[first == b]
        n = len(Xb)
        if n == 0:
            out.append(-1.0)
            continue
        S = Xb.sum(axis=0).astype(float)
        out.append(float((S * (n - S)).sum() * (2.0 / B) / n))
    return np.array(out)


def agent_scores(frame, i, p, tau, bins=3072):
    """The paper's scores of agent i's five first moves, by brute force."""
    X, first = masks(futures(frame, i, p, tau), p, bins)
    return var_scores(X, first)
