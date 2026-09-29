// The FSM decision of one 2D agent in JavaScript, for the viewer's bird panel, plus the trajectory replay and run
// analysis the viewer needs. Two evaluators:
//   evaluateAgent  exact, mirroring the CPU simulator (src/fsm2d/cpu_evaluator.hpp): DAG template, collision
//                  pruning by path counting, pseudo-angle visual states, W_b = (1 / (pi n_b)) integral S_b (n_b - S_b)
//   buildFastTemplate / prepareAgent / paintRange / combineFast
//                  the GPU kernel's algorithm (src/fsm2d/gpu_kernel.cuh): positions, units, bitmasks of 3072 bins and
//                  interval edges; split over web workers by the viewer
// Works as a browser global (FSMCore), inside web workers (FSMCore.workerSource()) and as a Node module.
function FSMCoreModule(root) {
  "use strict";
  const NM = 5;
  const TWO_PI = 2 * Math.PI;

  function buildTemplate(p) {
    const T = p.tau, H = 2 * T + 1, dth = (p.dtheta_deg * Math.PI) / 180;
    const h = [], ox = [], oy = [], lstart = [0, 1], eparent = [], echild = [], emove = [];
    const estart = new Array(T + 2).fill(0);
    const addNode = (key) => {
      let x = 0, y = 0;
      for (let slot = 0; slot < H; slot++) {
        const m = key[1 + 2 * slot], s = key[2 + 2 * slot];
        if (m === 0) continue;
        const dist = (m * p.v0 + s * p.dv) * p.dt, ang = (slot - T) * dth;
        x += dist * Math.cos(ang);
        y += dist * Math.sin(ang);
      }
      h.push(key[0]); ox.push(x); oy.push(y);
      return h.length - 1;
    };
    let curKeys = [new Array(1 + 2 * H).fill(0)];
    let curIds = [addNode(curKeys[0])];
    for (let k = 1; k <= T; k++) {
      estart[k] = eparent.length;
      const index = new Map(), nextKeys = [], nextIds = [];
      for (let q = 0; q < curKeys.length; q++) {
        for (let m = 0; m < NM; m++) {
          const ck = curKeys[q].slice();
          let hh = ck[0], ds = 0;
          if (m === 0) hh -= 1; else if (m === 1) hh += 1; else ds = m - 3;
          ck[0] = hh; ck[1 + 2 * (hh + T)] += 1; ck[2 + 2 * (hh + T)] += ds;
          const sk = ck.join(",");
          let cid = index.get(sk);
          if (cid === undefined) {
            cid = addNode(ck); index.set(sk, cid); nextKeys.push(ck); nextIds.push(cid);
          }
          eparent.push(curIds[q]); echild.push(cid); emove.push(m);
        }
      }
      lstart.push(h.length);
      curKeys = nextKeys; curIds = nextIds;
    }
    estart[T + 1] = eparent.length;
    const level = new Uint8Array(h.length);
    for (let k = 1; k <= T; k++) for (let n = lstart[k]; n < lstart[k + 1]; n++) level[n] = k;
    return { tau: T, dth, h: Int8Array.from(h), ox: Float64Array.from(ox), oy: Float64Array.from(oy), lstart,
      eparent: Int32Array.from(eparent), echild: Int32Array.from(echild), emove: Int8Array.from(emove), estart, level };
  }

  function pseudo(x, y) {
    const a = x / (Math.abs(x) + Math.abs(y));
    return y >= 0 ? 1 - a : 3 + a;
  }
  function pseudoToAngle(p) {
    if (p <= 0) return 0;
    if (p >= 4) return TWO_PI;
    if (p <= 2) { const a = 1 - p; return Math.atan2(1 - Math.abs(a), a); }
    const a = p - 3;
    return Math.atan2(-(1 - Math.abs(a)), a) + TWO_PI;
  }

  // frame: {x, y, th} typed arrays of length N. Returns everything the viewer draws for agent i.
  function evaluateAgent(p, t, frame, i, bins) {
    bins = bins || 360;
    const N = frame.x.length, T = t.tau, nn = t.h.length, R = p.R, R2 = R * R, D2 = 4 * R2;
    const asin = p.half_angle === "asin";
    // Ballistic futures of everyone.
    const QX = new Float64Array((T + 1) * N), QY = new Float64Array((T + 1) * N);
    for (let j = 0; j < N; j++) {
      const cx = Math.cos(frame.th[j]) * p.v0 * p.dt, cy = Math.sin(frame.th[j]) * p.v0 * p.dt;
      for (let k = 0; k <= T; k++) { QX[k * N + j] = frame.x[j] + k * cx; QY[k * N + j] = frame.y[j] + k * cy; }
    }
    const x0 = frame.x[i], y0 = frame.y[i], th0 = frame.th[i], c0 = Math.cos(th0), s0 = Math.sin(th0);
    const ch = new Float64Array(2 * T + 1), sh = new Float64Array(2 * T + 1);
    for (let hh = -T; hh <= T; hh++) { ch[hh + T] = Math.cos(th0 + hh * t.dth); sh[hh + T] = Math.sin(th0 + hh * t.dth); }
    const px = new Float64Array(nn), py = new Float64Array(nn);
    for (let n = 0; n < nn; n++) { px[n] = x0 + c0 * t.ox[n] - s0 * t.oy[n]; py[n] = y0 + s0 * t.ox[n] + c0 * t.oy[n]; }
    // Collision candidates per level (bbox test).
    const cand = [[]];
    for (let k = 1; k <= T; k++) {
      let bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity;
      for (let n = t.lstart[k - 1]; n < t.lstart[k + 1]; n++) {
        if (px[n] < bx0) bx0 = px[n]; if (px[n] > bx1) bx1 = px[n];
        if (py[n] < by0) by0 = py[n]; if (py[n] > by1) by1 = py[n];
      }
      bx0 -= 2 * R; bx1 += 2 * R; by0 -= 2 * R; by1 += 2 * R;
      const cl = [];
      for (let j = 0; j < N; j++) {
        if (j === i) continue;
        const ax = QX[(k - 1) * N + j], bx = QX[k * N + j], ay = QY[(k - 1) * N + j], by = QY[k * N + j];
        if (Math.max(ax, bx) < bx0 || Math.min(ax, bx) > bx1) continue;
        if (Math.max(ay, by) < by0 || Math.min(ay, by) > by1) continue;
        cl.push(j);
      }
      cand.push(cl);
    }
    const collides = (k, p0x, p0y, p1x, p1y) => {
      for (const j of cand[k]) {
        const r1x = p1x - QX[k * N + j], r1y = p1y - QY[k * N + j];
        if (r1x * r1x + r1y * r1y < D2) return true;
        const r0x = p0x - QX[(k - 1) * N + j], r0y = p0y - QY[(k - 1) * N + j];
        if (r0x * r0x + r0y * r0y < D2) continue;
        const dx = r1x - r0x, dy = r1y - r0y, dd = dx * dx + dy * dy;
        if (dd <= 0) continue;
        const s = -(r0x * dx + r0y * dy) / dd;
        if (s > 0 && s < 1) { const mx = r0x + s * dx, my = r0y + s * dy; if (mx * mx + my * my < D2) return true; }
      }
      return false;
    };
    // Path multiplicities per branch.
    const mult = new Int32Array(nn * NM);
    for (let e = t.estart[1]; e < t.estart[2]; e++) {
      const c = t.echild[e];
      if (!collides(1, x0, y0, px[c], py[c])) mult[c * NM + t.emove[e]] = 1;
    }
    for (let k = 2; k <= T; k++) {
      for (let e = t.estart[k]; e < t.estart[k + 1]; e++) {
        const a = t.eparent[e], c = t.echild[e];
        let any = 0;
        for (let b = 0; b < NM; b++) any |= mult[a * NM + b];
        if (!any) continue;
        if (collides(k, px[a], py[a], px[c], py[c])) continue;
        for (let b = 0; b < NM; b++) mult[c * NM + b] += mult[a * NM + b];
      }
    }
    // Visual states of reachable nodes -> exact bin coverage per branch and the sweep score.
    const lo = new Float64Array(2 * N + 2), hi = new Float64Array(2 * N + 2), ord = new Int32Array(2 * N + 2);
    const nb = new Float64Array(NM), pb = new Float64Array(NM * bins), bits = new Uint8Array(nn);
    const evTh = [], evNode = [], evSign = [];
    const bw = TWO_PI / bins;
    for (let k = 1; k <= T; k++) {
      const qo = k * N;
      for (let n = t.lstart[k]; n < t.lstart[k + 1]; n++) {
        let bt = 0;
        for (let b = 0; b < NM; b++) if (mult[n * NM + b] > 0) { bt |= 1 << b; nb[b] += mult[n * NM + b]; }
        bits[n] = bt;
        if (!bt) continue;
        const hh = t.h[n] + T, c = ch[hh], s = sh[hh], X = px[n], Y = py[n];
        let m = 0, full = false;
        for (let j = 0; j < N; j++) {
          if (j === i) continue;
          const rx = QX[qo + j] - X, ry = QY[qo + j] - Y;
          const u = c * rx + s * ry, v = -s * rx + c * ry, d2 = u * u + v * v;
          let a;
          if (asin) { if (d2 <= R2) { full = true; break; } a = Math.sqrt(d2 - R2); }
          else { if (d2 <= 0) continue; a = Math.sqrt(d2); }
          const l = pseudo(u * a + v * R, -u * R + v * a), r = pseudo(u * a - v * R, u * R + v * a);
          if (l <= r) { lo[m] = l; hi[m] = r; m++; }
          else { lo[m] = l; hi[m] = 4; m++; lo[m] = 0; hi[m] = r; m++; }
        }
        const merged = [];
        if (full) merged.push(0, TWO_PI);
        else if (m > 0) {
          for (let q = 0; q < m; q++) ord[q] = q;
          const o = Array.prototype.slice.call(ord, 0, m).sort((A, B) => lo[A] - lo[B]);
          let clo = lo[o[0]], chi = hi[o[0]];
          for (let q = 1; q < m; q++) {
            const L = lo[o[q]], H = hi[o[q]];
            if (L <= chi) { if (H > chi) chi = H; }
            else { merged.push(pseudoToAngle(clo), pseudoToAngle(chi)); clo = L; chi = H; }
          }
          merged.push(pseudoToAngle(clo), pseudoToAngle(chi));
        }
        for (let q = 0; q < merged.length; q += 2) {
          const a1 = merged[q], a2 = merged[q + 1];
          evTh.push(a1, a2); evNode.push(n, n); evSign.push(1, -1);
          const i0 = Math.min(bins - 1, Math.floor(a1 / bw)), i1 = Math.min(bins - 1, Math.floor(a2 / bw));
          for (let kk = i0; kk <= i1; kk++) {
            const cov = (Math.min(a2, (kk + 1) * bw) - Math.max(a1, kk * bw)) / bw;
            if (cov <= 0) continue;
            for (let b = 0; b < NM; b++) pb[b * bins + kk] += cov * mult[n * NM + b];
          }
        }
      }
    }
    for (let b = 0; b < NM; b++) for (let kk = 0; kk < bins; kk++) pb[b * bins + kk] = nb[b] > 0 ? pb[b * bins + kk] / nb[b] : 0;
    // Exact sweep score.
    const eo = evTh.map((_, q) => q).sort((A, B) => evTh[A] - evTh[B]);
    const S = new Float64Array(NM), acc = new Float64Array(NM);
    let prev = 0;
    for (const q of eo) {
      const d = evTh[q] - prev;
      if (d > 0) { for (let b = 0; b < NM; b++) acc[b] += S[b] * (nb[b] - S[b]) * d; prev = evTh[q]; }
      const n = evNode[q];
      for (let b = 0; b < NM; b++) S[b] += evSign[q] * mult[n * NM + b];
    }
    const score = new Float64Array(NM);
    for (let b = 0; b < NM; b++) {
      if (nb[b] <= 0) score[b] = -1;
      else score[b] = acc[b] / (Math.PI * nb[b]);
    }
    return { px, py, bits, level: t.level, p: pb, score, nb, bins };
  }

  // Current visual field of agent i (merged intervals relative to its heading), for the inner ring.
  function visualField(p, frame, i) {
    const N = frame.x.length, R = p.R, c = Math.cos(frame.th[i]), s = Math.sin(frame.th[i]);
    const iv = [];
    for (let j = 0; j < N; j++) {
      if (j === i) continue;
      const rx = frame.x[j] - frame.x[i], ry = frame.y[j] - frame.y[i];
      const u = c * rx + s * ry, v = -s * rx + c * ry, d = Math.hypot(u, v);
      if (d <= 0) continue;
      let ang = Math.atan2(v, u); if (ang < 0) ang += TWO_PI;
      const half = p.half_angle === "asin" ? Math.asin(Math.min(1, R / d)) : Math.atan(R / d);
      const a = ang - half, b = ang + half;
      if (a < 0) { iv.push([a + TWO_PI, TWO_PI], [0, b]); }
      else if (b > TWO_PI) { iv.push([a, TWO_PI], [0, b - TWO_PI]); }
      else iv.push([a, b]);
    }
    iv.sort((A, B) => A[0] - B[0]);
    const out = [];
    for (const [a, b] of iv) {
      if (out.length && a <= out[out.length - 1][1]) out[out.length - 1][1] = Math.max(out[out.length - 1][1], b);
      else out.push([a, b]);
    }
    return out;
  }

  // Connected components of agents within `link` of each other (grid hash + union-find).
  function clusterLabels(traj, k, N, link) {
    const parent = new Int32Array(N);
    for (let j = 0; j < N; j++) parent[j] = j;
    const find = (a) => { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
    const cells = new Map(), o = k * N * 3, L2 = link * link;
    const key = (cx, cy) => (cx + 1048576) * 2097152 + (cy + 1048576);
    for (let j = 0; j < N; j++) {
      const cx = Math.floor(traj[o + 3 * j] / link), cy = Math.floor(traj[o + 3 * j + 1] / link);
      for (let ddx = -1; ddx <= 1; ddx++) for (let ddy = -1; ddy <= 1; ddy++) {
        const list = cells.get(key(cx + ddx, cy + ddy));
        if (!list) continue;
        for (const q of list) {
          const dx = traj[o + 3 * j] - traj[o + 3 * q], dy = traj[o + 3 * j + 1] - traj[o + 3 * q + 1];
          if (dx * dx + dy * dy <= L2) { const a = find(j), b = find(q); if (a !== b) parent[a] = b; }
        }
      }
      const kk = key(cx, cy);
      if (!cells.has(kk)) cells.set(kk, []);
      cells.get(kk).push(j);
    }
    const lab = new Int32Array(N), ids = new Map();
    for (let j = 0; j < N; j++) { const r = find(j); if (!ids.has(r)) ids.set(r, ids.size); lab[j] = ids.get(r); }
    return { lab, n: ids.size };
  }

  function smooth1d(a, sigma) {
    if (sigma <= 0) return Float64Array.from(a);
    const r = Math.ceil(4 * sigma), w = [];
    let tot = 0;
    for (let q = -r; q <= r; q++) { const v = Math.exp((-0.5 * q * q) / (sigma * sigma)); w.push(v); tot += v; }
    const n = a.length, out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let q = -r; q <= r; q++) s += w[q + r] * a[Math.min(n - 1, Math.max(0, i + q))];
      out[i] = s / tot;
    }
    return out;
  }

  // Positions from the initial state and the moves, as the simulators integrate them (double precision):
  // theta += dh dtheta, v = v0 + ds dv, x += v dt cos(theta), y += v dt sin(theta). Paper moves (dh, ds).
  // stride > 1 keeps one frame in `stride` (frames 0, stride, 2 stride, ...), for runs too large to hold every frame;
  // every step is still integrated. Returns [ceil(T / stride), N, 3].
  function replay2d(P, init, mv, T, N, stride = 1) {
    const M = [[-1, 0], [1, 0], [0, -1], [0, 0], [0, 1]], dth = (P.dtheta_deg * Math.PI) / 180;
    const out = new Float32Array(Math.ceil(T / stride) * N * 3), x = new Float64Array(N), y = new Float64Array(N), th = new Float64Array(N);
    for (let j = 0; j < N; j++) { x[j] = init[4 * j]; y[j] = init[4 * j + 1]; th[j] = init[4 * j + 2]; }
    for (let k = 0; k < T; k++) {
      if (k > 0) {
        const mo = (k - 1) * N;
        for (let j = 0; j < N; j++) {
          const m = M[mv[mo + j]];
          th[j] += m[0] * dth;
          const v = P.v0 + m[1] * P.dv;
          x[j] += v * P.dt * Math.cos(th[j]);
          y[j] += v * P.dt * Math.sin(th[j]);
        }
      }
      if (k % stride) continue;
      const q = (k / stride) * N * 3;
      for (let j = 0; j < N; j++) { out[q + 3 * j] = x[j]; out[q + 3 * j + 1] = y[j]; out[q + 3 * j + 2] = th[j]; }
    }
    return out;
  }
  // The moves taken at the kept frames of a strided replay: [ceil(T / stride) - 1, N].
  function strideMoves(mv, T, N, stride) {
    if (stride === 1) return mv;
    const K = Math.ceil(T / stride), out = new Uint8Array(Math.max(0, K - 1) * N);
    for (let q = 0; q < K - 1; q++) out.set(mv.subarray(q * stride * N, (q * stride + 1) * N), q * N);
    return out;
  }
  // analyseRun's result as one Float32Array (viewer/tools/analyse.js precomputes it for bundled runs, so that the viewer
  // need not cluster every frame on load): [cam 4T | swarmCam 4T | order T | meanHead T | fragments T | focal, flockMostly].
  function packAnalysis(a, T) {
    const f = new Float32Array(11 * T + 2);
    f.set(a.cam, 0); f.set(a.swarmCam, 4 * T); f.set(a.order, 8 * T); f.set(a.meanHead, 9 * T);
    for (let k = 0; k < T; k++) f[10 * T + k] = a.fragments[k];
    f[11 * T] = a.focal; f[11 * T + 1] = a.flockMostly ? 1 : 0;
    return f;
  }
  function unpackAnalysis(f, T) {
    if (f.length !== 11 * T + 2) return null;  // an older layout: the viewer recomputes it
    return { cam: f.subarray(0, 4 * T), swarmCam: f.subarray(4 * T, 8 * T), order: f.subarray(8 * T, 9 * T),
             meanHead: f.subarray(9 * T, 10 * T), fragments: Int32Array.from(f.subarray(10 * T, 11 * T)), focal: f[11 * T],
             flockMostly: f[11 * T + 1] > 0 };
  }

  // Per-run derived data for the viewer: camera track following the largest (tracked) fragment,
  // order parameter, fragment counts and a default focal agent.
  // Median nearest-neighbour distance over a few frames (grid hash; neighbours further than `cell` are ignored).
  function medianNN(traj, T, N, cell) {
    const ds = [];
    for (const k of [0.25, 0.5, 0.75, 1].map((f) => Math.min(T - 1, Math.floor(f * (T - 1))))) {
      const o = k * N * 3, cells = new Map(), key = (cx, cy) => (cx + 1048576) * 2097152 + (cy + 1048576);
      for (let j = 0; j < N; j++) {
        const kk = key(Math.floor(traj[o + 3 * j] / cell), Math.floor(traj[o + 3 * j + 1] / cell));
        if (!cells.has(kk)) cells.set(kk, []);
        cells.get(kk).push(j);
      }
      for (let j = 0; j < N; j += Math.max(1, Math.floor(N / 400))) {
        const cx = Math.floor(traj[o + 3 * j] / cell), cy = Math.floor(traj[o + 3 * j + 1] / cell);
        let best = Infinity;
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++)
          for (const q of cells.get(key(cx + dx, cy + dy)) || []) {
            if (q === j) continue;
            best = Math.min(best, Math.hypot(traj[o + 3 * j] - traj[o + 3 * q], traj[o + 3 * j + 1] - traj[o + 3 * q + 1]));
          }
        if (best < Infinity) ds.push(best);
      }
    }
    ds.sort((a, b) => a - b);
    return ds.length ? ds[ds.length >> 1] : 0;
  }

  // Fragments use link distance max(20, 3 x median nn distance): large FSM flocks equilibrate dilute, and a fixed
  // 20-unit link would split a visibly coherent flock (python/summarize_runs.py reports the same measure).
  function analyseRun(traj, T, N, opts) {
    const link = (opts && opts.link) || Math.max(20, 3 * medianNN(traj, T, N, 40)), sm = (opts && opts.smooth) || 8, switchRatio = 1.5;
    const order = new Float32Array(T), fragments = new Int32Array(T), followed = new Int32Array(T);
    const meanHead = new Float32Array(T);
    const cx = new Float64Array(T), cy = new Float64Array(T), rad = new Float64Array(T), head = new Float64Array(T);
    const alwaysIn = new Uint8Array(N).fill(1), distSum = new Float64Array(N);
    // whole-swarm track (for fragmented runs, where the largest fragment keeps changing): median centre and the
    // 90th-percentile distance of a fixed sample of birds
    const stride = Math.max(1, Math.ceil(N / 2000)), ns = Math.ceil(N / stride), xs = new Float64Array(ns), ys = new Float64Array(ns);
    const wx = new Float64Array(T), wy = new Float64Array(T), wr = new Float64Array(T), wh = new Float64Array(T);
    const median = (a) => { const b = Float64Array.from(a).sort(); return b[b.length >> 1]; };
    let prev = null;
    for (let k = 0; k < T; k++) {
      const o = k * N * 3;
      {
        let q = 0;
        for (let j = 0; j < N; j += stride, q++) { xs[q] = traj[o + 3 * j]; ys[q] = traj[o + 3 * j + 1]; }
        wx[k] = median(xs); wy[k] = median(ys);
        for (q = 0; q < ns; q++) xs[q] = Math.hypot(xs[q] - wx[k], ys[q] - wy[k]);
        const b = xs.slice().sort();
        wr[k] = b[Math.min(ns - 1, Math.floor(0.9 * (ns - 1)))];
      }
      let sx = 0, sy = 0;
      for (let j = 0; j < N; j++) { sx += Math.cos(traj[o + 3 * j + 2]); sy += Math.sin(traj[o + 3 * j + 2]); }
      order[k] = Math.hypot(sx, sy) / N; meanHead[k] = Math.atan2(sy, sx); wh[k] = meanHead[k];
      const { lab, n } = clusterLabels(traj, k, N, link);
      fragments[k] = n;
      const sizes = new Int32Array(n);
      for (let j = 0; j < N; j++) sizes[lab[j]]++;
      let big = 0;
      for (let c = 1; c < n; c++) if (sizes[c] > sizes[big]) big = c;
      let c = big;
      if (prev) {
        const ov = new Int32Array(n);
        for (let j = 0; j < N; j++) if (prev[j]) ov[lab[j]]++;
        c = 0;
        for (let q = 1; q < n; q++) if (ov[q] > ov[c]) c = q;
        if (sizes[big] > switchRatio * sizes[c]) c = big;
      }
      const mask = new Uint8Array(N);
      let mx = 0, my = 0, hx = 0, hy = 0, m = 0;
      for (let j = 0; j < N; j++) if (lab[j] === c) {
        mask[j] = 1; m++;
        mx += traj[o + 3 * j]; my += traj[o + 3 * j + 1];
        hx += Math.cos(traj[o + 3 * j + 2]); hy += Math.sin(traj[o + 3 * j + 2]);
      }
      mx /= m; my /= m;
      followed[k] = m; cx[k] = mx; cy[k] = my; head[k] = Math.atan2(hy, hx);
      const d = [];
      for (let j = 0; j < N; j++) {
        const dd = Math.hypot(traj[o + 3 * j] - mx, traj[o + 3 * j + 1] - my);
        distSum[j] += dd;
        if (mask[j]) d.push(dd); else alwaysIn[j] = 0;
      }
      d.sort((a, b) => a - b);
      rad[k] = d[Math.min(d.length - 1, Math.floor(0.97 * (d.length - 1)))];
      prev = mask;
    }
    for (const h of [head, wh])
      for (let k = 1; k < T; k++) {  // unwrap heading
        let dh = h[k] - h[k - 1];
        while (dh > Math.PI) { h[k] -= TWO_PI; dh -= TWO_PI; }
        while (dh < -Math.PI) { h[k] += TWO_PI; dh += TWO_PI; }
      }
    const scx = smooth1d(cx, sm), scy = smooth1d(cy, sm), srad = smooth1d(rad, sm * 3), shead = smooth1d(head, sm * 2);
    const cam = new Float32Array(T * 4);
    for (let k = 0; k < T; k++) { cam[4 * k] = scx[k]; cam[4 * k + 1] = scy[k]; cam[4 * k + 2] = srad[k]; cam[4 * k + 3] = shead[k]; }
    // the swarm camera is smoothed much more (a slow drift and turn, whatever the fragments do)
    const swx = smooth1d(wx, 30), swy = smooth1d(wy, 30), swr = smooth1d(wr, 60), swh = smooth1d(wh, 80);
    const swarmCam = new Float32Array(T * 4);
    for (let k = 0; k < T; k++) { swarmCam[4 * k] = swx[k]; swarmCam[4 * k + 1] = swy[k]; swarmCam[4 * k + 2] = swr[k]; swarmCam[4 * k + 3] = swh[k]; }
    // "auto" follows the largest fragment when it is most of the flock most of the time and its track is steady
    // (99th percentile of the camera's speed within 2.5x the swarm camera's), else the whole swarm
    const share = Float64Array.from(followed, (m) => m / N).sort();
    const p99speed = (c) => {
      const v = new Float64Array(Math.max(1, T - 1));
      for (let k = 1; k < T; k++) v[k - 1] = Math.hypot(c[4 * k] - c[4 * k - 4], c[4 * k + 1] - c[4 * k - 3]);
      v.sort();
      return v[Math.floor(0.99 * (v.length - 1))];
    };
    const flockMostly = share[share.length >> 1] >= 0.5 && p99speed(cam) <= 2.5 * Math.max(1, p99speed(swarmCam));
    let cands = [];
    for (let j = 0; j < N; j++) if (alwaysIn[j]) cands.push(j);
    if (!cands.length) cands = Array.from({ length: N }, (_, j) => j);
    cands.sort((a, b) => distSum[a] - distSum[b]);
    const focal = cands[Math.floor(cands.length / 2)];
    return { cam, swarmCam, flockMostly, order, fragments, followed, meanHead, focal, link };
  }


  // ------------------------------------------------------------------------------------------------
  // Fast path, mirroring src/fsm2d/gpu_kernel.cuh: future positions share one agent-frame bitmask of
  // B angular bins (nodes within `merge` of each other share a position, as in the run), nodes with the same
  // heading at a position form a unit, and each unit adds +-(its path counts) at its mask's interval edges,
  // rotated by its heading, into a difference array of S_b. Work is split over positions, so several workers
  // can each take a slice (paintRange) and the main thread sums them (combineFast).
  // ------------------------------------------------------------------------------------------------
  const FAST_BINS = 3072;
  // Same fitted polynomials as the GPU kernel: atan2 max error 1.7e-6 rad, atan(x <= 1/2) 5e-7 rad (bin = 2e-3 rad).
  function fastAtan2(y, x) {
    const ax = Math.abs(x), ay = Math.abs(y), mx = ax > ay ? ax : ay, mn = ax > ay ? ay : ax;
    const t = mn / mx, s = t * t;
    let r = ((((-1.171559314e-02 * s + 5.263565089e-02) * s - 1.164120617e-01) * s + 1.935322990e-01) * s - 3.326208634e-01) * s + 9.999770681e-01;
    r *= t;
    if (ay > ax) r = 1.57079632679 - r;
    if (x < 0) r = 3.14159265359 - r;
    return y < 0 ? -r : r;
  }
  function halfAtan(x) {
    if (x > 0.5) return Math.atan(x);
    const s = x * x;
    return (((-9.578907725e-02 * s + 1.920534081e-01) * s - 3.328121439e-01) * s + 9.999906105e-01) * x;
  }

  function buildFastTemplate(p) {
    const t = buildTemplate(p);
    const T = t.tau, nn = t.h.length, B = FAST_BINS;
    const merge = p.merge === undefined || p.merge === null ? 0 : p.merge;
    const eps = merge > 0 ? merge : 1e-6;
    const binsPerDth = (t.dth * B) / TWO_PI;
    // positions: per level, nodes snapped to an eps grid; projection point = mean of the merged nodes
    const posX = [], posY = [], posLevel = [], posU0 = [], unitH = [], unitSh = [], unitN0 = [], unitNodes = [];
    for (let k = 1; k <= T; k++) {
      const map = new Map(), lists = [];
      for (let n = t.lstart[k]; n < t.lstart[k + 1]; n++) {
        const key = Math.round(t.ox[n] / eps) + "," + Math.round(t.oy[n] / eps);
        let q = map.get(key);
        if (q === undefined) { q = lists.length; map.set(key, q); lists.push([]); }
        lists[q].push(n);
      }
      for (const nodes of lists) {
        let x = 0, y = 0;
        for (const n of nodes) { x += t.ox[n]; y += t.oy[n]; }
        posX.push(x / nodes.length); posY.push(y / nodes.length); posLevel.push(k); posU0.push(unitH.length);
        nodes.sort((a, b) => t.h[a] - t.h[b]);
        for (let z = 0; z < nodes.length; z++) {
          if (z === 0 || t.h[nodes[z]] !== t.h[nodes[z - 1]]) {
            unitH.push(t.h[nodes[z]]);
            unitSh.push((((Math.round(t.h[nodes[z]] * binsPerDth)) % B) + B) % B);
            unitN0.push(unitNodes.length);
          }
          unitNodes.push(nodes[z]);
        }
      }
    }
    posU0.push(unitH.length); unitN0.push(unitNodes.length);
    // collision bboxes per level (nodes at levels k-1..k), in the agent frame
    const lbb = new Float64Array(4 * (T + 1));
    for (let k = 1; k <= T; k++) {
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (let n = t.lstart[k - 1]; n < t.lstart[k + 1]; n++) {
        x0 = Math.min(x0, t.ox[n]); x1 = Math.max(x1, t.ox[n]); y0 = Math.min(y0, t.oy[n]); y1 = Math.max(y1, t.oy[n]);
      }
      lbb.set([x0, x1, y0, y1], 4 * k);
    }
    return Object.assign(t, {
      B, merge, P: posX.length, posX: Float64Array.from(posX), posY: Float64Array.from(posY),
      posLevel: Int8Array.from(posLevel), posU0: Int32Array.from(posU0), unitSh: Int32Array.from(unitSh),
      unitN0: Int32Array.from(unitN0), unitNodes: Int32Array.from(unitNodes), lbb, nodeCount: nn,
    });
  }

  // Agent-frame futures of everyone else (level-major, [T+1][N-1]) and collision-free path counts per branch.
  function prepareAgent(p, t, frame, i) {
    const N = frame.x.length, T = t.tau, M = N - 1, R = p.R, D2 = 4 * R * R;
    const x0 = frame.x[i], y0 = frame.y[i], c0 = Math.cos(frame.th[i]), s0 = Math.sin(frame.th[i]);
    const rx = new Float64Array((T + 1) * M), ry = new Float64Array((T + 1) * M);
    for (let j = 0, q = 0; j < N; j++) {
      if (j === i) continue;
      const ux = Math.cos(frame.th[j]) * p.v0 * p.dt, uy = Math.sin(frame.th[j]) * p.v0 * p.dt;
      for (let k = 0; k <= T; k++) {
        const dx = frame.x[j] + k * ux - x0, dy = frame.y[j] + k * uy - y0;
        rx[k * M + q] = c0 * dx + s0 * dy; ry[k * M + q] = -s0 * dx + c0 * dy;
      }
      q++;
    }
    const mult = new Int32Array(t.nodeCount * NM);
    // Collision candidates per level on a flat grid keyed by segment midpoints: an edge (length <= (v0+dv)dt) can
    // only come within 2R of a candidate segment (length v0 dt) whose midpoint is within half of each plus 2R of
    // the edge's midpoint, so one cell ring around the edge midpoint covers it when the cell is at least that big.
    const m2 = 2 * R, C = 0.5 * (p.v0 + p.dv) * p.dt + 0.5 * p.v0 * p.dt + m2;
    const cand = new Int32Array(M);
    for (let k = 1; k <= T; k++) {
      const bb = t.lbb, gx0 = bb[4 * k] - m2 - C, gy0 = bb[4 * k + 2] - m2 - C;
      const GX = Math.max(1, Math.ceil((bb[4 * k + 1] - bb[4 * k] + 2 * m2 + 2 * C) / C)), GY = Math.max(1, Math.ceil((bb[4 * k + 3] - bb[4 * k + 2] + 2 * m2 + 2 * C) / C));
      const cnt = new Int32Array(GX * GY + 1), cellOf = new Int32Array(M).fill(-1);
      for (let q = 0; q < M; q++) {
        const ax = rx[(k - 1) * M + q], ay = ry[(k - 1) * M + q], bx = rx[k * M + q], by = ry[k * M + q];
        if (Math.max(ax, bx) < bb[4 * k] - m2 || Math.min(ax, bx) > bb[4 * k + 1] + m2) continue;
        if (Math.max(ay, by) < bb[4 * k + 2] - m2 || Math.min(ay, by) > bb[4 * k + 3] + m2) continue;
        const cx = Math.min(GX - 1, Math.max(0, Math.floor((0.5 * (ax + bx) - gx0) / C)));
        const cy = Math.min(GY - 1, Math.max(0, Math.floor((0.5 * (ay + by) - gy0) / C)));
        cellOf[q] = cy * GX + cx; cnt[cellOf[q] + 1]++;
      }
      for (let c = 0; c < GX * GY; c++) cnt[c + 1] += cnt[c];
      const fill = cnt.slice(0, GX * GY);
      for (let q = 0; q < M; q++) if (cellOf[q] >= 0) cand[fill[cellOf[q]]++] = q;
      for (let e = t.estart[k]; e < t.estart[k + 1]; e++) {
        const a = t.eparent[e], c = t.echild[e];
        if (k > 1) { let any = false; for (let b = 0; b < NM; b++) if (mult[a * NM + b]) { any = true; break; } if (!any) continue; }
        const p0x = t.ox[a], p0y = t.oy[a], p1x = t.ox[c], p1y = t.oy[c];
        const ex = Math.floor((0.5 * (p0x + p1x) - gx0) / C), ey = Math.floor((0.5 * (p0y + p1y) - gy0) / C);
        let hit = false;
        for (let cy = Math.max(0, ey - 1); cy <= Math.min(GY - 1, ey + 1) && !hit; cy++)
          for (let cx = Math.max(0, ex - 1); cx <= Math.min(GX - 1, ex + 1) && !hit; cx++) {
            const cc = cy * GX + cx;
            for (let z = cnt[cc]; z < cnt[cc + 1]; z++) {
              const q = cand[z];
              const r1x = p1x - rx[k * M + q], r1y = p1y - ry[k * M + q];
              if (r1x * r1x + r1y * r1y < D2) { hit = true; break; }
              const r0x = p0x - rx[(k - 1) * M + q], r0y = p0y - ry[(k - 1) * M + q];
              if (r0x * r0x + r0y * r0y < D2) continue;
              const ddx = r1x - r0x, ddy = r1y - r0y, dd = ddx * ddx + ddy * ddy;
              if (dd <= 0) continue;
              const s = -(r0x * ddx + r0y * ddy) / dd;
              if (s > 0 && s < 1) { const mx = r0x + s * ddx, my = r0y + s * ddy; if (mx * mx + my * my < D2) { hit = true; break; } }
            }
          }
        if (hit) continue;
        if (k === 1) mult[c * NM + t.emove[e]] += 1;
        else for (let b = 0; b < NM; b++) mult[c * NM + b] += mult[a * NM + b];
      }
    }
    return { rx, ry, M, mult };
  }

  // Paint positions [p0, p1) and scatter their units' edges into a difference array D[NM][B]; also returns the
  // part of S_b at body bin 0 carried by these units (base).
  function paintRange(p, t, prep, p0, p1) {
    const B = t.B, W = B / 32, R = p.R, asin = p.half_angle === "asin", invbw = B / TWO_PI;
    const { rx, ry, M, mult } = prep;
    const D = new Int32Array(NM * B), base = new Float64Array(NM), mask = new Uint32Array(W);
    const um = new Int32Array(NM), edges = new Int32Array(2 * B);
    for (let pi = p0; pi < p1; pi++) {
      const u0 = t.posU0[pi], u1 = t.posU0[pi + 1];
      let reach = false;
      for (let u = u0; u < u1 && !reach; u++)
        for (let z = t.unitN0[u]; z < t.unitN0[u + 1] && !reach; z++)
          for (let b = 0; b < NM; b++) if (mult[t.unitNodes[z] * NM + b]) { reach = true; break; }
      if (!reach) continue;
      const X = t.posX[pi], Y = t.posY[pi], o = t.posLevel[pi] * M;
      mask.fill(0);
      for (let q = 0; q < M; q++) {
        const dx = rx[o + q] - X, dy = ry[o + q] - Y, d2 = dx * dx + dy * dy;
        if (d2 <= 0) continue;
        let half;
        if (asin) { if (d2 <= R * R) { mask.fill(0xffffffff); break; } half = Math.asin(R / Math.sqrt(d2)); }
        else half = halfAtan(R / Math.sqrt(d2));
        const c = fastAtan2(dy, dx) * invbw + B, h = half * invbw;
        const k0 = Math.round(c - h), k1 = Math.round(c + h - 1);
        if (k1 < k0) continue;
        for (let kb = k0; kb <= k1;) {
          const kk = kb % B, w = kk >>> 5, bit = kk & 31, run = Math.min(32 - bit, k1 - kb + 1);
          mask[w] |= run === 32 ? 0xffffffff : (((1 << run) - 1) << bit);
          kb += run;
        }
      }
      let ne = 0;  // edges: entry = bin * 2 + (1 if an interval ends here)
      for (let w = 0; w < W; w++) {
        const m = mask[w], prev = mask[w === 0 ? W - 1 : w - 1];
        let ed = (m ^ ((m << 1) | (prev >>> 31))) >>> 0;
        while (ed) {
          const bit = 31 - Math.clz32(ed & -ed);
          ed = (ed & (ed - 1)) >>> 0;
          edges[ne++] = ((w * 32 + bit) << 1) | (((m >>> bit) & 1) ^ 1);
        }
      }
      for (let u = u0; u < u1; u++) {
        um.fill(0);
        let any = false;
        for (let z = t.unitN0[u]; z < t.unitN0[u + 1]; z++)
          for (let b = 0; b < NM; b++) { const v = mult[t.unitNodes[z] * NM + b]; um[b] += v; if (v) any = true; }
        if (!any) continue;
        const sh = t.unitSh[u];
        if ((mask[sh >>> 5] >>> (sh & 31)) & 1) for (let b = 0; b < NM; b++) base[b] += um[b];
        for (let e = 0; e < ne; e++) {
          let bin = (edges[e] >> 1) - sh;
          if (bin < 0) bin += B;
          if (bin === 0) continue;
          const sgn = edges[e] & 1 ? -1 : 1;
          for (let b = 0; b < NM; b++) if (um[b]) D[b * B + bin] += sgn * um[b];
        }
      }
    }
    return { D, base };
  }

  // Sum the slices, then S_b = base_b + prefix sum of D, W_b = (1/(pi n_b)) sum S (n_b - S) dtheta, and the
  // per-branch occupancy p_b(theta) = S_b / n_b averaged onto `bins` display bins. Returns what the viewer draws.
  function combineFast(p, t, frame, i, prep, parts, bins) {
    bins = bins || 384;
    const B = t.B, T = t.tau, nn = t.nodeCount, mult = prep.mult;
    const D = new Float64Array(NM * B), base = new Float64Array(NM), nb = new Float64Array(NM);
    for (const pt of parts) { for (let q = 0; q < NM * B; q++) D[q] += pt.D[q]; for (let b = 0; b < NM; b++) base[b] += pt.base[b]; }
    const bits = new Uint8Array(nn);
    for (let n = 1; n < nn; n++) {
      let bt = 0;
      for (let b = 0; b < NM; b++) { const v = mult[n * NM + b]; if (v > 0) { bt |= 1 << b; nb[b] += v; } }
      bits[n] = bt;
    }
    const score = new Float64Array(NM), pb = new Float64Array(NM * bins), per = B / bins;
    for (let b = 0; b < NM; b++) {
      let S = base[b], acc = 0;
      for (let q = 0; q < B; q++) {
        if (q > 0) S += D[b * B + q];
        acc += S * (nb[b] - S);
        if (nb[b] > 0) pb[b * bins + Math.floor(q / per)] += S / nb[b] / per;
      }
      score[b] = nb[b] > 0 ? (acc * (2 / B)) / nb[b] : -1;
    }
    const x0 = frame.x[i], y0 = frame.y[i], c0 = Math.cos(frame.th[i]), s0 = Math.sin(frame.th[i]);
    const px = new Float64Array(nn), py = new Float64Array(nn);
    for (let n = 0; n < nn; n++) { px[n] = x0 + c0 * t.ox[n] - s0 * t.oy[n]; py[n] = y0 + s0 * t.ox[n] + c0 * t.oy[n]; }
    return { px, py, bits, level: t.level, p: pb, score, nb, bins };
  }

  // Single-threaded convenience (Node checks, fallback)
  function evaluateFast(p, t, frame, i, bins) {
    const prep = prepareAgent(p, t, frame, i);
    return combineFast(p, t, frame, i, prep, [paintRange(p, t, prep, 0, t.P)], bins);
  }

  const api = { buildTemplate, evaluateAgent, visualField, analyseRun, replay2d, strideMoves, packAnalysis, unpackAnalysis, NM, FAST_BINS, buildFastTemplate, prepareAgent,
    paintRange, combineFast, evaluateFast };
  // Source that defines FSMCore in a worker (the viewer builds its workers from a Blob, which also works from file://).
  api.workerSource = () => "(" + FSMCoreModule.toString() + ")(self);";
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.FSMCore = api;
  return api;
}
FSMCoreModule(typeof self !== "undefined" ? self : this);
