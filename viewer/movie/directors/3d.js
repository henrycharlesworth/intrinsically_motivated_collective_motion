// Movie director for the 3D viewer (3d.html), run by movie.py after the shot file, which sets window.__SHOT:
//   speed   [[fromStep, stepsPerSecond], ...]   playback speed from each step on, ramped over `ramp` seconds (1.5)
//   view    [[step, azimuth, elevation, zoom, focus], ...]   the camera, interpolated smoothly between keys:
//             azimuth  radians around the flock, relative to "behind it" (pi/2: beside it, pi: in front); the
//                      reference follows the flock's heading with a lag
//             elevation radians above the horizontal; zoom 1 frames the flock
//             focus    (optional; needs `groups`) -1: aim at group A, 0: the whole swarm, +1: group B
//   groups  {at: step}   split the agents into two groups along the swarm's long axis at that step (for flocks that
//                        break in two and drift apart, so the camera can follow one half)
//   set     {colour, trail, size, predScale}   viewer settings
//   text    false        hide the title and the numbers of the HUD (the minimap and graphs stay)
//   badge   false        show the current playback speed
//   hold    0            seconds on the last frame
(() => {
  const C = window.__SHOT, FPS = C.fps || 30, RAMP = C.ramp || 1.5, HOLD = C.hold || 0, LAG = C.headingLag || 0.04;
  const K = C.view, ks = K.map((k) => k[0]);
  const fOff = monotoneCubic(ks, K.map((k) => k[1])), fEl = monotoneCubic(ks, K.map((k) => k[2]));
  const fZ = monotoneCubic(ks, K.map((k) => k[3])), fF = monotoneCubic(ks, K.map((k) => k[4] || 0));
  const target = (t) => { let v = C.speed[0][1]; for (const [s, sp] of C.speed) if (t >= s) v = sp; return v; };
  const times = [], speeds = [];
  let hdg = null, last = -2, grp = null;
  window.__D = {
    ready: () => !!(window.FSM3D && FSM3D.info.T > 0),
    setup() {
      const hp = document.getElementById("hidepanel");
      if (hp) hp.click();
      const st = document.createElement("style");
      st.textContent = `.transport,.hint,.runpick{display:none!important}
        #spd{position:absolute;left:32px;bottom:28px;z-index:20;display:flex;align-items:baseline;gap:10px;padding:10px 18px 11px;border-radius:10px;
          background:rgba(9,13,26,.78);border:1px solid rgba(131,144,173,.28);font:500 38px "IBM Plex Mono",ui-monospace,monospace;color:#eef4ff}
        #spd .u{font-size:21px;color:#9fb0d0} #spd .i{font-size:24px;color:#4cc9f0} #spd[hidden]{display:none!important}` +
        (C.text === false ? ".masthead,#hud>span,#hud>b{display:none!important}" : "");
      document.head.appendChild(st);
      const spd = document.createElement("div");
      spd.id = "spd";
      spd.hidden = !C.badge;
      document.body.appendChild(spd);
      if (C.set) FSM3D.set(C.set);
      const T = FSM3D.info.T - 1, maxRate = Math.max(...C.speed.map((s) => s[1])) / RAMP;
      let t = C.speed[0][0], v = target(t);
      while (t < T) {
        times.push(t); speeds.push(v);
        const dv = target(t) - v;
        v += Math.sign(dv) * Math.min(Math.abs(dv), maxRate / FPS);
        t += v / FPS;
      }
      times.push(T); speeds.push(0);
      for (let i = 0; i < HOLD * FPS; i++) { times.push(T); speeds.push(0); }
      this.n = times.length;
      return { T: FSM3D.info.T, N: FSM3D.info.N, seconds: +(this.n / FPS).toFixed(1) };
    },
    n: 0,
    frame(i) {
      const t = times[i];
      const d = FSM3D.dir(t), h = Math.atan2(-d[1], -d[0]);  // behind the mean heading
      if (hdg === null || i !== last + 1) hdg = h;
      else { const dh = Math.atan2(Math.sin(h - hdg), Math.cos(h - hdg)); hdg += dh * LAG; }
      FSM3D.view(hdg + fOff(t), fEl(t), fZ(t));
      if (C.groups) {
        if (!grp) grp = splitGroups(C.groups.at);
        const [cA, rA] = groupStats(t, false), [cB, rB] = groupStats(t, true), f = fF(t);
        const all = [0, 1, 2].map((q) => (cA[q] + cB[q]) / 2);
        const dAll = 8 * (rA + rB) / 2 + 0.8 * Math.hypot(cA[0] - cB[0], cA[1] - cB[1], cA[2] - cB[2]);
        const [cG, rG] = f < 0 ? [cA, rA] : [cB, rB], w = Math.abs(f);
        FSM3D.aim([0, 1, 2].map((q) => all[q] + (cG[q] - all[q]) * w).concat([dAll + (8 * rG - dAll) * w]));
      }
      FSM3D.seek(t, i !== last + 1);  // snap the camera on a jump (stills)
      last = i;
      document.getElementById("spd").innerHTML = speeds[i] > 0 ? `<span class="i">▶</span>${Math.round(speeds[i])}<span class="u">steps/s</span>` : "";
    },
  };

  function splitGroups(t0) {  // true / false per agent: which side of the swarm's principal plane at step t0
    const N = FSM3D.info.N, P = [];
    for (let j = 0; j < N; j++) P.push(FSM3D.pos(t0, j));
    const c = [0, 1, 2].map((q) => P.reduce((a, p) => a + p[q], 0) / N);
    let v = [1, 0.3, 0.1];
    for (let it = 0; it < 30; it++) {  // power iteration on the covariance
      const w = [0, 0, 0];
      for (const p of P) {
        const e = [p[0] - c[0], p[1] - c[1], p[2] - c[2]], dv = e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
        for (let q = 0; q < 3; q++) w[q] += dv * e[q];
      }
      const n = Math.hypot(...w);
      v = w.map((x) => x / n);
    }
    return P.map((p) => (p[0] - c[0]) * v[0] + (p[1] - c[1]) * v[1] + (p[2] - c[2]) * v[2] > 0);
  }
  function groupStats(t, side) {  // centre and rms radius of one group
    const N = FSM3D.info.N, c = [0, 0, 0], P = [];
    for (let j = 0; j < N; j++) if (grp[j] === side) { const p = FSM3D.pos(t, j); P.push(p); for (let q = 0; q < 3; q++) c[q] += p[q]; }
    for (let q = 0; q < 3; q++) c[q] /= P.length;
    let r = 0;
    for (const p of P) r += (p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 + (p[2] - c[2]) ** 2;
    return [c, Math.sqrt(r / P.length)];
  }
  // Monotone cubic (Fritsch-Carlson) interpolation through (xs, ys): smooth, and never overshoots the keys.
  function monotoneCubic(xs, ys) {
    const m = xs.length, d = [], s = new Array(m).fill(0);
    for (let i = 0; i < m - 1; i++) d.push((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]));
    s[0] = d[0]; s[m - 1] = d[m - 2];
    for (let i = 1; i < m - 1; i++)
      if (d[i - 1] * d[i] > 0) {
        const w1 = 2 * (xs[i + 1] - xs[i]) + (xs[i] - xs[i - 1]), w2 = (xs[i + 1] - xs[i]) + 2 * (xs[i] - xs[i - 1]);
        s[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]);
      }
    return (x) => {
      if (m === 1 || x <= xs[0]) return ys[0];
      if (x >= xs[m - 1]) return ys[m - 1];
      let i = 0;
      while (i < m - 2 && x > xs[i + 1]) i++;
      const h = xs[i + 1] - xs[i], u = (x - xs[i]) / h;
      return (2 * u ** 3 - 3 * u ** 2 + 1) * ys[i] + (u ** 3 - 2 * u ** 2 + u) * h * s[i] + (-2 * u ** 3 + 3 * u ** 2) * ys[i + 1] + (u ** 3 - u ** 2) * h * s[i + 1];
    };
  }
})();
