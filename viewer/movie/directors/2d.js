// Movie director for the 2D viewer (index.html), run by movie.py after the shot file, which sets window.__SHOT:
//   speed  [[fromStep, stepsPerSecond], ...]  playback speed from each step on, ramped over `ramp` seconds (1.5)
//   view   [[step, zoom, track], ...]         zoom (1: the whole flock) and track (0..1: glide the view centre from the
//                                             flock onto the tracked bird), interpolated smoothly between keys
//   track  {at, pct, dir}                     the tracked bird: at step `at`, the one at radius percentile `pct`
//                                             from the centre, nearest the direction `dir`
//   badge  true                               show the current playback speed (steps/s); no other text is drawn
//   set    {trail, colour, camera, turn}      viewer settings (defaults: 10-step trails, whole-swarm camera)
//   hold   2                                  seconds on the last frame
(() => {
  const C = window.__SHOT, FPS = C.fps || 30, RAMP = C.ramp || 1.5, HOLD = C.hold === undefined ? 2 : C.hold;
  const mono = monotoneCubic;
  const V = C.view, vx = V.map((k) => k[0]);
  const fLZ = mono(vx, V.map((k) => Math.log(k[1]))), fTr = mono(vx, V.map((k) => k[2] || 0));
  const target = (t) => { let v = C.speed[0][1]; for (const [s, sp] of C.speed) if (t >= s) v = sp; return v; };
  let bird = -1, E = 1;
  const times = [], speeds = [];
  window.__D = {
    ready: () => !!(window.FSM2D && FSM2D.info.ready),
    setup() {
      const st = document.createElement("style");
      st.textContent = `.transport,.hint,.runpick,.masthead,#gridlabel{display:none!important}
        #spd{position:absolute;left:32px;bottom:28px;z-index:20;display:flex;align-items:baseline;gap:10px;padding:10px 18px 11px;border-radius:10px;
          background:rgba(9,13,26,.78);border:1px solid rgba(131,144,173,.28);font:500 38px "IBM Plex Mono",ui-monospace,monospace;color:#eef4ff;font-variant-numeric:tabular-nums}
        #spd .u{font-size:21px;color:#9fb0d0;letter-spacing:.04em} #spd .i{font-size:24px;color:#4cc9f0} #spd[hidden]{display:none!important}`;
      document.head.appendChild(st);
      const spd = document.createElement("div");
      spd.id = "spd";
      spd.hidden = C.badge === false;
      document.body.appendChild(spd);
      E = FSM2D.info.every || 1;  // steps per stored frame (large runs keep one frame in several)
      FSM2D.set(Object.assign({ camera: "swarm", turn: false, trail: 10 }, C.set || {}));
      if (C.track) {
        const T0 = C.track.at / E, c = FSM2D.cam(T0), N = FSM2D.info.N, rs = [];
        for (let j = 0; j < N; j++) { const p = FSM2D.pos(T0, j); rs.push([Math.hypot(p[0] - c.x, p[1] - c.y), j, p]); }
        rs.sort((a, b) => a[0] - b[0]);
        const r0 = rs[Math.floor(C.track.pct * (N - 1))][0], [ux, uy] = C.track.dir || [1, 0];
        let best = -Infinity;
        for (const [r, j, p] of rs)
          if (Math.abs(r - r0) < 0.08 * r0) { const sc = ((p[0] - c.x) * ux + (p[1] - c.y) * uy) / r; if (sc > best) { best = sc; bird = j; } }
      }
      // integrate the speed schedule: the speed moves towards the current target at |dv| / RAMP per second
      const T = (FSM2D.info.T - 1) * E, maxRate = Math.max(...C.speed.map((s) => s[1])) / RAMP;
      let t = C.speed[0][0], v = target(t);
      while (t < T) {
        times.push(t); speeds.push(v);
        const dv = target(t) - v;
        v += Math.sign(dv) * Math.min(Math.abs(dv), maxRate / FPS);
        t += v / FPS;
      }
      for (let i = 0; i < HOLD * FPS; i++) { times.push(T); speeds.push(0); }
      this.n = times.length;
      return { T: FSM2D.info.T, N: FSM2D.info.N, every: E, bird, seconds: +(this.n / FPS).toFixed(1) };
    },
    n: 0,
    frame(i) {
      const t = times[i], w = fTr(t);
      let u = 0, v = 0;
      if (bird >= 0 && w > 0) { const c = FSM2D.cam(t / E), p = FSM2D.pos(t / E, bird); u = w * (p[0] - c.x); v = w * (p[1] - c.y); }
      FSM2D.view(Math.exp(fLZ(t)), u, v);
      FSM2D.seek(t / E);
      document.getElementById("spd").innerHTML = speeds[i] > 0
        ? `<span class="i">▶</span>${Math.round(speeds[i])}<span class="u">steps/s</span>`
        : `<span class="i">■</span>end<span class="u">step ${(FSM2D.info.T - 1) * E}</span>`;
    },
  };

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
