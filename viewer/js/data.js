// Runs bundled with the viewer (written by scripts/build_viewer.py).
//
// data/runs.js calls FSMData.list([...]) with one entry per run: {id, dims, label}. data/<id>.js calls
// FSMData.register(id, {meta, arrays}), where each array is {type: "u8" | "f32" | "f64", deflate: bool,
// chunks: [base64 strings]}. They are plain scripts, not fetched files, so the viewer works when index.html is opened
// straight from disk, with no web server. A run's file is only loaded when the run is opened.
//
// A run is stored as its initial state and every agent's move at every step ("replay" format): the viewer
// re-integrates the trajectory, which is deterministic, so this is exact and ~40x smaller than the positions.
//
// Served by viewer/serve.py, the viewer also lists the runs of a folder (runs/index.json, entries with src:
// "server") and reads their simulator files directly: PREFIX.json, the first frame of PREFIX.traj, PREFIX.moves
// (+ .pred, .seen), into the same form as a bundled run.
(function (root) {
  "use strict";
  const runs = [], loaded = new Map(), waiting = new Map();
  const MAX_FRAMES_BYTES = 1.2e9;  // as the exporter: a 2D run's frames are held as float32 [T, N, 3]

  function list(entries) { runs.push(...entries); }

  function register(id, run) {
    loaded.set(id, run);
    const w = waiting.get(id);
    if (w) { waiting.delete(id); w(run); }
  }

  function loadScript(src) {
    return new Promise((resolve, reject) => {
      const s = document.createElement("script");
      s.src = src;
      s.onload = resolve;
      s.onerror = () => reject(new Error("cannot load " + src));
      document.head.appendChild(s);
    });
  }

  function b64ToBytes(txt) {
    const bin = atob(txt), out = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }

  async function decode(a) {
    const parts = a.chunks.map(b64ToBytes);
    let buf;
    if (a.deflate) {
      const stream = new Blob(parts).stream().pipeThrough(new DecompressionStream("deflate-raw"));
      buf = await new Response(stream).arrayBuffer();
    } else {
      buf = await new Blob(parts).arrayBuffer();
    }
    return a.type === "f64" ? new Float64Array(buf) : a.type === "f32" ? new Float32Array(buf) : new Uint8Array(buf);
  }

  // The server's runs, when there is one (resolves once the menu is complete).
  const ready = (typeof location !== "undefined" && location.protocol.startsWith("http")
    ? fetch("runs/index.json", { cache: "no-store" }).then((r) => (r.ok ? r.json() : [])).then((l) => { runs.push(...l); }, () => {})
    : Promise.resolve());

  async function fetchRunFile(base, ext, range) {
    const r = await fetch(base + "." + ext, { cache: "no-store", headers: range ? { Range: "bytes=" + range } : {} });
    if (!r.ok) throw new Error("cannot read " + base + "." + ext + " (" + r.status + ")");
    return r.arrayBuffer();
  }

  async function loadServer(e) {
    const base = "runs/" + e.id.split("/").map(encodeURIComponent).join("/"), has = (x) => e.files.includes(x);
    const params = JSON.parse(new TextDecoder().decode(await fetchRunFile(base, "json")));
    const N = params.N, dims = params.dims === 3 ? 3 : 2, w = dims === 3 ? 6 : 4, keep = dims === 3 ? 5 : 3;
    const head = await fetch(base + ".traj", { method: "HEAD", cache: "no-store" });
    let T = Math.floor(+head.headers.get("Content-Length") / (8 * N * w));
    const moves = has("moves") ? new Uint8Array(await fetchRunFile(base, "moves")) : null;
    if (moves) T = Math.min(T, Math.floor(moves.length / N) + 1);
    if (T < 2) throw new Error(e.id + " has no complete step yet");
    const meta = { params, T, N, dims }, out = { meta };
    if (!moves || (dims === 2 && params.noise > 0)) {  // the moves do not determine the trajectory: read it all
      const raw = new Float64Array(await fetchRunFile(base, "traj", "0-" + (T * N * w * 8 - 1)));
      const tr = new Float32Array(T * N * keep);
      for (let q = 0; q < T * N; q++) for (let c = 0; c < keep; c++) tr[keep * q + c] = raw[w * q + c];
      out.traj = tr;
      out.moves = moves ? moves.subarray(0, (T - 1) * N) : new Uint8Array((T - 1) * N).fill(dims === 3 ? 5 : 3);
    } else {
      out.init = new Float64Array(await fetchRunFile(base, "traj", "0-" + (N * w * 8 - 1)));
      out.moves = moves.subarray(0, (T - 1) * N);
      const stride = Math.ceil((T * N * 12) / MAX_FRAMES_BYTES);
      if (dims === 2 && stride > 1) meta.stride = stride;
    }
    if (dims === 3 && has("pred")) {
      const n = (params.predator && params.predator.n) || 1;
      const pr = new Float64Array(await fetchRunFile(base, "pred"));
      out.pred = Float32Array.from(pr.subarray(0, Math.min(pr.length, T * n * 7)));
      meta.pred = true; meta.pred_n = n;
    }
    if (dims === 3 && has("seen")) {
      out.seen = new Uint8Array(await fetchRunFile(base, "seen")).subarray(0, (T - 1) * N);
      meta.seen = true;
    }
    return out;
  }

  // {meta, init, moves, [traj, analysis, pred, seen]}: meta as written by the exporter, arrays decoded.
  async function load(id) {
    const e = runs.find((r) => r.id === id && r.src === "server");
    if (e) return loadServer(e);
    if (!loaded.has(id)) {
      const got = new Promise((resolve) => waiting.set(id, resolve));
      await loadScript("data/" + encodeURIComponent(id) + ".js");
      await got;
    }
    const run = loaded.get(id), out = { meta: run.meta };
    for (const [name, a] of Object.entries(run.arrays)) out[name] = await decode(a);
    return out;
  }

  root.FSMData = { list, register, load, ready, runs: () => runs.slice() };
})(typeof self !== "undefined" ? self : this);
