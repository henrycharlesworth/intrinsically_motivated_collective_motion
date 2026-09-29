// Runs bundled with the viewer (written by scripts/build_viewer.py).
//
// data/runs.js calls FSMData.list([...]) with one entry per run: {id, dims, label}. data/<id>.js calls
// FSMData.register(id, {meta, arrays}), where each array is {type: "u8" | "f32" | "f64", deflate: bool,
// chunks: [base64 strings]}. They are plain scripts, not fetched files, so the viewer works when index.html is opened
// straight from disk, with no web server. A run's file is only loaded when the run is opened.
//
// A run is stored as its initial state and every agent's move at every step ("replay" format): the viewer
// re-integrates the trajectory, which is deterministic, so this is exact and ~40x smaller than the positions.
(function (root) {
  "use strict";
  const runs = [], loaded = new Map(), waiting = new Map();

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

  // {meta, init, moves, [analysis, pred, seen]}: meta as written by the exporter, arrays decoded.
  async function load(id) {
    if (!loaded.has(id)) {
      const got = new Promise((resolve) => waiting.set(id, resolve));
      await loadScript("data/" + encodeURIComponent(id) + ".js");
      await got;
    }
    const run = loaded.get(id), out = { meta: run.meta };
    for (const [name, a] of Object.entries(run.arrays)) out[name] = await decode(a);
    return out;
  }

  root.FSMData = { list, register, load, runs: () => runs.slice() };
})(typeof self !== "undefined" ? self : this);
