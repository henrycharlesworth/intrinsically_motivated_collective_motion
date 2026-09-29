// The 2D viewer's per-frame analysis of a run (camera tracks, order, fragments, default bird: fsm-core.js
// analyseRun), computed ahead of time so that opening a large run in the viewer is fast (clustering every frame takes
// ~30 s at N = 10000, 5000 steps). Called by the exporter (python/fsmtools/viewer_export.py):
//   node viewer/tools/analyse.js DIR     # DIR holds init.f64, moves.u8, meta.json {params, T, N, stride}; writes analysis.f32
const fs = require("fs");
const path = require("path");
const FSM = require(path.join(__dirname, "..", "js", "fsm-core.js"));

const dir = process.argv[2];
const m = JSON.parse(fs.readFileSync(path.join(dir, "meta.json"), "utf8"));
const read = (f) => { const b = fs.readFileSync(path.join(dir, f)); return b.buffer.slice(b.byteOffset, b.byteOffset + b.length); };
const st = m.stride || 1, T = Math.ceil(m.T / st);
const traj = FSM.replay2d(m.params, new Float64Array(read("init.f64")), new Uint8Array(read("moves.u8")), m.T, m.N, st);
const out = FSM.packAnalysis(FSM.analyseRun(traj, T, m.N), T);
fs.writeFileSync(path.join(dir, "analysis.f32"), Buffer.from(out.buffer));
