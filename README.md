# Intrinsically motivated collective motion

Simulation code for Charlesworth & Turner, *Intrinsically motivated collective motion*, PNAS 116, 15362 (2019)
([doi:10.1073/pnas.1822069116](https://doi.org/10.1073/pnas.1822069116)), rewritten in 2026 for large flocks, long
horizons and 3D.

Each agent ("bird") looks τ steps ahead: it imagines every sequence of moves it could make, assumes the others fly
straight on, discards the futures that collide, and takes the first move whose surviving futures have the most
varied visual states. This is Future State Maximisation (FSM); flocking emerges from it alone.

| | |
|---|---|
| `src/fsm2d/` | 2D simulator: exact CPU version (`bin/fsm2d`, OpenMP) and CUDA version (`bin/fsm2d_gpu`) |
| `src/fsm3d/` | 3D simulator, CUDA (`bin/fsm3d_gpu`), with optional predators |
| `src/common/` | Command line, random streams, run files and resume, shared by the simulators |
| `python/fsmtools/` | Loading runs, flock measures, brute-force reference models (2D and 3D), viewer export |
| `scripts/` | Parameter sweeps on many GPUs, run summaries, the empowerment check, viewer bundles |
| `viewer/` | Browser viewer (2D and 3D): a local server for a folder of runs, or standalone bundles; and a video tool |
| `tests/` | Test suite: exactness, GPU against CPU and reference, reproducibility |
| [`docs/NOTES.md`](docs/NOTES.md) | The model, the derivation of the score and its link to empowerment, 3D, and results |
| [`docs/OPTIMIZATIONS.md`](docs/OPTIMIZATIONS.md) | How the code is made fast, from the original implementation to multi-GPU |

## Build

```sh
make                          # bin/fsm2d: needs a C++17 compiler with OpenMP (g++ or clang)
make gpu CUDA_ARCH=sm_89      # bin/fsm2d_gpu, bin/fsm3d_gpu: needs the CUDA toolkit; sm_89 = RTX 40xx / L4,
                              #   sm_90 = H100, sm_120 = RTX 50xx / RTX PRO 6000 (default)
pip install -r requirements.txt && pip install -e .   # fsmtools (python/) and the tools' dependencies
make test                     # builds what it can and runs tests/ (the GPU tests skip without a GPU)
```

3D runs need an NVIDIA GPU; on a CPU, `python -m fsmtools.reference3d simulate --N 60 --tau 3` runs tiny 3D flocks
(the brute-force reference: N ≲ 100, τ ≤ 3).

## Run

```sh
# 2D, N = 500 birds, lookahead tau = 6, 1000 steps from a random start in a 250 x 250 box
bin/fsm2d     --N 500 --tau 6 --box 250 --steps 1000 --out runs/n500_t6        # CPU, exact (~1 s/step on 24 cores)
bin/fsm2d_gpu --N 500 --tau 6 --box 250 --steps 1000 --out runs/n500_t6        # GPU (~20 ms/step)

# 3D, N = 1000, tau = 6, on a 192 x 96 visual-field grid, merging futures within 1.5 of each other
bin/fsm3d_gpu --N 1000 --tau 6 --box 110 --grid 192x96 --merge 1.5 --steps 2000 --out runs/n1000_t6_3d
# ... with two predators (one from behind, one diving from above), from step 100
bin/fsm3d_gpu --N 1000 --tau 6 --box 110 --grid 192x96 --merge 1.5 --steps 500 --out runs/pred \
    --predator com --pred-start 100 --pred-dist 150 --pred-speed 20 --pred-grow 4 --pred2 intercept,0,70

# on several GPUs (results are identical to one GPU), and continuing or extending an interrupted run
bin/fsm2d_gpu --N 50000 --tau 6 --box 2500 --steps 10000 --gpus all --out runs/n50000_t6
bin/fsm2d_gpu --N 50000 --tau 6 --box 2500 --steps 20000 --gpus all --out runs/n50000_t6 --resume

# a sweep over N and tau, several runs at a time on 8 GPUs (resumable)
python3 scripts/run_matrix.py --dims 2 --N 500,1000,2500,5000,10000 --tau 6,7,8,9,10 --gpus 0-7 --out runs/matrix
python3 scripts/summarize.py runs/matrix/n*_t6
```

Main flags (all simulators): `--N --tau --steps --seed --box --heading-sd --init FILE --out PREFIX --resume`, and the
model parameters `--v0 --dv --dtheta --R --dt --half-angle atan|asin` (2D also `--noise`; 3D `--phi-max`). GPU:
`--device D`, `--gpus 0,1,..|all`; 2D `--bins` (angular resolution) and `--merge EPS` (automatic: exact for τ ≤ 7,
1.5 above); 3D `--grid NLONxNLAT` and `--merge EPS`. Predators (3D): `--predator line|com|nearest|intercept`,
`--pred-start --pred-dist --pred-from AZ,EL --pred-speed --pred-radius --pred-grow --pred-turn --pred2 MODE,AZ,EL`,
and `--pred-vision on` (only birds that see a predator react; `--pred-blind --pred-memory`). The header of each
`main` file lists every flag; `docs/NOTES.md` explains the model behind them.

**Output**, written step by step (a run can be watched or resumed while it goes; `fsmtools.load_run` reads it):
`PREFIX.json` (parameters), `PREFIX.traj` (float64 frames: x, y, θ, v in 2D; x, y, z, θ, φ, v in 3D), `PREFIX.moves`
(int8, the move each bird took), `PREFIX.scores` (float32, each bird's score for every first move), and for
predators `PREFIX.pred` and `PREFIX.seen`. All randomness is indexed by (seed, step, bird), so runs are
reproducible bit for bit, whatever the number of threads or GPUs, and across interruptions.

## View

Save runs anywhere under `runs/` (as in the examples above: `--out runs/<name>`), then start the viewer:

```sh
python viewer/serve.py                  # http://localhost:8000; needs only Python 3, no GPU
python viewer/serve.py --runs ~/results --port 8080   # another folder of runs
```

The **Experiment** menu of `index.html` (2D) and `3d.html` (3D) lists every run in the folder, including its
subfolders, by path. A run that is still going can be watched: reload the page to pick up the new steps. Space
plays, ← → step, [ ] change speed, drag pans (2D) or orbits (3D). Click a bird to select it. In 2D, **Bird view**
shows that bird's futures and the score of each first move, recomputed in the browser. The 3D controls are
**Follow**, **Shape** (the flock's ellipsoid) and, for predator runs, **Slice** and colouring by who has seen the
predator. `viewer/README.txt` lists all the controls.

The viewer also works with no server at all, opened straight from disk (`scripts/build_viewer.py` bundles chosen
runs, and can zip them with the viewer to send), and `viewer/movie/movie.py` renders videos from it.

## Original release (2019)

The 2019 code (discrete visual states, `main.cpp`) is in the git history of `master`. Interactive explanations of the
[visual state](https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/visState.html) and
the [future-state tree search](https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/treeSearch.html).
