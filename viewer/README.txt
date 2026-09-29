FSM flock viewer
================

Runs of the Future State Maximisation (FSM) flocking model, replayed in a web browser. No installation, server or
GPU is needed.

Open
----
Unzip the folder anywhere (keep its contents together), then double-click

    index.html   the 2D runs
    3d.html      the 3D runs, including the predator experiments

Use a recent Chrome, Edge, Firefox or Safari. The "Experiment" menu at the top left lists every run, 2D and 3D; picking
a run of the other kind switches page. A run is stored as its starting positions plus each bird's moves, and the
viewer rebuilds the trajectories when it loads: small runs open at once, the largest take a few seconds.

Controls (both pages)
---------------------
    Space           play / pause
    <- ->           one step back / forward (Shift: 10 steps)
    [ ]             slower / faster
    time slider     jump to any step
    scroll          zoom
    click a bird    select it and show it in the side panel

2D (index.html)
    drag            pan
    Bird view       the selected bird's futures: for each of its possible first moves, the score W that FSM
                    maximises, and a ring around the bird showing which directions of its view vary between its
                    futures (the panel explains the picture). Pause first: it is computed in the browser.
    camera          follow the largest flock or the whole swarm
    trail, colour   length of the tails; colour by heading or by speed

3D (3d.html)
    drag            orbit the camera
    F / Follow      the camera rides behind the flock
    Orbit           slow automatic orbit
    Shape           the flock's principal-axes ellipsoid
    Slice           (predator runs) show only a thin slab through the predator, to see the hole it bores through
                    the flock
    Bird view       the selected bird's whole visual sphere, unrolled
    colour          by move, elevation or depth; in runs where birds must see a predator to react, by whether each
                    bird sees it, remembers it, or does not know about it

The runs
--------
Each menu entry gives the number of birds N, the lookahead tau (how many steps ahead each bird plans), the number of
steps, and for the predator runs what happens. The largest runs (tens of thousands of birds) need a few GB of memory
and play every second or fourth step.
