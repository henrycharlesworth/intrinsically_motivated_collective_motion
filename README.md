# Code for "Intrinsically Motivated Collective Motion"
This is example code for the collective motion model I've been working on during my PhD. Running it as it is will run using the default parameters for 1000 timesteps.

## Example Result (click image below to watch video):
<a href="http://www.henrycharlesworth.com/fileStorage/N=50_treeSearch.mp4"><img src="http://www.henrycharlesworth.com/fileStorage/typicalResult.png" width="400"/></a>

## Visualizing the Model:
I made some apps in D3.js to help understand what's going on.
### Visual State Definition (click image below to open app):
Each agent (or "bird") has a visual state which is defined in terms of the <i>projection</i> of the other birds onto its angular field of view. This is then divided up into a discrete number if sensors which read either a 1 or a 0 depending on whether or not they are more than half filled.

<a href="https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/visState.html"><img src="https://www.henrycharlesworth.com/fileStorage/visStatePreview.png" width="500" /> </a>
### Future State Maximization (click image below to open app):
Each bird has five decisions available to it at each time step. It can continue to move in its current direction at v0, v0+dv or v0-dv, or it can move at v0 but change its orientation by +/- dTheta. To make this decision the agent models a <i>tree of possible futures</i> for each move that it currently has available to it whilst (by default) modelling the other agents as moving in a straight line (although we also experiment with different future models). Each of these leads to a branch of different possible future states that each has a visual state (a vector of 1s and 0s) associated with it and the agent chooses the move which leads to the <i>largest number of unique future visual states</i>. Each agent chooses how to update itself in exactly the same way.

<a href="https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/treeSearch.html"><img src="https://www.henrycharlesworth.com/fileStorage/treeSearchPreview.png" width="650" /> </a>

## Training a Neural Network to Mimic This Model
<a href="https://www.henrycharlesworth.com/fileStorage/sideBySide_fullAndDeepHeuristic2.mp4"> <img src="https://www.henrycharlesworth.com/fileStorage/heuristicvsfull.png" width="800" /> </a>
