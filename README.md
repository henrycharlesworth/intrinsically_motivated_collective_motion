# Code for "Intrinsically Motivated Collective Motion"
Update: paper has now been published in PNAS: https://www.pnas.org/content/early/2019/07/15/1822069116. Note that the code here is only for the discrete sensor version of the model, using the ballistic heuristic for the future trajectories of the other agents.

This is example code for the collective motion model I've been working on during my PhD. Running it as it is will run using the default parameters for 1000 timesteps.

## Example Results (click images below to watch video):
<a href="http://www.henrycharlesworth.com/fileStorage/N=50_treeSearch.mp4"><img src="http://www.henrycharlesworth.com/fileStorage/typicalResult.png" width="400"/></a>
<a href="https://www.henrycharlesworth.com/fileStorage/largeContinuousFlock.mp4"><img src="https://www.henrycharlesworth.com/fileStorage/largerFlockThumbnail.png" width="400"/></a>

## Visualizing the Model:
I made some apps in D3.js to help understand what's going on.
### Visual State Definition (click image below to open app):
Each agent (or "bird") has a visual state which is defined in terms of the <i>projection</i> of the other birds onto its angular field of view. This is then divided up into a discrete number if sensors which read either a 1 or a 0 depending on whether or not they are more than half filled.

<a href="https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/visState.html"><img src="https://www.henrycharlesworth.com/fileStorage/visStatePreview.png" width="500" /> </a>
### Future State Maximization (click image below to open app):
Each bird has five decisions available to it at each time step. It can continue to move in its current direction at v0, v0+dv or v0-dv, or it can move at v0 but change its orientation by +/- dTheta. To make this decision the agent models a <i>tree of possible futures</i> for each move that it currently has available to it whilst (by default) modelling the other agents as moving in a straight line (although we also experiment with different future models). Each of these leads to a branch of different possible future states that each has a visual state (a vector of 1s and 0s) associated with it and the agent chooses the move which leads to the <i>largest number of unique future visual states</i>. Each agent chooses how to update itself in exactly the same way.

<a href="https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/treeSearch.html"><img src="https://www.henrycharlesworth.com/fileStorage/treeSearchPreview.png" width="650" /> </a>

## Training a Neural Network to Mimic This Model
The full model requires modelling a large number of future trajectories for each bird at each timestep. This is clearly very computationally intense (especially the further into the future we look) and not at all realistic for a biological organism. As such, we decided to see if it would be possible to train a neural network to give qualitatively similar results to the full model but using only the <i>current visual state</i> (and it turns out to be necessary to include the previous state too). We do this by running the full model and recording the visual states when each decision is made, as well as recording each decision and using this as a target for the neural network to fit. A comparison with the full model is shown below:
<a href="https://www.henrycharlesworth.com/fileStorage/sideBySide_fullAndDeepHeuristic2.mp4"> <img src="https://www.henrycharlesworth.com/fileStorage/heuristicvsfull.png" width="800" /> </a>

and we see it works very nicely (albeit not quite perfectly)!
