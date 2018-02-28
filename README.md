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

<a href="https://henrycharlesworth.com/IntrinsicallyMotivatedCollectiveMotion/Visualizations/treeSearch.html"><img src="https://www.henrycharlesworth.com/fileStorage/treeSearchPreview.png" width="650" /> </a>
