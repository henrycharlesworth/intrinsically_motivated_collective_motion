// 3D, two predators (N = 1000, tau = 6, 500 steps, predators from step 100): one constant speed from the random start,
// the camera circling the flock through the splits and rejoins. No text, no trails.
window.__SHOT = {
  speed: [[0, 5]],
  view: [
    // step  azimuth elevation zoom
    [0,   3.0,  0.35, 1.1], [20,  2.65, 0.28, 1.15], [40,  2.3,  0.2,  1.25], [72,  1.6,  0.3,  0.95],
    [95,  0.7,  0.22, 1.05], [108, 0.25, 0.12, 1.2], [124, 0.75, 0.35, 0.95], [150, 1.55, 0.55, 0.8],
    [185, 2.25, 0.35, 0.85], [212, 2.85, 0.18, 0.9], [245, 3.45, 0.05, 0.95], [262, 3.9, -0.18, 0.9],
    [295, 4.4, -0.25, 1.0], [318, 4.9,  0.1,  0.95], [350, 5.4,  0.75, 1.15], [400, 5.9,  1.05, 1.2],
    [445, 6.4,  0.55, 1.0], [480, 6.9,  0.28, 0.85], [499, 7.3,  0.32, 1.35],
  ],
  set: { colour: "move", trail: 0, size: 1.5, predScale: 1.8 },
  text: false,
};
