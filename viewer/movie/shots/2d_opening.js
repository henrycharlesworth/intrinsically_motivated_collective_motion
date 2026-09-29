// 2D, 5000 steps from a random start (the N = 2500 and N = 10000, tau = 6 videos): 15 steps/s through the opening
// (square -> crescent -> disc -> starburst) and a close-up that follows one bird near the flock's edge, then 50 steps/s.
window.__SHOT = {
  speed: [[0, 15], [660, 50]],
  track: { at: 500, pct: 0.3, dir: [1, 0] },
  view: [[0, 1, 0], [440, 1, 0], [490, 3.5, 1], [550, 4, 1], [610, 1.3, 0.2], [660, 1, 0], [4700, 1, 0], [5000, 0.9, 0]],
};
