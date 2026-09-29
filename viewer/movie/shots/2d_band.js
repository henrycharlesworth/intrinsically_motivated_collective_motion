// 2D, 5000 steps (the N = 10000, tau = 10 video): 15 steps/s through the opening and a close-up of the long band the
// flock turns into, 50 steps/s otherwise.
window.__SHOT = {
  speed: [[0, 15], [400, 50], [1950, 15], [2230, 50]],
  track: { at: 2050, pct: 0.12, dir: [1, 1] },
  view: [[0, 1, 0], [1900, 1, 0], [1990, 1.6, 0.4], [2030, 3.5, 1], [2090, 4, 1], [2160, 1.3, 0.2], [2220, 1, 0], [4700, 1, 0], [5000, 0.9, 0]],
};
