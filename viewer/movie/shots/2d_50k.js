// 2D, N = 50000, 10000 steps: slow through the opening and two close-ups (the edge of the young disc, the inside of
// the mature flock), faster through the slowly changing later flock. No trails: at this density single birds read
// better as dots.
window.__SHOT = {
  speed: [[0, 15], [660, 50], [2500, 100], [5900, 15], [6200, 100]],
  track: { at: 500, pct: 0.3, dir: [1, 0] },
  view: [[0, 1, 0], [440, 1, 0], [490, 3.5, 1], [550, 4, 1], [610, 1.3, 0.2], [660, 1, 0],
         [5850, 1, 0], [5930, 4, 0], [6150, 4, 0], [6230, 1, 0], [9700, 1, 0], [10000, 0.9, 0]],
  set: { trail: 0 },
};
