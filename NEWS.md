
<!-- NEWS.md is generated from NEWS.Rmd. Please edit NEWS.Rmd file -->

# QGA 1.1

- Introduced parallelization (optional, parameter `parallel`, default FALSE)

- Added the quantum interference crossover operator (parameters
  `interference_flag`, default FALSE, and `interference_rate`, default 0.1)

- Added `interference_type = "amplitude"` (with strength `interference_lambda`): fitness-weighted interference that sums the qubit amplitudes of the population, constructively for individuals above the average fitness and destructively for those below

- Fixed a bug in the rotation gate: when an individual was at least as good as the best one, the case x = 1, b = 0 was never handled and the case x = 0, b = 1 was rotated twice in opposite directions

- Added two parameters to the QGA function:

  - stop_limit (value to stop the iterations if the fitness is higher
    than a given value)
  - stop_iters (number of iterations without improvement of fitness
    before stopping)
