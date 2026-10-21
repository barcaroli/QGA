
<!-- NEWS.md is generated from NEWS.Rmd. Please edit NEWS.Rmd file -->

# QGA 1.1

- Introduced parallelization (optional, parameter `parallel`, default FALSE)

- Added the quantum interference crossover operator (parameters
  `interference_flag`, default FALSE, and `interference_rate`, default 0.1)

- Added `interference_type = "amplitude"` (with strength `interference_lambda`): fitness-weighted interference that sums the qubit amplitudes of the population, constructively for individuals above the average fitness and destructively for those below

- Added the vignette "Avoiding local optima", illustrating the parameters that counteract premature convergence and how to use them

- Added three options against premature convergence (all off by default): the H-epsilon gate (`epsilon`), which keeps the probability of each bit value in [epsilon, 1 - epsilon]; a partial restart towards equal probabilities after `restart_iters` generations without improvement (`restart_strength`); personal attractors with local and global migration (`attractor_type = "personal"`, `migration_group_size`, `migration_local_period`, `migration_global_period`)

- Added `restart_theta` (default FALSE): at each restart the rotation angle is brought back to `thetainit`; the vignette "Avoiding local optima" now gives guidelines for choosing `restart_iters` and `restart_strength`

- Added the function `parameters_hint()`, which suggests values for `thetainit`, `thetaend`, `restart_iters`, `restart_strength`, `attractor_type` and other parameters given the number of generations, the genome and the kind of problem

- Fixed a bug in the rotation gate: when an individual was at least as good as the best one, the case x = 1, b = 0 was never handled and the case x = 0, b = 1 was rotated twice in opposite directions

- Added two parameters to the QGA function:

  - stop_limit (value to stop the iterations if the fitness is higher
    than a given value)
  - stop_iters (number of iterations without improvement of fitness
    before stopping)
