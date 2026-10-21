#' Quantum Genetic Algorithm
#'
#' @description
#' Runs a Quantum Genetic Algorithm (QGA) to optimize a user-defined
#' fitness function over discrete decision variables. Internally it
#' initializes a quantum population, iterates measurement, rotation,
#' and optional mutation and interference, repairs invalid genes, evaluates
#' fitness (in parallel if `parallel = TRUE`), and tracks best/average fitness
#' across generations.
#'
#' @details
#' Solutions are encoded as bitstrings and decoded to integer values in
#' the range 1..`nvalues_sol` for each of the `Genome` genes. The user must
#' supply a fitness function `eval_fitness(solution_integers, eval_func_inputs)`
#' that returns a numeric score to maximize.
#'
#' @param popsize Integer. Population size.
#' @param generation_max Integer. Maximum number of generations to run.
#' @param nvalues_sol Integer. Number of allowed discrete values per gene (1..nvalues_sol).
#' @param Genome Integer. Number of genes (decision variables).
#' @param thetainit Numeric (radians). Initial rotation step for the update rule.
#' @param thetaend Numeric (radians). Final rotation step (linearly decayed to this value).
#' @param pop_mutation_rate_init Numeric in [0,1]. Initial per-individual mutation probability (default `1/(popsize+1)` when `mutation_flag`).
#' @param pop_mutation_rate_end Numeric in [0,1]. Final per-individual mutation probability (default `1/(popsize+1)` when `mutation_flag`).
#' @param mutation_rate_init Numeric in [0,1]. Initial per-bit mutation probability (default `1/(Genome+1)` when `mutation_flag`).
#' @param mutation_rate_end Numeric in [0,1]. Final per-bit mutation probability (default `1/(Genome+1)` when `mutation_flag`).
#' @param mutation_flag Logical. Whether to apply the mutation operator.
#' @param plotting Logical. If TRUE, plots fitness history over generations.
#' @param verbose Logical. If TRUE, prints per-generation summary to console.
#' @param progress Logical. If TRUE, shows a text progress bar.
#' @param eval_fitness Function. User-provided function with signature
#'   `function(solution_integers, eval_func_inputs)` returning a numeric fitness value.
#' @param eval_func_inputs Any. Second argument passed through to `eval_fitness`.
#' @param stop_limit Numeric. Early-stop threshold; stops when best fitness >= `stop_limit`.
#' @param stop_iters Integer or NULL. If set, stops when there is no improvement
#'   over `stop_iters` generations.
#' @param interference_flag Logical. Whether to apply the interference operator
#'   after mutation (see `interference_type`). Default FALSE.
#' @param interference_rate Numeric in [0,1]. Probability of applying the
#'   interference operator to the population at each iteration. Default 0.1.
#' @param interference_type Character. Type of interference operator:
#'   `"diagonal"` (default) copies the amplitudes of each qubit from another
#'   individual following a diagonal scheme (Narayanan & Moore, 1996), without
#'   using fitness; `"amplitude"` adds to each qubit the fitness-weighted sum
#'   of the amplitudes of the whole population, so that amplitudes of
#'   individuals above the average fitness interfere constructively and those
#'   below the average destructively. See `interference()`.
#' @param interference_lambda Numeric >= 0. Strength of the amplitude
#'   interference (used only with `interference_type = "amplitude"`). Default 0.5.
#' @param epsilon Numeric in \[0, 0.5). H-epsilon gate (Han & Kim, 2004): the
#'   probability of each bit value is kept in \[`epsilon`, 1 - `epsilon`\], so
#'   that qubits never collapse to a deterministic state. Default 0 (no gate).
#' @param restart_iters Integer or NULL. If set, when the best fitness does not
#'   improve for `restart_iters` generations, the population (except the best
#'   individual) is partially brought back towards equal probabilities.
#'   Should be smaller than `stop_iters`. Default NULL (no restart).
#' @param restart_strength Numeric in \[0, 1\]. Fraction of the way each qubit
#'   angle is moved towards pi/4 at a restart (1 = full reset). Default 0.5.
#' @param restart_theta Logical. If TRUE, at each restart the rotation angle is
#'   brought back to `thetainit` and decays again linearly to `thetaend` over
#'   the remaining generations. Default FALSE.
#' @param attractor_type Character. What each individual is rotated towards:
#'   `"generation"` (default), the best individual of the current generation;
#'   `"personal"`, its own attractor, i.e. the best solution it has found so
#'   far, shared through local and global migration (Han & Kim, 2002).
#' @param migration_group_size Integer. Size of the groups of consecutive
#'   individuals for local migration (`attractor_type = "personal"`). Default 2.
#' @param migration_local_period Integer or NULL. Every this many generations,
#'   the best attractor of each group is copied to the group. Default 1.
#' @param migration_global_period Integer or NULL. Every this many generations,
#'   the best attractor is copied to all individuals. Default 50.
#' @param parallel Logical. If TRUE, fitness evaluation is distributed over
#'   `detectCores() - 1` workers via `foreach`/`doParallel`. Default FALSE
#'   (sequential evaluation, no cluster is created).
#'
#' @return
#' A list with two elements:
#' - `best_solution_integers`: Integer vector of length `Genome` with values in 1..`nvalues_sol`.
#' - `fitness_history_df`: Data frame with columns `generation`, `fitness_average`, `fitness_best`.
#'
#' @examples
#' # Minimal example (toy fitness: sum of gene values)
#' set.seed(1)
#' f <- function(sol, data) sum(sol)
#' 
#' # This call is illustrative; tune popsize/generation_max for real problems
#' # out <- QGA(popsize = 10, generation_max = 5,
#' #            nvalues_sol = 4, Genome = 3,
#' #            eval_fitness = f, eval_func_inputs = NULL,
#' #            plotting = FALSE, verbose = FALSE, progress = FALSE)
#'
#' @seealso `rotation()`, `mutation()`, `measure()`, `repair()`
#' @export
QGA <- function(popsize = 20,  
                generation_max = 200,
                nvalues_sol,
                Genome,
                thetainit = 3.1415926535 * 0.05,
                thetaend  = 3.1415926535 * 0.025,
                pop_mutation_rate_init = NULL,
                pop_mutation_rate_end  = NULL,
                mutation_rate_init     = NULL,
                mutation_rate_end      = NULL,
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose  = TRUE,
                progress = TRUE,
                eval_fitness,
                eval_func_inputs,
                stop_limit = NULL,
                stop_iters = NULL,
                interference_flag = FALSE,
                interference_rate = 0.1,
                interference_type = c("diagonal", "amplitude"),
                interference_lambda = 0.5,
                epsilon = 0,
                restart_iters = NULL,
                restart_strength = 0.5,
                restart_theta = FALSE,
                attractor_type = c("generation", "personal"),
                migration_group_size = 2,
                migration_local_period = 1,
                migration_global_period = 50,
                parallel = FALSE) {
  
  # ---- Parallel deps ----
  if (parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
    if (!requireNamespace("foreach",    quietly = TRUE)) install.packages("foreach")
    library(doParallel)
    library(foreach)
  }
  
  # ---- Checks & defaults ----
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome))      stop("Genome parameter value missing!")
  if (is.null(pop_mutation_rate_init) && mutation_flag) pop_mutation_rate_init <- 1/(popsize+1)
  if (is.null(pop_mutation_rate_end)  && mutation_flag) pop_mutation_rate_end  <- 1/(popsize+1)
  if (is.null(mutation_rate_init)     && mutation_flag) mutation_rate_init     <- 1/(Genome+1)
  if (is.null(mutation_rate_end)      && mutation_flag) mutation_rate_end      <- 1/(Genome+1)
  if (is.null(stop_limit)) stop_limit <- Inf
  interference_type <- match.arg(interference_type)
  attractor_type    <- match.arg(attractor_type)
  if (epsilon < 0 || epsilon >= 0.5) stop("epsilon must be in [0, 0.5)")
  
  # ---- Genome encoding ----
  n <- 0L
  while (nvalues_sol > 2^n) n <- n + 1L
  geneLength   <- n
  genomeLength <- Genome * geneLength
  if (geneLength <= 0L || genomeLength <= 0L) stop("Invalid geneLength/genomeLength.")
  
  # ---- Helper: decode bitstring -> integer solution (1..nvalues_sol) ----
  decode_solution <- function(bits, geneLength, Genome) {
    out <- integer(Genome)
    idx <- 1L
    for (x in seq_len(Genome)) {
      acc <- 0L
      for (y in seq_len(geneLength)) {
        b <- bits[idx]
        if (b != 0L && b != 1L) b <- as.integer(b != 0)
        acc <- acc + b * 2L^(geneLength - y)
        idx <- idx + 1L
      }
      out[x] <- acc + 1L
    }
    out
  }
  
  # ---- Parallel backend (only if parallel = TRUE) ----
  if (parallel) {
    numCores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
    
    # Export una tantum del decoder ai worker (evita i warning)
    parallel::clusterExport(cl, varlist = c("decode_solution"),
                            envir = environment())
    
    # Fallback operator: seriale se un solo worker
    `%op%` <- if (foreach::getDoParWorkers() > 1) `%dopar%` else `%do%`
  }
  
  # ---- Helper: evaluate fitness of the whole population ----
  eval_population <- function(chromosome) {
    if (parallel) {
      fitness_list <- foreach::foreach(i = 1:popsize, .combine = c,
                                       .export = c("decode_solution", "geneLength", "Genome",
                                                   "eval_fitness", "eval_func_inputs")) %op% {
        bits <- chromosome[i, ]
        sol  <- decode_solution(bits, geneLength, Genome)
        eval_fitness(sol, eval_func_inputs)
      }
    } else {
      fitness_list <- lapply(seq_len(popsize), function(i) {
        sol <- decode_solution(chromosome[i, ], geneLength, Genome)
        eval_fitness(sol, eval_func_inputs)
      })
    }
    unlist(fitness_list, use.names = FALSE)
  }
  
  # ---- Working vars ----
  q_alphabeta      <- array(0.0, c(genomeLength, 2L, popsize))
  work_q_alphabeta <- array(0.0, c(genomeLength, 2L, popsize))
  chromosome       <- array(0L,  c(popsize,     genomeLength))
  
  # Hadamard gate
  h <- array(c(1 / sqrt(2.0), 1 / sqrt(2.0),
               1 / sqrt(2.0), -1 / sqrt(2.0)), c(2L, 2L))
  # Rotation gate
  rot <- array(0.0, c(2L, 2L))
  
  # ---- Output history ----
  res <- data.frame(
    generation      = seq_len(generation_max + 1L),
    fitness_average = numeric(generation_max + 1L),
    fitness_best    = numeric(generation_max + 1L)
  )
  
  # ---- Best tracking (disambiguated) ----
  best_idx_by_gen <- integer(generation_max + 1L) # index of best per generation
  best_solution   <- integer(genomeLength)        # best chromosome (bitstring) of current gen
  solution_best   <- integer(genomeLength)        # global best chromosome (bitstring)
  fitness_best    <- -Inf
  generation      <- 1L
  
  # ---- Init population ----
  theta <- thetainit
  q_alphabeta <- generate_pop(popsize, genomeLength, q_alphabeta, rot, theta, h, array(c(1,0), c(2L,1L)))
  chromosome  <- measure(popsize, genomeLength, q_alphabeta, chromosome)
  chromosome  <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
  
  # ---- Initial fitness (decode -> evaluate) ----
  fitness         <- eval_population(chromosome)
  fitness_max     <- max(fitness)
  fitness_average <- mean(fitness)
  
  best_index                  <- which.max(fitness)
  best_idx_by_gen[generation] <- best_index
  best_solution               <- chromosome[best_index, ]
  solution_max                <- best_solution
  
  if (fitness_max > fitness_best) {
    fitness_best  <- fitness_max
    solution_best <- solution_max
  }
  last_improvement <- generation
  n_restarts       <- 0L
  theta_start      <- 0L   # generation from which theta decays (moved by restart_theta)

  # personal attractors: initially each individual's own solution
  if (attractor_type == "personal") {
    attractor         <- chromosome
    attractor_fitness <- fitness
  }
  
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation]    <- fitness_best
  if (plotting) plot_Output(res[seq_len(generation), ])
  if (verbose)  cat("\n", generation, ",", fitness_average, ",", fitness_max)
  
  # ---- Loop ----
  if (progress) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  iter        <- 0L
  old_fitness <- -Inf
  
  while (generation <= generation_max &&
         stop_limit  > fitness_max &&
         (is.null(stop_iters) || (!res$fitness_best[iter + 1L] - old_fitness == 0))) {
    
    iter <- iter + 1L
    if (!is.null(stop_iters) && iter > stop_iters) {
      old_fitness <- res$fitness_best[iter - stop_iters]
    }
    
    if (progress) setTxtProgressBar(pb, generation)
    
    # linear decay of theta
    theta <- thetainit - ((thetainit - thetaend) / (generation_max - theta_start)) *
      (generation - theta_start)
    if (theta < 0) theta <- 0
    
    # attractors: best individual of the current generation for everybody
    if (attractor_type == "generation") {
      b <- best_idx_by_gen[generation]
      attractor         <- matrix(chromosome[b, ], popsize, genomeLength, byrow = TRUE)
      attractor_fitness <- rep(fitness[b], popsize)
    }
    
    q_alphabeta <- rotation(chromosome,
                            attractor,
                            attractor_fitness,
                            genomeLength,
                            solution_best,      # best bitstring so far
                            q_alphabeta,
                            work_q_alphabeta,
                            popsize,
                            fitness, 
                            theta)
    
    generation <- generation + 1L
    
    # mutation rates
    pop_mutation_rate <- pop_mutation_rate_init - ((pop_mutation_rate_init - pop_mutation_rate_end) / generation_max) * generation
    mutation_rate     <- mutation_rate_init     - ((mutation_rate_init     - mutation_rate_end)     / generation_max) * generation
    
    if (mutation_flag) {
      q_alphabeta <- mutation(pop_mutation_rate,
                              mutation_rate,
                              popsize,
                              chromosome,
                              solution_best,
                              q_alphabeta,
                              work_q_alphabeta,
                              genomeLength)
    }
    
    if (interference_flag) {
      q_alphabeta <- interference(interference_rate,
                                  popsize,
                                  chromosome,
                                  solution_best,
                                  q_alphabeta,
                                  genomeLength,
                                  interference_type,
                                  fitness,
                                  interference_lambda)
    }
    
    if (epsilon > 0) q_alphabeta <- h_epsilon(q_alphabeta, epsilon)
    
    chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
    chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
    
    # fitness (decode -> evaluate)
    fitness         <- eval_population(chromosome)
    fitness_max     <- max(fitness)
    fitness_average <- mean(fitness)
    
    best_index                  <- which.max(fitness)
    best_idx_by_gen[generation] <- best_index
    best_solution               <- chromosome[best_index, ]
    solution_max                <- best_solution
    
    if (fitness_max > fitness_best) {
      fitness_best     <- fitness_max
      solution_best    <- solution_max
      last_improvement <- generation
    }
    
    # personal attractors: update with improved solutions, then migrate
    if (attractor_type == "personal") {
      improved <- fitness > attractor_fitness
      attractor[improved, ]       <- chromosome[improved, ]
      attractor_fitness[improved] <- fitness[improved]
      mig <- migration(attractor, attractor_fitness, generation,
                       migration_group_size, migration_local_period,
                       migration_global_period)
      attractor         <- mig$attractor
      attractor_fitness <- mig$attractor_fitness
    }
    
    # partial restart after restart_iters generations without improvement
    if (!is.null(restart_iters) && generation - last_improvement >= restart_iters) {
      q_alphabeta      <- partial_restart(q_alphabeta, restart_strength,
                                          chromosome, solution_best)
      last_improvement <- generation
      n_restarts       <- n_restarts + 1L
      if (restart_theta && generation < generation_max) theta_start <- generation
      if (verbose) cat("\n *** Partial restart at generation", generation)
    }
    
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation]    <- fitness_best
    if (plotting) plot_Output(res[seq_len(generation), ])
    if (verbose)  cat("\n", generation, ",", fitness_average, ",", fitness_best)
  }
  
  if (progress) close(pb)
  cat("\n *** Best fitness: ", fitness_best)
  
  # history plot safe-guard
  last_row <- max(1L, min(nrow(res), ifelse(is.null(iter), 1L, iter)))
  plot_Output(res[seq_len(last_row), ])
  
  # ---- Decode final best solution (bitstring -> integers 1..nvalues_sol) ----
  solution <- decode_solution(solution_best, geneLength, Genome)
  
  out <- list(solution, res[seq_len(last_row), ])
  return(out)
}
