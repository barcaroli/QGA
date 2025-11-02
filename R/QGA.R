#' Quantum Genetic Algorithm
#'
#' @description
#' Runs a Quantum Genetic Algorithm (QGA) to optimize a user-defined
#' fitness function over discrete decision variables. Internally it
#' initializes a quantum population, iterates measurement, rotation,
#' and optional mutation, repairs invalid genes, evaluates fitness in
#' parallel when possible, and tracks best/average fitness across generations.
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
                stop_iters = NULL) {
  
  # ---- Parallel deps ----
  if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
  if (!requireNamespace("foreach",    quietly = TRUE)) install.packages("foreach")
  library(doParallel)
  library(foreach)
  
  # ---- Checks & defaults ----
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome))      stop("Genome parameter value missing!")
  if (is.null(pop_mutation_rate_init) && mutation_flag) pop_mutation_rate_init <- 1/(popsize+1)
  if (is.null(pop_mutation_rate_end)  && mutation_flag) pop_mutation_rate_end  <- 1/(popsize+1)
  if (is.null(mutation_rate_init)     && mutation_flag) mutation_rate_init     <- 1/(Genome+1)
  if (is.null(mutation_rate_end)      && mutation_flag) mutation_rate_end      <- 1/(Genome+1)
  if (is.null(stop_limit)) stop_limit <- Inf
  
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
  
  # ---- Parallel backend ----
  numCores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  # Export una tantum del decoder ai worker (evita i warning)
  parallel::clusterExport(cl, varlist = c("decode_solution"),
                          envir = environment())
  
  # Fallback operator: seriale se un solo worker
  `%op%` <- if (foreach::getDoParWorkers() > 1) `%dopar%` else `%do%`
  
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
  fitness_list <- foreach(i = 1:popsize, .combine = c) %op% {
    bits <- chromosome[i, ]
    sol  <- decode_solution(bits, geneLength, Genome)
    eval_fitness(sol, eval_func_inputs)
  }
  fitness         <- unlist(fitness_list, use.names = FALSE)
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
    theta <- thetainit - ((thetainit - thetaend) / generation_max) * generation
    if (theta < 0) theta <- 0
    
    # rotation expects best indices history + current generation
    q_alphabeta <- rotation(chromosome,
                            best_idx_by_gen,    # pass indices, not bitstrings
                            generation,
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
    
    chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
    chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
    
    # fitness (decode -> evaluate)
    fitness_list <- foreach(i = 1:popsize, .combine = c) %op% {
      bits <- chromosome[i, ]
      sol  <- decode_solution(bits, geneLength, Genome)
      eval_fitness(sol, eval_func_inputs)
    }
    fitness         <- unlist(fitness_list, use.names = FALSE)
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
