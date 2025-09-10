#' Quantum Genetic Algorithm
#'
#' @description 
#' Main function to execute a Quantum Genetic Algorithm
#'
#' @details
#' This function is the 'engine', which performs the quantum genetic algorithm calling
#' the function for the evaluation of the fitness that is specific for the particular
#' problem to be optimized.
#'
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration (default 20)
#' @param generation_max the number of iterations to be performed (default 200)
#' @param Genome the length of the genome (or chromosome), representing a possible solution 
#' @param nvalues_sol the number of possible integer values contained in each element (gene) of the solution 
#' @param thetainit angle (radians) for the rotation gate at the beginning (default pi*0.05)
#' @param thetaend  angle (radians) for the rotation gate at the end (default pi*0.025)
#' @param pop_mutation_rate_init initial mutation rate for X-Pauli gate at population level (default 1/(popsize+1))
#' @param pop_mutation_rate_end  final   mutation rate for X-Pauli gate at population level (default 1/(popsize+1))
#' @param mutation_rate_init initial mutation rate for X-Pauli at gene level (default 1/(Genome+1))
#' @param mutation_rate_end  final   mutation rate for X-Pauli at gene level (default 1/(Genome+1))
#' @param mutation_flag apply mutation gate or not (default TRUE)
#' @param plotting plot during iterations
#' @param verbose print fitness during iterations
#' @param progress show progress bar
#' @param eval_fitness function to evaluate the fitness of a solution
#' @param eval_func_inputs list of inputs required by eval_fitness
#' @param stop_limit stop if best fitness exceeds this value (default Inf)
#' @param stop_iters stop after this many iterations without improvement (default NULL = disabled)
#'
#' @export
#' @return list(best_solution_as_integers, fitness_history_dataframe)
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
  
  # ---- DEPENDENCIES (parallel) ----
  if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
  if (!requireNamespace("foreach", quietly = TRUE))    install.packages("foreach")
  library(doParallel)
  library(foreach)
  
  # ---- CHECKS & DEFAULTS ----
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome))      stop("Genome parameter value missing!")
  
  if (is.null(pop_mutation_rate_init) && mutation_flag) pop_mutation_rate_init <- 1/(popsize+1)
  if (is.null(pop_mutation_rate_end)  && mutation_flag) pop_mutation_rate_end  <- 1/(popsize+1)
  if (is.null(mutation_rate_init)     && mutation_flag) mutation_rate_init     <- 1/(Genome+1)
  if (is.null(mutation_rate_end)      && mutation_flag) mutation_rate_end      <- 1/(Genome+1)
  
  if (is.null(stop_limit)) stop_limit <- Inf
  
  # ---- GENOME ENCODING ----
  n <- 0
  while (nvalues_sol > 2^n) n <- n + 1
  geneLength   <- n
  genomeLength <- Genome * geneLength
  
  # ---- PARALLEL BACKEND ----
  numCores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  # Fallback operator: serial if only one worker
  `%op%` <- if (foreach::getDoParWorkers() > 1) `%dopar%` else `%do%`
  
  # ---- WORKING VARIABLES ----
  # (alcune variabili non sono usate direttamente qui ma lasciate per completezza del flusso)
  qubit_0 <- array(c(1, 0), c(2, 1))
  qubit_1 <- array(c(0, 1), c(2, 1))
  
  q_alphabeta      <- array(0.0, c(genomeLength, 2, popsize))
  work_q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  chromosome       <- array(0L,  c(popsize, genomeLength))
  
  # Hadamard gate
  h <- array(c(1 / sqrt(2.0), 1 / sqrt(2.0),
               1 / sqrt(2.0), -1 / sqrt(2.0)), c(2, 2))
  # Rotation Q-gate
  rot <- array(0.0, c(2, 2))
  
  # ---- OUTPUT HISTORY ----
  res <- data.frame(
    generation      = seq_len(generation_max + 1),
    fitness_average = numeric(generation_max + 1),
    fitness_best    = numeric(generation_max + 1)
  )
  
  # ---- STATE (best tracking, disambiguated) ----
  best_idx_by_gen <- integer(generation_max + 1)   # indice dell'individuo best per generazione
  best_solution   <- integer(genomeLength)         # cromosoma best alla generazione corrente
  solution_best   <- integer(genomeLength)         # miglior cromosoma globale
  fitness_best    <- -Inf
  generation      <- 1
  
  # ---- INITIALIZE POPULATION ----
  theta <- thetainit
  q_alphabeta <- generate_pop(popsize,
                              genomeLength,
                              q_alphabeta,
                              rot,
                              theta,
                              h,
                              qubit_0)
  chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
  chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
  
  # ---- FITNESS EVALUATION (PARALLEL/serial) ----
  fitness_wrapper <- function(sol) eval_fitness(sol, eval_func_inputs)
  
  fitness_list <- foreach(i = 1:popsize, .combine = c,
                          .export = c("fitness_wrapper")) %op% {
                            fitness_wrapper(chromosome[i, ])
                          }
  fitness        <- unlist(fitness_list, use.names = FALSE)
  fitness_max    <- max(fitness)
  fitness_average<- mean(fitness)
  
  best_index                 <- which.max(fitness)
  best_idx_by_gen[generation]<- best_index
  best_solution              <- chromosome[best_index, ]
  solution_max               <- best_solution
  
  if (fitness_max > fitness_best) {
    fitness_best  <- fitness_max
    solution_best <- solution_max
  }
  
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation]    <- fitness_best
  if (plotting) plot_Output(res[seq_len(generation), ])
  if (verbose)  cat("\n", generation, ",", fitness_average, ",", fitness_max)
  
  # ---- LOOP ----
  if (progress) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  iter        <- 0
  old_fitness <- -Inf
  
  while (generation <= generation_max &&
         stop_limit  > fitness_max &&
         (is.null(stop_iters) || (!res$fitness_best[iter + 1] - old_fitness == 0))) {
    
    iter <- iter + 1
    if (!is.null(stop_iters) && iter > stop_iters) {
      old_fitness <- res$fitness_best[iter - stop_iters]
    }
    
    if (progress) setTxtProgressBar(pb, generation)
    
    # aggiorna theta (decadimento lineare)
    theta <- thetainit - ((thetainit - thetaend) / generation_max) * generation
    if (theta < 0) theta <- 0
    
    # NOTA: rotation si aspetta l'INDICE del best per generazione (non il cromosoma):
    q_alphabeta <- rotation(chromosome,
                            best_idx_by_gen,     # <--- passiamo il vettore di INDICI
                            generation,
                            genomeLength,
                            solution_best,       # miglior cromosoma globale
                            q_alphabeta,
                            work_q_alphabeta,
                            popsize,
                            fitness,
                            theta)
    
    generation <- generation + 1
    
    # tassi di mutazione dinamici
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
    
    # fitness (parallel/serial)
    fitness_list <- foreach(i = 1:popsize, .combine = c,
                            .export = c("fitness_wrapper")) %op% {
                              fitness_wrapper(chromosome[i, ])
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
  plot_Output(res[seq_len(max(1, iter)), ])
  
  # decode soluzione_best (da bit a interi 1..nvalues_sol)
  solution1 <- array(solution_best, c(geneLength, Genome))
  solution  <- integer(Genome)
  for (x in seq_len(Genome)) {
    acc <- 0L
    for (y in seq_len(geneLength)) {
      acc <- acc + solution1[y, x] * 2^(geneLength - y)
    }
    solution[x] <- acc
  }
  solution <- solution + 1L
  
  out <- list(solution, res[seq_len(max(1, iter)), ])
  return(out)
}
