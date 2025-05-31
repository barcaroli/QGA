#' Quantum Genetic Algorithm
#'
#' @description 
#' 
#' Main function to execute a Quantum Genetic Algorithm
#' 
#' @details
#' 
#' This function is the 'engine', which performs the quantum genetic algorithm calling
#' the function for the evaluation of the fitness that is specific for the particulare
#' problem to be optmized.
#' 
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration
#' (default is 20)
#' @param generation_max the number of iterations to be performed
#' (default is 200)
#' @param Genome the length of the genome (or chromosome), representing a possible solution 
#' @param nvalues_sol the number of possible integer values contained in each element (gene) of the solution 
#' @param thetainit the angle (expressed in radiants) to be used when applying the rotation gate
#' when starting the iterations 
#' (default is pi * 0.05, where pi = 3.1415926535)
#' @param thetaend the angle (expressed in radiants) to be used when applying the rotation gate 
#' at the end of the iterations
#' (default is pi * 0.025, where pi = 3.1415926535)
#' @param pop_mutation_rate_init initial mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param pop_mutation_rate_end final mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param mutation_rate_init initial mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome  (default is 1/(Genome+1)))
#' @param mutation_rate_end final mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome (default is 1/(Genome+1))
#' @param mutation_flag flag indicating if the mutation gate is to be applied or not (default is TRUE)
#' @param plotting flag indicating plotting during iterations
#' @param verbose flag indicating printing fitness during iterations
#' @param progress flag indicating progress bar during iterations
#' @param eval_fitness name of the function that will be used to evaluate the fitness of each solution
#' @param eval_func_inputs specific inputs required by the eval_fitness function
#' @param stop_limit value to stop the iterations if the fitness is higher than a given value
#' @param stop_iters number of iterations without improvement of fitness before stopping 
#' 
#' 
#' @export
#' 
#' @return A numeric vector (positive integers) giving the best solution obtained by the QGA
#' 
#' @examples 
#' #----------------------------------------
#' # Fitness evaluation for Knapsack Problem
#' #----------------------------------------
#' KnapsackProblem <- function(solution,
#'                             eval_func_inputs) {
#'   solution <- solution - 1
#'   items <- eval_func_inputs[[1]]
#'   maxweight <- eval_func_inputs[[2]]
#'   tot_items <- sum(solution)
#'   # Penalization
#'   if (sum(items$weight[solution]) > maxweight) {
#'     tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
#'   }
#'   return(tot_items)
#' }
#' #----------------------------------------
#' # Prepare data for fitness evaluation
#' items <- as.data.frame(list(Item = paste0("item",c(1:300)),
#'                             weight = rep(NA,300)))
#' set.seed(1234)
#' items$weight <- rnorm(300,mean=50,sd=20)
#' hist(items$weight)
#' sum(items$weight)
#' maxweight = sum(items$weight) / 2
#' maxweight
#' #----------------------
#' # Perform optimization
#' popsize = 20
#' Genome = nrow(items)
#' solutionQGA <- QGA(popsize = 20,
#'                 generation_max = 500,
#'                 nvalues_sol = 2,
#'                 Genome = nrow(items),
#'                 thetainit = 3.1415926535 * 0.05,
#'                 thetaend = 3.1415926535 * 0.025,
#'                 pop_mutation_rate_init = 1/(popsize + 1),
#'                 pop_mutation_rate_end = 1/(popsize + 1),
#'                 mutation_rate_init = 1,
#'                 mutation_rate_end = 1,
#'                 mutation_flag = TRUE,
#'                 plotting = FALSE,
#'                 verbose = FALSE,
#'                 progress = FALSE,
#'                 eval_fitness = KnapsackProblem,
#'                 eval_func_inputs = list(items,
#'                                         maxweight),
#'                 stop_iters = NULL)
#' #----------------------
#' # Analyze results
#' solution <- solutionQGA[[1]]
#' solution <- solution - 1
#' sum(solution)
#' sum(items$weight[solution])
#' maxweight
#' 
 
QGA <- function(popsize = 20,  
                generation_max = 200,
                nvalues_sol,
                Genome,
                thetainit = pi * 0.05,
                thetaend = pi * 0.025,
                pop_mutation_rate_init = NULL,
                pop_mutation_rate_end = NULL,
                mutation_rate_init = NULL,
                mutation_rate_end = NULL,
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose = TRUE,
                progress = TRUE,
                eval_fitness,
                eval_func_inputs,
                stop_limit = NULL,
                stop_iters = NULL,
                use_parallel = TRUE) {
  # Parallel setup (optional)
  if (use_parallel && popsize > 4) {
    n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({parallel::stopCluster(cl); doParallel::registerDoSEQ()}, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }
  
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome)) stop("Genome parameter value missing!")
  
  if (is.null(pop_mutation_rate_init) & mutation_flag) pop_mutation_rate_init <- 1/(popsize+1)
  if (is.null(pop_mutation_rate_end) & mutation_flag) pop_mutation_rate_end <- 1/(popsize+1)
  if (is.null(mutation_rate_init) & mutation_flag) mutation_rate_init <- 1/(Genome+1)
  if (is.null(mutation_rate_end) & mutation_flag) mutation_rate_end <- 1/(Genome+1)
  
  n <- 0
  while (nvalues_sol > 2^n) n <- n + 1
  geneLength <- n
  genomeLength <- Genome * geneLength
  
  qubit_0 <- array(c(1, 0), c(2, 1))
  h <- array(c(1 / sqrt(2), 1 / sqrt(2), 1 / sqrt(2), -1 / sqrt(2)), c(2, 2))
  
  # Initialize result storage
  res <- data.frame(
    generation = seq_len(generation_max + 1),
    fitness_average = rep(0, generation_max + 1),
    fitness_best = rep(0, generation_max + 1)
  )
  
  fitness_best <- -Inf
  solution_best <- rep(0, genomeLength)
  generation <- 1
  theta <- thetainit
  
  q_alphabeta <- generate_pop(popsize, genomeLength, array(0.0, c(genomeLength, 2, popsize)), 
                              array(0.0, c(2, 2)), theta, h, qubit_0)
  chromosome <- measure(popsize, genomeLength, q_alphabeta, array(0, c(popsize, genomeLength)))
  chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
  
  # Fitness evaluation (parallel if enabled)
  fitness <- foreach(i = 1:popsize, .combine = c) %dopar% {
    eval_fitness(chromosome[i, ], eval_func_inputs)
  }
  fitness_max <- max(fitness)
  fitness_average <- mean(fitness)
  best_idx <- which.max(fitness)
  best_chromosome <- chromosome[best_idx, ]
  solution_max <- best_chromosome
  
  if (fitness_max > fitness_best) {
    fitness_best <- fitness_max
    solution_best <- solution_max
  }
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation] <- fitness_best
  if (plotting) plot_Output(res[seq_len(generation), ])
  if (verbose) cat("\n", generation, ",", fitness_average, ",", fitness_max)
  
  if (progress) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  if (is.null(stop_limit)) stop_limit <- Inf
  iter <- 0
  old_fitness <- -Inf
  
  # Main evolution loop
  while (
    generation <= generation_max &&
    stop_limit > fitness_max &&
    (is.null(stop_iters) || iter < stop_iters || (res$fitness_best[generation] - old_fitness != 0))
  ) {
    iter <- iter + 1
    if (!is.null(stop_iters) && iter > stop_iters) {
      old_fitness <- res$fitness_best[generation - stop_iters]
    }
    if (progress) setTxtProgressBar(pb, generation)
    
    theta <- max(0, thetainit - ((thetainit - thetaend) / generation_max) * generation)
    q_alphabeta <- rotation(chromosome, best_chromosome, generation, genomeLength,
                            solution_best, q_alphabeta, array(0.0, c(genomeLength, 2, popsize)), popsize, fitness, theta)
    generation <- generation + 1
    pop_mutation_rate <- pop_mutation_rate_init - ((pop_mutation_rate_init - pop_mutation_rate_end) / generation_max) * generation
    mutation_rate <- mutation_rate_init - ((mutation_rate_init - mutation_rate_end) / generation_max) * generation
    if (mutation_flag) {
      q_alphabeta <- mutation(pop_mutation_rate, mutation_rate, popsize, chromosome,
                              solution_best, q_alphabeta, array(0.0, c(genomeLength, 2, popsize)), genomeLength)
    }
    chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
    chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
    
    # Fitness evaluation (parallel if enabled)
    fitness <- foreach(i = 1:popsize, .combine = c) %dopar% {
      eval_fitness(chromosome[i, ], eval_func_inputs)
    }
    fitness_max <- max(fitness)
    fitness_average <- mean(fitness)
    best_idx <- which.max(fitness)
    best_chromosome <- chromosome[best_idx, ]
    solution_max <- best_chromosome
    
    if (fitness_max > fitness_best) {
      fitness_best <- fitness_max
      solution_best <- solution_max
    }
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation] <- fitness_best
    if (plotting) plot_Output(res[seq_len(generation), ])
    if (verbose) cat("\n", generation, ",", fitness_average, ",", fitness_best)
  }
  if (progress) close(pb)
  cat("\n *** Best fitness: ", fitness_best)
  plot_Output(res[seq_len(generation), ])
  
  # Vectorized solution decoding
  powers <- 2^((geneLength - 1):0)
  solution1 <- matrix(solution_best, nrow = geneLength)
  solution <- as.integer(colSums(solution1 * powers) + 1)
  
  out <- list(solution, res[seq_len(generation), ])
  return(out)
}
