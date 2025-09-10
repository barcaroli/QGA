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
                thetainit = 3.1415926535 * 0.05,
                thetaend = 3.1415926535 * 0.025,
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
                stop_iters = NULL) {
  # ---- LIBRERIE PARALLEL ----
  if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
  if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
  library(doParallel)
  library(foreach)

  numCores <- parallel::detectCores() - 1
  cl <- makeCluster(numCores)
  doParallel::registerDoParallel(cl)

  # ---- CHECK E DEFAULT ----
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome)) stop("Genome parameter value missing!")
  if (is.null(pop_mutation_rate_init) & mutation_flag == TRUE) pop_mutation_rate_init = 1/(popsize+1)
  if (is.null(pop_mutation_rate_end) & mutation_flag == TRUE) pop_mutation_rate_end = 1/(popsize+1)
  if (is.null(mutation_rate_init) & mutation_flag == TRUE) mutation_rate_init = 1/(Genome+1)
  if (is.null(mutation_rate_end) & mutation_flag == TRUE) mutation_rate_end = 1/(Genome+1)

  # ---- CODIFICA GENOMA ----
  n = 0
  while (nvalues_sol > 2^n) {
    n = n+1
  }
  geneLength = n 
  genomeLength <- Genome * geneLength 

  # ---- VARIABILI DI LAVORO ----
  qubit_0 <- array(c(1, 0), c(2, 1))
  qubit_1 <- array(c(0, 1), c(2, 1))

  fitness <- array(0.0, c(1, popsize))
  q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  work_q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  chromosome <- array(0, c(popsize, genomeLength))
  best_chromosome <- array(0.0, generation_max)

  # Hadamard gate
  h <- array(c(1 / sqrt(2.0), 1 / sqrt(2.0), 1 / sqrt(2.0), -1 / sqrt(2.0)), c(2, 2))

  # Rotation Q-gate
  rot <- array(0.0, c(2, 2))

  # ---- OUTPUT ----
  res <- NULL
  res$generation <- c(1:(generation_max + 1))
  res$fitness_average <- rep(0, (generation_max + 1))
  res$fitness_best <- rep(0, (generation_max + 1))
  res <- as.data.frame(res)

  fitness_best <- -999999
  solution_best <- rep(0, genomeLength)
  generation <- 1

  # ---- INIZIALIZZAZIONE ----
  theta <- thetainit
  q_alphabeta <- generate_pop(popsize,
                              genomeLength,
                              q_alphabeta,
                              rot,
                              theta,
                              h,
                              qubit_0)
  chromosome <- measure(popsize,
                        genomeLength,
                        q_alphabeta,
                        chromosome)
  chromosome <- repair(popsize,
                       chromosome,
                       geneLength,
                       genomeLength,
                       nvalues_sol,
                       Genome)
  a <-evaluate(chromosome,
                    best_chromosome,
                    popsize,
                    Genome,
                    geneLength,
                    nvalues_sol,
                    generation,
                    eval_fitness,
                    eval_func_inputs)
  fitness <- a$fitness
  fitness_max <- a$fitness_max
  fitness_average <- a$fitness_average
  best_chromosome <- a$best_chromosome
  solution_max <- a$solution_max
  if (fitness_max > fitness_best) {
    fitness_best <- fitness_max
    solution_best <- solution_max
  }
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation] <- fitness_best
  if (plotting == TRUE) plot_Output(res[c(1:generation), ])
  if (verbose == TRUE) cat("\n", generation, ",", fitness_average, ",", fitness_max)

  if (progress == TRUE) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  if (is.null(stop_limit)) stop_limit <- Inf
  iter <- 0
  old_fitness <- -Inf

  while (generation <= generation_max & 
         stop_limit > fitness_max &
         (is.null(stop_iters) | (!res$fitness_best[iter+1] - old_fitness == 0))) {
    iter <- iter + 1
    if (!is.null(stop_iters)) {
      if (iter > stop_iters) {
        old_fitness <- res$fitness_best[iter-stop_iters]
      }
    } 
    if (progress == TRUE) setTxtProgressBar(pb, generation)
    theta <- thetainit - ((thetainit - thetaend) / generation_max) * generation
    if (theta < 0) theta <- 0
    q_alphabeta <- rotation(chromosome,
                            best_chromosome,
                            generation,
                            genomeLength,
                            solution_best,
                            q_alphabeta,
                            work_q_alphabeta,
                            popsize,
                            fitness, 
                            theta)
    generation <- generation + 1
    pop_mutation_rate = pop_mutation_rate_init - ((pop_mutation_rate_init - pop_mutation_rate_end) / generation_max) * generation
    mutation_rate = mutation_rate_init - ((mutation_rate_init - mutation_rate_end) / generation_max) * generation
    if (mutation_flag == TRUE) {
      q_alphabeta <- mutation(pop_mutation_rate, 
                              mutation_rate,
                              popsize,
                              chromosome,
                              solution_best,
                              q_alphabeta,
                              work_q_alphabeta,
                              genomeLength)      
    }
    chromosome <- measure(popsize,
                          genomeLength,
                          q_alphabeta,
                          chromosome)
    chromosome <- repair(popsize,
                         chromosome,
                         geneLength,
                         genomeLength,
                         nvalues_sol,
                         Genome)
    # --- VALUTAZIONE PARALLELA FITNESS ---
    fitness_list <- foreach(i = 1:popsize, .combine = c, .packages = c()) %dopar% {
      eval_fitness(chromosome[i,], eval_func_inputs)
    }
    fitness <- unlist(fitness_list)
    fitness_max <- max(fitness)
    fitness_average <- mean(fitness)
    best_index <- which.max(fitness)
    best_chromosome <- chromosome[best_index,]
    solution_max <- chromosome[best_index,]

    if (fitness_max > fitness_best) {
      fitness_best <- fitness_max
      solution_best <- solution_max
    }
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation] <- fitness_best
    if (plotting == TRUE) plot_Output(res[c(1:generation), ])
    if (verbose == TRUE) cat("\n", generation, ",", fitness_average, ",", fitness_best)
  }
  if (progress == TRUE) close(pb)
  cat("\n *** Best fitness: ",fitness_best)
  plot_Output(res[c(1:iter), ])
  solution1 <- array(solution_best,c(geneLength,Genome))
  solution <- c(rep(0,Genome))
  for (x in c(1:Genome)) {
    for (y in c(1:geneLength)) {
      solution[x] <- solution[x] + solution1[y,x]*2^(geneLength - y) 
    }
  }
  solution <- solution + 1
  out <- list(solution,res[c(1:iter), ])

  stopCluster(cl) # <-- IMPORTANTE, chiude il cluster parallelo!
  return(out)
}
