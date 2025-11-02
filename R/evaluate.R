#---------------------------
# EVALUATION OF THE SOLUTION                   
#---------------------------

#' Evaluate Population Fitness
#'
#' @description
#' Decodes each chromosome into integer gene values and evaluates the
#' user-provided fitness function. Tracks best and average fitness and
#' updates the `best_chromosome` index for the current generation.
#'
#' @param chromosome Integer matrix `[popsize, geneLength*Genome]` of 0/1 bits.
#' @param best_chromosome Integer vector; best individual index per generation (updated in place for `generation`).
#' @param popsize Integer. Population size.
#' @param Genome Integer. Number of genes.
#' @param geneLength Integer. Number of bits per gene.
#' @param nvalues_sol Integer. Allowed values per gene (not used internally).
#' @param generation Integer. Current generation index (1-based).
#' @param eval_fitness Function. `function(solution_integers, eval_func_inputs)` returning numeric fitness.
#' @param eval_func_inputs Any. Second argument passed to `eval_fitness`.
#'
#' @return A list with elements:
#' - `fitness`: numeric vector length `popsize`.
#' - `fitness_max`: numeric scalar, best fitness.
#' - `fitness_average`: numeric scalar, mean fitness.
#' - `best_chromosome`: updated integer vector of best indices.
#' - `solution_max`: integer vector bitstring of the best chromosome.
#' @keywords internal
evaluate <- function(chromosome,
                     best_chromosome,
                     popsize,
                     Genome,
                     geneLength,
                     nvalues_sol,
                     generation,
                     eval_fitness,
                     eval_func_inputs){
  fitness <- array(0.0, c(1, popsize))
  fitness_total <- 0
  fitness_average <- -99999999
  fitness_max <- -999999999
  the_best_chromosome <- 0
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(geneLength,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:geneLength)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(geneLength - y) 
      }
    }
    solution <- solution + 1
    fitness[i] <- eval_fitness(solution,eval_func_inputs)
    fitness_total <- fitness_total + fitness[i]
    if (fitness[i] >= fitness_max) {
      fitness_max <- fitness[i]
      the_best_chromosome <- i
      solution_max <- chromosome[i, ]
    }
  }
  fitness_average <- fitness_total / popsize
  best_chromosome[generation] <- the_best_chromosome
  return(list(
    fitness = fitness,
    fitness_max = fitness_max,
    fitness_average = fitness_average,
    best_chromosome = best_chromosome,
    solution_max = solution_max)
  )
}
