#---------------------------
# EVALUATION OF THE SOLUTION                   
#---------------------------
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