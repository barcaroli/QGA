#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the Knapsack Problem
#---------------------------------------------

#----------------------
library(QGA)
#----------------------

#----------------------------------------
# Fitness evaluation for Knapsack Problem
#----------------------------------------
KnapsackProblem <- function(solution,
                            eval_func_inputs) {
  solution <- solution - 1
  items <- eval_func_inputs[[1]]
  maxweight <- eval_func_inputs[[2]]
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(tot_items)
}


#----------------------------------------
# Prepare data for fitness evaluation
items <- as.data.frame(list(Item = paste0("item",c(1:300)),
                            weight = rep(NA,300)))
set.seed(1234)
items$weight <- rnorm(300,mean=50,sd=20)
hist(items$weight)
sum(items$weight)
maxweight = sum(items$weight) / 2
maxweight


#----------------------
# Perform optimization
popsize = 20
Genome = nrow(items)
solution <- QGA(popsize = 20,
                generation_max = 500,
                nvalues_sol = 2,
                Genome = nrow(items),
                thetainit = 3.1415926535 * 0.05,
                thetaend = 3.1415926535 * 0.025,
                pop_mutation_rate_init = 1/(popsize + 1),
                pop_mutation_rate_end = 1/(popsize + 1),
                mutation_rate_init = 1,
                mutation_rate_end = 1,
                mutation_flag = TRUE,
                plotting = FALSE,
                verbose = FALSE,
                eval_fitness = KnapsackProblem,
                eval_func_inputs = list(items,maxweight))

#----------------------
# Analyze results
solution <- solution - 1
sum(solution)
sum(items$weight[solution])
maxweight

#-----------------------------------------
# Compare with classical genetic algorithm
library(genalg)
evaluate <- function(solution) {
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(-tot_items)
}
rbga.results <- rbga.bin(size=nrow(items),
         suggestions=NULL,
         popSize=20, 
         iters=1000, 
         elitism=NA, 
         zeroToOneRatio=round(maxweight*10/sum(items$weight)),
         evalFunc=evaluate)
plot(rbga.results)
filter = rbga.results$evaluations == min(rbga.results$evaluations)
bestObjectCount = sum(rep(1, rbga.results$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = rbga.results$population[filter, ][1, 
  ]
} else {
  bestSolution = rbga.results$population[filter, ]
}
sum(bestSolution)
sum(items$weight[bestSolution])
maxweight

