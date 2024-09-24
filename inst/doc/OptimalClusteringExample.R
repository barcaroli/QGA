#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# for optimal multivariate clustering
#---------------------------------------------

#----------------------
library(QGA)
#-------------------------------------------
# Fitness evaluation for best stratification
#-------------------------------------------
clustering <- function(solution, eval_func_inputs) {
  maxvalue <- 5
  penalfactor <- 2
  df <- eval_func_inputs[[1]]
  vars <- eval_func_inputs[[2]]
  fitness <- 0
  for (v in vars) {
    cv <- tapply(df[,v],solution,FUN=sd) / tapply(df[,v],solution,FUN=mean)
    cv <- ifelse(is.na(cv),maxvalue,cv)
    fitness <- fitness + sum(cv)
  }
  # Penalisation on unbalanced clusters
  b <- table(solution)/nrow(df)
  fitness <- fitness + penalfactor * (sum(abs(b - c(rep(1/(length(b)),length(b))))))
  return(-fitness)
}


#-------------------------------------------
# Prepare data for fitness evaluation
data(iris)
vars <- colnames(iris)[1:4]
nclust = 3

# Check value of the fitness with Species
fitness <- 0
for (v in vars) {
  fitness <- fitness + sum(tapply(iris[,v],iris$Species,FUN=sd) / tapply(iris[,v],iris$Species,FUN=mean))
  cat("\n Var: ", v, "  sum of cvs: ", sum(tapply(iris[,v],iris$Species,FUN=sd) / tapply(iris[,v],iris$Species,FUN=mean)))
}
fitness
# [1] 1.627784

#----------------------
# Perform optimization
popsize = 20
Genome = nrow(iris)
set.seed(1234)
solutionQGA <- QGA(popsize,
                generation_max = 1000,
                nvalues_sol = nclust,
                Genome,
                thetainit = 3.1415926535 * 0.1,
                thetaend = 3.1415926535 * 0.05,
                pop_mutation_rate_init = 1/(popsize + 1),
                pop_mutation_rate_end = 1/(popsize + 1),
                mutation_rate_init = 1/(Genome + 1),
                mutation_rate_end = 1/(Genome + 1),
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose = FALSE,
                eval_fitness = clustering,
                eval_func_inputs = list(iris, vars),
                stop_iters = 300)
#----------------------
# Analyze results
solution <- solutionQGA[[1]]
table(solution)
fitness <- 0
df <- iris
for (v in vars) {
  cv <- tapply(df[,v],solution,FUN=sd) / tapply(df[,v],solution,FUN=mean)
  cv <- ifelse(is.na(cv),maxvalue,cv)
  fitness <- fitness + sum(cv)
}
fitness
# [1] 1.624756
iris$stratum <- solution
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      0  0 50
# versicolor 45  5  0
# virginica   5 45  0


#-------------------------------
# Comparison with genalg
library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  maxvalue <- 5
  penalfactor <- 2
  fitness <- 0
  for (v in vars) {
    cv <- tapply(iris[,v],solution,FUN=sd) / tapply(iris[,v],solution,FUN=mean)
    cv <- ifelse(is.na(cv),maxvalue,cv)
    fitness <- fitness + sum(cv)
  }
  b <- table(solution)/nrow(iris)
  # Penalisation on unbalanced clusters
  fitness <- fitness + penalfactor * (sum(abs(b - c(rep(1/(length(b)),length(b))))))
  return(fitness)
}
set.seed(1234)
solution_genalg <- rbga(stringMin=c(rep(1,nrow(iris))), 
                   stringMax=c(rep(nclust,nrow(iris))),
                   popSize=20, 
                   iters=2500, 
                   elitism=NA, 
                   mutationChance = 5/(nrow(iris)+1),
                   evalFunc=evaluate)
plot(solution_genalg)
filter = solution_genalg$evaluations == min(solution_genalg$evaluations)
bestObjectCount = sum(rep(1, solution_genalg$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = solution_genalg$population[filter, ][1, 
  ]
} else {
  bestSolution = solution_genalg$population[filter, ]
}
bestSolution <- round(bestSolution)
fitness <- 0
for (v in vars) {
  cv <- tapply(iris[,v],bestSolution,FUN=sd) / tapply(iris[,v],bestSolution,FUN=mean)
  cv <- ifelse(is.na(cv),maxvalue,cv)
  fitness <- fitness + sum(cv)
}
fitness
# [1] 1.697082
iris$stratum <- bestSolution
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      0 50  0
# versicolor  9  0 41
# virginica  41  0  9
