#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the optimization of a sampling frame
#---------------------------------------------

#----------------------
require(SamplingStrata)
library(QGA)
#-------------------------------------------
# Fitness evaluation for best stratification
#-------------------------------------------
BestStratification <- function(solution,
                               eval_func_inputs) {
  frame <- eval_func_inputs[[1]]
  cv <- eval_func_inputs[[2]]
  strata = SamplingStrata::aggrStrata2(dataset=frame,
                                       model=NULL,
                                       vett=solution,
                                       dominio=1)
  fitness <- -sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
  return(fitness)
}

#-------------------------------------------
# Prepare data for fitness evaluation
data(iris)
iris$id <- c(1:nrow(iris))
iris$dom <- 1
frame <- buildFrameDF(
  df = iris,
  id = "id",
  domainvalue = "dom",
  X = c("id"),
  Y = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
)
cv <- as.data.frame(list(
  DOM = "DOM1",
  CV1 = 0.03,
  CV2 = 0.03,
  CV3 = 0.03,
  CV4 = 0.03,
  domainvalue = 1
))
nstrat = 3

#----------------------
# Set parameters
# popsize = 20
# generation_max = 1000
# nvalues_sol = nstrat
# Genome = nrow(iris)
# thetainit = 3.1415926535 * 0.05
# thetaend = 3.1415926535 * 0.025
# pop_mutation_rate_init = 1/(popsize + 1)
# pop_mutation_rate_end = 1/(popsize + 1)
# mutation_rate_init = 1/(Genome + 1)
# mutation_rate_end = 1/(Genome + 1)
# mutation_flag = TRUE
# eval_fitness = BestStratification
# eval_func_inputs = list(frame, cv)

#----------------------
# Perform optimization
popsize = 20
Genome = nrow(iris)
set.seed(1234)
solutionQGA <- QGA(popsize,
                generation_max = 2000,
                nvalues_sol = nstrat,
                Genome,
                thetainit = 3.1415926535 * 0.15,
                thetaend = 3.1415926535 * 0.0125,
                pop_mutation_rate_init = 1/(popsize + 1),
                pop_mutation_rate_end = 1/(popsize + 1),
                mutation_rate_init = 1/(Genome + 1),
                mutation_rate_end = 1/(Genome + 1),
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose = FALSE,
                eval_fitness = BestStratification,
                eval_func_inputs = list(frame, cv))
#----------------------
# Analyze results
solutionQGA <- solutionQGA[[1]]
table(solutionQGA)
strata <- aggrStrata2(dataset = frame, 
                      vett = solutionQGA, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
# 1] 19.75488
iris$stratum <- solutionQGA
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      0 50  0
# versicolor 48  0  2
# virginica   4  0 46

#-------------------------------
# Comparison with SamplingStrata
set.seed(1234)
solution_SamplingStrata <- optimStrata(method = "atomic",
                  framesamp = frame,
                  nStrata = nstrat,
                  errors = cv,
                  pops = popsize,
                  minnumstr = 1,
                  iter = 2000)
sum(solution_SamplingStrata$aggr_strata$SOLUZ)
# [1] 19.75488
iris$stratum <- solution_SamplingStrata$framenew$LABEL
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      6  5 39
# versicolor 40  3  7
# virginica   4 41  5
#-------------------------------
# Comparison with genalg
library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  strata = SamplingStrata::aggrStrata2(dataset=frame,
                                       model=NULL,
                                       vett=solution,
                                       dominio=1)
  fitness <- sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
  return(fitness)
}
solution_genalg <- rbga(stringMin=c(rep(1,nrow(iris))), 
                   stringMax=c(rep(nstrat,nrow(iris))),
                   popSize=20, 
                   iters=2000, 
                   elitism=NA, 
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
strata <- aggrStrata2(dataset = frame, 
                      vett = bestSolution, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
# [1] 19.75488
iris$stratum <- bestSolution
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      0 50  0
# versicolor 48  0  2
# virginica   4  0 46
save.image("run_best_stratification.RData")
