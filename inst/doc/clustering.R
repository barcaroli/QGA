#---------------------
# Clustering with QGA
#---------------------

library(QGA)


## -----------------------------------------------------------------------------------------------
clustering <- function(solution, eval_func_inputs) {
  maxvalue <- 5
  penalfactor <- 2
  df <- eval_func_inputs[[1]]
  vars <- eval_func_inputs[[2]]
  # Fitness function
  fitness <- 0
  for (v in vars) {
    cv <- tapply(df[,v],solution,FUN=sd) / tapply(df[,v],solution,FUN=mean)
    cv <- ifelse(is.na(cv),maxvalue,cv)
    fitness <- fitness + sum(cv)
  }
  # Penalization on unbalanced clusters
  b <- table(solution)/nrow(df)
  fitness <- fitness + penalfactor * (sum(abs(b - c(rep(1/(length(b)),length(b))))))
  return(-fitness)
}


## -----------------------------------------------------------------------------------------------
data(iris)
vars <- colnames(iris)[1:4]
vars


## -------------------------------------------------------------------------------------
nclust = 3
popsize = 20
Genome = nrow(iris)
set.seed(1234)
solutionQGA <- QGA(popsize,
                generation_max = 1500,
                nvalues_sol = nclust,
                Genome,
                thetainit = 3.1415926535 * 0.1,
                thetaend = 3.1415926535 * 0.05,
                pop_mutation_rate_init = 1/(popsize + 1),
                pop_mutation_rate_end = 1/(popsize + 1),
                mutation_rate_init = 1/(Genome + 1),
                mutation_rate_end = 1/(Genome + 1),
                mutation_flag = TRUE,
                plotting = FALSE,
                verbose = FALSE,
                progress = FALSE,
                eval_fitness = clustering,
                eval_func_inputs = list(iris, vars))


## ----eval=TRUE, echo=FALSE, include=FALSE-------------------------------------------------------
load("clustering.RData")


## -----------------------------------------------------------------------------------------------
QGA:::plot_Output(CLUSTsolution[[2]])


## -----------------------------------------------------------------------------------------------
solution <- CLUSTsolution[[1]]
table(solution)


## -----------------------------------------------------------------------------------------------
iris$cluster <- solution
xtabs( ~ Species + cluster, data=iris)

