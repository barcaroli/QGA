## ----setup, include = FALSE---------------------------------------------------------------------
options(width = 999)
knitr::opts_chunk$set(fig.width=6, fig.height=4,
  collapse = TRUE,
  comment = "#>"
)
# if (!require(bookdown)) install.packages("bookdown", dependencies=TRUE)
# library(bookdown)
library(QGA)


## -----------------------------------------------------------------------------------------------
KnapsackProblem <- function(solution,eval_func_inputs) {
  solution <- solution - 1
  items <- eval_func_inputs[[1]]
  maxweight <- eval_func_inputs[[2]]
  # Fitness
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(tot_items)
}


## -----------------------------------------------------------------------------------------------
items <- as.data.frame(list(Item = paste0("item",c(1:500)),
                            weight = rep(NA,500)))
set.seed(1234)
items$weight <- rnorm(500,mean=200,sd=80)
head(items)



## -----------------------------------------------------------------------------------------------
sum(items$weight)
maxweight = sum(items$weight) / 5
maxweight


## -----------------------------------------------------------------------------------------------
popsize = 20
generation_max = 500
nvalues_sol = 2
Genome = nrow(items)
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.025
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome+1)
mutation_rate_end = 2/(Genome+1)
mutation_flag = TRUE
plotting = FALSE
verbose = FALSE
progress = FALSE
eval_fitness = KnapsackProblem
eval_func_inputs = list(items,maxweight)


## ----eval=FALSE---------------------------------------------------------------------------------
## set.seed(1234)
## knapsackSolution  <- QGA(popsize,
##                 generation_max,
##                 nvalues_sol,
##                 Genome,
##                 thetainit,
##                 thetaend,
##                 pop_mutation_rate_init,
##                 pop_mutation_rate_end,
##                 mutation_rate_init,
##                 mutation_rate_end,
##                 mutation_flag,
##                 plotting,
##                 verbose,
##                 progress,
##                 eval_fitness,
##                 eval_func_inputs)


## ----eval=TRUE, echo=FALSE, include=FALSE-------------------------------------------------------
load("knapsackSolution.RData")


## -----------------------------------------------------------------------------------------------
QGA:::plot_Output(knapsackSolution [[2]])


## -----------------------------------------------------------------------------------------------
best <- knapsackSolution[[1]] - 1
sum(best)


## -----------------------------------------------------------------------------------------------
sum(items$weight[best])
maxweight

