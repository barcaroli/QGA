## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------
f <- function(x) {
  abs((x - 5) / (2 + sin(x)))
}
x <- seq(0,15,0.01)
y <- f(x)
plot(x,y)


## ----message=FALSE, warning=FALSE-------------------------------------------------------
# Analytic solution
if (!require(numDeriv)) install.packages("numDeriv", dependencies=TRUE)
library(numDeriv)
g <- function(x) {
  (x - 5) / (2 + sin(x))
}
f_prime <- function(x) {
  sgn <- sign(g(x))
  g_deriv <- grad(g, x)
  sgn * g_deriv
}
if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}
library(rootSolve)
root1 <- uniroot(f_prime, c(5.1, 12))$root
root2 <- uniroot(f_prime, c(12.1, 15))$root
# cat("Root 1: x ≈", root1, "\n")
# cat("Root 2: x ≈", root2, "\n")
f_root1 <- f(root1)
# f_root2 <- f(root2)
cat("f(", root1, ") ≈", f_root1, "\n")
# cat("f(", root2, ") ≈", f_root2, "\n")


## ---------------------------------------------------------------------------------------
functionMax <- function(solution, eval_func_inputs) {
  solution <- solution - 1
  # translate from binary to decimal value
  x=0
  for (j in c(1:Genome)) {
    x=x+solution[j]*2^(Genome-j)
  }
  # Normalize the obtained x into the interval 0-15
  x=x*16/(2^Genome-1)-1
  # replaces the value of x in the function f(x)
  y= abs((x-5)/(2+sin(x)))
  fitness=y
  # Store the current solution
  X <<- c(X,x)
  Y <<- c(Y,y)
  # cat("\nSolution=",solution,"  x=",x,"  y=",y)
  return(fitness)
}


## ---------------------------------------------------------------------------------------
popsize = 20
generation_max = 200
Genome = 16
nvalues_sol = 2
thetainit = 3.1415926535 * 0.15
thetaend = 3.1415926535 * 0.015
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(popsize + 1)
mutation_rate_end = 1/(popsize + 1)
mutation_flag = TRUE


## ---------------------------------------------------------------------------------------
library(QGA)
X <- NULL
Y <- NULL
set.seed(1234)
solutionQGA <- QGA(popsize,
                   generation_max,
                   nvalues_sol,
                   Genome=16,
                   thetainit,
                   thetaend,
                   pop_mutation_rate_init,
                   pop_mutation_rate_end,
                   mutation_rate_init,
                   mutation_rate_end,
                   mutation_flag = TRUE,
                   plotting = FALSE,
                   verbose = FALSE,
                   progress = FALSE,
                   eval_fitness = functionMax,
                   eval_func_inputs = list(X,Y))


## ---------------------------------------------------------------------------------------
solution <- solutionQGA[[1]]
solution <- solution - 1
solution
x=0
for (j in c(1:Genome)) {
  x=x+solution[j]*2^(Genome-j)
  # cat("\n",j," solution[j]:",solution[j],"  2^solution[Genome-j-1]:",2^(Genome-j))
}
x=x*16/(2^Genome-1)-1
y= abs((x-5)/(2+sin(x)))
cat("Function maximum: ",y," at x=",x)


## ---------------------------------------------------------------------------------------
cat("f(", root1, ") ≈", f_root1, "\n")


## ---------------------------------------------------------------------------------------
plot(X,Y,type="h",col="red")
title("f(x)=abs((x-5)/(2+sin(x)))")

