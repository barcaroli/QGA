# FUNCTION MAXIMUM
# (same example than in Lahoz-Beltra at
# https://github.com/ResearchCodesHub/QuantumGeneticAlgorithms/blob/master/QGA.py )
# Let f(x)=abs(x-5/2+sin(x)) be a function that takes   #
# values in the range 0<=x<=15. Within this range f(x)  #
# has a maximum value at x=11 (binary is equal to 1011) #

functionMax <- function(solution, eval_func_inputs) {
  solution <- solution - 1
  # translate from binary to decimal value
  x=0
  for (j in c(1:Genome)) {
    x=x+solution[j]*2^(Genome-j)
    # cat("\n",j," solution[j]:",solution[j],"  2^solution[Genome-j-1]:",2^(Genome-j))
  }
  x=x-1
  # replaces the value of x in the function f(x)
  y= abs((x-5)/(2+sin(x)))
  # the fitness value is calculated below:
  # (Note that in this example is multiplied
  # by a scale value, e.g. 100)
  fitness=y*100
  cat("\nSolution=",solution,"  x=",x,"  y=",y)
  return(fitness)
}
#----------------------
# Perform optimization
popsize = 20
Genome = 4
set.seed(1234)
solutionQGA <- QGA(popsize = 5,
                   generation_max = 50,
                   nvalues_sol = 2,
                   Genome,
                   thetainit = 3.1415926535 * 0.05,
                   thetaend = 3.1415926535 * 0.025,
                   pop_mutation_rate_init = 1/(popsize + 1),
                   pop_mutation_rate_end = 1/(popsize + 1),
                   mutation_rate_init = 1/(popsize + 1),
                   mutation_rate_end = 1/(popsize + 1),
                   mutation_flag = TRUE,
                   plotting = TRUE,
                   verbose = FALSE,
                   progress = TRUE,
                   eval_fitness = functionMax,
                   eval_func_inputs = NULL)
#----------------------
# Analyze results
solution <- solutionQGA[[1]]
solution <- solution - 1
solution
x=0
for (j in c(1:Genome)) {
  x=x+solution[j]*2^(Genome-j)
  # cat("\n",j," solution[j]:",solution[j],"  2^solution[Genome-j-1]:",2^(Genome-j))
}
y= abs((x-5)/(2+sin(x)))
cat("Function maximum: ",y," at x=",x)
