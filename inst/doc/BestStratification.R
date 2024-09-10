#--------------------------------------------------
# Best stratification of a sampling frame with QGA
#--------------------------------------------------

if (!require(SamplingStrata)) install.packages("SamplingStrata")
library(SamplingStrata)
library(QGA)


## -----------------------------------------------------------------------------------------------
# Sampling frame
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
head(frame)

## -----------------------------------------------------------------------------------------------
# Precision constraints
cv <- as.data.frame(list(
  DOM = "DOM1",
  CV1 = 0.03,
  CV2 = 0.03,
  CV3 = 0.03,
  CV4 = 0.03,
  domainvalue = 1
))


## -----------------------------------------------------------------------------------------------
nstrat = 3


## -----------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------
set.seed(1234)
solutionQGA <- QGA(popsize=20,
                generation_max = 2000,
                nvalues_sol = nstrat,
                Genome = nrow(iris),
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


## -----------------------------------------------------------------------------------------------
QGA:::plot_Output(solutionQGA[[2]])


## -----------------------------------------------------------------------------------------------
solution_QGA <- solutionQGA[[1]]
strata <- aggrStrata2(dataset = frame, 
                      vett = solution_QGA, 
                      dominio = 1)


## ----message=F, warning=FALSE-------------------------------------------------------------------
sum(bethel(strata, cv, realAllocation = TRUE))


## -----------------------------------------------------------------------------------------------
iris$stratum <- solution_QGA
table(iris$Species, iris$stratum)


## ----eval=FALSE---------------------------------------------------------------------------------
set.seed(1234)
solution_SamplingStrata <- optimStrata(method = "atomic",
                  framesamp = frame,
                  nStrata = nstrat,
                  errors = cv,
                  pops = popsize,
                  minnumstr = 1,
                  iter = 2000)


## -----------------------------------------------------------------------------------------------
sum(solution_SamplingStrata$aggr_strata$SOLUZ)


## -----------------------------------------------------------------------------------------------
iris$stratum <- solution_SamplingStrata$framenew$LABEL
table(iris$Species, iris$stratum)

