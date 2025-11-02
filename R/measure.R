#---------------------------
# MEASUREMENT                     
#---------------------------

#' Measure Quantum Chromosomes
#'
#' @description
#' Collapses quantum amplitudes to classical 0/1 bits for each individual
#' and position according to the probability mass in `q_alphabeta`.
#'
#' @param popsize Integer. Number of individuals.
#' @param genomeLength Integer. Total number of bits per chromosome.
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` with amplitudes.
#' @param chromosome Integer matrix `[popsize, genomeLength]` to fill with 0/1.
#'
#' @return The updated `chromosome` matrix of 0/1 values.
#' @keywords internal
measure <- function(popsize,
                    genomeLength,
                    q_alphabeta,
                    chromosome) {
  for (i in (1:popsize)) {
    for (j in (1:genomeLength)) {
      p_alpha <- runif(1)
      if (p_alpha <= 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 0
      if (p_alpha > 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 1
    }
  }
  return(chromosome)
}
