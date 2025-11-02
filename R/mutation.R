#----------
# MUTATION                   
#----------

#' Quantum Mutation Operator
#'
#' @description
#' With probability `pop_mutation_rate` per individual and `mutation_rate`
#' per bit, swaps the alpha/beta amplitudes to introduce diversity for
#' individuals that differ from the current global best.
#'
#' @param pop_mutation_rate Numeric in [0,1]. Per-individual mutation probability.
#' @param mutation_rate Numeric in [0,1]. Per-bit mutation probability.
#' @param popsize Integer. Population size.
#' @param chromosome Integer matrix `[popsize, genomeLength]` with 0/1 bits.
#' @param solution_best Integer vector (bitstring) of the global best.
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` with amplitudes.
#' @param work_q_alphabeta Numeric workspace array like `q_alphabeta`.
#' @param genomeLength Integer. Total number of bits per chromosome.
#'
#' @return Updated `q_alphabeta` amplitude array.
#' @keywords internal
mutation <- function(pop_mutation_rate, 
                     mutation_rate,
                     popsize,
                     chromosome,
                     solution_best,
                     q_alphabeta,
                     work_q_alphabeta,
                     genomeLength) {
  work_q_alphabeta <- q_alphabeta
  for (i in c(1:popsize)) {
    if (sum(chromosome[i, ] != solution_best) != 0) {
      rnd1 <- runif(1)
      if (rnd1 < pop_mutation_rate) {
        for (j in c(1:genomeLength)) {
          rnd2 <- runif(1)
          if (rnd2 < mutation_rate) {
            work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 2, i]
            work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 1, i]
          }
          if (rnd2 >= mutation_rate) {
            work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 1, i]
            work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 2, i]
          }
        }
      }
    }
  }
  q_alphabeta <- work_q_alphabeta
  return(q_alphabeta)
}
