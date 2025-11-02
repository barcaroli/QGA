#---------------------------
# POPULATION INITIALIZATION                     
#---------------------------

#' Initialize Quantum Population
#'
#' @description
#' Initializes the amplitudes for each qubit of every individual by
#' applying a random rotation to the Hadamard-initialized |0⟩ state.
#'
#' @param popsize Integer. Number of individuals in the population.
#' @param genomeLength Integer. Total number of bits per chromosome.
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` to fill with amplitudes.
#' @param rot Numeric 2x2 matrix used as rotation workspace.
#' @param theta Numeric. Input angle (overwritten internally by random angles).
#' @param h Numeric 2x2 Hadamard matrix.
#' @param qubit_0 Numeric length-2 vector representing the |0⟩ state.
#'
#' @return The updated `q_alphabeta` array of shape `[genomeLength, 2, popsize]`.
#' @keywords internal
generate_pop <- function(popsize,
                         genomeLength,
                         q_alphabeta,
                         rot,
                         theta,
                         h,
                         qubit_0) {
  for (i in c(1:popsize)) {
    for (j in c(1:genomeLength)) {
      theta <- runif(1) * 360
      theta <- pi*theta
      rot[1, 1] <- cos(theta)
      rot[1, 2] <- -sin(theta)
      rot[2, 1] <- sin(theta)
      rot[2, 2] <- cos(theta)
      q_alphabeta[j, 1, i] <- rot[1, 1] * h[1, 1] * qubit_0[1] + rot[1, 2] * h[1, 2] * qubit_0[2]
      q_alphabeta[j, 2, i] <- rot[2, 1] * h[2, 1] * qubit_0[1] + rot[2, 2] * h[2, 2] * qubit_0[2]
    }
  }
  return(q_alphabeta)
}
