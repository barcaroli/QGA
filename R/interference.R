#--------------
# INTERFERENCE
#--------------

#' Quantum interference operator
#'
#' @description
#' Applied to the whole population with probability `interference_rate`.
#' Two variants are available:
#'
#' - `"diagonal"`: quantum interference crossover (Narayanan & Moore, 1996).
#'   The amplitudes of the j-th qubit of individual i are copied from
#'   individual ((i + j - 2) mod popsize) + 1, so that each new individual is
#'   a "diagonal" combination of the whole population. Fitness is not used.
#'
#' - `"amplitude"`: fitness-weighted amplitude interference. Amplitudes are
#'   first brought to a canonical form (alpha, beta >= 0, which leaves the
#'   measurement probabilities unchanged), then each qubit receives the
#'   signed sum of the amplitudes of the whole population:
#'   \deqn{\psi_i \leftarrow \psi_i + \lambda \sum_k w_k \psi_k, \quad
#'   w_k = (f_k - \bar{f}) / \sum_k |f_k - \bar{f}|}
#'   Individuals above the average fitness contribute with positive sign
#'   (constructive interference), those below with negative sign
#'   (destructive interference). Amplitudes that would become negative are
#'   set to 0 (complete cancellation), and the qubit is renormalized to its
#'   original norm; if both amplitudes cancel, it is reset to the
#'   equal-probability state.
#'
#' Individuals whose measured chromosome equals the best solution are left
#' unchanged (as in rotation and mutation).
#'
#' @keywords internal
interference <- function(interference_rate,
                         popsize,
                         chromosome,
                         solution_best,
                         q_alphabeta,
                         genomeLength,
                         interference_type = "diagonal",
                         fitness = NULL,
                         interference_lambda = 0.5) {
  if (runif(1) >= interference_rate) return(q_alphabeta)
  if (interference_type == "diagonal") {
    work_q_alphabeta <- q_alphabeta
    for (j in c(1:genomeLength)) {
      src <- ((c(1:popsize) + j - 2) %% popsize) + 1
      work_q_alphabeta[j, , ] <- q_alphabeta[j, , src]
    }
  }
  if (interference_type == "amplitude") {
    dev <- fitness - mean(fitness)
    if (sum(abs(dev)) == 0) return(q_alphabeta)
    w <- dev / sum(abs(dev))
    alpha <- abs(q_alphabeta[, 1, ])                 # genomeLength x popsize
    beta  <- abs(q_alphabeta[, 2, ])
    norm  <- sqrt(alpha^2 + beta^2)
    # signed sum of the amplitudes of the whole population (one per qubit)
    field_alpha <- as.vector(alpha %*% w)
    field_beta  <- as.vector(beta  %*% w)
    new_alpha <- pmax(alpha + interference_lambda * field_alpha, 0)
    new_beta  <- pmax(beta  + interference_lambda * field_beta,  0)
    new_norm  <- sqrt(new_alpha^2 + new_beta^2)
    zero <- new_norm == 0
    new_alpha[zero] <- 1
    new_beta[zero]  <- 1
    new_norm[zero]  <- sqrt(2)
    work_q_alphabeta <- q_alphabeta
    work_q_alphabeta[, 1, ] <- new_alpha / new_norm * norm
    work_q_alphabeta[, 2, ] <- new_beta  / new_norm * norm
  }
  for (i in c(1:popsize)) {
    if (sum(chromosome[i, ] != solution_best) == 0) {
      work_q_alphabeta[, , i] <- q_alphabeta[, , i]
    }
  }
  return(work_q_alphabeta)
}
