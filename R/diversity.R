#-----------------------------------------
# DIVERSITY: H-EPSILON GATE, RESTART,
# PERSONAL ATTRACTORS AND MIGRATION
#-----------------------------------------

#' H-epsilon gate
#'
#' @description
#' Prevents qubits from collapsing to a deterministic state (Han & Kim, 2004).
#' The probability of measuring 0 of every qubit is kept in
#' \[`epsilon`, 1 - `epsilon`\]: qubits outside this interval are moved to its
#' boundary, keeping the signs of the amplitudes and the norm of the state.
#'
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` with amplitudes.
#' @param epsilon Numeric in \[0, 0.5). Minimum probability of each bit value.
#'
#' @return Updated `q_alphabeta` amplitude array.
#' @keywords internal
h_epsilon <- function(q_alphabeta, epsilon) {
  alpha <- q_alphabeta[, 1, ]
  beta  <- q_alphabeta[, 2, ]
  r2 <- alpha^2 + beta^2
  p0 <- alpha^2 / r2
  p0_new <- pmin(pmax(p0, epsilon), 1 - epsilon)
  s_alpha <- ifelse(alpha < 0, -1, 1)
  s_beta  <- ifelse(beta  < 0, -1, 1)
  q_alphabeta[, 1, ] <- s_alpha * sqrt(p0_new * r2)
  q_alphabeta[, 2, ] <- s_beta  * sqrt((1 - p0_new) * r2)
  return(q_alphabeta)
}

#' Partial restart
#'
#' @description
#' Brings the population partially back towards the equal-probability state:
#' the angle of every qubit is moved towards pi/4 by the fraction `strength`
#' (0 = no change, 1 = full reset to equal probabilities), keeping the signs
#' of the amplitudes and the norm of the state. Individuals whose measured
#' chromosome equals the best solution are left unchanged.
#'
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` with amplitudes.
#' @param strength Numeric in \[0, 1\].
#' @param chromosome Integer matrix `[popsize, genomeLength]` of 0/1 bits.
#' @param solution_best Integer vector, the global best bitstring.
#'
#' @return Updated `q_alphabeta` amplitude array.
#' @keywords internal
partial_restart <- function(q_alphabeta, strength, chromosome, solution_best) {
  work_q_alphabeta <- q_alphabeta
  alpha <- q_alphabeta[, 1, ]
  beta  <- q_alphabeta[, 2, ]
  r   <- sqrt(alpha^2 + beta^2)
  phi <- atan2(abs(beta), abs(alpha))
  phi <- phi + strength * (pi / 4 - phi)
  work_q_alphabeta[, 1, ] <- ifelse(alpha < 0, -1, 1) * r * cos(phi)
  work_q_alphabeta[, 2, ] <- ifelse(beta  < 0, -1, 1) * r * sin(phi)
  for (i in seq_len(dim(q_alphabeta)[3])) {
    if (sum(chromosome[i, ] != solution_best) == 0) {
      work_q_alphabeta[, , i] <- q_alphabeta[, , i]
    }
  }
  return(work_q_alphabeta)
}

#' Migration of personal attractors
#'
#' @description
#' Personal attractors (Han & Kim, 2002): every individual is rotated towards
#' its own attractor, i.e. the best solution it has found so far. Every
#' `global_period` generations the best attractor is copied to all individuals
#' (global migration); otherwise, every `local_period` generations the best
#' attractor of each group of `group_size` consecutive individuals is copied
#' to the other members of the group (local migration).
#'
#' @param attractor Integer matrix `[popsize, genomeLength]` of attractor bitstrings.
#' @param attractor_fitness Numeric vector of length `popsize`.
#' @param generation Integer. Current generation.
#' @param group_size Integer. Size of the local migration groups.
#' @param local_period Integer or NULL. Period of local migration.
#' @param global_period Integer or NULL. Period of global migration.
#'
#' @return A list with the updated `attractor` and `attractor_fitness`.
#' @keywords internal
migration <- function(attractor, attractor_fitness, generation,
                      group_size, local_period, global_period) {
  popsize <- nrow(attractor)
  if (!is.null(global_period) && generation %% global_period == 0) {
    groups <- rep(1L, popsize)
  } else if (!is.null(local_period) && generation %% local_period == 0) {
    groups <- ceiling(seq_len(popsize) / group_size)
  } else {
    return(list(attractor = attractor, attractor_fitness = attractor_fitness))
  }
  for (g in unique(groups)) {
    members <- which(groups == g)
    b <- members[which.max(attractor_fitness[members])]
    attractor[members, ] <- matrix(attractor[b, ], length(members),
                                   ncol(attractor), byrow = TRUE)
    attractor_fitness[members] <- attractor_fitness[b]
  }
  return(list(attractor = attractor, attractor_fitness = attractor_fitness))
}
