#--------------
# ROTATION                   
#--------------

#' Apply Rotation Update
#'
#' @description
#' Implements the Han–Kim lookup-table rule to rotate qubit amplitudes
#' towards the bit values of high-fitness solutions, based on the current
#' generation's best individual and the global best.
#'
#' @param chromosome Integer matrix `[popsize, genomeLength]` of 0/1 bits.
#' @param best_chromosome Integer vector; index of best individual per generation.
#' @param generation Integer. Current generation index (1-based).
#' @param genomeLength Integer. Total number of bits per chromosome.
#' @param solution_best Integer vector of length `genomeLength` representing the global best bitstring.
#' @param q_alphabeta Numeric array `[genomeLength, 2, popsize]` with amplitudes.
#' @param work_q_alphabeta Numeric workspace array like `q_alphabeta`.
#' @param popsize Integer. Population size.
#' @param fitness Numeric vector of length `popsize` with fitness values.
#' @param theta Numeric (radians). Rotation step size.
#'
#' @return Updated `q_alphabeta` amplitude array.
#' @keywords internal
rotation <- function(chromosome,
                     best_chromosome,
                     generation,
                     genomeLength,
                     solution_best,
                     q_alphabeta,
                     work_q_alphabeta,
                     popsize,
                     fitness, 
                     theta) {
  rot <- array(0, c(2, 2))
  for (i in c(1:popsize)) {
    if (sum(chromosome[i, ] != solution_best) != 0) {
      for (j in c(1:genomeLength)) {
        #----------------------
        # Han-Kim lookup table 
        #----------------------
        # f(x) > f(b) FALSE
        if (fitness[i] < fitness[best_chromosome[generation]]) {
          # x = 0 b = 1
          if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              s = 1
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              s = -1
            }
            if (q_alphabeta[j, 1, i] == 0) {
              s = 0
            }
            if (q_alphabeta[j, 2, i] == 0) {
              if (runif(1) < 0.5) s <- 1 else s <- -1
            }
            rot[1, 1] <- cos(s * theta)
            rot[1, 2] <- -sin(s * theta)
            rot[2, 1] <- sin(s * theta)
            rot[2, 2] <- cos(s * theta)
            work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
            work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
            q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
            q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          }
          # x = 1 b = 0
          if (chromosome[i, j] == 1 & chromosome[best_chromosome[generation], j] == 0) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              s = -1
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              s = 1
            }
            if (q_alphabeta[j, 1, i] == 0) {
              if (runif(1) < 0.5) s <- 1 else s <- -1
            }
            if (q_alphabeta[j, 2, i] == 0) {
              s = 0
            }
            rot[1, 1] <- cos(s * theta)
            rot[1, 2] <- -sin(s * theta)
            rot[2, 1] <- sin(s * theta)
            rot[2, 2] <- cos(s * theta)
            work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
            work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
            q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
            q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          }
        }
        # f(x) > f(b) TRUE
        if (fitness[i] >= fitness[best_chromosome[generation]]) {
          # x = 0 b = 1
          if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              s = -1
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              s = 1
            }
            if (q_alphabeta[j, 1, i] == 0) {
              if (runif(1) < 0.5) s <- 1 else s <- -1
            }
            if (q_alphabeta[j, 2, i] == 0) {
              s = 0
            }
            rot[1, 1] <- cos(s * theta)
            rot[1, 2] <- -sin(s * theta)
            rot[2, 1] <- sin(s * theta)
            rot[2, 2] <- cos(s * theta)
            work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
            work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
            q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
            q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          }
          # x = 1 b = 0
          if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              s = 1
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              s = -1
            }
            if (q_alphabeta[j, 1, i] == 0) {
              s = 0
            }
            if (q_alphabeta[j, 2, i] == 0) {
              if (runif(1) < 0.5) s <- 1 else s <- -1
            }
            rot[1, 1] <- cos(s * theta)
            rot[1, 2] <- -sin(s * theta)
            rot[2, 1] <- sin(s * theta)
            rot[2, 2] <- cos(s * theta)
            work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
            work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
            q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
            q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          }
          # x = 1 b = 1
          if (chromosome[i, j] == 1 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              s = 1
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              s = -1
            }
            if (q_alphabeta[j, 1, i] == 0) {
              s = 0
            }
            if (q_alphabeta[j, 2, i] == 0) {
              if (runif(1) < 0.5) s <- 1 else s <- -1
            }
            rot[1, 1] <- cos(s * theta)
            rot[1, 2] <- -sin(s * theta)
            rot[2, 1] <- sin(s * theta)
            rot[2, 2] <- cos(s * theta)
            work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
            work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
            q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
            q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          }
        }
      }
    }
  }
  return(q_alphabeta)
}
