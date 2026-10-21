#-----------------------------------------
# PARAMETERS HINT
#-----------------------------------------

#' Suggested parameter values for QGA
#'
#' @description
#' Suggests values for the parameters of [QGA()] that control convergence
#' (`thetainit`, `thetaend`, `restart_iters`, `restart_strength`,
#' `attractor_type`, ...) given the number of generations, the size of the
#' genome and the kind of problem. The suggestions follow the guidelines of
#' the vignette "Avoiding local optima":
#'
#' - a qubit always rotated in the same direction saturates after about
#'   (pi/4)/theta generations (the *saturation time*); the angles are chosen
#'   so that the saturation time is about `generation_max / 40` when restarts
#'   are used, and `generation_max / 10` otherwise, never faster than the
#'   default angles (about 6.7 generations);
#' - `restart_iters` is a multiple of the saturation time (1.5 times for
#'   deceptive problems, 3 times when the problem is unknown, 6 times for
#'   problems with many local optima), not a fraction of `generation_max`;
#' - `restart_strength = 0.5`; `restart_theta` and `epsilon` are left off.
#'
#' These are starting points, derived from benchmarks on two test problems:
#' always compare a few configurations on several seeds.
#'
#' @param generation_max Integer. Number of generations of the run.
#' @param Genome Integer. Number of genes (decision variables).
#' @param nvalues_sol Integer. Number of values of each gene.
#' @param popsize Integer. Population size. Default 20.
#' @param problem_type Character. Kind of problem: `"unknown"` (default),
#'   `"deceptive"` (good partial solutions lead away from the optimum),
#'   `"multimodal"` (many local optima scattered in the search space),
#'   `"unimodal"` (a single optimum, no traps).
#' @param stop_iters Logical. If TRUE, a value for `stop_iters` compatible with
#'   the restarts is suggested too. Default FALSE.
#' @param eval_seconds Numeric or NULL. Average time (seconds) of one fitness
#'   evaluation, used to estimate the duration of the run. Default NULL.
#' @param verbose Logical. If TRUE (default), the suggestions are printed with
#'   their motivation.
#'
#' @return Invisibly, a named list of parameter values that can be passed to
#'   [QGA()], e.g. `do.call(QGA, c(list(...), hint))`.
#'
#' @examples
#' hint <- parameters_hint(generation_max = 300, Genome = 48, nvalues_sol = 2,
#'                         problem_type = "deceptive")
#' # out <- do.call(QGA, c(list(popsize = 20, generation_max = 300,
#' #                            nvalues_sol = 2, Genome = 48,
#' #                            eval_fitness = f, eval_func_inputs = NULL),
#' #                       hint))
#'
#' @seealso [QGA()], `vignette("local_optima", package = "QGA")`
#' @export
parameters_hint <- function(generation_max,
                            Genome,
                            nvalues_sol,
                            popsize = 20,
                            problem_type = c("unknown", "deceptive",
                                             "multimodal", "unimodal"),
                            stop_iters = FALSE,
                            eval_seconds = NULL,
                            verbose = TRUE) {
  problem_type <- match.arg(problem_type)
  if (generation_max < 1 || Genome < 1 || nvalues_sol < 2)
    stop("generation_max and Genome must be >= 1, nvalues_sol >= 2")

  # genome encoding, as in QGA()
  geneLength   <- ceiling(log2(nvalues_sol))
  genomeLength <- Genome * geneLength

  # ---- Rotation angles from the target saturation time ----
  use_restart <- problem_type != "unimodal"
  t_sat <- if (use_restart) generation_max / 40 else generation_max / 10
  t_sat <- min(max(t_sat, 20 / 3), 100)   # never faster than the default angles
  theta_mean <- (pi / 4) / t_sat        # average angle over the run
  thetainit  <- theta_mean * 4 / 3      # thetainit = 2 * thetaend, as the defaults
  thetaend   <- theta_mean * 2 / 3

  # ---- Restart ----
  restart_mult  <- c(deceptive = 1.5, unknown = 3, multimodal = 6)
  restart_iters <- if (use_restart) max(5L, as.integer(round(restart_mult[[problem_type]] * t_sat))) else NULL

  # ---- Attractors ----
  attractor_type <- if (problem_type %in% c("unknown", "multimodal")) "personal" else "generation"

  hint <- list(thetainit = thetainit,
               thetaend = thetaend,
               restart_iters = restart_iters,
               restart_strength = 0.5,
               restart_theta = FALSE,
               epsilon = 0,
               attractor_type = attractor_type)
  if (attractor_type == "personal")
    hint$migration_global_period <- as.integer(min(50, max(5, round(generation_max / 6))))
  if (stop_iters)
    hint$stop_iters <- as.integer(if (use_restart) 3 * restart_iters else round(generation_max / 4))

  # ---- Warnings ----
  warn <- character(0)
  if (problem_type %in% c("deceptive", "unknown") && generation_max < 10 * restart_iters)
    warn <- c(warn, sprintf("generation_max allows few restarts: consider at least %d generations.",
                            10L * restart_iters))
  if (generation_max < 3 * t_sat + 50)
    warn <- c(warn, "generation_max is small: the population may not converge before the end.")
  if (popsize < 10)
    warn <- c(warn, "popsize < 10: the search may depend too much on the first measurements.")

  if (verbose) {
    cat("\n QGA parameters hint")
    cat("\n -------------------")
    cat(sprintf("\n problem_type: %s; genome of %d genes x %d bits = %d qubits; %d generations",
                problem_type, Genome, geneLength, genomeLength, generation_max))
    cat(sprintf("\n saturation time (pi/4)/theta: about %.1f generations", t_sat))
    cat(sprintf("\n\n thetainit        = %.4f   (= %.4f * pi)", thetainit, thetainit / pi))
    cat(sprintf("\n thetaend         = %.4f   (= %.4f * pi)", thetaend, thetaend / pi))
    if (use_restart) {
      cat(sprintf("\n restart_iters    = %d   (%.1f x saturation time)", restart_iters,
                  restart_mult[[problem_type]]))
      cat("\n restart_strength = 0.5")
    } else {
      cat("\n restart_iters    = NULL   (no traps expected: no restart)")
    }
    cat("\n restart_theta    = FALSE  (experimental)")
    cat("\n epsilon          = 0")
    cat(sprintf("\n attractor_type   = \"%s\"", attractor_type))
    if (!is.null(hint$migration_global_period))
      cat(sprintf("\n migration_global_period = %d", hint$migration_global_period))
    if (!is.null(hint$stop_iters))
      cat(sprintf("\n stop_iters       = %d", hint$stop_iters))
    n_eval <- popsize * (generation_max + 1)
    cat(sprintf("\n\n fitness evaluations: up to %d", n_eval))
    if (!is.null(eval_seconds)) {
      cat(sprintf(" (about %.1f minutes)", n_eval * eval_seconds / 60))
      if (eval_seconds > 0.05) cat("\n evaluation is expensive: consider parallel = TRUE")
    }
    if (length(warn) > 0) cat(paste0("\n Warning: ", warn), sep = "")
    cat("\n These are starting points: compare a few configurations on several seeds.\n")
  }
  if (length(warn) > 0 && !verbose) for (w in warn) warning(w, call. = FALSE)
  invisible(hint)
}
