QGA_par <- function(
    popsize = 20,
    generation_max = 200,
    nvalues_sol,
    Genome,
    thetainit = pi * 0.05,
    thetaend  = pi * 0.025,
    pop_mutation_rate_init = NULL,
    pop_mutation_rate_end  = NULL,
    mutation_rate_init     = NULL,
    mutation_rate_end      = NULL,
    mutation_flag = TRUE,
    plotting = TRUE,
    verbose  = TRUE,
    progress = TRUE,
    eval_fitness,
    eval_func_inputs,
    stop_limit = NULL,
    patience   = NULL,                          # early-stopping
    cores      = max(1L, parallel::detectCores() - 1L),
    chunking   = TRUE,                           # valutazione a blocchi
    dt_threads = 1L,                             # thread data.table dentro i worker
    pkgs       = c("data.table"),                # pacchetti da caricare nei worker
    clamp_sol  = TRUE                            # forza i geni nel range [1..nvalues_sol]
) {
  
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome))      stop("Genome parameter value missing!")
  if (is.null(stop_limit))  stop_limit <- Inf
  
  if (is.null(pop_mutation_rate_init) && mutation_flag) pop_mutation_rate_init <- 1/(popsize + 1)
  if (is.null(pop_mutation_rate_end)  && mutation_flag) pop_mutation_rate_end  <- 1/(popsize + 1)
  if (is.null(mutation_rate_init)     && mutation_flag) mutation_rate_init     <- 1/(Genome + 1)
  if (is.null(mutation_rate_end)      && mutation_flag) mutation_rate_end      <- 1/(Genome + 1)
  
  # geneLength = ceil(log2(nvalues_sol))
  n <- 0L; while (nvalues_sol > 2^n) n <- n + 1L
  geneLength   <- n
  genomeLength <- Genome * geneLength
  if (geneLength <= 0L || genomeLength <= 0L) stop("Invalid geneLength/genomeLength.")
  
  # -- util: decodifica cromosoma -> soluzione [1..nvalues_sol]
  decode_solution <- function(bits, geneLength, Genome, nvalues_sol, clamp_sol) {
    out <- integer(Genome)
    idx <- 1L
    for (x in seq_len(Genome)) {
      acc <- 0L
      for (y in seq_len(geneLength)) {
        b <- bits[idx]; if (b != 0L && b != 1L) b <- as.integer(b != 0)
        acc <- acc + b * 2L^(geneLength - y)
        idx <- idx + 1L
      }
      val <- acc + 1L
      if (clamp_sol) {
        if (val < 1L)           val <- 1L
        if (val > nvalues_sol)  val <- ((val - 1L) %% nvalues_sol) + 1L
      }
      out[x] <- val
    }
    out
  }
  
  # -- parallel setup ----------------------------------------------------------
  use_parallel <- isTRUE(cores > 1L)
  if (use_parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' required.")
    if (!requireNamespace("foreach", quietly = TRUE))    stop("Package 'foreach' required.")
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
    
    # carica pacchetti e setta threads nei worker
    if (length(pkgs)) {
      parallel::clusterExport(cl, varlist = c("pkgs"), envir = environment())
      parallel::clusterEvalQ(cl, {
        invisible(lapply(pkgs, function(p) { suppressWarnings(require(p, character.only=TRUE, quietly=TRUE)); NULL }))
        NULL
      })
    }
    if (!is.null(dt_threads) && is.finite(dt_threads)) {
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("data.table", quietly = TRUE)) {
          try(data.table::setDTthreads(as.integer(dt_threads)), silent = TRUE)
        }
        NULL
      })
    }
    
    # esporta oggetti/stato invarianti
    parallel::clusterExport(
      cl,
      varlist = c("decode_solution","eval_fitness","eval_func_inputs",
                  "geneLength","Genome","nvalues_sol","clamp_sol"),
      envir = environment()
    )
    
    `%op%` <- `%dopar%`
    nworkers <- foreach::getDoParWorkers()
  } else {
    `%op%` <- `%do%`
    nworkers <- 1L
  }
  
  # split helper
  split_chunks <- function(n, k) {
    if (!chunking || k <= 1L) return(list(seq_len(n)))
    parallel::splitIndices(n, min(k, n))
  }
  
  # -- stato -------------------------------------------------------------------
  q_alphabeta      <- array(0, c(genomeLength, 2L, popsize))
  work_q_alphabeta <- array(0, c(genomeLength, 2L, popsize))
  chromosome       <- array(0L, c(popsize, genomeLength))
  h        <- array(c(1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)), c(2L, 2L))
  rot      <- array(0, c(2L, 2L))
  qubit_0  <- array(c(1, 0), c(2L, 1L))
  
  res <- data.frame(
    generation = seq_len(generation_max + 1L),
    fitness_average = numeric(generation_max + 1L),
    fitness_best    = numeric(generation_max + 1L)
  )
  best_idx_by_gen <- integer(generation_max + 1L)
  
  # -- init popolazione --------------------------------------------------------
  generation <- 1L
  theta <- thetainit
  
  q_alphabeta <- QGA:::generate_pop(popsize, genomeLength, q_alphabeta, rot, theta, h, qubit_0)
  chromosome  <- QGA:::measure(popsize, genomeLength, q_alphabeta, chromosome)
  chromosome  <- QGA:::repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
  
  # -- valutazione popolazione (con chunking + diagnostica) --------------------
  eval_population <- function(chromosome) {
    idx_chunks <- split_chunks(popsize, nworkers)
    vals <- foreach::foreach(ch = idx_chunks, .combine = c) %op% {
      v <- numeric(length(ch))
      j <- 1L
      for (i in ch) {
        bits <- chromosome[i, ]
        sol  <- decode_solution(bits, geneLength, Genome, nvalues_sol, clamp_sol)
        # diagnostica "strict mode"
        if (!clamp_sol) {
          if (any(sol < 1L | sol > nvalues_sol)) {
            rng <- range(sol)
            stop(sprintf("Decoded solution out-of-range per individuo %d: min=%d max=%d (nvalues_sol=%d)",
                         i, rng[1], rng[2], nvalues_sol))
          }
        }
        # protezione aggiuntiva con messaggio sorgente
        v[j] <- tryCatch(
          eval_fitness(sol, eval_func_inputs),
          error = function(e) {
            stop(sprintf("Errore in eval_fitness per individuo %d: %s\nSol=%s",
                         i, conditionMessage(e), paste(sol, collapse=",")))
          }
        )
        j <- j + 1L
      }
      v
    }
    unlist(vals, use.names = FALSE)
  }
  
  # -- fitness iniziale --------------------------------------------------------
  fitness <- eval_population(chromosome)
  fitness_max     <- max(fitness)
  fitness_average <- mean(fitness)
  best_index      <- which.max(fitness)
  best_idx_by_gen[generation] <- best_index
  best_chromosome <- chromosome[best_index, ]
  solution_best   <- best_chromosome
  fitness_best    <- fitness_max
  
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation]    <- fitness_best
  
  if (plotting) QGA:::plot_Output(res[seq_len(generation), ])
  if (verbose)  cat("\n", generation, ",", fitness_average, ",", fitness_max)
  if (progress) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  
  stall <- 0L
  
  # -- loop evolutivo ----------------------------------------------------------
  while (generation <= generation_max && fitness_max < stop_limit) {
    if (progress) setTxtProgressBar(pb, generation)
    
    theta <- thetainit - ((thetainit - thetaend) / generation_max) * generation
    if (theta < 0) theta <- 0
    
    q_alphabeta <- QGA:::rotation(
      chromosome,
      best_idx_by_gen,     # atteso da rotation: indici best per generazione
      generation,
      genomeLength,
      solution_best,
      q_alphabeta,
      work_q_alphabeta,
      popsize,
      fitness,
      theta
    )
    
    generation <- generation + 1L
    
    if (mutation_flag) {
      pop_mutation_rate <- pop_mutation_rate_init - ((pop_mutation_rate_init - pop_mutation_rate_end) / generation_max) * generation
      mutation_rate     <- mutation_rate_init     - ((mutation_rate_init     - mutation_rate_end)     / generation_max) * generation
      q_alphabeta <- QGA:::mutation(
        pop_mutation_rate, mutation_rate, popsize, chromosome,
        solution_best, q_alphabeta, work_q_alphabeta, genomeLength
      )
    }
    
    chromosome <- QGA:::measure(popsize, genomeLength, q_alphabeta, chromosome)
    chromosome <- QGA:::repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
    
    fitness <- eval_population(chromosome)
    fitness_max     <- max(fitness)
    fitness_average <- mean(fitness)
    best_index      <- which.max(fitness)
    best_idx_by_gen[generation] <- best_index
    best_chromosome <- chromosome[best_index, ]
    
    if (fitness_max > fitness_best) {
      fitness_best  <- fitness_max
      solution_best <- best_chromosome
      stall <- 0L
    } else {
      stall <- stall + 1L
    }
    
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation]    <- fitness_best
    
    if (plotting) QGA:::plot_Output(res[seq_len(generation), ])
    if (verbose)  cat("\n", generation, ",", fitness_average, ",", fitness_best)
    
    if (!is.null(patience) && patience > 0L && stall >= patience) {
      if (verbose) cat("\nEarly stopping: nessun miglioramento per", patience, "generazioni.")
      break
    }
  }
  
  if (progress) close(pb)
  cat("\n *** Best fitness: ", fitness_best)
  
  last_row <- max(1L, min(nrow(res), generation))
  if (plotting) QGA:::plot_Output(res[seq_len(last_row), ])
  
  # decode best globale
  bits <- solution_best
  solution <- integer(Genome)
  idx <- 1L
  for (x in seq_len(Genome)) {
    acc <- 0L
    for (y in seq_len(geneLength)) {
      b <- bits[idx]; if (b != 0L && b != 1L) b <- as.integer(b != 0)
      acc <- acc + b * 2L^(geneLength - y)
      idx <- idx + 1L
    }
    val <- acc + 1L
    if (clamp_sol) {
      if (val < 1L)          val <- 1L
      if (val > nvalues_sol) val <- ((val - 1L) %% nvalues_sol) + 1L
    }
    solution[x] <- val
  }
  
  list(solution = solution, history = res[seq_len(last_row), ])
}
