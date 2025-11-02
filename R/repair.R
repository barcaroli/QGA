#' Repair Decoded Genes
#'
#' @description
#' Ensures that each decoded gene value falls within 1..`nvalues_sol` by
#' folding out-of-range values back into the valid domain and re-encoding
#' the chromosome to bits.
#'
#' @param popsize Integer. Number of individuals.
#' @param chromosome Integer matrix `[popsize, genomeLength]` with 0/1 bits.
#' @param geneLength Integer. Number of bits per gene.
#' @param genomeLength Integer. Total number of bits (`Genome * geneLength`).
#' @param nvalues_sol Integer. Allowed value count per gene (1..`nvalues_sol`).
#' @param Genome Integer. Number of genes.
#'
#' @return The repaired `chromosome` bit matrix.
#' @keywords internal
repair <- function(popsize,
                   chromosome,
                   geneLength,
                   genomeLength,
                   nvalues_sol,
                   Genome) {
  diff = 2^geneLength - nvalues_sol
  acceptable_values <- c(1:(2^geneLength - diff))
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(geneLength,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:geneLength)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(geneLength - y) 
      }
    }
    solution <- solution + 1
    solution[order(solution)]
    table(solution)
    if (max(solution) > nvalues_sol) { 
      solution[!(solution %in% acceptable_values)] <- solution[!(solution %in% acceptable_values)] - diff
    }
    a = array(c(1:genomeLength),c(geneLength,Genome))
    for (x in c(1:Genome)) {
      y1 = a[1,x]
      y2 = a[geneLength,x]
      chromosome[i,c(y1:y2)] <- as.binary(solution[x]-1,n=geneLength)
    }
  }  
  return(chromosome)
}
