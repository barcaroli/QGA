#---------------------------
# FROM DECIMAL TO BINARY                    
#---------------------------  

#' Integer To Binary Vector
#'
#' @description
#' Converts a non-negative integer to a length-`n` binary vector using
#' big-endian ordering (index 1 is the most-significant bit).
#'
#' @param number Integer. The value to convert to binary.
#' @param n Integer. The number of bits to produce.
#'
#' @return Integer vector of length `n` containing 0/1 digits.
#' @examples
#' as.binary(5, 4)  # 0101
#' @keywords internal
as.binary <- function(number,n) {
  bin <- rep(NA,n)
  i = n
  for (i in c(n:1)) {
    digit <- number %% 2
    number <- floor(number / 2)
    bin[i] <- digit
  }
  return(bin)
}
