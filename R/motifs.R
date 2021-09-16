#' Random motifs generator
#' 
#' Function generates random motif based on a given alphabet.
#' The maximum range of a motif equals `n + d`. 
#' @param alphabet elements used to generate a motif
#' @param n maximum number of alphabet elements
#' @param d number of possible gaps
#' @param motifProbs alphabet elements' probabilites
#' @return motif built on a given alphabet
#' @export
#' @examples
#' generate_motif(1:4, n = 2, d = 0)
#' generate_motif(c("a", "b", "c"), n = 6, d = 1)
#' generate_motif(1:4, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))
generate_motif <- function(alphabet, n, d, motifProbs = NULL) {
  
  checkmate::assert_number(n)
  checkmate::assert_number(d)
  checkmate::assert(length(alphabet) == length(unique(alphabet)),
         sum(motifProbs) == 1)
  
  if (!is.null(motifProbs)){
    checkmate::assertNumeric(motifProbs)
    checkmate::assert(length(alphabet) == length(motifProbs))
  }
  
  # generate contiguous motif 
  contiguous_motif <- sample(alphabet, sample(2:n, 1), replace = TRUE, prob = motifProbs)
  motif <- contiguous_motif
  
  if (d > 0) {
    # generate vector of gaps
    gaps <- expand.grid(list(0:d)[rep(1, length(contiguous_motif) - 1)])
    possibleGaps <- gaps[apply(gaps, 1, sum) <= d, , drop = FALSE]
    gap <- possibleGaps[sample(1:nrow(possibleGaps), 1), , drop = FALSE]
    
    # merge motif with gaps
    motif <- vector()
    for (i in 1:(length(contiguous_motif) -1)) {
      motif <- c(motif, contiguous_motif[i], rep("_", gap[1, i]))
    }
    motif <- c(motif, contiguous_motif[length(contiguous_motif)])
    
  }
  
  motif  
}
