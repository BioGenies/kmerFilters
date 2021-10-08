#' generates sequence of elements from alphabet with replacement
#' @param alphabet elements used to build sequence
#' @param len length of generated sample sequence
#' @param seqProbs alphabet probabilites for sequences
#' @return randomly generated sequence
#' @export
#' @examples
#' generate_sequence(5, 1L:4)
#' generate_sequence(10, c("a", "b", "c"))
#' generate_sequence(10, c("a", "b", "c"), c(0.6, 0.2, 0.2))
generate_sequence <- function(len, alphabet, seqProbs = NULL){
  sample(alphabet, size = len, replace = TRUE, prob = seqProbs)
}

#' function generates sequences (both positive & negative)
#' @param n_seq number of sequences to be generated
#' @param len sequence length
#' @param alphabet elements used to build sequence
#' @param motifs_list list of injected motifs
#' @param n_motifs number of motifs injected to each positive sequence
#' @param fraction of positive sequences
#' @param seqProbs alphabet probabilites for sequences
#' @return generated sequences
#' @export
#' @examples
#' n_seq <- 20
#' len <- 10
#' alph <- 1L:4
#' motifs <- generate_motifs(alph, 2, 3, 0)
#' generate_sequence_data(n_seq, len, alph, motifs, 1)
#' generate_sequence_data(n_seq, len, alph, motifs, 1, fraction = 0.8)
#' generate_sequence_data(n_seq, len, alph, motifs, 1, seqProbs = c(0.7, 0.1, 0.1, 0.1))
generate_sequence_data <- function(n_seq,
                               len,
                               alphabet,
                               motifs_list,
                               n_motifs,
                               fraction = 0.5,
                               seqProbs = NULL,
                               sequence_length = 10) {
  
  n_pos <- round(fraction*n_seq, 0)
  
  list_of_motifs <- list()
  list_of_masks <- list()
  
  target <- logical(n_seq)
  target[1:n_pos] <- TRUE
  sequences <- matrix(nrow = n_seq, ncol = len)
  
  for (i in 1:n_pos) {
    
    correct_motifs <- FALSE
    motifs <- motifs_list[sample(1:length(motifs_list), n_motifs)]
    
    while (!validate_motifs(motifs, sequence_length)) {
      motifs <- motifs_list[sample(1:length(motifs_list), n_motifs)]
    } 
    
    new_seq <- add_motifs(motifs, generate_sequence(len, alphabet))
    list_of_motifs[[i]] <- attr(new_seq, "motifs")
    list_of_masks[[i]] <- attr(new_seq, "masks")
    sequences[i, ] <- new_seq
  }
  for (i in 1:(n_seq - n_pos)) {
    sequences[n_pos + i, ] <- generate_sequence(len, alphabet, seqProbs)
  }
  attr(sequences, "motifs") <- list_of_motifs
  attr(sequences, "masks") <- list_of_masks
  attr(sequences, "target") <- target
  sequences
}

#' wrapper for seqR counters
#' @importFrom seqR count_kmers count_multimers
#' @export
count_ngrams <- function(sequences, alphabet, n = 4, d = 6) {
  
  if (n == 1) {
    
    test_res <- count_kmers(sequences, 1, alphabet, with_kmer_counts = FALSE)
    
  } else {
    # element & gaps' positions for count_multigrams
    ns <- c()
    ds <- c()
    for (i in 1:(n-1)) {
      ds_ <- expand.grid(list(0:d)[rep(1, i)])
      ds_ <- ds_[apply(ds_, 1, sum) <= d, , drop = FALSE]
      ns <- c(ns, rep(i+1, nrow(ds_)))
      ds <- c(ds, split(ds_, 1:nrow(ds_)))
    }
    ds <- lapply(ds, unlist)
    
    test_res <- count_multimers(sequences,
                                c(1, ns),
                                alphabet,
                                kmer_gaps_list = c(list(c()), ds),
                                with_kmer_counts = FALSE)
    
  }
  
  test_res
}

#' function counts n-grams in given sequences
#' @param n_seq number of sequences to be generated
#' @param l_seq sequence length
#' @param alphabet elements used to build sequence
#' @param motifs_list list of injected motifs
#' @param n_motifs number of motifs injected to each positive sequence
#' @param fraction TODO: add fraction: of positive sequences / change approach
#' @param seqProbs alphabet probabilites for sequences
#' @param n maximum number of alphabet elements in n-gram
#' @param d maximum number of gaps in n-gram
#' @return generated sequences
#' @export
#' @examples
#' n_seq <- 20
#' len <- 1200
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 2, 3, 0)
#' results <- generate_sequences(n_seq, len, alph, motifs, 1)
#' results <- generate_sequences(n_seq, len, alph, motifs, 1, seqProbs = c(0.7, 0.1, 0.1, 0.1))
#' results
#' attributes(results)
generate_kmer_data <- function(n_seq,
                               l_seq,
                               alphabet,
                               motifs_list,
                               n_motifs,
                               fraction = 0.5,
                               seqProbs = NULL,
                               n = 4,
                               d = 6,
                               sequence_length = 10) {
  # generate sequence data
  test_dat <- generate_sequence_data(n_seq, 
                                     l_seq, 
                                     alphabet, 
                                     motifs_list, 
                                     n_motifs, 
                                     fraction, 
                                     seqProbs,
                                     sequence_length)
  
  test_res <- count_ngrams(test_dat, alphabet, n, d)
  
  attr(test_res, "sequences") <- matrix(test_dat, nrow = nrow(test_dat), ncol = ncol(test_dat))
  attr(test_res, "motifs") <- attr(test_dat, "motifs")
  attr(test_res, "masks") <- attr(test_dat, "masks")
  attr(test_res, "target") <- attr(test_dat, "target")
  test_res
}