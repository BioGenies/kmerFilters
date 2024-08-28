
#' Simple motifs check
#'
#' This function check if at least one of provided motifs is contained in the
#' provided sequence.
#'
#' @inheritParams generate_kmer_data
#' @param sequence a vector of characters
#'
#' @return randomly generated sequences
#'
#' @examples
#' alph <- 1:4
#' motifs <- generate_motifs(alph, 3, 3, 3, 2)
#' sequence <- sample(alph,  100, replace = TRUE)
#' contains_motif(sequence, motifs)
#'
#' @export

contains_motif <- function(sequence, motifs) {

    sequence <- paste0(sequence, collapse = "")

    for(i in 1:length(motifs)) {
        ith_motif <- motifs[[i]]
        ith_motif[ith_motif == "_"] <- "."
        ith_motif <- paste0(ith_motif, collapse = "")

        find_motifs <- regexec(ith_motif, sequence)[[1]]

        if(find_motifs != -1)
            return(TRUE)
    }
    FALSE
}


#' Sampling from alphabet
#'
#' This function generates sequence of elements from alphabet with replacement
#'
#' @inheritParams generate_kmer_data
#'
#' @return randomly generated sequences
#'
#' @examples
#' set.seed(2)
#' alph <- 1:4
#' motifs <- generate_motifs(alph, 3, 3, 3, 2)
#' generate_negative_sequence(5, alph, motifs)
#' generate_negative_sequence(10, c("a", "b", "c"), motifs)
#' generate_negative_sequence(10, c("a", "b", "c"), c(0.6, 0.2, 0.2), motifs)
#'
#' @export

generate_negative_sequence <- function(sequence_length,
                                       alphabet,
                                       motifs,
                                       seqProbs = NULL,
                                       attempts = 100){

    candidate <- sample(alphabet, size = sequence_length,
                        replace = TRUE, prob = seqProbs)

    if(contains_motif(candidate, motifs)) {

        mask <- rep(FALSE, sequence_length)

        for(i in 1:length(motifs)) {
            ith_motif <- motifs[[i]]
            ith_motif[ith_motif == "_"] <- "."
            ith_motif <- paste0(ith_motif, collapse = "")

            found_motifs <- gregexpr(ith_motif, paste0(candidate, collapse = ""))[[1]]

            if(length(found_motifs) == 1 && found_motifs == -1) {
                next
            } else {
                motifs_ids <- as.vector(
                    sapply(found_motifs, function(x) x:(x + nchar(ith_motif) - 1))
                )
                mask[motifs_ids] <- TRUE
            }
        }

        attempt <- 1

        while (contains_motif(candidate, motifs) & attempt <= attempts) {
            attempt <- attempt + 1
            candidate[mask] <- sample(candidate[mask])
        }

        if(attempt == attempts)
            stop(paste0("It was impossible to generate negative sequence under
                 provided assumptions within ", attempts, " attempts."))
    }

    candidate
}

#' Sequences generation
#'
#' This function generates sequences (both positive & negative).
#'
#' @inheritParams generate_kmer_data
#'
#' @return generated sequences with motifs injected
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 10
#' alph <- 1L:4
#' motifs <- generate_motifs(alph, 3, 3, 3, 2)
#' generate_sequence_data(n_seq, sequence_length, alph, motifs, 3)
#'
#' @export

generate_sequence_data <- function(n_seq,
                                   sequence_length,
                                   alphabet,
                                   motifs,
                                   n_injections,
                                   fraction = 0.5,
                                   seqProbs = NULL
) {

    n_pos <- round(fraction * n_seq, 0)

    list_of_masks <- list()

    max_injection <- n_injections
    n_injections <- sample(1:n_injections, n_pos, replace = TRUE)

    target <- logical(n_seq)
    target[1:n_pos] <- TRUE
    motifs_map <- matrix(0, nrow = n_seq, ncol = length(motifs))
    sequences <- matrix(nrow = n_seq, ncol = sequence_length)

    for (i in 1:n_pos) {

        motifs_ids <- sample(1:length(motifs), n_injections[i])
        motifs_map[i, motifs_ids] <- 1
        selected_motifs <- motifs[motifs_ids]

        new_seq <- add_motifs(selected_motifs,
                              generate_negative_sequence(sequence_length,
                                                         alphabet, seqProbs))
        list_of_masks[[i]] <- attr(new_seq, "masks")
        sequences[i, ] <- new_seq
    }

    for (i in 1:(n_seq - n_pos)) {
        sequences[n_pos + i, ] <- generate_negative_sequence(sequence_length,
                                                             alphabet,
                                                             seqProbs)
    }
    attr(sequences, "max_injection") <- max_injection
    attr(sequences, "motifs_map") <- motifs_map
    attr(sequences, "masks") <- list_of_masks
    attr(sequences, "target") <- target
    attr(sequences, "motifs_set") <- motifs
    sequences
}

#' Counting kmers
#'
#' This function is a wrapper for seqR counters
#'
#' @importFrom seqR count_kmers count_multimers
#'
#' @inheritParams generate_motifs
#' @param sequences input data for k-mer counting
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 10
#' alph <- letters[1:6]
#' motifs <- generate_motifs(alph, 3, 3, 3, 2)
#' seq_data <- generate_sequence_data(n_seq, sequence_length, alph, motifs, 1)
#' count_seq_kmers(seq_data, alph)
#'
#' @export

count_seq_kmers <- function(sequences, alphabet, n = 4, d = 6) {

    sequences <- apply(sequences, 1, function(x) paste(x, collapse=""))

    if (n == 1) {
        test_res <- count_kmers(sequences, 1, alphabet,
                                with_kmer_counts = FALSE)
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

#' K-mer data generation
#'
#' This function generates sequences with provided motifs and constructs k-mer
#' representation table based on this data.
#'
#' @param n_seq number of sequences to be generated
#' @param sequence_length sequence length
#' @param alphabet elements used to build sequence
#' @param motifs list of motifs
#' @param n_injections maximal number of motifs injected to each positive
#' sequence (from 1 to \code{n_injections} will be injected)
#' @param fraction fraction of positive sequences
#' @param seqProbs alphabet probabilities for sequences
#' @param n maximum number of alphabet elements in n-gram
#' @param d maximum number of gaps in n-gram
#'
#' @return generated sequences
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 3, 3, 4, 6)
#' results <- generate_kmer_data(n_seq, sequence_length, alph, motifs, 1)
#'
#' @export

generate_kmer_data <- function(n_seq,
                               sequence_length,
                               alphabet,
                               motifs,
                               n_injections,
                               fraction = 0.5,
                               seqProbs = NULL,
                               n = 4,
                               d = 6
) {
    # generate sequence data
    test_dat <- generate_sequence_data(n_seq,
                                       sequence_length,
                                       alphabet,
                                       motifs,
                                       n_injections,
                                       fraction,
                                       seqProbs)

    test_res <- count_seq_kmers(test_dat, alphabet, n, d)

    test_res@Dimnames[[2]][1:length(alphabet)] <-
        paste0(test_res@Dimnames[[2]][1:length(alphabet)], "_0")

    attr(test_res, "sequences") <- matrix(test_dat,
                                          nrow = nrow(test_dat),
                                          ncol = ncol(test_dat))
    attr(test_res, "max_injection") <- attr(test_dat, "max_injection")
    attr(test_res, "motifs_set") <- attr(test_dat, "motifs_set")
    attr(test_res, "motifs_map") <- attr(test_dat, "motifs_map")
    attr(test_res, "masks") <- attr(test_dat, "masks")
    attr(test_res, "target") <- attr(test_dat, "target")
    test_res
}

