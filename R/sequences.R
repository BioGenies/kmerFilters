#' Sampling from alphabet
#'
#' This function generates sequence of elements from alphabet with replacement
#'
#' @inheritParams generate_kmer_data
#'
#' @return randomly generated sequences
#'
#' @examples
#' generate_sequence(5, 1L:4)
#' generate_sequence(10, c("a", "b", "c"))
#' generate_sequence(10, c("a", "b", "c"), c(0.6, 0.2, 0.2))
#'
#' @export

generate_sequence <- function(sequence_length, alphabet, seqProbs = NULL){
    sample(alphabet, size = sequence_length, replace = TRUE, prob = seqProbs)
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
#' generate_sequence_data(n_seq, sequence_length, alph, motifs, 1)
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

    list_of_motifs <- list()
    list_of_masks <- list()
    list_of_motifs_ids <- list()

    max_injection <- n_injections
    n_injections <- sample(1:n_injections, n_pos, replace = TRUE)

    target <- logical(n_seq)
    target[1:n_pos] <- TRUE
    sequences <- matrix(nrow = n_seq, ncol = sequence_length)

    for (i in 1:n_pos) {

        motifs_ids <- sample(1:length(motifs), n_injections[i])
        selected_motifs <- motifs[motifs_ids]
        new_seq <- add_motifs(selected_motifs,
                              generate_sequence(sequence_length,
                                                alphabet, seqProbs))
        list_of_motifs_ids[[i]] <- motifs_ids
        list_of_motifs[[i]] <- attr(new_seq, "motifs")
        list_of_masks[[i]] <- attr(new_seq, "masks")
        sequences[i, ] <- new_seq
    }

    for (i in 1:(n_seq - n_pos)) {
        sequences[n_pos + i, ] <- generate_sequence(sequence_length,
                                                    alphabet,
                                                    seqProbs)
    }
    attr(sequences, "max_injection") <- max_injection
    attr(sequences, "motifs_ids") <- list_of_motifs_ids
    attr(sequences, "motifs") <- list_of_motifs
    attr(sequences, "masks") <- list_of_masks
    attr(sequences, "target") <- target
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

    attr(test_res, "sequences") <- matrix(test_dat,
                                          nrow = nrow(test_dat),
                                          ncol = ncol(test_dat))
    attr(test_res, "max_injection") <- attr(test_dat, "max_injection")
    attr(test_res, "motifs_set") <- motifs
    attr(test_res, "motifs_ids") <- attr(test_dat, "motifs_ids")
    attr(test_res, "motifs") <- attr(test_dat, "motifs")
    attr(test_res, "masks") <- attr(test_dat, "masks")
    attr(test_res, "target") <- attr(test_dat, "target")
    test_res
}


#' Interaction model noise
#'
#' This function samples target variable according to the binomial model with
#' interactions
#'
#' @importFrom stats runif
#'
#' @param kmer_dat output of \code{\link{generate_kmer_data}}
#' @param probs an increasing vector of probabilities of success corresponding
#' to concurrent occurrence of n motifs in a sequence, where n denotes a number
#' between \code{1} and \code{n_injections}. This vector should have exactly
#' \code{n_injections} of elements from 0-1 interval. For example, when
#' \code{n_injections} equals \code{2}, the vector of probabilities should have
#' two elements, for example, a vector \code{c(0.7, 0.8)} means that we assume
#' the probability of success equal to \code{0.7} when one motif occurs in a
#' sequence and \code{0.8} when two motifs occur. The default value is
#' \code{NULL} meaning that the probabilities will be calculated (See details).
#' @param zero_prob a single value denoting the probability of success in the
#' case when no motifs occur in the sequence. Default to \code{0.1}.
#'
#' @return a binary vector of target variable sampled based on interaction model
#' and provided/calculated probabilities.
#'
#' @details
#' This function assumes the following interaction binomial model:
#'
#' \eqn{g(EY) = w_0 + w_1 (X_{m_1} + X_{m_2} + \ldots + X_{m_k}) +
#' w_2 \left(\sum_{i = 1}^{k-1}\sum_{j = i + 1}^{k} X_{m_i}X_{m_j}\right) +
#' \ldots + w_m X_{m_1}\ldots X_{m_m}}
#'
#' In the case when \code{probs} is \code{NULL} we calculate the probabilities
#' based on the formula \eqn{ exp(x_i)/(1 + exp(x_i))} where xi denotes the
#' number of motifs in ith sequence.
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' results <- generate_kmer_data(n_seq, sequence_length, alph,
#'                               motifs, n_injections = 4)
#' get_target_interactions(results)
#'
#' @export

get_target_interactions <- function(kmer_dat,
                                    probs = NULL,
                                    zero_prob = 0.1) {

    if(length(probs) != attr(kmer_dat, "max_injection") & !is.null(probs))
        stop("The length of prob vector should equal the max_injection number!")
    if(!is.null(probs) & !all(probs <= 1 & probs >= 0))
        stop("The provided probabilities should be less than 1 and greated than 0.")
    if(is.unsorted(probs))
        stop("The vector of probabilities should be increasing.")

    motifs_counts <- lengths(attr(kmer_dat, "motifs"))
    target <- attr(kmer_dat, "target")
    target_probs <- target

    target_probs[target == 0] <- zero_prob

    if(is.null(probs))
        target_probs[target != 0] <- exp(motifs_counts)/(1 + exp(motifs_counts))
    else
        target_probs[target != 0] <- probs[motifs_counts]

    rbinom_vec(target_probs)
}

#' Additive model noise
#'
#' This function samples target variable according to the binomial model with
#' additive impact
#'
#' @param kmer_dat output of \code{\link{generate_kmer_data}}
#' @param weights a vector of weights of motifs' impact on the outcome. The
#' length of \code{weights} should be the same as the number of motifs provided
#' during sequences generation (it is the \code{motifs} parameter in the
#' \code{\link{generate_kmer_data}} function). If \code{weights} parameter is
#' \code{NULL}, then weights will be sampled from the uniform distribution on
#' 0-1 interval. The probability of success for target sampling will be
#' calculated based on the formula provided in details section. Default to
#' \code{NULL}.
#' @param zero_weight a single value denoting the weight of no-motifs case. If
#' \code{NULL}, then we sample the weight from the uniform distribution on the
#' [-2, -1] interval. Default to \code{NULL}.
#'
#' @return a binary vector of target variable sampled based on additive model.
#'
#' @details
#' This function assumes the following additive binomial model:
#'
#' \eqn{g(EY) = w_0 + w_1 X_{m_1} + w_2 X_{m_2} + \ldots + w_m X_{m_m}}
#'
#' where \eqn{w_1, \ldots, w_m} are weights related to motifs.
#'
#' In the case when \code{weights} is \code{NULL} we calculate the probabilities
#' based on the formula
#'
#' \eqn{ exp(1 + x_i)/(1 + exp(1 + x_i))} where xi denotes the sum of weights of
#' motifs occurring in ith sequence.
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' results <- generate_kmer_data(n_seq, sequence_length, alph,
#'                               motifs, n_injections = 4)
#' get_target_additive(results)
#'
#' @export

get_target_additive <- function(kmer_dat,
                                weights = NULL,
                                zero_weight = NULL) {

    if(length(weights) != length(attr(kmer_dat, "motifs_set")) & !is.null(weights))
        stop("The length of weights vector should equal number of motifs!")

    target <- attr(kmer_dat, "target")
    ids <- attr(kmer_dat, "motifs_ids")
    motifs_set <- attr(kmer_dat, "motifs_set")

    if(is.null(weights))
        weights <- runif(length(motifs_set))
    if(is.null(zero_weight))
        zero_weight <- runif(1, -2, -1)

    target_weights <- sapply(ids, function(ith_motifs) {
        sum(weights[ith_motifs])
    })

    target[target] <- target_weights
    target[!target] <- zero_weight

    probs <- exp(target)/(1 + exp(target))

    rbinom_vec(probs)
}


#' Binomial sampling
#'
#' This function samples from binomial distribution using a vector of
#' probabilities.
#'
#' @importFrom purrr map
#' @importFrom stats rbinom
#'
#' @param probs a vector of probabilities
#'
#' @return a binary vector of binomial observations corresponding to provided
#' probabilities.
#'
#' @examples
#' probs <- runif(100)
#' rbinom_vec(probs)
#'
#' @export

rbinom_vec <- function(probs) {
    if(any(probs >= 1 | probs <= 0))
        stop("Provided probabilities should be greater or equal to 0 and less
             or equal to 1!")

    as.integer(unlist(
        purrr::map(probs, function(ith_prob) rbinom(1, 1, prob = ith_prob)))
    )
}



