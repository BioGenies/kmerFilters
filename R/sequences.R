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
                              generate_sequence(sequence_length,
                                                alphabet, seqProbs))
        list_of_masks[[i]] <- attr(new_seq, "masks")
        sequences[i, ] <- new_seq
    }

    for (i in 1:(n_seq - n_pos)) {
        sequences[n_pos + i, ] <- generate_sequence(sequence_length,
                                                    alphabet,
                                                    seqProbs)
    }
    attr(sequences, "max_injection") <- max_injection
    attr(sequences, "motifs_map") <- motifs_map
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
    attr(test_res, "motifs_map") <- attr(test_dat, "motifs_map")
    attr(test_res, "masks") <- attr(test_dat, "masks")
    attr(test_res, "target") <- attr(test_dat, "target")
    test_res
}


#' Logistic regression with interactions
#'
#' This function samples target variable according to the logiistic model with
#' interactions
#'
#' @importFrom stats runif
#'
#' @param kmer_dat output of \code{\link{generate_kmer_data}}
#' @param zero_weight a single value denoting the weight of no-motifs case. If
#' \code{NULL}, then we sample the weight from the uniform distribution on the
#' [-2, -1] interval. Default to \code{NULL}.
#'
#' @return a binary vector of target variable sampled based on interaction model
#' and provided/calculated probabilities.
#'
#' @details
#' approach is based on logistic regression with interactions indicating that
#' the effect of one predictor depends on the value of another predictor. Let's
#' define maximum number of motifs per sequence \eqn{k = \max\lbrace k_i, i = 1,
#' \ldots, n\rbrace}. Let \eqn{w_{1}, \ldots, w_{k}} denote weights of single
#' effects. Namely:
#'
#' \eqn{g(EY) = w_0 + \sum_{i = 1}^{k} w_{i} X_{m_i} +
#' \left(\sum_{i = 1}^{k-1}\sum_{j = i + 1}^{k} w_{ij} X_{m_i}X_{m_j}\right) +
#' \ldots + w_{1\ldots k} X_{m_1}\ldots X_{m_k}}
#'
#' In the case when \code{probs} is \code{NULL} we calculate the probabilities
#' based on the formula \eqn{ exp(x_i)/(1 + exp(x_i))} where \eqn{x_i} denotes the
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
                                    zero_weight = NULL) {

    target <- attr(kmer_dat, "target")
    n_motifs <- length(attr(kmer_dat, "motifs_set"))
    motifs_map <- attr(kmer_dat, "motifs_map")

    # be careful, the following code is too dummy (but working though)
    dummy_df <- as.data.frame(motifs_map)
    n_seq <- nrow(dummy_df)
    dummier_names <- paste0(colnames(dummy_df), collapse = "*")

    formula <- paste0(" runif(", n_seq, ") ~ ", dummier_names)
    model <- lm(formula, data = dummy_df)
    interactions_matrix <- model.matrix(model)

    weights <- runif(ncol(interactions_matrix))

    if(is.null(zero_weight))
        zero_weight <- runif(1, -2, -1)

    target_weights <- interactions_matrix %*% weights

    target[target] <- target_weights[target]
    target[!target] <- zero_weight

    probs <- exp(target)/(1 + exp(target))

    rbinom_vec(probs)
}

#' Logistic regression response
#'
#' This function samples target variable according to the logistic model with
#' additive impact
#'
#' @inheritParams get_target_interactions
#'
#' @param weights a vector of weights of motifs' impact on the outcome. The
#' length of \code{weights} should be the same as the number of motifs provided
#' during sequences generation (it is the \code{motifs} parameter in the
#' \code{\link{generate_kmer_data}} function). If \code{weights} parameter is
#' \code{NULL}, then weights will be sampled from the uniform distribution on
#' 0-1 interval. The probability of success for target sampling will be
#' calculated based on the formula provided in details section. Default to
#' \code{NULL}.
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
    motifs_map <- attr(kmer_dat, "motifs_map")
    motifs_set <- attr(kmer_dat, "motifs_set")

    if(is.null(weights))
        weights <- runif(length(motifs_set))
    if(is.null(zero_weight))
        zero_weight <- runif(1, -2, -1)

    target_weights <- motifs_map %*% weights

    target[target] <- target_weights[target]
    target[!target] <- zero_weight

    probs <- exp(target)/(1 + exp(target))

    rbinom_vec(probs)
}



#' Logic regression noise
#'
#' This function samples target variable according to the logic regression
#' model (assuming that the occurrence of certain combinations of motifs affects
#' the feature). In the case of logical models, simulating a binary variable
#' involves defining logical conditions that determine the variable's value
#' based on motifs, e.g., the binary variable takes the value 1 if certain l
#' ogical criteria are met, and 0 if these criteria are not met.
#'
#' @inheritParams get_target_interactions
#'
#' @param weights a vector of weights of considered logic expression based on
#' available motifs. The length of \code{weights} should be the same as the
#' provided number  of expressions to use \code{n_exp}. If \code{weights}
#' parameter is \code{NULL}, then weights will be sampled from the uniform
#' distribution on 0-1 interval. The probability of success for target sampling
#' will be calculated based on the formula provided in details section. Default
#' to \code{NULL}.
#'
#' @param random a logical. Indicating whether expressions have to be generated
#' randomly. Default to \code{TRUE}.
#'
#' @param n_exp number of random logic expressions to create. It is used only
#' when \code{random} equals \code{TRUE}.
#'
#' @param max_exp_depth a maximum number of motifs used in a logic expression.
#' Default to 3.
#'
#' @param expressions a matrix of binary variables corresponding to custom
#' logic expressions. It's dimension should be related to the length of
#' \code{weights} vector if it's provided. Default to \code{NULL}.
#'
#' @details
#' Here, we consider new variables, \eqn{L_1, \ldots, L_l} where each of them
#' is a logic expression based on a subset of motifs \eqn{m_1, \ldots, m_m}. For
#' example,
#'
#' \eqn{L_1(m_1, m_2, m_3) = (X_{m_1} \land X_{m_2}) \lor X_{m_3}.}
#'
#' Each variable \eqn{L_i} obtains its own weight in the model. Our model is
#' following:
#'
#' \eqn{g(EY) = w_0 + \sum_{i = 1}^{l} w_i L_i.}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' results <- generate_kmer_data(n_seq, sequence_length, alph,
#'                               motifs, n_injections = 4)
#' get_target_logic(results)
#'
#' @export


get_target_logic <- function(kmer_dat,
                             random = TRUE,
                             zero_weight = NULL,
                             weights = NULL,
                             n_exp = NULL,
                             max_exp_depth = NULL,
                             expressions = NULL) {

    motifs_set <- attr(kmer_dat, "motifs_set")
    motifs_map <- attr(kmer_dat, "motifs_map")
    target <- attr(kmer_dat, "target")

    if((!is.null(max_exp_depth) && max_exp_depth < 2) | length(motifs_set) < 2)
        stop("You need at least 2 motifs to create a logic expression.")

    if(is.null(n_exp))
        n_exp <- min(length(motifs_set) - 1, 3)
    if(is.null(max_exp_depth))
        max_exp_depth <- min(length(motifs_set) - 1, 3)
    if(is.null(weights))
        weights <- runif(n_exp)
    if(is.null(zero_weight))
        zero_weight <- runif(1, -2, -1)

    if(!is.null(expressions)) {
        if(length(weights) != ncol(expressions))
            stop("You have to provide weight for each column
                 from expressions matrix!")
        else
            target_weights <- expressions %*% weights
    } else {

        expressions_motifs <- sapply(1:n_exp, function(ith_exp) {
            depth <- sample(2:max_exp_depth, 1)
            sample_motifs_ids <- sample(1:length(motifs_set), depth,
                                        replace = FALSE)
            operator <- sample(c("&", "|"), depth - 1, replace = TRUE)

            expr <- paste0(paste0(" motifs_map[, %i] ",
                                  sep = operator,
                                  collapse = ""), " motifs_map[, %i] ")
            expr_params <- c(list(expr), as.list(sample_motifs_ids))

            as.numeric(eval(parse(text = do.call(sprintf, expr_params))))
        })
    }

    target_weights <- expressions_motifs %*% weights

    target[target] <- target_weights[target]
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

