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
#' @param binary logical, indicating whether the produced target variable should
#' be binary or continuous.
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
                                    zero_weight = NULL,
                                    binary = TRUE) {

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

    if(binary) {
        probs <- exp(target)/(1 + exp(target))
        rbinom_vec(probs)
    } else {
        target
    }
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
                                zero_weight = NULL,
                                binary = TRUE) {

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

    if(binary) {
        probs <- exp(target)/(1 + exp(target))
        rbinom_vec(probs)
    } else {
        target
    }
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
#' @param random a logical. Indicating whether expressions have to be generated
#' randomly. Default to \code{TRUE}.
#' @param n_exp number of random logic expressions to create. It is used only
#' when \code{random} equals \code{TRUE}.
#' @param max_exp_depth a maximum number of motifs used in a logic expression.
#' Default to 3.
#' @param expressions a matrix of binary variables corresponding to custom
#' logic expressions. You can create them based on motifs. It's dimension should
#' be related to the length of \code{weights} vector if it's provided. Default
#' to \code{NULL}. If \code{NULL}, random logic expressions will be created.
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
                             expressions = NULL,
                             binary = TRUE) {

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
        weights <- runif(ifelse(is.null(expressions), n_exp, ncol(expressions)))
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

        target_weights <- expressions_motifs %*% weights
    }

    target[target] <- target_weights[target]
    target[!target] <- zero_weight

    if(binary) {
        probs <- exp(target)/(1 + exp(target))
        rbinom_vec(probs)
    } else {
        target
    }
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
