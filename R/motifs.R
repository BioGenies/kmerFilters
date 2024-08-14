#' Random motifs generator
#'
#' This function generates random motif based on a given alphabet.
#' The maximum range of a motif equals `n + d`.
#'
#' @param alphabet elements used to generate a motif
#' @param n maximum number of alphabet elements
#' @param d number of possible gaps
#' @param motifProbs alphabet elements' probabilites
#'
#' @return motif built on a given alphabet
#'
#' @examples
#' generate_motif(1:4, n = 2, d = 0)
#' generate_motif(c("a", "b", "c"), n = 6, d = 1)
#' generate_motif(1:4, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))
#'
#' @export

generate_motif <- function(alphabet, n, d, motifProbs = NULL) {

    checkmate::assert_number(n)
    checkmate::assert_number(d)
    checkmate::assert(length(alphabet) == length(unique(alphabet)))

    if (!is.null(motifProbs)){
        checkmate::assert(sum(motifProbs) == 1)
        checkmate::assertNumeric(motifProbs)
        checkmate::assert(length(alphabet) == length(motifProbs))
    }

    # generate contiguous motif
    contiguous_motif <- sample(alphabet, sample(2:n, 1), replace = TRUE,
                               prob = motifProbs)
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

#' Motif to a sequence injection
#'
#' This function injects motifs to a sequence
#'
#' @param motifs list of motifs to be injected
#' @param sequence vector of alphabet elements
#'
#' @return list(sequence, motifs, masks)
#'
#' @examples
#' # simple injection
#' add_motifs(list(c(1, "_", 1), c(1, 1)), c(2, 2, 3, 4))
#' # little bit more interesting
#' alph <- as.character(1L:4)
#' motifs <- generate_motifs(alph, 2, 2, n = 4, d = 6)
#' example_sequence <- sample(alph, size = 10, replace = TRUE)
#' add_motifs(motifs, example_sequence)
#'
#' @export

add_motifs <- function(motifs, sequence) {
    sequence_len <- length(sequence)

    #create grid of possible motifs' positions
    maximum_motifs_positions <- lapply(motifs, function(x)
        seq(sequence_len - length(x) + 1))
    motifs_grid <- expand.grid(maximum_motifs_positions)
    motifs_grid <- motifs_grid[sample(1:nrow(motifs_grid)), , drop = FALSE]

    if(any(motifs_grid < 1)) {
        stop("Some of positions are invalid!")
    }

    for (i in 1:nrow(motifs_grid)) {
        list_of_masks <- list()
        injected_sequence <- sequence
        injected_positions <- logical(length(sequence))

        for (j in 1:ncol(motifs_grid)) {
            mask <- rep(FALSE, sequence_len)
            new_injected_sequence <- injected_sequence
            motif <- motifs[[j]]
            ids <- 0:(length(motif) - 1)
            ids <- ids[motif != "_"] + motifs_grid[i, j]
            mask[ids] <- TRUE
            new_injected_sequence[ids] <- motif[motif != "_"]

            if (j == 1) {
                injected_sequence <- new_injected_sequence
                injected_positions <- mask
            } else {
                if (all(injected_sequence[injected_positions] ==
                        new_injected_sequence[injected_positions])){
                    injected_sequence <- new_injected_sequence
                    injected_positions <- (injected_positions | mask)
                } else {
                    break
                }
            }

            list_of_masks[[j]] <- mask

            if (j == ncol(motifs_grid)){
                attr(injected_sequence, "motifs") <- motifs
                attr(injected_sequence, "masks") <- list_of_masks
                return(injected_sequence)
            }
        }
    }
    stop("Given motifs cannot be injected to a sequence!")
}

#' Validate if given set of motifs can occur in a sequence at the same time
#'
#' This function validates if given motifs can be injected to a sequence of
#' given length
#'
#' @param motifs list of motifs we are checking
#' @param sequence_length length of sequence we want to inject
#' @return logical value if such injection is possible
#'
#' @examples
#' set.seed(42)
#' motifs <- generate_motifs(1:4, n_motifs = 2, n_injections = 2, n = 3, d = 3)
#' validate_motifs(motifs, 7)
#' validate_motifs(motifs, 9)
#'
#' @export

validate_motifs <- function(motifs, sequence_length) {
    result <- tryCatch(add_motifs(motifs, rep("*", sequence_length)),
                       error = function(dummy) FALSE)
    ifelse(class(result) == "character", TRUE, FALSE)
}

#' Function generates list of motifs
#'
#' This function generates multiple motifs from alphabet
#'
#' @importFrom utils combn
#'
#' @inheritParams generate_motif
#' @param n_motifs number of motifs to generate
#' @param n_injections number of injections (for validation purposes:
#' checks if each subset of motifs of size `n_injections` can be injected
#' to a sequence of length `sequence_length`)
#' @param validate if true, returns a set of motifs that can be injected
#'  to a sequence of length 10
#' @param sequence_length length of a sequence that must contain all motifs
#'
#' @return list of generated motifs
#'
#' @examples
#' generate_motifs(1:4, 5, 3, n = 6, d = 6)
#' generate_motifs(1:4, 5, 3, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))
#'
#' @export

generate_motifs <- function(alphabet,
                            n_motifs,
                            n_injections,
                            n,
                            d,
                            motifProbs = NULL,
                            validate = TRUE,
                            sequence_length = 10) {

    if (!validate) {
        motifs <- lapply(1L:n_motifs, function(dummy)
            generate_motif(alphabet, n, d, motifProbs))
    }
    else {
        validated <- FALSE
        while (!validated) {

            motifs <- lapply(1L:n_motifs, function(dummy)
                generate_motif(alphabet, n, d, motifProbs))

            combns <- t(combn(1:n_motifs, n_injections))

            possible_injections <- lapply(1:nrow(combns), function(i) {
                validate_motifs(motifs[as.numeric(combns[i, ])],
                                sequence_length)
            })

            if (all(unlist(possible_injections))) {
                validated <- TRUE
            }
        }
    }

    motifs
}
