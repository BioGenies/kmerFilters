
#' QuiPT - Quick Permutation Test
#'
#' This function is a wrapper over the QuiPT test.
#'
#' @importFrom biogram test_features
#'
#' @param target a numeric response variable
#' @param kmers a matrix of kmers with named columns or an object obtained via
#' \code{\link{generate_kmer_data}} function.
#' @param ... other arguments for \code{\link[biogram]{test_features}} function.
#'
#' @return a numeric vector of named p-values corresponding to k-mers in the
#' feature space.
#'
#' @details
#' This function uses \code{\link[biogram]{test_features}}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_quipt(target, kmers)
#'
#' @export

filter_quipt <- function(target, kmers, ...) {

    pvals <- biogram::test_features(target, kmers, ...)

    pvals
}


#' Chi-squared independence test
#'
#' This function is a wrapper over chi-squared independence test for k-mer data.
#'
#' @importFrom stats chisq.test
#'
#' @inheritParams filter_quipt
#'
#' @return a numeric vector of named p-values corresponding to k-mers in the
#' feature space.
#'
#' @details
#' This function uses \code{\link[stats]{chisq.test}}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_chisq(target, kmers)
#'
#' @export

filter_chisq <- function(target, kmers) {

    pvals <- unlist(lapply(1:ncol(kmers), function(i) {

        ith_kmer <- as.vector(kmers[, i])
        suppressWarnings({
            if(length(unique(ith_kmer)) == 1) pval <- 1
            else pval <- chisq.test(ith_kmer, target, B = 200)$p.value
        })
        pval
    }))

    names(pvals) <- colnames(kmers)

    pvals
}


#' Fast Correlation Based Filter
#'
#' This function is a wrapper over Fast Correlation Based Filter for k-mer data.
#'
#' @importFrom FCBF fcbf
#'
#' @inheritParams filter_quipt
#' @param thresh a threshold for symmetrical uncertainty between a k-mer and
#' a target variable. Default to \code{0.25}.
#'
#' @return a character vector of names of selected kmers
#'
#' @details
#' This function uses \code{\link[FCBF]{fcbf}}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_fcbf(target, kmers)
#'
#' @export


filter_fcbf <- function(target, kmers, thresh = 0.25) {

    features <- apply(kmers, 1, as.factor)
    kmers_names <- rownames(features)

    suppressMessages({
        fcbf_results <- fcbf(features, target, minimum_su = thresh)
    })

    rownames(fcbf_results)
}


#' Entropy-Based Feature Selection Algorithms
#'
#' This function is a wrapper over Entropy-Based Feature Selection Algorithms
#' for k-mer data.
#'
#' @importFrom FSelectorRcpp information_gain
#'
#' @inheritParams filter_quipt
#'
#' @param method a character name of a filter type. One of "infogain",
#' "gainratio" or "symuncert". For more details see
#' \code{\link[FSelectorRcpp]{information_gain}}.
#'
#' @return a numeric vector of named p-values corresponding to k-mers in the
#' feature space.
#'
#' @details
#' This function uses \code{\link[FSelectorRcpp]{information_gain}}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_ig(target, kmers, "infogain", 0.01)
#'
#' @export

filter_ig <- function(target, kmers, method) {

    method <- match.arg(method, c("infogain", "gainratio", "symuncert"))

    res <- information_gain(x = kmers, y = target, discIntegers = FALSE,
                            type = "gainratio")

    scores <- res[["importance"]]

    names(scores) <- res[["attributes"]]

    scores
}


#' Entropy-Based Feature Selection Algorithms
#'
#' This function is a wrapper over Entropy-Based Feature Selection Algorithms
#' for k-mer data.
#'
#' @importFrom praznik MIM
#' @importFrom praznik MRMR
#' @importFrom praznik JMI
#' @importFrom praznik JMIM
#' @importFrom praznik DISR
#' @importFrom praznik NJMIM
#' @importFrom praznik CMIM
#' @importFrom praznik CMI
#'
#' @inheritParams filter_quipt
#'
#' @param method a character name of a filter type. One of "MIM", "MRMR", "JMI",
#' "JMIM", "DISR", "NJMIM", "CMIM", "CMI". See details for more information.
#'
#' @return a character vector of names of selected kmers
#'
#' @details
#' This function uses
#' - \code{\link[praznik]{MIM}}
#' - \code{\link[praznik]{MRMR}}
#' - \code{\link[praznik]{JMI}}
#' - \code{\link[praznik]{JMIM}}
#' - \code{\link[praznik]{DISR}}
#' - \code{\link[praznik]{NJMIM}}
#' - \code{\link[praznik]{CMIM}}
#' - \code{\link[praznik]{CMI}}
#'
#' @examples
#' n_seq <- 20
#' sequence_length <- 20
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_praznik(target, kmers, "infogain", 0.01)
#'
#' @export

filter_praznik <- function(target, kmers, method) {

    praznik_methods <-
        c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")

    method <- match.arg(method, praznik_methods)

    res <-  get(method)(X = as.data.frame(as.matrix(kmers)), Y = target,
                        k = ncol(kmers))
    res[["score"]]
}



#' Stepwise information criteria based filtering
#'
#' This function is filtering method based on Information Criteria
#'
#' @importFrom bigstep stepwise
#' @importFrom bigstep prepare_data
#' @importFrom bigstep aic
#' @importFrom bigstep maic
#' @importFrom bigstep maic2
#' @importFrom bigstep bic
#' @importFrom bigstep mbic
#' @importFrom bigstep mbic2
#' @importFrom Matrix colSums
#'
#' @inheritParams filter_quipt
#'
#' @param ic character name of information criterium. One of "aic", "maic",
#' "bic", "maic2", "mbic", "mbic2". See bigstep package for more information.
#' @param reduce a numeric value from (0, 1) interval. Denotes significance
#' level for preliminary reduction before execution of stepwise procedure.
#' Default to 0.2.
#' @param attach_correlated a logical value indicating whether the highly
#' correlated k-mers should be chosen.
#' @param threshold a numeric threshold from 0 to 1 denoting a threshold for
#' correlation coefficient when \code{attach_correlated} is TRUE. Ignored when
#' \code{attach_correlated} is FALSE.
#'
#' @return a character vector of names of selected kmers
#'
#' @details
#' This function uses bigstep package.
#'
#' @examples
#' n_seq <- 10
#' sequence_length <- 10
#' alph <- letters[1:20]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_ic(target, kmers, "mbic2", 0.2)
#'
#' @export

filter_ic <- function(target, kmers, ic = "mbic2", reduce = 0.15,
                      attach_correlated = TRUE, threshold = 0.9) {

    ic <- match.arg(ic, c("aic", "maic", "bic", "maic2", "mbic", "mbic2"))

    dat <- prepare_data(target, as.matrix(kmers))
    dat <- reduce_matrix(dat, minpv = reduce)

    res_ic <- stepwise(dat, crit = get(ic))[["model"]]

    chosen_kmers <- res_ic

    if(attach_correlated) {
        n <- nrow(kmers)
        candidates_kmers <- kmers[, dat[["candidates"]]]
        candidates_sums <- Matrix::colSums(candidates_kmers)

        for(y_name in res_ic) {
            xy <- Matrix::colSums(candidates_kmers * candidates_kmers[, y_name])
            correlations <- sparse_cor(y_name, xy, candidates_sums, n)
            chosen_kmers <- c(chosen_kmers, names(which(correlations >= threshold)))
        }
    }

    unique(chosen_kmers)
}

#' Sparse correlation
#'
#' This function calculates correlation between one k-mer and a matrix of k-mers
#' under sparsity assumption
#'
#' @param y_name a character name of reference k-mer
#' @param xy a vector of sums of xy products (between y and x(s))
#' @param candidates_sums a vector of sums of k-mers
#' @param n a numeric number of observations
#'
#' @return a numeric vector of correlations between y and x's.
#'

sparse_cor <- function(y_name, xy, candidates_sums, n) {
    nom <- n * xy - candidates_sums * candidates_sums[y_name]
    denom <- sqrt(n * candidates_sums - candidates_sums^2) *
        sqrt(n * candidates_sums[y_name] - candidates_sums[y_name]^2)

    nom / denom
}





#' Get a list of filtering methods
#'
#' @return a character vector of names of filtering functions
#'
#' @examples
#' list_filters()
#'
#' @export

list_filters <- function()
    ls("package:kmerFilters")[grep(pattern = "filter_",
                                   x = ls("package:kmerFilters"))]


