
#' QuiPT - Quick Permutation Test
#'
#' This function is a wrapper over the QuiPT test.
#'
#' @importFrom biogram test_features
#'
#' @param target a numeric response variable
#' @param kmers a matrix of kmers with named columns or an object obtained via
#' \code{\link{generate_kmer_data}} function.
#' @param significance_level a number from 0-1 interval denoting significance
#' level for testing.
#' @param ... other arguments for \code{\link[biogram]{test_features}} function.
#'
#' @return a character vector of names of selected kmers
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

filter_quipt <- function(target, kmers, significance_level = 0.05, ...) {

    pvals <- biogram::test_features(target, kmers, ...)

    names(pvals[which(pvals < significance_level)])
}


#' Chi-squared independence test
#'
#' This function is a wrapper over chi-squared independence test for k-mer data.
#'
#' @importFrom stats chisq.test
#'
#' @inheritParams filter_quipt
#'
#' @return a character vector of names of selected kmers
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

filter_chisq <- function(target, kmers, significance_level = 0.05) {

    pvals <- unlist(lapply(1:ncol(kmers[, 1:10]), function(i) {

        ith_kmer <- as.vector(kmers[, i])
        suppressWarnings({
            if(length(unique(ith_kmer)) == 1) pval <- 1
            else pval <- chisq.test(ith_kmer, target, B = 200)$p.value
        })
        pval
    }))

    colnames(kmers)[which(pvals < significance_level)]
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
#' @param thresh a numeric threshold for variable selection.
#'
#' @return a character vector of names of selected kmers
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

filter_ig <- function(target, kmers, method, thresh) {

    method <- match.arg(method, c("infogain", "gainratio", "symuncert"))

    res <- information_gain(x = kmers, y = target, discIntegers = FALSE,
                            type = "gainratio")

    res[["attributes"]][res[["importance"]] < thresh]
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
#' @param thresh a threshold for corresponfing score.
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

filter_praznik <- function(target, kmers, method, thresh) {

    praznik_methods <-
        c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")

    method <- match.arg(method, praznik_methods)

    res <-  get(method)(X = as.data.frame(as.matrix(kmers)), Y = target,
                        k = ncol(kmers))
    names(res[["score"]])[res[["score"]] < thresh]
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
#'
#' @inheritParams filter_quipt
#'
#' @param ic character name of information criterium. One of "aic", "maic",
#' "bic", "maic2", "mbic", "mbic2". See bigstep package for more information.
#'
#' @param reduce a numeric value from (0, 1) interval. Denotes significance
#' level for preliminary reduction before execution of stepwise procedure.
#' Default to 0.2.
#'
#' @return a character vector of names of selected kmers
#'
#' @details
#' This function uses bigstep package.
#'
#' @examples
#' n_seq <- 200
#' sequence_length <- 200
#' alph <- letters[1:20]
#' motifs <- generate_motifs(alph, 4, 4, 4, 6)
#' kmers <- generate_kmer_data(n_seq, sequence_length, alph,
#'                             motifs, n_injections = 4)
#' target <- get_target_additive(kmers)
#' filter_ic(target, kmers, "mbic2", 0.2)
#'
#' @export

filter_ic <- function(target, kmers, ic = "mbic2", reduce = 0.1) {

    ic <- match.arg(ic, c("aic", "maic", "bic", "maic2", "mbic", "mbic2"))

    dat <- prepare_data(target, as.matrix(kmers))
    dat <- reduce_matrix(dat, minpv = reduce)
    candidates_ids <- dat[["candidates"]]

    res_ic <- stepwise(dat, crit = get(ic))

}


