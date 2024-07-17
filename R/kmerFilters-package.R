#' kmerFilters
#'
#' @aliases kmerFilters
#'
#' @import checkmate
#' @importFrom stats lm model.matrix
#'
#' @name kmerFilters
#'

if(getRversion() >= "2.15.1")
    utils::globalVariables(c('TARGET.GROUP', 'ACP', 'AMP', 'ID', 'SEQUENCE',
                             'NA_frac'))
