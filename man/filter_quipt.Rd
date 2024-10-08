% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{filter_quipt}
\alias{filter_quipt}
\title{QuiPT - Quick Permutation Test}
\usage{
filter_quipt(target, kmers, ...)
}
\arguments{
\item{target}{a numeric response variable}

\item{kmers}{a matrix of kmers with named columns or an object obtained via
\code{\link{generate_kmer_data}} function.}

\item{...}{other arguments for \code{\link[biogram]{test_features}} function.}
}
\value{
a numeric vector of named p-values corresponding to k-mers in the
feature space.
}
\description{
This function is a wrapper over the QuiPT test.
}
\details{
This function uses \code{\link[biogram]{test_features}}
}
\examples{
n_seq <- 20
sequence_length <- 20
alph <- letters[1:4]
motifs <- generate_motifs(alph, 4, 4, 4, 6)
kmers <- generate_kmer_data(n_seq, sequence_length, alph,
                            motifs, n_injections = 4)
target <- get_target_additive(kmers)
filter_quipt(target, kmers)

}
