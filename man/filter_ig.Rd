% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{filter_ig}
\alias{filter_ig}
\title{Entropy-Based Feature Selection Algorithms}
\usage{
filter_ig(target, kmers, method)
}
\arguments{
\item{target}{a numeric response variable}

\item{kmers}{a matrix of kmers with named columns or an object obtained via
\code{\link{generate_kmer_data}} function.}

\item{method}{a character name of a filter type. One of "infogain",
"gainratio" or "symuncert". For more details see
\code{\link[FSelectorRcpp]{information_gain}}.}
}
\value{
a numeric vector of named p-values corresponding to k-mers in the
feature space.
}
\description{
This function is a wrapper over Entropy-Based Feature Selection Algorithms
for k-mer data.
}
\details{
This function uses \code{\link[FSelectorRcpp]{information_gain}}
}
\examples{
n_seq <- 20
sequence_length <- 20
alph <- letters[1:4]
motifs <- generate_motifs(alph, 4, 4, 4, 6)
kmers <- generate_kmer_data(n_seq, sequence_length, alph,
                            motifs, n_injections = 4)
target <- get_target_additive(kmers)
filter_ig(target, kmers, "infogain", 0.01)

}
