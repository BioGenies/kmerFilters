% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{filter_fcbf}
\alias{filter_fcbf}
\title{Fast Correlation Based Filter}
\usage{
filter_fcbf(target, kmers, thresh = 0.25)
}
\arguments{
\item{target}{a numeric response variable}

\item{kmers}{a matrix of kmers with named columns or an object obtained via
\code{\link{generate_kmer_data}} function.}

\item{thresh}{a threshold for symmetrical uncertainty between a k-mer and
a target variable. Default to \code{0.25}.}
}
\value{
a character vector of names of selected kmers
}
\description{
This function is a wrapper over Fast Correlation Based Filter for k-mer data.
}
\details{
This function uses \code{\link[FCBF]{fcbf}}
}
\examples{
n_seq <- 20
sequence_length <- 20
alph <- letters[1:4]
motifs <- generate_motifs(alph, 4, 4, 4, 6)
kmers <- generate_kmer_data(n_seq, sequence_length, alph,
                            motifs, n_injections = 4)
target <- get_target_additive(kmers)
filter_fcbf(target, kmers)

}
