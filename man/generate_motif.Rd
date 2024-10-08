% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifs.R
\name{generate_motif}
\alias{generate_motif}
\title{Random motifs generator}
\usage{
generate_motif(alphabet, n, d, motifProbs = NULL)
}
\arguments{
\item{alphabet}{elements used to generate a motif}

\item{n}{maximum number of alphabet elements}

\item{d}{number of possible gaps}

\item{motifProbs}{alphabet elements' probabilites}
}
\value{
motif built on a given alphabet
}
\description{
This function generates random motif based on a given alphabet.
The maximum range of a motif equals `n + d`.
}
\examples{
generate_motif(1:4, n = 2, d = 0)
generate_motif(c("a", "b", "c"), n = 6, d = 1)
generate_motif(1:4, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))

}
