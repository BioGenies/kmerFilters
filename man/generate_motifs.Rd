% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifs.R
\name{generate_motifs}
\alias{generate_motifs}
\title{Function generates list of motifs}
\usage{
generate_motifs(
  alphabet,
  n_motifs,
  n_injections,
  n,
  d,
  motifProbs = NULL,
  validate = TRUE,
  sequence_length = 10
)
}
\arguments{
\item{alphabet}{elements used to generate a motif}

\item{n_motifs}{number of motifs to generate}

\item{n_injections}{number of injections (for validation purposes:
checks if each subset of motifs of size `n_injections` can be injected
to a sequence of length `sequence_length`)}

\item{n}{maximum number of alphabet elements}

\item{d}{number of possible gaps}

\item{motifProbs}{alphabet elements' probabilites}

\item{validate}{if true, returns a set of motifs that can be injected
to a sequence of length 10}

\item{sequence_length}{length of a sequence that must contain all motifs}
}
\value{
list of generated motifs
}
\description{
This function generates multiple motifs from alphabet
}
\examples{
generate_motifs(1:4, 5, 3, n = 6, d = 6)
generate_motifs(1:4, 5, 3, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))

}
