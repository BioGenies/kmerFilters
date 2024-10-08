% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_target.R
\name{get_target_logic}
\alias{get_target_logic}
\title{Logic regression noise}
\usage{
get_target_logic(
  kmer_dat,
  random = TRUE,
  zero_weight = NULL,
  weights = NULL,
  n_exp = NULL,
  max_exp_depth = NULL,
  expressions = NULL,
  binary = TRUE
)
}
\arguments{
\item{kmer_dat}{output of \code{\link{generate_kmer_data}}}

\item{random}{a logical. Indicating whether expressions have to be generated
randomly. Default to \code{TRUE}.}

\item{zero_weight}{a single value denoting the weight of no-motifs case. If
\code{NULL}, then we sample the weight from the uniform distribution on the
[-2, -1] interval. Default to \code{NULL}.}

\item{weights}{a vector of weights of considered logic expression based on
available motifs. The length of \code{weights} should be the same as the
provided number  of expressions to use \code{n_exp}. If \code{weights}
parameter is \code{NULL}, then weights will be sampled from the uniform
distribution on 0-1 interval. The probability of success for target sampling
will be calculated based on the formula provided in details section. Default
to \code{NULL}.}

\item{n_exp}{number of random logic expressions to create. It is used only
when \code{random} equals \code{TRUE}.}

\item{max_exp_depth}{a maximum number of motifs used in a logic expression.
Default to 3.}

\item{expressions}{a matrix of binary variables corresponding to custom
logic expressions. You can create them based on motifs. It's dimension should
be related to the length of \code{weights} vector if it's provided. Default
to \code{NULL}. If \code{NULL}, random logic expressions will be created.}

\item{binary}{logical, indicating whether the produced target variable should
be binary or continuous.}
}
\description{
This function samples target variable according to the logic regression
model (assuming that the occurrence of certain combinations of motifs affects
the feature). In the case of logical models, simulating a binary variable
involves defining logical conditions that determine the variable's value
based on motifs, e.g., the binary variable takes the value 1 if certain l
ogical criteria are met, and 0 if these criteria are not met.
}
\details{
Here, we consider new variables, \eqn{L_1, \ldots, L_l} where each of them
is a logic expression based on a subset of motifs \eqn{m_1, \ldots, m_m}. For
example,

\eqn{L_1(m_1, m_2, m_3) = (X_{m_1} \land X_{m_2}) \lor X_{m_3}.}

Each variable \eqn{L_i} obtains its own weight in the model. Our model is
following:

\eqn{g(EY) = w_0 + \sum_{i = 1}^{l} w_i L_i.}
}
\examples{
n_seq <- 20
sequence_length <- 20
alph <- letters[1:4]
motifs <- generate_motifs(alph, 4, 4, 4, 6)
results <- generate_kmer_data(n_seq, sequence_length, alph,
                              motifs, n_injections = 4)
get_target_logic(results)

}
