library(kmerFilters)
library(testthat)

context("Motif generation")

test_that("generate single motif", {

  set.seed(42)

  alphabet <- letters[1:4]
  n <- 4
  d <- 2

  motif <- generate_motif(alphabet, n, d)
  expect_equal(motif, c("a", "a"))

  #weighted alphabet elements
  weights <- runif(4)
  weights <- weights / sum(weights)

  motif <- generate_motif(alphabet, n, d, weights)
  expect_equal(motif, c("b", "_", "_", "a"))
})


test_that("checkmate assertions", {

  set.seed(42)

  alphabet <- letters[1:4]
  n <- 4
  d <- 2

  weights <- runif(4)
  weights <- weights / sum(weights)

  expect_error(generate_motif(alphabet, "n", d))
  expect_error(generate_motif(alphabet, n, "d"))
  expect_error(generate_motif(c(alphabet, alphabet), n, d))
  expect_error(generate_motif(alphabet, n, d, weights*2))
  expect_error(generate_motif(alphabet, n, d, c("w", "e", "i", "g")))
  expect_error(generate_motif(alphabet, n, d,
                              c("w", "e", "i", "g", "h", "t", "s")))

  # check if this setup works when called properly
  motif <- generate_motif(alphabet, n , d, weights)
  expect_equal(motif, c("a", "d", "b"))

  })

test_that("generate list of motifs", {

  set.seed(42)

  alphabet <- letters[1:4]
  n <- 4
  d <- 2
  n_motifs <- 3
  n_injections <- 3

  weights <- runif(4)
  weights <- weights / sum(weights)

  motifs <- generate_motifs(alphabet, n_motifs, n_injections, n, d, weights)

  expect_equal(motifs,
               list(c("a", "d", "b"),
                    c("d", "_", "c", "_", "b", "a"),
                    c("b", "_", "_", "a", "a")))

  motifs <- generate_motifs(alphabet, 1, 1, n, d, weights, validate = FALSE)

  expect_equal(motifs, list(c("c", "_", "_", "a", "a", "a")))
})

test_that("motif injection", {

  set.seed(1)

  sequence <- rep(letters[1:4], each=2)
  motifs <- generate_motifs(letters[1:4], 2, 2, 4, 2)
  injected_sequence <- add_motifs(motifs, sequence)
  expect_equal(injected_sequence,
               structure(c("a", "c", "b", "d", "c", "c", "d", "d"),
                         motifs = list(c("d", "c"),
                                       c("a", "c", "_", "_", "c")),
                         masks = list(c(FALSE, FALSE, FALSE, TRUE,
                                        TRUE, FALSE, FALSE, FALSE),
                                      c(TRUE, TRUE, FALSE, FALSE,
                                        TRUE, FALSE, FALSE, FALSE))))
})

test_that("set of motifs cannot be injected", {

  set.seed(1)

  expect_error(add_motifs(generate_motifs(letters[1:4], 5, 5, 3, 1),
                          sequence))
})

