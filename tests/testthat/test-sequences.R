library(kmerFilters)
library(testthat)
context("Sequence generation")

test_that("generate single sequence", {

  set.seed(42)
  seq <- generate_sequence(5, 1L:4)
  expect_equal(seq, c(1L, 1L, 1L, 1L, 2L))
})

test_that("sequence probabilites", {
  seqProbs <- runif(6)
  seqProbs <- seqProbs / sum(seqProbs)

  long_sequence <- generate_sequence(10e6, 1:6, seqProbs)
  seqFreqs <- c(table(long_sequence) / length(long_sequence))
  names(seqFreqs) <- NULL
  expect_equal(seqProbs, c(seqFreqs), tolerance = 10e-4)
})

test_that("generate sequence data", {

  set.seed(42)

  n_sequences <- 20
  sequence_length <- 10
  alphabet <- letters[1:6]
  n_motifs <- 3
  n_injections <- 2
  fraction = 0.8

  motifs <- generate_motifs(alphabet, n_motifs, n_injections, n = 4, d = 6)
  seq_data <- generate_sequence_data(n_sequences,
                                     sequence_length,
                                     alphabet,
                                     motifs,
                                     n_injections,
                                     fraction)

  # correct fraction of positive sequences
  target_counts <- table(attr(seq_data, "target"))
  expect_equivalent(fraction, target_counts["TRUE"] / n_sequences)

  # all positive sequences contain `n_injections` motifs
  expect_equal(n_injections, unique(unlist(lapply(attr(seq_data, "motifs"),
                                                  length))))

  expect_equal(alphabet, sort(unique(c(seq_data))))
  expect_equal(n_sequences, nrow(seq_data))
  expect_equal(sequence_length, ncol(seq_data))
})

test_that("k-mer counts", {

  set.seed(42)

  n_seq <- 20
  sequence_length <- 10
  alph <- letters[1:6]
  motifs <- generate_motifs(alph, 3, 3, 3, 2)
  seq_data <- generate_sequence_data(n_seq, sequence_length, alph, motifs, 1)

  kmers <- count_seq_kmers(seq_data, alph)

  expect_equal(dim(kmers), c(20L, 6027L))

  kmers <- count_seq_kmers(seq_data, alph, n = 1, d = 0)
  expect_equal(ncol(kmers), length(alph))
})

test_that("k-mer data", {

  set.seed(42)

  n_seq <- 20
  sequence_length <- 10
  alph <- letters[1:6]
  n_motifs <- 3
  n_injections <- 3
  n <- 4
  d <- 6
  fraction = 0.8

  motifs <- generate_motifs(alph, n_motifs, n_injections, n, d)
  kmer_data <- generate_kmer_data(n_seq, sequence_length, alph, motifs,
                                  n_injections, fraction)


  expect_equal(dim(kmer_data), c(20L, 2616L))
  expect_equal(n_seq, nrow(kmer_data))

  # correct fraction of positive sequences
  target_counts <- table(attr(kmer_data, "target"))
  expect_equivalent(fraction, target_counts["TRUE"] / n_seq)

  # all positive sequences contain `n_injections` motifs
  expect_equal(n_injections, unique(unlist(lapply(attr(kmer_data, "motifs"),
                                                  length))))
})
