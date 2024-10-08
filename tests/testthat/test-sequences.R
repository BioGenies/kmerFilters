library(kmerFilters)
library(testthat)

context("Sequence generation")


test_that("generate single sequence", {
    set.seed(42)

    motifs <- generate_motifs(1:4, 3, 1, n = 4, d = 6)
    sequence <- generate_negative_sequence(5, 1L:4, motifs)
    expect_equal(sequence, c(3L, 2L, 1L, 2L, 4L))
})

test_that("sequence probabilites", {
    set.seed(42)
    seqProbs <- runif(6)
    seqProbs <- seqProbs / sum(seqProbs)

    motifs <- generate_motifs(1:4, 3, 1, n = 4, d = 6)

    long_sequence <- generate_negative_sequence(10e5, 1:6, motifs, seqProbs)
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

    # all positive sequences contain maximum of`n_injections` motifs
    expect_true(
        all(unique(unlist(lengths(attr(seq_data, "motifs")))) <= n_injections)
    )

    expect_equal(alphabet, sort(unique(c(seq_data))))
    expect_equal(n_sequences, nrow(seq_data))
})

test_that("k-mer counts", {

    set.seed(42)

    n_seq <- 20
    sequence_length <- 10
    alph <- letters[1:6]
    motifs <- generate_motifs(alph, 3, 3, 3, 2)
    seq_data <- generate_sequence_data(n_seq, sequence_length, alph, motifs, 1)

    kmers <- count_seq_kmers(seq_data, alph)

    expect_equal(dim(kmers), c(20L, 6201L))

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

    expect_equal(dim(kmer_data), c(20L, 5022L))
    expect_equal(n_seq, nrow(kmer_data))

    # correct fraction of positive sequences
    target_counts <- table(attr(kmer_data, "target"))
    expect_equivalent(fraction, target_counts["TRUE"] / n_seq)

    # all positive sequences contain `n_injections` motifs
    expect_true(
        all(unique(unlist(lengths(attr(kmer_data, "motifs")))) <= n_injections)
    )
})


test_that("Interaction model works", {

    set.seed(1)

    n_seq <- 20
    sequence_length <- 20
    alph <- letters[1:4]
    motifs <- generate_motifs(alph, 4, 4, 4, 6)
    results <- generate_kmer_data(n_seq, sequence_length, alph,
                                  motifs, n_injections = 3)
    expect_identical(
        c(1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L,
          0L, 0L, 1L, 0L, 0L), get_target_interactions(results)
    )

    expect_equal(
        c(2.3877928876318, 4.19229120924138, 2.12166230170988, 0.988599377218634,
          4.47515172231942, 1.58585822791792, 1.17262414400466, 0.266789596760646,
          4.73356275283732, 2.3877928876318, -1.03832712979056, -1.03832712979056,
          -1.03832712979056, -1.03832712979056, -1.03832712979056, -1.03832712979056,
          -1.03832712979056, -1.03832712979056, -1.03832712979056, -1.03832712979056
        ), get_target_interactions(results, binary = FALSE)
    )
})



test_that("Additive model works", {

    set.seed(1)

    n_seq <- 20
    sequence_length <- 20
    alph <- letters[1:4]
    motifs <- generate_motifs(alph, 4, 4, 4, 6)
    results <- generate_kmer_data(n_seq, sequence_length, alph,
                                  motifs, n_injections = 3)
    expect_identical(
        c(0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L,
          0L, 1L, 0L, 0L, 0L), get_target_additive(results)
    )

    expect_equal(
        c(0.675145752727985, 1.19026312511414, 1.00235199206509),
        get_target_additive(results, binary = FALSE)[1:3]
    )

    expect_error(get_target_additive(results, weights = 0.1),
                 "The length of weights vector should equal number of motifs!")
})




test_that("Logic model works", {

    set.seed(1)

    n_seq <- 20
    sequence_length <- 20
    alph <- letters[1:4]
    motifs <- generate_motifs(alph, 4, 4, 4, 6)
    results <- generate_kmer_data(n_seq, sequence_length, alph,
                                  motifs, n_injections = 3)
    expect_identical(
        c(1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L), get_target_logic(results)
    )

    expect_equal(
        c(0.583639514865354, 1.4445967883803, 1.4445967883803, 0.583639514865354
        ), get_target_logic(results, binary = FALSE)[1:4]
    )

    expressions = matrix(rbinom(n_seq*2, 1, .5), nrow = n_seq)

    expect_error(
        get_target_logic(results, expressions = expressions, weights = 1),
        "You have to provide weight for each column"
    )

    expect_identical(
        c(0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L),
        get_target_logic(results, expressions = expressions, weights = c(1, 2))
    )

    motifs <- generate_motifs(alph, 1, 1, 4, 6)
    results <- generate_kmer_data(n_seq, sequence_length, alph,
                                  motifs, n_injections = 1)

    expect_error(get_target_logic(results), "You need at least 2 motifs")
})


test_that("rbinom_vec works", {

    set.seed(1)
    probs <- runif(10)

    expect_identical(rbinom_vec(probs),
                     c(0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
    expect_error(rbinom_vec(2),"Provided probabilities should be greater")
})


