
library(kmerFilters)
library(testthat)


test_that("list filtering methods works", {
    expect_equal(list_filters(),
                 c("filter_chisq", "filter_fcbf", "filter_ic", "filter_ig",
                   "filter_praznik", "filter_quipt"))
})



n_seq <- 50
sequence_length <- 10
alph <- letters[1:4]
motifs <- generate_motifs(alph, 4, 4, 4, 6)
kmers <- generate_kmer_data(n_seq, sequence_length, alph,
                            motifs, n_injections = 4)
target <- get_target_additive(kmers)


test_that("quipt filter works", {
    res <- filter_quipt(target, kmers)

    expect_equal(res[1:10],
                 c(a_0 = 0.433460286708944, d_0 = 0.537216985854405,
                   b_0 = 0.862413139242882, d.c_0 = 1, d.b_0 = 0.904876795162511,
                   d.d_0 = 0.0668944277862377, b.a_0 = 0.904876795162511,
                   a.d_0 = 0.741734233008649, a.a_0 = 0.904876795162511,
                   c.b_0 = 0.904876795162511))
})


test_that("chisq filter works", {
    res <- filter_chisq(target, kmers)

    expect_equal(res[1:10],
                 c(a_0 = 0.359163595525953, d_0 = 0.619081488098256,
                   b_0 = 0.772148036951171, c_0 = 1, d.c_0 = 0.999999999999999,
                   d.b_0 = 0.693010969821401,  d.d_0 = 0.00137293131428169,
                   b.a_0 = 0.504838137702313, a.d_0 = 0.113732669593295,
                   a.a_0 = 0.438514651784921))
})


test_that("fcbf filter works", {
    res <- filter_fcbf(target, kmers)

    expect_equal(res, c("d.d.c_0.4", "d.d.b_0.1"))
})


test_that("FSelectorRcpp filter works", {
    res <- filter_ig(target, kmers, method = "infogain")

    expect_equal(res[1:5], c(a_0 = 0.0926487792836429, d_0 = NaN, b_0 = NaN,
                             c_0 = NaN, d.a_0 = 0.0318668917281004))
})


test_that("FSelectorRcpp filter works", {
    res <- filter_ig(target, kmers, method = "infogain")

    expect_equal(res[1:5], c(a_0 = 0.0926487792836429, d_0 = NaN, b_0 = NaN,
                             c_0 = NaN, d.a_0 = 0.0318668917281004))
})


test_that("praznik filter works", {
    res <- filter_praznik(target, kmers, method = "MIM")

    expect_equal(res[1:5], c(d.d.c_1.4 = 0.309448173049744,
                             b.c.a.d_0.0.1 = 0.309448173049744,
                             b.c.b.c_0.5.1 = 0.309448173049744,
                             b.c.d.b_0.2.2 = 0.309448173049744,
                             b.c.d_0.2 = 0.29557289101432))
})
