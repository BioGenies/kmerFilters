
library(kmerFilters)
library(testthat)


test_that("list filtering methods works", {
    expect_equal(list_filters(),
                 c("filter_chisq", "filter_fcbf", "filter_ic", "filter_ig",
                   "filter_praznik", "filter_quipt"))
})

set.seed(12)

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
                 c(d_0 = 0.75732779982674, b_0 = 0.728316101771109,
                   d.c_0 = 1, d.a_0 = 1, a.b_0 = 0.728316101771109, d.d_0 = 1,
                   b.a_0 = 1, b.c_0 = 1, a.d_0 = 0.728316101771109, c.d_0 = 1))
})


test_that("chisq filter works", {
    res <- filter_chisq(target, kmers)

    expect_equal(res[1:10],
                 c(a_0 = 1, d_0 = 0.885992079772807, b_0 = 0.137259222349999,
                   c_0 = 1, d.c_0 = 0.710251149725874, d.a_0 = 0.470784349666462,
                   a.b_0 = 0.056129794661529, d.d_0 = 0.793891166336373,
                   b.a_0 = 0.514976937835016, b.c_0 = 1))
})


test_that("fcbf filter works", {
    res <- filter_fcbf(target, kmers, thresh = 0.1)

    expect_equal(res, c("c.d.c_1.3", "a.d.d.c_0.2.1", "b.d.a_2.1", "d.a.d_4.1", "b.a_6",
                        "a.b.b_0.1", "d.a.b_0.1", "c.a.a_2.2", "d.c.b_1.4", "d.c.c_3.0",
                        "c.d.c_3.1", "b.c.c_0.0", "d.b.d_0.1", "c.a.d.c_1.0.3", "b.a.b_2.0",
                        "a.a.b_3.0", "a.c.a_2.1", "c.d.d_0.2", "a.b.d_0.3", "c.c.d_0.3",
                        "a.b.a_1.3", "d.c.d_2.3", "d.c.a_1.4", "d.d.d.d_0.1.1"))
})


test_that("FSelectorRcpp filter works", {
    res <- filter_ig(target, kmers, method = "infogain")

    expect_equal(res[1:5], c(a_0 = NaN, d_0 = 0.0243373426620957,
                             b_0 = 0.0545796620922217, c_0 = NaN,
                             d.c_0 = 0.00758121297699884))
})


test_that("FSelectorRcpp filter works", {
    res <- filter_ig(target, kmers, method = "infogain")

    expect_equal(res[1:5], c(a_0 = NaN, d_0 = 0.0243373426620957,
                             b_0 = 0.0545796620922217, c_0 = NaN,
                             d.c_0 = 0.00758121297699884))
})


test_that("praznik filter works", {
    res <- filter_praznik(target, kmers, method = "MIM")

    expect_equal(res[1:5], c(c.d.c_1.3 = 0.127675835930895,
                             a.d.d.c_0.2.1 = 0.127675835930895,
                             d.d.a_2.2 = 0.111515430511026,
                             c.a.d.d_0.0.0 = 0.111515430511026,
                             c.c_5 = 0.101033251383654))
})
