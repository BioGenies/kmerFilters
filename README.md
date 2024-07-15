
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kmerFilters : R Package for k-mer Data Simulation and Filtering

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/jakubkala/kmerFilters/branch/main/graph/badge.svg)](https://codecov.io/gh/jakubkala/kmerFilters?branch=main)
<!-- badges: end -->

## Overview

The kmerFilters is an R package providing tools for simulating k-mer
data and benchmarking various filtering techniques. It is designed to
assist researchers and practitioners in evaluating and comparing
different approaches for preprocessing sequence data, particularly for
applications such as protein function prediction.

## Features

Our package provides fast and efficient generation of synthetic k-mer
datasets with customizable parameters, including sequence length, number
of sequences, k-mer size and motifs impact. The latest embraces
modelling the response variable in discrete or continuous way depending
on few biological cases, such as

- additive impact of each motif,

- additive impact of interactions of motifs,

- additive impact of custom interactions (logical expressions using
  motifs).

### How to install

You can install the latest version from GitHub:

    # install.packages("devtools")
    devtools::install_github("BioGenies/kmerFilters")

### How to use

You can use the `generate_motifs` function for motifs generation. For
example:

``` r
library(kmerFilters)

alph <- letters[1:4]
n_injections <- 4

motifs <- generate_motifs(alphabet = alph, 
                          n_motifs = 4, 
                          n_injections = 4, 
                          n = 4, 
                          d = 6)

motifs
#> [[1]]
#> [1] "d" "_" "_" "_" "_" "_" "d"
#> 
#> [[2]]
#> [1] "b" "_" "_" "_" "_" "_" "a"
#> 
#> [[3]]
#>  [1] "b" "_" "_" "d" "_" "_" "_" "_" "d" "b"
#> 
#> [[4]]
#> [1] "c" "_" "d" "_" "c"
```

Using simulated motifs we can simulate positive and negative sequences
and consrtuct a k-mer feature space:

``` r
results <- generate_kmer_data(n_seq = 20, 
                              sequence_length = 20, 
                              alphabet = alph,
                              motifs = motifs, 
                              n_injections = 4)
results
#> 20 x 15169 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 33 column names 'a', 'd', 'b' ... ]]
#>                                                                               
#>  [1,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 1 ......
#>  [2,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . 1 1 1 1 . 1 1 1 1 1 . 1 1 1 . . 1 1 ......
#>  [3,] 1 1 1 1 1 1 1 1 . . 1 1 1 1 . 1 . 1 . . 1 1 1 . 1 . . 1 1 1 1 1 . ......
#>  [4,] 1 1 1 1 1 1 1 1 1 . 1 . 1 . 1 1 . 1 1 1 1 . 1 1 1 1 1 . . 1 1 1 1 ......
#>  [5,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . 1 1 1 1 . . 1 1 1 . . 1 1 . ......
#>  [6,] 1 1 1 1 1 . 1 1 1 . 1 1 1 . 1 1 1 1 . . 1 . 1 1 . . 1 1 1 1 1 1 1 ......
#>  [7,] 1 1 1 1 1 1 1 1 1 . 1 . 1 1 . 1 . . 1 . 1 1 1 . 1 . . 1 . 1 1 . 1 ......
#>  [8,] 1 1 1 1 1 . 1 1 1 1 1 1 . 1 . 1 . 1 1 1 1 1 1 1 1 . 1 . . 1 1 1 1 ......
#>  [9,] 1 1 1 1 . . 1 1 1 1 . . 1 . 1 . . 1 1 1 . . 1 1 . 1 1 1 1 1 . . 1 ......
#> [10,] 1 1 1 1 1 . 1 1 1 . 1 1 1 1 . 1 1 . . . 1 1 1 1 . 1 1 . . . 1 1 . ......
#> [11,] 1 1 1 1 1 1 1 1 1 . 1 . 1 1 1 1 1 . . 1 1 1 1 . . . 1 1 1 1 1 1 1 ......
#> [12,] 1 1 1 1 . . 1 1 1 1 1 1 1 1 . 1 . 1 1 1 . 1 1 1 1 1 1 1 1 1 1 . 1 ......
#> [13,] 1 1 1 1 1 1 . 1 . . 1 1 . . 1 1 1 1 1 . 1 . . 1 1 . 1 . 1 . 1 1 1 ......
#> [14,] 1 1 1 1 1 1 . 1 . 1 . 1 1 1 1 . . 1 1 . 1 1 . 1 . 1 . 1 . 1 1 1 1 ......
#> [15,] 1 1 1 1 1 1 . 1 1 . . 1 1 . 1 . 1 1 1 1 . . . 1 1 1 1 1 1 1 1 . 1 ......
#> [16,] 1 1 1 1 . 1 1 1 . 1 1 . 1 . 1 . 1 1 1 1 . 1 . 1 1 1 1 1 1 1 . 1 1 ......
#> [17,] 1 1 1 1 . 1 . 1 1 1 . 1 1 . . 1 1 1 1 1 1 . 1 1 1 . . 1 1 . . . 1 ......
#> [18,] 1 1 1 1 1 . 1 1 . . 1 . 1 1 1 . 1 1 . . 1 . 1 . . 1 1 . 1 1 . 1 1 ......
#> [19,] 1 1 1 1 1 1 . . 1 1 1 1 . 1 1 1 . 1 . . 1 1 1 . . . . . 1 . 1 . . ......
#> [20,] 1 1 1 1 1 . 1 1 1 1 1 1 . . 1 1 1 1 1 . 1 1 1 . 1 . . 1 1 1 1 . . ......
#> 
#>  .....suppressing 15136 columns in show(); maybe adjust 'options(max.print= *, width = *)'
#>  ..............................
```

Using obtained data you can choose how to generate a response variable.
We provide three functions for that:

- `get_target_additive()`,
- `get_target_interactions()`,
- `get_target_logic()`.

For example:

``` r
get_target_additive(results)
#>  [1] 1 1 0 1 1 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0
```

### Contributing

Contributions to kmerFilters are welcome! If you have suggestions for
new features, improvements, or bug fixes, please submit an issue or pull
request on GitHub.
