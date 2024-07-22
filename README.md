
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kmerFilters : R Package for k-mer Data Simulation and Filtering

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/BioGenies/kmerFilters/branch/main/graph/badge.svg)](https://app.codecov.io/gh/BioGenies/kmerFilters?branch=main)
<!-- badges: end -->

## Overview

k-mers (n-grams) refer to k-length substrings derived from longer
sequences, which can be continuous (representing a block of subsequent
residues) or discontinuous (where the wildcards representing the gaps
between residues are allowed). In biological applications, k-mers are
used for various purposes, including genome assembly, sequence
alignment, motif discovery, variant calling, phylogenetics, and CRISPR
target identification. Due to the vast number of variables introduced by
k-mer notation, we need tools to filter the variables.

In this package we provide tools for simulating k-mer data and
benchmarking various filtering techniques. It is designed to assist
researchers and practitioners in evaluating and comparing different
approaches for preprocessing sequence data, particularly for
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

To learn more, please visit: <https://biogenies.info/kmerFilters/>

Specific resources:

- [Simulation
  description](https://biogenies.info/kmerFilters/articles/simulation_description.html)

## How to install

You can install the latest version from GitHub:

    # install.packages("devtools")
    devtools::install_github("BioGenies/kmerFilters")

## Examples of usage

### How to simulate data

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
#> [1] "d" "_" "_" "_" "_" "_" "b"
#> 
#> [[2]]
#> [1] "c" "_" "_" "d" "c"
#> 
#> [[3]]
#> [1] "b" "_" "c" "_" "_" "c"
#> 
#> [[4]]
#> [1] "c" "_" "_" "b" "_" "_" "c"
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
#> 20 x 14940 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 33 column names 'a', 'd', 'b' ... ]]
#>                                                                               
#>  [1,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 . ......
#>  [2,] 1 1 1 1 . 1 1 1 1 1 1 1 1 1 1 1 1 . . . 1 1 . 1 1 1 1 1 . 1 . 1 1 ......
#>  [3,] 1 1 1 1 . 1 . 1 . . 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . 1 1 1 1 . 1 1 1 ......
#>  [4,] 1 1 1 1 . 1 1 1 1 1 . 1 1 . 1 1 1 1 1 1 . 1 1 1 . 1 . 1 1 1 1 1 1 ......
#>  [5,] 1 1 1 1 1 1 1 . 1 1 1 1 1 1 . 1 1 1 . . 1 1 1 . 1 1 . 1 1 1 1 1 . ......
#>  [6,] 1 1 1 1 1 . 1 1 1 1 1 1 . 1 . 1 . 1 1 . 1 1 . 1 . . . 1 1 1 1 . 1 ......
#>  [7,] 1 1 1 1 1 1 . 1 1 . 1 1 1 1 1 1 . 1 1 1 . 1 1 . 1 1 . 1 1 1 1 . . ......
#>  [8,] 1 1 1 1 1 . 1 . 1 1 1 1 . 1 . 1 1 . . . 1 . . 1 1 . 1 1 . 1 1 1 . ......
#>  [9,] 1 1 1 1 1 1 . 1 1 1 1 1 1 1 . 1 . . 1 1 1 . . 1 1 . 1 1 1 1 1 1 1 ......
#> [10,] 1 1 1 1 1 1 1 1 . 1 1 1 1 . . . 1 1 . . 1 1 1 . 1 1 1 1 . 1 1 1 . ......
#> [11,] 1 1 1 1 . . . 1 . 1 . 1 1 1 1 1 1 1 1 . 1 1 1 1 1 1 1 1 1 . 1 . . ......
#> [12,] 1 1 1 1 1 1 1 . . 1 1 1 1 1 1 1 1 . 1 . 1 1 . . 1 1 1 1 . 1 1 1 . ......
#> [13,] 1 1 1 1 1 1 1 . 1 1 . 1 . 1 1 . 1 1 1 1 1 1 . . . 1 . 1 1 1 1 1 1 ......
#> [14,] 1 1 1 1 1 . 1 1 1 . 1 1 1 1 1 1 . . 1 1 1 1 1 1 . 1 1 1 1 1 . 1 . ......
#> [15,] 1 1 . 1 1 1 . 1 . 1 1 . 1 . . . 1 . 1 1 1 . 1 . 1 . 1 . 1 1 . 1 1 ......
#> [16,] 1 1 1 1 . 1 1 . 1 1 . 1 . 1 . 1 1 1 1 . 1 . 1 . . 1 . 1 . 1 1 . . ......
#> [17,] 1 1 1 1 1 1 1 1 1 1 1 1 . 1 1 . . 1 1 . 1 1 1 1 1 1 . . 1 1 1 1 1 ......
#> [18,] 1 1 1 1 1 1 1 1 1 1 . 1 1 . 1 . . 1 1 1 1 1 . 1 . 1 . 1 1 1 . 1 1 ......
#> [19,] 1 1 1 1 . . 1 1 . 1 . 1 1 1 1 1 1 1 1 1 . 1 . 1 . . . 1 1 1 1 . 1 ......
#> [20,] 1 1 1 1 1 . 1 1 1 . 1 1 1 1 1 1 . 1 1 1 1 1 . . . 1 1 1 1 1 1 . 1 ......
#> 
#>  .....suppressing 14907 columns in show(); maybe adjust 'options(max.print= *, width = *)'
#>  ..............................
```

Using obtained data you can choose how to generate a response variable.
We provide three functions for that:

- `get_target_additive()`,
- `get_target_interactions()`,
- `get_target_logic()`.

For example, the following code:

``` r
get_target_additive(results)
#>  [1] 0 1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0
```

creates a binary response variable based on the logistic regression
model assumptions.

### How to filter

## Contributing

Contributions to kmerFilters are welcome! If you have suggestions for
new features, improvements, or bug fixes, please submit an issue or pull
request on GitHub.
