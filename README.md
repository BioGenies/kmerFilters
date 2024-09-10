
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
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

alph <- letters[1:4]
n_injections <- 4

motifs <- generate_motifs(alphabet = alph, 
                          n_motifs = 4, 
                          n_injections = 4, 
                          n = 4, 
                          d = 6)

motifs
#> [[1]]
#> [1] "b" "_" "_" "d"
#> 
#> [[2]]
#> [1] "c" "_" "c" "a" "c"
#> 
#> [[3]]
#> [1] "a" "_" "d"
#> 
#> [[4]]
#> [1] "c" "b" "_" "_" "_" "_" "b"
```

Using simulated motifs we can simulate positive and negative sequences
and consrtuct a k-mer feature space:

``` r
results <- generate_kmer_data(n_seq = 200, 
                              sequence_length = 20, 
                              alphabet = alph,
                              motifs = motifs, 
                              n_injections = 4)
results[1:10, ]
#> 10 x 23406 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 33 column names 'a_0', 'd_0', 'b_0' ... ]]
#>                                                                               
#>  [1,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 1 ......
#>  [2,] 1 1 1 1 . . 1 . . 1 . . 1 1 1 1 1 1 1 . . 1 1 1 . 1 . 1 1 . . 1 1 ......
#>  [3,] 1 1 1 1 1 1 . 1 . 1 1 1 . . 1 . 1 1 1 1 1 1 . 1 1 . 1 1 . 1 1 . 1 ......
#>  [4,] 1 1 1 1 1 1 . 1 . 1 1 1 . 1 1 . 1 1 1 . 1 . 1 1 1 . . 1 1 1 1 1 . ......
#>  [5,] 1 1 1 1 1 . 1 1 . 1 . 1 1 1 1 1 1 1 1 1 . . . 1 1 . 1 1 1 1 . 1 1 ......
#>  [6,] 1 1 1 1 . . 1 . . 1 1 . 1 1 1 1 1 . 1 1 1 1 1 1 . 1 . 1 1 . . 1 1 ......
#>  [7,] 1 1 1 1 1 . . 1 1 1 . 1 . 1 1 1 . 1 1 . 1 1 . 1 1 1 . 1 . 1 1 . 1 ......
#>  [8,] 1 1 1 1 . 1 . . . 1 1 . 1 1 1 . 1 . 1 . . . 1 1 1 1 . 1 1 1 1 1 . ......
#>  [9,] 1 1 1 1 1 1 1 . 1 1 1 1 . 1 1 1 1 . 1 . . . 1 1 1 . 1 1 1 1 1 1 . ......
#> [10,] 1 1 1 1 1 1 . . . 1 1 1 . 1 1 1 1 1 1 . 1 . 1 1 1 1 1 1 1 1 1 1 1 ......
#> 
#>  .....suppressing 23373 columns in show(); maybe adjust options(max.print=, width=)
#>  ..............................
```

Using obtained data you can choose how to generate a response variable.
We provide three functions for that:

- `get_target_additive()`,
- `get_target_interactions()`,
- `get_target_logic()`.

For example, the following code:

``` r
target <- get_target_additive(results)
target
#>   [1] 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 0 1 1 0 1 0 1 1 1 1 1 1 1
#>  [38] 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 1 1 0 1 0 1 0 0 0 1 1 0 1 1 0
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1
#> [112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
#> [149] 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
#> [186] 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1
```

creates a binary response variable based on the logistic regression
model assumptions.

### How to filter

We have implemented many filtering methods. You can list them easily
using the function

``` r
list_filters()
#> [1] "filter_chisq"   "filter_fcbf"    "filter_ic"      "filter_ig"     
#> [5] "filter_praznik" "filter_quipt"
```

Using k-mer space and target variable you can filter k-mers as follows:

``` r
filter_quipt(target, results, significance_level = 0.05)
#>  [1] "a.c_1"         "b.d_2"         "b.b_4"         "c.b_5"        
#>  [5] "c.a.c_0.0"     "c.c.a_1.0"     "c.c.c_1.1"     "b.c.c_3.1"    
#>  [9] "c.b.b_0.4"     "c.c.a.c_1.0.0" "d.c.c.a_0.1.0" "c.b.b.b_0.3.0"
#> [13] "c.c.c.c_1.1.1" "b.c.c.c_1.1.1" "d.c.c.c_3.1.1" "b.c.b.b_0.2.1"
#> [17] "b.c.b.b_0.3.1" "c.b.c.b_0.0.3"
```

## Contributing

Contributions to kmerFilters are welcome! If you have suggestions for
new features, improvements, or bug fixes, please submit an issue or pull
request on GitHub.
