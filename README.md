# Tailor
An R package for clustering via Gaussian Mixture Modeling, optimized for heavy-tailed distributions in flow cytometry data.

## Installation
The vignette may take a few minutes to build; if it causes any trouble, just set `build_vignettes = FALSE`.

`devtools::install_github("matei-ionita/Tailor", build_vignettes = TRUE)`

## Summary
Due to increasing size and complexity of flow cytometry datasets, there has been a lot of recent interest in unsupervised 
clustering approaches for distinguishing phenotypically distinct cell populations in the data. One particular challenge
are unimodal distributions (i.e. single peak) which contain one or more small populations hiding in the tail of a larger
one. The Tailor package uses Gaussian mixtures to model such distributions. There are two main hurdles to overcome:

* Mixture modeling is computationally expensive: O(N\*k\*d^2) per iteration, where N is the number of datapoints, k
the number of mixture components, and d the number of features.

* The Expectation-Maximization (EM) algorithm, commonly used to find an iterative solution to the mixture modeling problem,
optimizes a non-convex function. Therefore, the result depends strongly on the initialization, even more so as d increases.

Tailor addresses both of these inconveniences with the help of a preliminary binning procedure. First, the data is subset into
bins of phenotypically similar data points. Then certain bin representatives are chosen, some of which serve as a weighted
subsample of the original data, and some of which are used to provide a sensible initialization for the EM
algorithm.

Our weighted subsampling scheme ensures that small populations are not lost in the subsampling step. However, it also
requires some modifications of the standard EM algorithm, to take into account the weights, and
to correct for a systematic underestimation of variance due to subsampling.

After the modified EM algorithm is run, some of the resulting mixture components are deemed to be similar enough and merged 
into the same categorical cluster. Tailor then outputs information about the phenotypes of these categorical clusters.
Given new data, Tailor can quickly predict the cluster assignment for each new datapoint, based on the model learned previously.

## Warning
Tailor is a work in progress. Let me know if you find bugs or inconsistencies.
