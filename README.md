# elastic

This repository has theory and code for the Elastic Multidimensional Scaling Method
proposed by Victor McGee in 1965. We have implemented the metric (ratio) and
non-metric (ordinal) versions in both R and in C (where the C version has
a driver in  R). On moderate size examples the C version is 20 times as
fast as the R version in the metric case and 100 times as fast in the non-metric
case.

The C version in a shared library called by the R program in smacofSSElasticC.R, 
the R version is smacofSSElasticR.R. Currently the programs assume square symmetric data.
There are R utilities for data manupulation and for plotting.