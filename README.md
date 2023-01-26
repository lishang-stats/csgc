
<!-- README.md is generated from README.Rmd. Please edit that file -->

# csgc

<!-- badges: start -->
<!-- badges: end -->

The goal of csgc is to perform statistical analysis of network data
based on centred subgraph counts.

## Installation

You can install the development version of csgc from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lishang-stats/csgc")
```

## Example

Assume an undirected network has $n$ vertices, with $A=(A_{ij})_{nn}$ be
the adjacency matrix and `A_{ii}=0` for all $i$. Suppose the parameter
matrix is $K=(k_{rs})_{k\times k}$. Let $z=(z_1,\dots,z_n)$ be the label
vector, where $z_i=r$ means the $i$-th vertex is in group $r$. Consider
the Bernoulli-type Stochastic Block Model, i.e.
$$A_{ij}=A_{ji}\overset{i.i.d.}{\sim} \text{Be}(k_{z_iz_j})$$ for
$1\leq i,j\leq n, i\ne j$. We first generate matrix $A$ and $P$ from
$K$.

``` r
library(csgc)
#> 
#> Attaching package: 'csgc'
#> The following object is masked from 'package:base':
#> 
#>     trace
## basic example code
```
