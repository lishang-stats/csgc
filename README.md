
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

Assume an undirected unweighted network has $n$ vertices, with $A$ be
the adjacency matrix and $P$ be the corresponding probability matrix.
The parameter matrix $K$ (a $k$ by $k$ matrix) presents the connection
probability within/between each block. Let $z$ be the label vector with
length $n$, where each entry can be from $1$ to $k$. Consider the
Bernoulli-type Stochastic Block Model
![equation](https://latex.codecogs.com/png.image?\dpi%7B200%7DA_%7Bij%7D=A_%7Bji%7D\overset%7Bi.i.d.%7D%7B\sim%7D%20\mathop%7B\mathrm%7BBe%7D%7D(k_%7Bz_iz_j%7D),%20\quad\text%7Bfor%20$1\leq%20i,j\leq%20n,%20i\ne%20j$%7D)
. We first generate matrix $A$ and $P$ from $K$.

``` r
set.seed(123)
library(csgc)
#> 
#> Attaching package: 'csgc'
#> The following object is masked from 'package:base':
#> 
#>     trace
k = 4
n = 200
K = matrix(0,k,k)
K[upper.tri(K,diag=T)] = runif(k*(k+1)/2)
K[lower.tri(K,diag=T)] = t(K)[lower.tri(K,diag=T)]
z = rep(1:k,each=n/k)
# mat = gen_adj_sbm(K,z)
# A = mat$A
# P = mat$P
## basic example code
```
