
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

## csgc: SBM

Assume an undirected unweighted network has $n$ vertices, with $A$ be
the adjacency matrix and $P$ be the corresponding probability matrix.
The parameter matrix $K$ (a $k$ by $k$ matrix) presents the connection
probability within/between each block. Let $z$ be the label vector with
length $n$, where each entry can be from $1$ to $k$. Consider the
Bernoulli-type Stochastic Block Model
![equation](https://latex.codecogs.com/svg.latex?A_%7Bij%7D=A_%7Bji%7D\overset%7Bi.i.d.%7D%7B\sim%7D%20\mathop%7B\mathrm%7BBe%7D%7D\left(k_%7Bz_iz_j%7D\right))
for $1\leq i,j\leq n, i\ne j$. We first generate matrix $A$ and $P$ from
$K$ and $z$.
<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

<!-- ```{r gen_sbm} -->
<!-- set.seed(123) -->
<!-- library(csgc) -->
<!-- k = 4 -->
<!-- n = 100 -->
<!-- K = matrix(c(0.8, 0.5, 0.1, 0.1, -->
<!--              0.5, 0.1, 0.1, 0.1, -->
<!--              0.1, 0.1, 0.8, 0.5, -->
<!--              0.1, 0.1, 0.5, 0.1),4,4) -->
<!-- z = rep(1:k,each=n/k) -->
<!-- mat = gen_adj_sbm(K,z) -->
<!-- A = mat$A -->
<!-- P = mat$P -->
<!-- ``` -->
<!-- Then, we calculate the maximum likelihood estimated $K$ and $P$ from $A$ and true labels $z$. -->
<!-- ```{r sbm_mle} -->
<!-- mat2 = sbm_mle(A,z) -->
<!-- Khat = mat2$K -->
<!-- Phat = mat2$P -->
<!-- ``` -->
<!-- We can also get the estimated labels given $A$ and number of communities $k$ using spectral clustering method, and compare the estimated labels with true labels. -->
<!-- ```{r spectral} -->
<!-- zhat = spectral_sbm(A,4) -->
<!-- ccrate = ccr(z,zhat) -->
<!-- ccrate -->
<!-- ``` -->
<!-- If we use the estimated labels to further get the MLE for $P$, we can calculate the centred subgraph count statistics and check if these statistics are deviate from 0 or not. -->
<!-- ```{r csgc} -->
<!-- Phat2 = sbm_mle(A,zhat)$P -->
<!-- t = csgc(A,Phat2,"binomial")$t -->
<!-- t -->
<!-- ``` -->
<!-- Sum of squares of these statistics should approximately follow a chi-square distribution if SBM fits the network data. Suppose we run the above simulation multiple times, we can get a batch of $p$-values, which are supposed to follow a uniform distribution. We can calculate the Wasserstein distance between these $p$-values and a uniform distribution on [0,1]. -->
<!-- ```{r dw} -->
<!-- p = runif(100) -->
<!-- wasserstein_uniform(p) -->
<!-- ``` -->
<!-- If these csgc statistics are deviate from 0, we can use our csgc greedy algorithm to adjust label estimations, and see if correct classification rate increases. -->
<!-- ```{r greedy} -->
<!-- out = csgc_greedy(A,zhat,parallel=T) -->
<!-- zout = out$zout -->
<!-- chisqout = out$chisqout -->
<!-- statsout = out$statsout -->
<!-- chisqout -->
<!-- statsout -->
<!-- ccr(z,zout) -->
<!-- ``` -->
<!-- ## csgc: DCSBM -->
<!-- Assume further that $\theta_i$ is the degree parameter for vertex $i$. Consider the Poisson-type Degree Corrected Stochastic Block Model -->
<!-- ![equation](https://latex.codecogs.com/svg.latex?A_{ij}=A_{ji}\overset{i.i.d.}{\sim}%20\mathop{\mathrm{Po}}\left(\theta_i%20\theta_j%20k_{z_iz_j}\right)) -->
<!-- for $1\leq i,j\leq n, i\ne j$. We can generate matrix $A$ and $P$ from $K$ and $\Theta$ and $z$. -->
<!-- ```{r gen_dcsbm} -->
<!-- Theta = runif(n,.8,1) -->
<!-- dmat = gen_adj_dcsbm(K,Theta,z) -->
<!-- dA = dmat$A -->
<!-- dP = dmat$P -->
<!-- ``` -->
<!-- Though no closed form, numerical MLE for $K$ and $P$ can be obtained, assuming $A$ and true labels $z$ are known. -->
<!-- ```{r dcsbm_mle} -->
<!-- dmat2 = dcsbm_mle(dA,z) -->
<!-- dKhat = mat2$K -->
<!-- dPhat = mat2$P -->
<!-- ``` -->
<!-- Using spectral clustering method for DCSBM, we can also get the estimated labels and with true labels. -->
<!-- ```{r dcspectral} -->
<!-- dzhat = spectral_dcsbm(dA,4) -->
<!-- dccrate = ccr(z,dzhat) -->
<!-- dccrate -->
<!-- ``` -->
