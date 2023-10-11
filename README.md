
<!-- README.md is generated from README.Rmd. Please edit that file -->

# csgc

<!-- badges: start -->
<!-- badges: end -->

The goal of csgc is to perform statistical analysis of network data
based on centred subgraph counts.

## Installation

You can install the development version of csgc R package from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("lishang-stats/csgc")
```

See also "csgc.py" file for the Python version and test code, which
allows GPU acceleration.

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

``` r
set.seed(123)
library(csgc)
#> 
#> Attaching package: 'csgc'
#> The following object is masked from 'package:base':
#> 
#>     trace
k = 4
n = 100
K = matrix(c(0.8, 0.5, 0.1, 0.1,
             0.5, 0.1, 0.1, 0.1,
             0.1, 0.1, 0.8, 0.5,
             0.1, 0.1, 0.5, 0.1),4,4)
z = rep(1:k,each=n/k)
mat = gen_adj_sbm(K,z)
A = mat$A
P = mat$P
```

Then, we calculate the maximum likelihood estimated $K$ and $P$ from $A$
and true labels $z$.

``` r
mat2 = sbm_mle(A,z)
Khat = mat2$K
Phat = mat2$P
```

We can also get the estimated labels given $A$ and number of communities
$k$ using spectral clustering method, and compare the estimated labels
with true labels.

``` r
zhat = spectral_sbm(A,4)
ccrate = ccr(z,zhat)
ccrate
#> [1] 0.71
```

If we use the estimated labels to further get the MLE for $P$, we can
calculate the centred subgraph count statistics and check if these
statistics deviate from 0 or not.

``` r
Phat2 = sbm_mle(A,zhat)$P
t = csgc(A,Phat2,"bernoulli")$t
t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>     21.025806      2.647350     19.910428     -7.621849      1.819534 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>     31.237910      2.883050     16.087230     50.357134     28.512798
```

Sum of squares of these statistics should approximately follow a
chi-square distribution if SBM fits the network data. Suppose we run the
above simulation multiple times, we can get a batch of $p$-values, which
are supposed to follow a uniform distribution. We can calculate the
Wasserstein distance between these $p$-values and a uniform distribution
on \[0,1\].

``` r
p = numeric(100)
for (i in 1:100){
  A1 = gen_adj_sbm(K,z)$A
  t = csgc(A1,sbm_mle(A1,z)$P,"bernoulli")$t
  p[i] = pchisq(sum(t^2), 10, lower.tail = FALSE)
}
wasserstein_uniform(p)
#> [1] 0.1996285
```

If these csgc statistics deviate from 0, we can use our csgc greedy
algorithm to adjust label estimations, and see if correct classification
rate increases.

``` r
out = csgc_greedy(A,zhat,parallel=T)
#> chisq stats: 4486.72334570905, change location: 46
#> adjusted labels: 1111212122121211211111211212222122122222212212212334444444444343443444343443323433342323334332333333
#> chisq stats: 3794.46576351085, change location: 12
#> adjusted labels: 1111212122111211211111211212222122122222212212212334444444444343443444343443323433342323334332333333
#> chisq stats: 3143.59391850798, change location: 17
#> adjusted labels: 1111212122111211111111211212222122122222212212212334444444444343443444343443323433342323334332333333
#> chisq stats: 2607.28078700959, change location: 45
#> adjusted labels: 1111212122111211111111211212222122122222212222212334444444444343443444343443323433342323334332333333
#> chisq stats: 2145.2566815541, change location: 23
#> adjusted labels: 1111212122111211111111111212222122122222212222212334444444444343443444343443323433342323334332333333
#> chisq stats: 1722.4355753283, change location: 67
#> adjusted labels: 1111212122111211111111111212222122122222212222212334444444444343444444343443323433342323334332333333
#> chisq stats: 1433.63752250575, change location: 5
#> adjusted labels: 1111112122111211111111111212222122122222212222212334444444444343444444343443323433342323334332333333
#> chisq stats: 1187.77189810095, change location: 84
#> adjusted labels: 1111112122111211111111111212222122122222212222212334444444444343444444343443323433332323334332333333
#> chisq stats: 979.848831098498, change location: 27
#> adjusted labels: 1111112122111211111111111222222122122222212222212334444444444343444444343443323433332323334332333333
#> chisq stats: 796.193813607788, change location: 14
#> adjusted labels: 1111112122111111111111111222222122122222212222212334444444444343444444343443323433332323334332333333
#> chisq stats: 631.131864281572, change location: 7
#> adjusted labels: 1111111122111111111111111222222122122222212222212334444444444343444444343443323433332323334332333333
#> chisq stats: 497.972661133027, change location: 10
#> adjusted labels: 1111111121111111111111111222222122122222212222212334444444444343444444343443323433332323334332333333
#> chisq stats: 389.124723523015, change location: 91
#> adjusted labels: 1111111121111111111111111222222122122222212222212334444444444343444444343443323433332323333332333333
#> chisq stats: 309.464505727188, change location: 80
#> adjusted labels: 1111111121111111111111111222222122122222212222212334444444444343444444343443323333332323333332333333
#> chisq stats: 247.787921772089, change location: 42
#> adjusted labels: 1111111121111111111111111222222122122222222222212334444444444343444444343443323333332323333332333333
#> chisq stats: 198.380908702665, change location: 51
#> adjusted labels: 1111111121111111111111111222222122122222222222212344444444444343444444343443323333332323333332333333
#> chisq stats: 156.216304713406, change location: 64
#> adjusted labels: 1111111121111111111111111222222122122222222222212344444444444344444444343443323333332323333332333333
#> chisq stats: 120.185170752735, change location: 48
#> adjusted labels: 1111111121111111111111111222222122122222222222222344444444444344444444343443323333332323333332333333
#> chisq stats: 85.6620549432261, change location: 32
#> adjusted labels: 1111111121111111111111111222222222122222222222222344444444444344444444343443323333332323333332333333
#> chisq stats: 57.7334426869356, change location: 71
#> adjusted labels: 1111111121111111111111111222222222122222222222222344444444444344444444443443323333332323333332333333
#> chisq stats: 34.5689706679674, change location: 62
#> adjusted labels: 1111111121111111111111111222222222122222222222222344444444444444444444443443323333332323333332333333
#> chisq stats: 17.7029872605905, change location: 73
#> adjusted labels: 1111111121111111111111111222222222122222222222222344444444444444444444444443323333332323333332333333
#> chisq stats: 11.4538275370245, change location: 50
#> adjusted labels: 1111111121111111111111111222222222122222222222222244444444444444444444444443323333332323333332333333
#> chisq stats: 8.96372989223866, change location: 35
#> adjusted labels: 1111111121111111111111111222222222222222222222222244444444444444444444444443323333332323333332333333
#> chisq stats: 7.00403784223967, change location: 77
#> adjusted labels: 1111111121111111111111111222222222222222222222222244444444444444444444444443423333332323333332333333
#> chisq stats: 5.01317144962853, change location: 94
#> adjusted labels: 1111111121111111111111111222222222222222222222222244444444444444444444444443423333332323333333333333
#> chisq stats: 4.32975497604902, change location: 78
#> adjusted labels: 1111111121111111111111111222222222222222222222222244444444444444444444444443433333332323333333333333
#> chisq stats: 3.41388701672236, change location: 39
#> adjusted labels: 1111111121111111111111111222222222222232222222222244444444444444444444444443433333332323333333333333
#> chisq stats: 3.29350918373957, change location: 9
#> adjusted labels: 1111111111111111111111111222222222222232222222222244444444444444444444444443433333332323333333333333
#> chisq stats: 3.03554642320018, change location: 18
#> adjusted labels: 1111111111111111121111111222222222222232222222222244444444444444444444444443433333332323333333333333
#> chisq stats: 3.05305783065938, change location: 87
#> adjusted labels: 1111111111111111121111111222222222222232222222222244444444444444444444444443433333332333333333333333
zout = out$zout
chisqout = out$chisqout
statsout = out$statsout
chisqout
#> [1] 3.035546
statsout
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>    0.04279361   -0.29176806    1.35066437   -0.47974516    0.43144080 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>   -0.47381646   -0.15771418    0.52246022   -0.42474161    0.07245537
ccr(z,zout)
#> [1] 0.95
```

## csgc: DCSBM

Assume further that $\theta_i$ is the degree parameter for vertex $i$.
Consider the Poisson-type Degree Corrected Stochastic Block Model
![equation](https://latex.codecogs.com/svg.latex?A_%7Bij%7D=A_%7Bji%7D\overset%7Bi.i.d.%7D%7B\sim%7D%20\mathop%7B\mathrm%7BPo%7D%7D\left(\theta_i%20\theta_j%20k_%7Bz_iz_j%7D\right))
for $1\leq i,j\leq n, i\ne j$. We can generate matrix $A$ and $P$ from
$K$ and $\Theta$ and $z$.

``` r
Theta = runif(n,.2,1)
dmat = gen_adj_dcsbm(K,Theta,z)
dA = dmat$A
dP = dmat$P
```

Though no closed form, numerical MLE for $P$ can be obtained, assuming
$A$ and true labels $z$ are known.

``` r
dmat2 = dcsbm_mle(dA,z)
dPhat = dmat2$P
```

Using spectral clustering method for DCSBM, we can also get the
estimated labels and with true labels.

``` r
dzhat = spectral_dcsbm(dA,4)
dccrate = ccr(z,dzhat)
dccrate
#> [1] 0.49
```

## csgc: Compatible with ergm/blockmodels/greed package

Finally, we make the csgc function compatible with other packages so
that people can have easy access to these values. Before that, we will
need the following packages: “ergm”, “blockmodels” and “greed”.

``` r
# install.packages(c("ergm", "blockmodels", "greed", "igraphdata"))
library(ergm)
#> Warning: package 'ergm' was built under R version 4.1.3
#> Loading required package: network
#> Warning: package 'network' was built under R version 4.1.3
#> 
#> 'network' 1.18.1 (2023-01-24), part of the Statnet Project
#> * 'news(package="network")' for changes since last version
#> * 'citation("network")' for citation information
#> * 'https://statnet.org' for help, support, and other information
#> 
#> 'ergm' 4.4.0 (2023-01-26), part of the Statnet Project
#> * 'news(package="ergm")' for changes since last version
#> * 'citation("ergm")' for citation information
#> * 'https://statnet.org' for help, support, and other information
#> 'ergm' 4 is a major update that introduces some backwards-incompatible
#> changes. Please type 'news(package="ergm")' for a list of major
#> changes.
library(blockmodels)
#> Warning: package 'blockmodels' was built under R version 4.1.3
library(greed)
#> Warning: package 'greed' was built under R version 4.1.3
library(igraphdata)
#> Warning: package 'igraphdata' was built under R version 4.1.3
```

For Exponential Random Graph Model, we can extract adjacency matrix (A)
and predicted probability matrix (P), so that csgc statistics can be
calculated.

``` r
data(faux.dixon.high)
fit <- ergm(faux.dixon.high ~ edges + mutual)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
#> Starting Monte Carlo maximum likelihood estimation (MCMLE):
#> Iteration 1 of at most 60:
#> Warning: 'glpk' selected as the solver, but package 'Rglpk' is not available;
#> falling back to 'lpSolveAPI'. This should be fine unless the sample size and/or
#> the number of parameters is very big.
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0304.
#> Convergence test p-value: 0.3016. Not converged with 99% confidence; increasing sample size.
#> Iteration 2 of at most 60:
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0430.
#> Convergence test p-value: 0.0012. Converged with 99% confidence.
#> Finished MCMLE.
#> Evaluating log-likelihood at the estimate. Fitting the dyad-independent submodel...
#> Bridging between the dyad-independent submodel and the full model...
#> Setting up bridge sampling...
#> Using 16 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 .
#> Bridging finished.
#> This model was fit using MCMC.  To examine model diagnostics and check
#> for degeneracy, use the mcmc.diagnostics() function.
out = A_P_from_ergm(model=fit)
#> Please check that the model satisfies the conditional independence of the edges.
#> Warning in A_P_from_ergm(model = fit): Convert directed graph to undirected
#> graph.
csgc(out$A, out$P)$t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>      59.44233      58.13349      52.32496     172.34371     142.62009 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>     288.92686     264.16784     197.18525     494.36571     330.54601
```

For “blockmodels” object, we can easily extract adjacency matrix (A).
The “blockmodels” package will estimate the number of communities (k),
labels (z) and parameter matrix (K). From that, we can get the estimate
probability matrix (P).

``` r
npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
P<-matrix(runif(Q*Q),Q,Q)
M<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z)) ## adjacency matrix
fit <- BM_bernoulli("SBM",M, plotting='')
out = A_P_from_blockmodels(model=fit)
#> Warning in A_P_from_blockmodels(model = fit): Convert directed graph to
#> undirected graph.
csgc(out$A, out$P, out$modeltype)$t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>      84.53258      40.94152     162.60744     326.20236    1896.20998 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>     427.65261    1187.47921     581.77047     683.47687    9040.87610
```

We can also use “greed” package to do model fitting for SBM and DCSBM.

``` r
data(Books)
sbm = greed(Books$X, model = Sbm())
#> 
#> -- Fitting a guess SBM model --
#> 
#> i Initializing a population of 20 solutions.
#> i Generation 1 : best solution with an ICL of -1278 and 7 clusters.
#> i Generation 2 : best solution with an ICL of -1254 and 6 clusters.
#> i Generation 3 : best solution with an ICL of -1254 and 6 clusters.
#> -- Final clustering --
#> 
#> -- Clustering with a SBM model 5 clusters and an ICL of -1254
out1 = A_P_from_greed(data=Books$X, blockmodel=sbm)
csgc(out1$A, out1$P, out1$modeltype)$t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>     2.1604290    11.5804480    13.0885441    -3.5990661     2.7017219 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>    -4.5114451     7.7160612    15.5304911    -0.9049469    -1.3578917
```

``` r
library(igraphdata)
data(karate)
dcsbm = greed(karate, model= DcSbm())
#> Warning: Sbm model used with an igraph object. Vertex and nodes attributes were
#> removed, the clustering will only use the adjacency matrix.
#> 
#> -- Fitting a guess DCSBM model --
#> 
#> i Initializing a population of 20 solutions.
#> i Generation 1 : best solution with an ICL of -229 and 2 clusters.
#> i Generation 2 : best solution with an ICL of -227 and 1 clusters.
#> i Generation 3 : best solution with an ICL of -227 and 1 clusters.
#> -- Final clustering --
#> 
#> -- Clustering with a DCSBM model 1 clusters and an ICL of -227
out2 = A_P_from_greed(data=karate, blockmodel=dcsbm)
csgc(out2$A, out2$P, out2$modeltype)$t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>    -2.4185709    -2.2530913     2.4026465     1.9677440     0.6853484 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>    -1.7426575    -9.2650511    10.3745336     1.3081061     0.9861794
```
