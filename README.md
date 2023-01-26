
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
![equation](https://latex.codecogs.com/png.latex?A_%7Bij%7D=A_%7Bji%7D\overset%7Bi.i.d.%7D%7B\sim%7D%20\mathop%7B\mathrm%7BBe%7D%7D\left(k_%7Bz_iz_j%7D\right))
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
#> [1] 0.75
```

If we use the estimated labels to further get the MLE for $P$, we can
calculate the centred subgraph count statistics and check if these
statistics are deviate from 0 or not.

``` r
Phat2 = sbm_mle(A,zhat)$P
t = csgc(A,Phat2,"binomial")$t
t
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>     19.096477      2.711380     16.798103     -5.050725      3.795722 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>     26.266534      3.224632     13.393930     42.925778     31.525595
```

Sum of squares of these statistics should approximately follow a
chi-square distribution if SBM fits the network data. Suppose we run the
above simulation multiple times, we can get a batch of $p$-values, which
are supposed to follow a uniform distribution. We can calculate the
Wasserstein distance between these $p$-values and a uniform distribution
on \[0,1\].

``` r
p = runif(100)
wasserstein_uniform(p)
#> [1] 0.02510469
```

If these csgc statistics are deviate from 0, we can use our csgc greedy
algorithm to adjust label estimations, and see if correct classification
rate increases.

``` r
out = csgc_greedy(A,zhat,parallel=T)
#> chisq stats: 3477.41824600513, change location: 46
#> adjusted labels: 1111212111121211211111211212222122122222212213212344444444444343443444343443323433342323333332333333
#> chisq stats: 2911.83152192601, change location: 17
#> adjusted labels: 1111212111121211111111211212222122122222212213212344444444444343443444343443323433342323333332333333
#> chisq stats: 2380.32624922314, change location: 67
#> adjusted labels: 1111212111121211111111211212222122122222212213212344444444444343444444343443323433342323333332333333
#> chisq stats: 1914.58029015496, change location: 12
#> adjusted labels: 1111212111111211111111211212222122122222212213212344444444444343444444343443323433342323333332333333
#> chisq stats: 1525.01605625007, change location: 23
#> adjusted labels: 1111212111111211111111111212222122122222212213212344444444444343444444343443323433342323333332333333
#> chisq stats: 1201.77413872641, change location: 45
#> adjusted labels: 1111212111111211111111111212222122122222212223212344444444444343444444343443323433342323333332333333
#> chisq stats: 962.352174941613, change location: 84
#> adjusted labels: 1111212111111211111111111212222122122222212223212344444444444343444444343443323433332323333332333333
#> chisq stats: 761.504965556248, change location: 14
#> adjusted labels: 1111212111111111111111111212222122122222212223212344444444444343444444343443323433332323333332333333
#> chisq stats: 580.225188062817, change location: 5
#> adjusted labels: 1111112111111111111111111212222122122222212223212344444444444343444444343443323433332323333332333333
#> chisq stats: 423.602715956301, change location: 7
#> adjusted labels: 1111111111111111111111111212222122122222212223212344444444444343444444343443323433332323333332333333
#> chisq stats: 319.961580453473, change location: 27
#> adjusted labels: 1111111111111111111111111222222122122222212223212344444444444343444444343443323433332323333332333333
#> chisq stats: 238.523487885342, change location: 80
#> adjusted labels: 1111111111111111111111111222222122122222212223212344444444444343444444343443323333332323333332333333
#> chisq stats: 183.464997091376, change location: 64
#> adjusted labels: 1111111111111111111111111222222122122222212223212344444444444344444444343443323333332323333332333333
#> chisq stats: 136.063695659569, change location: 71
#> adjusted labels: 1111111111111111111111111222222122122222212223212344444444444344444444443443323333332323333332333333
#> chisq stats: 105.665342029417, change location: 62
#> adjusted labels: 1111111111111111111111111222222122122222212223212344444444444444444444443443323333332323333332333333
#> chisq stats: 79.3813784141062, change location: 48
#> adjusted labels: 1111111111111111111111111222222122122222212223222344444444444444444444443443323333332323333332333333
#> chisq stats: 51.9617173594888, change location: 32
#> adjusted labels: 1111111111111111111111111222222222122222212223222344444444444444444444443443323333332323333332333333
#> chisq stats: 31.2018515423755, change location: 73
#> adjusted labels: 1111111111111111111111111222222222122222212223222344444444444444444444444443323333332323333332333333
#> chisq stats: 22.5627662840749, change location: 42
#> adjusted labels: 1111111111111111111111111222222222122222222223222344444444444444444444444443323333332323333332333333
#> chisq stats: 13.3587628725007, change location: 50
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444443323333332323333332333333
#> chisq stats: 9.1939554248852, change location: 78
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444443333333332323333332333333
#> chisq stats: 5.90366384688834, change location: 74
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444343333333332323333332333333
#> chisq stats: 4.90820415913339, change location: 77
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444343433333332323333332333333
#> chisq stats: 3.96484549348372, change location: 94
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444343433333332323333333333333
#> chisq stats: 3.75748103193621, change location: 85
#> adjusted labels: 1111111111111111111111111222222222122222222223222244444444444444444444444343433333333323333333333333
#> chisq stats: 3.08379262648469, change location: 35
#> adjusted labels: 1111111111111111111111111222222222222222222223222244444444444444444444444343433333333323333333333333
#> chisq stats: 3.11013929518305, change location: 39
#> adjusted labels: 1111111111111111111111111222222222222232222223222244444444444444444444444343433333333323333333333333
zout = out$zout
chisqout = out$chisqout
statsout = out$statsout
chisqout
#> [1] 3.083793
statsout
#>     t.twostar    t.triangle   t.fourcycle   t.threepath   t.threestar 
#>    0.22174965   -0.08725584    1.08678563    0.51817374    0.56995710 
#> t.triangleapp t.twotriangle   t.fivecycle    t.fourpath    t.fourstar 
#>   -0.13082998   -0.43474782    0.40048793   -0.91447892    0.22307754
ccr(z,zout)
#> [1] 0.96
```

## csgc: DCSBM

Assume further that $\theta_i$ is the degree parameter for vertex $i$.
Consider the Poisson-type Degree Corrected Stochastic Block Model
![equation](https://latex.codecogs.com/png.latex?A_%7Bij%7D=A_%7Bji%7D\overset%7Bi.i.d.%7D%7B\sim%7D%20\mathop%7B\mathrm%7BPo%7D%7D\left(\theta_i%20\theta_j%20k_%7Bz_iz_j%7D\right))
for $1\leq i,j\leq n, i\ne j$. We can generate matrix $A$ and $P$ from
$K$ and $\Theta$ and $z$.

``` r
Theta = runif(n,.8,1)
dmat = gen_adj_dcsbm(K,Theta,z)
dA = dmat$A
dP = dmat$P
```

Though no closed form, numerical MLE for $K$ and $P$ can be obtained,
assuming $A$ and true labels $z$ are known.

``` r
dmat2 = dcsbm_mle(dA,z)
dKhat = mat2$K
dPhat = mat2$P
```

Using spectral clustering method for DCSBM, we can also get the
estimated labels and with true labels.

``` r
dzhat = spectral_dcsbm(dA,4)
dccrate = ccr(z,dzhat)
dccrate
#> [1] 0.51
```
