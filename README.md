
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BiDirectCausal

<!-- badges: start -->

<!-- badges: end -->

The goal of R package **BiDirectCausal** is to inferring possibly
bi-directional causal effects between two traits with GWAS summary data.

## Installation

You can install BiDirectCausal from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xue-hr/BiDirectCausal")
```

## Example

Here is an example showing how to apply different bi-direction methods
to simulated data. We first generate data using our main simulation
setup, which could be found in our draft.

``` r
n = 50000 # sample size
p.snp.X = 15 # number of SNPs in g_X
p.snp.Y = 10 # number of SNPs in g_Y
p.snp.U = 10 # number of SNPs in g_B
xi.range = 0.2 # generate non-zero xi's from uniform distribution, implying correlated pleiotropy
theta.XtoY = 0.1 # causal effect from X to Y
theta.YtoX = 0 # causal effect from Y to X
p.all = p.snp.X + p.snp.Y + p.snp.U

set.seed(1)
### generate alpha
alpha.effect = runif(p.snp.X,0.2,0.3)*(rbinom(p.snp.X,1,0.5)*2-1)

### generate beta
beta.effect = runif(p.snp.Y,0.2,0.3)*(rbinom(p.snp.Y,1,0.5)*2-1)

### generate effects of invalid IVs
gamma.effect = runif(p.snp.U,0.2,0.3)*(rbinom(p.snp.U,1,0.5)*2-1)
eta.effect = runif(p.snp.U,0.2,0.3)*(rbinom(p.snp.U,1,0.5)*2-1)
xi.effect = runif(p.snp.U,-xi.range,xi.range)


### Generate Individual Data
MAF = rep(0.3,p.all)
Z = matrix(rbinom(2*n*p.all,2,MAF),ncol = p.all,byrow = T)
U = Z%*%c(rep(0,p.snp.X),rep(0,p.snp.Y),xi.effect) + rnorm(2*n,0,sqrt(2))
error_X = rnorm(2*n,0,1)
error_Y = rnorm(2*n,0,1)

X = 
  (Z%*%c(alpha.effect,theta.YtoX*beta.effect,gamma.effect + theta.YtoX*eta.effect) + 
     (1+theta.YtoX)*U + error_X + theta.YtoX*error_Y)/(1-theta.XtoY*theta.YtoX)

Y = 
  (Z%*%c(theta.XtoY*alpha.effect,beta.effect,theta.XtoY*gamma.effect + eta.effect) + 
     (1+theta.XtoY)*U + theta.XtoY*error_X + error_Y)/(1-theta.XtoY*theta.YtoX)

### Generate two independent samples (Z1,X1,Y1) and (Z2,X2,Y2)
Z1 = Z[1:n,]
X1 = as.matrix(X[1:n])
Y1 = as.matrix(Y[1:n])

Z2 = Z[(n+1):(2*n),]
X2 = as.matrix(X[(n+1):(2*n)])
Y2 = as.matrix(Y[(n+1):(2*n)])

Z1 = scale(Z1,scale = F)
X1 = scale(X1,scale = F)
Y1 = scale(Y1,scale = F)

### Get summary statistics for X from marginal linear regression
z1Tz1 = t(Z1)%*%Z1
b_X = as.numeric(1/diag(z1Tz1)*(t(Z1)%*%X1))
rep_X1 = X1[,rep(1,p.all)]
se_X = 
  sqrt(colSums((rep_X1 - Z1%*%diag(b_X))^2)/(n-2)/diag(z1Tz1))

### Get summary statistics for Y from marginal linear regression
Z2 = scale(Z2,scale = F)
X2 = scale(X2,scale = F)
Y2 = scale(Y2,scale = F)

z2Tz2 = t(Z2)%*%Z2
b_Y = as.numeric(1/diag(z2Tz2)*(t(Z2)%*%Y2))
rep_Y2 = Y2[,rep(1,p.all)]
se_Y = 
  sqrt(colSums((rep_Y2 - Z2%*%diag(b_Y))^2)/(n-2)/diag(z2Tz2))
```

Now we have GWAS summary statistics `b_X` is the vector of estimated
effect sizes for trait **X**, `se_X` is the vector of standard errors of
`b_X`; `b_Y` is the vector of estimated effect sizes for trait **Y**,
`se_Y` is the vector of standard errors of `b_Y`. Now we apply our
bi-directional CDcML methods:

``` r
library(MRcML)
library(MRCD)
library(BiDirectCausal)
BiDirCDcML(b_X = b_X,
           b_Y = b_Y,
           se_X = se_X,
           se_Y = se_Y,
           n_X = n,
           n_Y = n,
           sig.cutoff = 0.05/35, 
           num_pert = 100,random_start = 0, random.seed = 1)
#> $XtoY.est.NoS
#> [1] 0.1042717
#> 
#> $XtoY.se.NoS
#> [1] 0.01349092
#> 
#> $YtoX.est.NoS
#> [1] -0.0130484
#> 
#> $YtoX.se.NoS
#> [1] 0.01733937
#> 
#> $XtoY.est.S
#> [1] 0.1042717
#> 
#> $XtoY.se.S
#> [1] 0.01349092
#> 
#> $YtoX.est.S
#> [1] -0.01304843
#> 
#> $YtoX.se.S
#> [1] 0.01733937
#> 
#> $XtoY.est.NoS.DP
#> [1] 0.1067298
#> 
#> $XtoY.se.NoS.DP
#> [1] 0.01261763
#> 
#> $YtoX.est.NoS.DP
#> [1] -0.01223717
#> 
#> $YtoX.se.NoS.DP
#> [1] 0.01609297
#> 
#> $XtoY.est.S.DP
#> [1] 0.1048932
#> 
#> $XtoY.se.S.DP
#> [1] 0.01390019
#> 
#> $YtoX.est.S.DP
#> [1] -0.0120089
#> 
#> $YtoX.se.S.DP
#> [1] 0.01674305
```

Apply bi-directional MRcML methods:

``` r
BiDirMRcML(b_X = b_X,
           b_Y = b_Y,
           se_X = se_X,
           se_Y = se_Y,
           n_X = n,
           n_Y = n,
           sig.cutoff = 0.05/35, 
           num_pert = 100,random_start = 0, random.seed = 1)
#> $XtoY.est.NoS
#> [1] 0.1089637
#> 
#> $XtoY.se.NoS
#> [1] 0.01410568
#> 
#> $YtoX.est.NoS
#> [1] -0.01242067
#> 
#> $YtoX.se.NoS
#> [1] 0.0165682
#> 
#> $XtoY.est.S
#> [1] 0.1089637
#> 
#> $XtoY.se.S
#> [1] 0.01410568
#> 
#> $YtoX.est.S
#> [1] -0.01242066
#> 
#> $YtoX.se.S
#> [1] 0.0165682
#> 
#> $XtoY.est.NoS.DP
#> [1] 0.1115354
#> 
#> $XtoY.se.NoS.DP
#> [1] 0.01319329
#> 
#> $YtoX.est.NoS.DP
#> [1] -0.01398398
#> 
#> $YtoX.se.NoS.DP
#> [1] 0.01642966
#> 
#> $XtoY.est.S.DP
#> [1] 0.1114202
#> 
#> $XtoY.se.S.DP
#> [1] 0.01250001
#> 
#> $YtoX.est.S.DP
#> [1] -0.0112045
#> 
#> $YtoX.se.S.DP
#> [1] 0.01736859
```

Apply bi-directional CD methods:

``` r
BiDirCDMethod(b_X = b_X,
              b_Y = b_Y,
              se_X = se_X,
              se_Y = se_Y,
              n_X = n,
              n_Y = n,
              sig.cutoff = 0.05/35, random.seed = 1)
#> $CDRatio.XtoY.est.NoS
#>          K 
#> 0.04681161 
#> 
#> $CDRatio.XtoY.se.NoS
#>      se(K) 
#> 0.01129413 
#> 
#> $CDEgger.XtoY.est.NoS
#>           K 
#> 0.006165783 
#> 
#> $CDEgger.XtoY.se.NoS
#>     se(K) 
#> 0.1114091 
#> 
#> $CDRatio.YtoX.est.NoS
#>            K 
#> -0.007732285 
#> 
#> $CDRatio.YtoX.se.NoS
#>      se(K) 
#> 0.01405544 
#> 
#> $CDEgger.YtoX.est.NoS
#>           K 
#> -0.06710133 
#> 
#> $CDEgger.YtoX.se.NoS
#>     se(K) 
#> 0.1868871 
#> 
#> $CDRatio.XtoY.est.S
#>          K 
#> 0.03687256 
#> 
#> $CDRatio.XtoY.se.S
#>      se(K) 
#> 0.01176274 
#> 
#> $CDEgger.XtoY.est.S
#>          K 
#> 0.01614943 
#> 
#> $CDEgger.XtoY.se.S
#>      se(K) 
#> 0.06214906 
#> 
#> $CDRatio.YtoX.est.S
#>           K 
#> 0.001098699 
#> 
#> $CDRatio.YtoX.se.S
#>      se(K) 
#> 0.01431168 
#> 
#> $CDEgger.YtoX.est.S
#>           K 
#> -0.03139753 
#> 
#> $CDEgger.YtoX.se.S
#>    se(K) 
#> 0.124538
```

With `est` and corresponding `se` we could get the p-value. For example,
for **X** to **Y**, CD-cML-DP-S gives estimate 0.1045348 and standard
error 0.01432697, we could calculate p-value as

``` r
pnorm(-abs(0.1045348/0.01432697))*2
#> [1] 2.956466e-13
```

## Contact

Our draft is under review, please contact Haoran Xue at
<xuexx268@umn.edu> for any questions or comments\!
