# Example 2


- [Joint structural equation model: Latent growth with strict invariance
  model](#joint-structural-equation-model-latent-growth-with-strict-invariance-model)
- [Measurement (Free Structural)
  Model](#measurement-free-structural-model)
  - [Partial Invariance](#partial-invariance)
- [Two-Stage Path Analysis](#two-stage-path-analysis)
  - [With Regression factor scores](#with-regression-factor-scores)
  - [With Bartlett factor scores](#with-bartlett-factor-scores)
- [Compare to Using Just Factor
  Scores](#compare-to-using-just-factor-scores)

``` r
library(lavaan)
```

    This is lavaan 0.6-17
    lavaan is FREE software! Please report any bugs.

``` r
# install.packages("remotes")
# remotes::install_github("Gengrui-Zhang/R2spa")
library(R2spa)
library(semlrtp)  # for likelihood ratio test; can be obtained from https://github.com/sfcheung/semlrtp
```

We use the example from [this
tutorial](https://quantdev.ssri.psu.edu/tutorials/growth-modeling-chapter-14-modeling-change-latent-variables-measured-continuous)
from Chapter 14 of the book *Growth Modeling* (Grimm, Ram & Estabrook,
2017) to demonstrate how to perform linear growth modeling with
two-stage path analysis (2S-PA)

``` r
# Load data
eclsk <- read.csv(
    "https://quantdev.ssri.psu.edu/sites/qdev/files/ECLS_Science.csv",
    header = TRUE
)
```

## Joint structural equation model: Latent growth with strict invariance model

In the
[tutorial](https://quantdev.ssri.psu.edu/tutorials/growth-modeling-chapter-14-modeling-change-latent-variables-measured-continuous),
the authors first performed longitudinal invariance testing and found
support for strict invariance. They moved on to fit a latent growth
model based on the strict invariance model, as shown below.

``` r
measurement_mod <- "
# factor loadings (with constraints)
eta1 =~ l1 * s_g3 + l2 * r_g3 + l3 * m_g3
eta2 =~ l1 * s_g5 + l2 * r_g5 + l3 * m_g5
eta3 =~ l1 * s_g8 + l2 * r_g8 + l3 * m_g8

# unique variances/covariances 
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# observed variable intercepts
s_g3 ~ i1 * 1
s_g5 ~ i1 * 1
s_g8 ~ i1 * 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"

jsem_growth_mod <- paste0(measurement_mod, "
# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

# latent basis model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

i ~~ vi * i + start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 0 * 1
s ~ 1

# standardize first time point
1 == psi + vi
")
jsem_growth_fit <- sem(jsem_growth_mod,
                       data = eclsk,
                       meanstructure = TRUE,
                       auto.fix.first = FALSE,
                       estimator = "ML",
                       missing = "listwise")
```

Intercept-only model

``` r
jsem_location_mod <- paste0(measurement_mod, "
# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

# location model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3

i ~~ vi * i + start(.8) * i

i ~ 0 * 1

# standardize first time point
1 == psi + vi
")
jsem_location_fit <- sem(jsem_location_mod,
                         data = eclsk,
                         meanstructure = TRUE,
                         auto.fix.first = FALSE,
                         estimator = "ML",
                         missing = "listwise")
```

## Measurement (Free Structural) Model

``` r
strict_mod <- paste0(measurement_mod, "
# factor variances
eta1 ~~ 1 * eta1 + eta2 + eta3
eta2 ~~ eta2 + eta3
eta3 ~~ eta3

# latent variable intercepts
eta1 ~ 0 * 1
eta2 ~ 1
eta3 ~ 1
")
strict_fit <- lavaan(strict_mod,
                     data = eclsk,
                     meanstructure = TRUE,
                     estimator = "ML",
                     missing = "listwise")
```

Interpretational Confounding

- The definition of the latent variables (i.e., the loadings) changes
  depending on the structural model

``` r
cbind(lavInspect(jsem_growth_fit, what = "est")$lambda[1:3, 1],
      lavInspect(jsem_location_fit, what = "est")$lambda[1:3, 1],
      lavInspect(strict_fit, what = "est")$lambda[1:3, 1]) |>
    `colnames<-`(c("latent basis", "location", "measurement only")) |>
    knitr::kable()
```

|      | latent basis | location | measurement only |
|:-----|-------------:|---------:|-----------------:|
| s_g3 |     14.87219 | 18.57194 |         14.82956 |
| r_g3 |     21.47327 | 28.18858 |         21.39170 |
| m_g3 |     20.20097 | 25.93159 |         20.10590 |

### Partial Invariance

Modification indices on intercept constraints

``` r
modindices(strict_fit, op = "~1", free.remove = FALSE)
```

    Warning in modindices(strict_fit, op = "~1", free.remove = FALSE): lavaan WARNING: the modindices() function ignores equality constraints;
              use lavTestScore() to assess the impact of releasing one 
              or multiple constraints

        lhs op rhs      mi    epc sepc.lv sepc.all sepc.nox
    28 s_g3 ~1      10.563  0.825   0.825    0.050    0.050
    29 s_g5 ~1     121.492 -2.738  -2.738   -0.166   -0.166
    30 s_g8 ~1      86.851  2.751   2.751    0.168    0.168
    31 r_g3 ~1       2.758  0.716   0.716    0.028    0.028
    32 r_g5 ~1       5.306  0.953   0.953    0.037    0.037
    33 r_g8 ~1      19.558 -2.076  -2.076   -0.081   -0.081
    34 m_g3 ~1      22.922 -1.693  -1.693   -0.071   -0.071
    35 m_g5 ~1      78.050  2.912   2.912    0.122    0.122
    36 m_g8 ~1      25.706 -1.938  -1.938   -0.082   -0.082
    44 eta2 ~1       0.000  0.000   0.000    0.000    0.000
    45 eta3 ~1       0.000  0.000   0.000    0.000    0.000

The measurement intercept for `s_g5` at Time 2 has a large MI, and that
for `s_g8` at Time 3 is big too. This suggests the science indicator may
not be invariant over time, so we free the intercepts for those.

``` r
measurement2_mod <- "
# factor loadings (with constraints)
eta1 =~ l1 * s_g3 + l2 * r_g3 + l3 * m_g3
eta2 =~ l1 * s_g5 + l2 * r_g5 + l3 * m_g5
eta3 =~ l1 * s_g8 + l2 * r_g8 + l3 * m_g8

# unique variances/covariances 
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# observed variable intercepts
s_g3 ~ 1
s_g5 ~ 1
s_g8 ~ 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"
pscalar_mod <- paste0(measurement2_mod, "
# factor variances
eta1 ~~ 1 * eta1 + eta2 + eta3
eta2 ~~ eta2 + eta3
eta3 ~~ eta3

# latent variable intercepts
eta1 ~ 0 * 1
eta2 ~ 1
eta3 ~ 1
")
pscalar_fit <- lavaan(pscalar_mod,
                      data = eclsk,
                      meanstructure = TRUE,
                      estimator = "ML",
                      missing = "listwise")
anova(strict_fit, pscalar_fit)
```


    Chi-Squared Difference Test

                Df   AIC   BIC  Chisq Chisq diff   RMSEA Df diff Pr(>Chisq)    
    pscalar_fit 27 62878 63007 205.67                                          
    strict_fit  29 63163 63283 494.84     289.17 0.40211       2  < 2.2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Second-order growth model with partial invariance

``` r
# Constrain one loading as done in Grimm et al. (2017)
measurement3_mod <- "
# factor loadings (with constraints)
eta1 =~ l1 * s_g3 + 22.603 * r_g3 + l3 * m_g3
eta2 =~ l1 * s_g5 + 22.603 * r_g5 + l3 * m_g5
eta3 =~ l1 * s_g8 + 22.603 * r_g8 + l3 * m_g8

# unique variances/covariances 
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# observed variable intercepts
s_g3 ~ 1
s_g5 ~ 1
s_g8 ~ 1
r_g3 ~ 130.453 * 1
r_g5 ~ 130.453 * 1
r_g8 ~ 130.453 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"
jsem_growth2_mod <- paste0(measurement3_mod, "
# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

# latent basis model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 1
s ~ 1
")
jsem_growth2_fit <- sem(jsem_growth2_mod,
                        data = eclsk,
                        meanstructure = TRUE,
                        auto.fix.first = FALSE,
                        estimator = "ML",
                        missing = "listwise")
```

## Two-Stage Path Analysis

### With Regression factor scores

Obtain factor scores using `R2spa::get_fs_lavaan()`, and the needed
input for `R2spa::tspa()`

``` r
fs_rs <- R2spa::get_fs_lavaan(pscalar_fit)
fs_loading <- attr(fs_rs, which = "fsL")
fs_errorvar <- attr(fs_rs, which = "fsT")
fs_intercept <- attr(fs_rs, which = "fsb")
```

Alternatively, can be done by hand

``` r
fs_rs <- lavPredict(pscalar_fit, acov = TRUE, fsm = TRUE)
# Obtain loadings, error covariances, and intercepts
pars <- lavInspect(pscalar_fit, what = "est")
psi <- pars$psi
alpha <- pars$alpha
fs_loading <- diag(3) - attr(fs_rs, which = "acov")[[1]] %*% solve(psi)
fs_names <- paste0("fs_", c("eta1", "eta2", "eta3"))
rownames(fs_loading) <- fs_names
fs_errorvar <- fs_loading %*% attr(fs_rs, which = "acov")[[1]]
dimnames(fs_errorvar) <- rep(list(fs_names), 2)
fs_intercept <- c(alpha - fs_loading %*% alpha)
names(fs_intercept) <- fs_names
colnames(fs_rs)[1:3] <- fs_names
```

Growth Model with 2S-PA

``` r
# Growth model
growth_mod <- "
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 1
s ~ 1
"
# Fit the growth model
tspa_growth_fit <- tspa(growth_mod, data = data.frame(fs_rs),
                        fsL = fs_loading,
                        fsT = fs_errorvar,
                        fsb = fs_intercept,
                        estimator = "ML")
summary(tspa_growth_fit)
```

    lavaan 0.6.17 ended normally after 28 iterations

      Estimator                                         ML
      Optimization method                           NLMINB
      Number of model parameters                         9
      Number of equality constraints                     2

      Number of observations                           888

    Model Test User Model:
                                                          
      Test statistic                                16.882
      Degrees of freedom                                 2
      P-value (Chi-square)                           0.000

    Parameter Estimates:

      Standard errors                             Standard
      Information                                 Expected
      Information saturated (h1) model          Structured

    Latent Variables:
                       Estimate  Std.Err  z-value  P(>|z|)
      eta1 =~                                             
        fs_eta1           0.648                           
        fs_eta2           0.219                           
        fs_eta3           0.109                           
      eta2 =~                                             
        fs_eta1           0.161                           
        fs_eta2           0.519                           
        fs_eta3           0.146                           
      eta3 =~                                             
        fs_eta1           0.126                           
        fs_eta2           0.197                           
        fs_eta3           0.671                           
      i =~                                                
        eta1              1.000                           
        eta2              1.000                           
        eta3              1.000                           
      s =~                                                
        eta1              0.000                           
        eta2              0.576    0.006   92.757    0.000
        eta3              1.000                           

    Covariances:
                       Estimate  Std.Err  z-value  P(>|z|)
     .fs_eta1 ~~                                          
       .fs_eta2           0.061                           
       .fs_eta3           0.050                           
     .fs_eta2 ~~                                          
       .fs_eta3           0.055                           
      i ~~                                                
        s                -0.076    0.021   -3.706    0.000

    Intercepts:
                       Estimate  Std.Err  z-value  P(>|z|)
       .fs_eta1          -0.412                           
       .fs_eta2           0.152                           
       .fs_eta3           0.458                           
        i                 0.001    0.035    0.022    0.982
        s                 1.874    0.018  101.860    0.000

    Variances:
                       Estimate  Std.Err  z-value  P(>|z|)
       .fs_eta1           0.072                           
       .fs_eta2           0.066                           
       .fs_eta3           0.070                           
       .eta1     (psi)    0.029    0.004    7.780    0.000
       .eta2     (psi)    0.029    0.004    7.780    0.000
       .eta3     (psi)    0.029    0.004    7.780    0.000
        i                 0.979    0.053   18.467    0.000
        s                 0.100    0.016    6.256    0.000

### With Bartlett factor scores

``` r
fs_bs <- R2spa::get_fs_lavaan(pscalar_fit, method = "Bartlett")
fs_loading <- attr(fs_bs, which = "fsL")
fs_errorvar <- attr(fs_bs, which = "fsT")
fs_intercept <- attr(fs_bs, which = "fsb")
```

Alternatively, can be done by hand

``` r
fs_bs <- lavPredict(pscalar_fit, acov = TRUE, fsm = TRUE,
                    method = "Bartlett")
# Obtain loadings, error covariances, and intercepts
pars <- lavInspect(pscalar_fit, what = "est")
psi <- pars$psi
alpha <- pars$alpha
fs_loading <- diag(3)  # property of Bartlett score
fs_names <- paste0("fs_", c("eta1", "eta2", "eta3"))
rownames(fs_loading) <- fs_names
fs_errorvar <- attr(fs_bs, which = "acov")[[1]]
dimnames(fs_errorvar) <- rep(list(fs_names), 2)
fs_intercept <- rep(0, 3)  # property of Bartlett score
names(fs_intercept) <- fs_names
colnames(fs_bs)[1:3] <- fs_names
```

Growth Model with 2S-PA

``` r
# Fit the growth model
tspa_growth2_fit <- tspa(growth_mod, data = data.frame(fs_bs),
                         fsL = fs_loading,
                         fsT = fs_errorvar,
                         fsb = fs_intercept,
                         estimator = "ML")
summary(tspa_growth2_fit)
```

    lavaan 0.6.17 ended normally after 28 iterations

      Estimator                                         ML
      Optimization method                           NLMINB
      Number of model parameters                         9
      Number of equality constraints                     2

      Number of observations                           888

    Model Test User Model:
                                                          
      Test statistic                                16.882
      Degrees of freedom                                 2
      P-value (Chi-square)                           0.000

    Parameter Estimates:

      Standard errors                             Standard
      Information                                 Expected
      Information saturated (h1) model          Structured

    Latent Variables:
                       Estimate  Std.Err  z-value  P(>|z|)
      eta1 =~                                             
        fs_eta1           1.000                           
        fs_eta2          -0.000                           
        fs_eta3           0.000                           
      eta2 =~                                             
        fs_eta1          -0.000                           
        fs_eta2           1.000                           
        fs_eta3           0.000                           
      eta3 =~                                             
        fs_eta1          -0.000                           
        fs_eta2          -0.000                           
        fs_eta3           1.000                           
      i =~                                                
        eta1              1.000                           
        eta2              1.000                           
        eta3              1.000                           
      s =~                                                
        eta1              0.000                           
        eta2              0.576    0.006   92.757    0.000
        eta3              1.000                           

    Covariances:
                       Estimate  Std.Err  z-value  P(>|z|)
     .fs_eta1 ~~                                          
       .fs_eta2           0.061                           
       .fs_eta3           0.038                           
     .fs_eta2 ~~                                          
       .fs_eta3           0.049                           
      i ~~                                                
        s                -0.076    0.021   -3.706    0.000

    Intercepts:
                       Estimate  Std.Err  z-value  P(>|z|)
       .fs_eta1           0.000                           
       .fs_eta2           0.000                           
       .fs_eta3           0.000                           
        i                 0.001    0.035    0.022    0.982
        s                 1.874    0.018  101.860    0.000

    Variances:
                       Estimate  Std.Err  z-value  P(>|z|)
       .fs_eta1           0.110                           
       .fs_eta2           0.110                           
       .fs_eta3           0.109                           
       .eta1     (psi)    0.029    0.004    7.780    0.000
       .eta2     (psi)    0.029    0.004    7.780    0.000
       .eta3     (psi)    0.029    0.004    7.780    0.000
        i                 0.979    0.053   18.467    0.000
        s                 0.100    0.016    6.256    0.000

## Compare to Using Just Factor Scores

``` r
fs_growth_mod <- "
i =~ 1 * fs_eta1 + 1 * fs_eta2 + 1 * fs_eta3
s =~ 0 * fs_eta1 + start(.5) * fs_eta2 + 1 * fs_eta3

# factor variances
fs_eta1 ~~ psi * fs_eta1
fs_eta2 ~~ psi * fs_eta2
fs_eta3 ~~ psi * fs_eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 1
s ~ 1
"
growth_fit <- growth(fs_growth_mod,
                     data = data.frame(fs_rs),
                     meanstructure = TRUE,
                     estimator = "ML",
                     missing = "listwise")
growth2_fit <- growth(fs_growth_mod,
                      data = data.frame(fs_bs),
                      meanstructure = TRUE,
                      estimator = "ML",
                      missing = "listwise")
```

The table below shows the difference in estimates and test statistics:

``` r
options(knitr.kable.NA = "")
rbind(as.data.frame(
        lrtp(jsem_growth2_fit, op = c("~~", "~1"))
    )[c(50, 47), c("est", "se", "Chisq")],
    as.data.frame(
        lrtp(tspa_growth_fit, op = c("~~", "~1"))
    )[c(32, 29), c("est", "se", "Chisq")],
    as.data.frame(
        lrtp(tspa_growth2_fit, op = c("~~", "~1"))
    )[c(32, 29), c("est", "se", "Chisq")],
    as.data.frame(
        lrtp(growth_fit, op = c("~~", "~1"))
    )[c(14, 11), c("est", "se", "Chisq")],
    as.data.frame(
        lrtp(growth2_fit, op = c("~~", "~1"))
    )[c(14, 11), c("est", "se", "Chisq")]) |>
    base::`[`(c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10), ) |>
    cbind(Parameter = c("Mean slope", rep("", 4),
                        "Var slope", rep("", 4)),
          Model = rep(c("JSEM", "2S-PA (Reg)", "2S-PA (Bart)", "FS (Reg)", "FS (Bart)"), 2)) |>
    base::`[`(, c(4:5, 1:3)) |>
    knitr::kable(col.names = c("Parameter", "Model", "Est", "SE", "$\\chi^2$"), row.names = FALSE, digits = 3)
```

| Parameter  | Model        |   Est |    SE | $\chi^2$ |
|:-----------|:-------------|------:|------:|---------:|
| Mean slope | JSEM         | 1.873 | 0.025 | 2223.513 |
|            | 2S-PA (Reg)  | 1.874 | 0.018 | 2271.428 |
|            | 2S-PA (Bart) | 1.874 | 0.018 | 2271.428 |
|            | FS (Reg)     | 1.874 | 0.010 | 3282.137 |
|            | FS (Bart)    | 1.874 | 0.019 | 2248.001 |
| Var slope  | JSEM         | 0.099 | 0.017 |          |
|            | 2S-PA (Reg)  | 0.100 | 0.016 |          |
|            | 2S-PA (Bart) | 0.100 | 0.016 |          |
|            | FS (Reg)     | 0.065 | 0.004 |          |
|            | FS (Bart)    | 0.141 | 0.016 |          |
