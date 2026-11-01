# rprev

`rprev` estimates disease prevalence from registry data using the Monte Carlo
framework of [Crouch et al. (2014)](https://doi.org/10.1016/j.canep.2014.02.005).
The method combines observed registry cases with simulated incident populations
and bootstrapped survival models.

This repository continues the development of
[`stulacy/rprev-dev`](https://github.com/stulacy/rprev-dev).

## What is new in version 1.1.0?

The original implementation estimates prevalence for a **single index date**.
Version 1.1.0 extends the framework to support **multiple index dates** within
one coherent simulation pipeline.

Key features:

- The `index` argument accepts either a single date or a `c(t1, ..., tK)` vector of dates. Input dates are parsed, deduplicated, and sorted chronologically.
- Within each Monte Carlo replicate, the incident population is generated once and then evaluated at all requested index dates, so that the status of a simulated case evolves consistently over time: for example, a case that is no longer alive at an earlier index date cannot reappear as alive at a later one.
- For each simulated case `i`, one shared latent threshold is sampled per replicate, `U_i ~ Unif(0,1)`, and alive indicators are set as `A_{ik} = 1{U_i <= S_i(t_k)}` for each index date `t_k`. This preserves the marginal target `P(A_{ik}=1)=S_i(t_k)` while enforcing cross-time coherence.
- Backward compatibility is retained for the single-index case (`K = 1`).

The output contains prevalence estimates and uncertainty summaries for every
combination of index date and estimation horizon.

## Thesis context

This codebase is part of an **Applied Mathematics MSc thesis** completed at [Óbuda University (Óbudai Egyetem)](https://nik.uni-obuda.hu/en/home-english/).

**Thesis title:**
**Extending the Simulation-Based Prevalence Estimation of Crouch to Continuous Temporal Modeling**

The thesis provides:

- a methodological background on prevalence estimation from incidence and survival,
- a literature-grounded treatment of the Crouch et al. simulation framework,
- a formal and implementation-level extension to the multi-index setting.

The central practical motivation is computational efficiency: one simulated population and one fitted survival structure can be reused across several index dates, rather than rerunning full simulation workflows independently for each date.

## Supplementary notebooks

Supplementary thesis materials are available in the GitHub
[`notebooks/` directory](https://github.com/lrado1/rprev-ext/tree/master/notebooks).
They cover:

- generation and export of synthetic registry-style data,
- consistency and accuracy checks for the multi-index implementation,
- runtime benchmarking against the reference single-index workflow.

The notebooks are excluded from the CRAN source package and remain available
through GitHub. Running them requires Jupyter with an R kernel such as
[IRkernel](https://irkernel.github.io/installation/).

## Installation

Note: the repository name is `rprev-ext`, but the package name remains `rprev`.

Install the released version from CRAN:

```r
install.packages("rprev")
```

Install the development version from GitHub with `devtools`:

```r
# install.packages("devtools")
devtools::install_github("lrado1/rprev-ext")
```

If you cloned the repository locally, you can either load it directly for development or install it from the local checkout:

```r
# install.packages("devtools")
devtools::load_all(".")   # load directly from the local repository
# or
devtools::install(".")    # install the package from the local repository
```

## Multi-index example

```r
library(rprev)
library(survival)
data(prevsim)

set.seed(123)
results <- prevalence(
  index = c("2010-01-01", "2011-01-01", "2012-01-01"),
  num_years_to_estimate = c(5, 10),
  data = prevsim,
  inc_formula = entrydate ~ sex,
  surv_formula = Surv(time, status) ~ age + sex,
  dist = "weibull",
  N_boot = 100,
  population_size = 1e6,
  death_column = "eventdate"
)

print(results)
summary(results)
```
