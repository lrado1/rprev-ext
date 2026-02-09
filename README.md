# rprev-ext (multi-index extension)

This repository is a **fork** of the development repository [`stulacy/rprev-dev`](https://github.com/stulacy/rprev-dev).  
The upstream project implements the Monte Carlo prevalence framework of [Crouch et al. (2014)](https://doi.org/10.1016/j.canep.2014.02.005), where prevalence is estimated from incidence and survival using simulated incident populations and bootstrapped survival models.

## What is new in this fork?

The original implementation estimates prevalence for a **single index date**.  
This fork extends the framework to support **multiple index dates** in one coherent simulation pipeline.

Key features:

- The `index` argument accepts either a single date or a vector `c(t1, ..., tK)`.
- Within each Monte Carlo / bootstrap replicate, the incident population is generated once and then evaluated at all index dates.
- **Methodological deviation from the original single-index sampling (analogue, not contradiction):** for each simulated case `i`, one shared latent threshold is sampled per replicate, `U_i ~ Unif(0,1)`, and alive indicators are set as `A_{ik} = 1{U_i <= S_i(t_k)}` for each index date `t_k`. This preserves the original marginal target `P(A_{ik}=1)=S_i(t_k)` while enforcing cross-time coherence (`A_{ik}` evolves monotonically with `S_i(t_k)`), avoiding trajectories that can be non-coherent under independent per-time Bernoulli draws.

The output is a timepoint-wise set of prevalence estimates with uncertainty summaries.

## Thesis context

This codebase is part of an **Applied Mathematics MSc thesis** completed at [Óbuda University (Óbudai Egyetem)](https://nik.uni-obuda.hu/en/home-english/).

**Thesis title:**  
**Extending the Simulation-Based Prevalence Estimation of Crouch to Continuous Temporal Modeling**

The thesis provides:

- a methodological background on prevalence estimation from incidence and survival,
- a literature-grounded treatment of the Crouch et al. simulation framework,
- a formal and implementation-level extension to the multi-index setting.

The central practical motivation is computational efficiency: one simulated population and one fitted survival structure can be reused across several index dates, rather than rerunning full simulation workflows independently for each date.

## Notebooks

Additional thesis materials (exploratory analyses, accuracy checks, runtime evaluations) are available in:

- `notebooks/`

## Installation

Install from GitHub with `devtools`:

```r
# install.packages("devtools")
devtools::install_github("lrado1/rprev-ext", ref = "master")
```

## Minimal Example

```r
library(rprev)
library(survival)
data(prevsim)

results <- prevalence(
  index = c("2010-01-01", "2011-01-01", "2012-01-01"),
  num_years_to_estimate = c(5, 10, 20),
  data = prevsim,
  inc_formula = entrydate ~ sex,
  surv_formula = Surv(time, status) ~ age + sex,
  dist = "weibull",
  population_size = 1e6,
  death_column = "eventdate"
)

print(results)
summary(results)
```
