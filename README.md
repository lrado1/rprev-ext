# rprev-ext (multi-index extension)

This repository is a **fork** of the development repository [`stulacy/rprev-dev`](https://github.com/stulacy/rprev-dev).  
The upstream project implements the Monte Carlo prevalence framework of [Crouch et al. (2014)](https://doi.org/10.1016/j.canep.2014.02.005), where prevalence is estimated from incidence and survival using simulated incident populations and bootstrapped survival models.

## What is new in this fork?

The original implementation estimates prevalence for a **single index date**.  
This fork extends the framework to support **multiple index dates** within one coherent simulation pipeline.

Key features:

- The `index` argument accepts either a single date or a `c(t1, ..., tK)` vector of dates.
- Within each Monte Carlo replicate, the incident population is generated once and then evaluated at all requested index dates, so that the status of a simulated case evolves consistently over time: for example, a case that is no longer alive at an earlier index date cannot reappear as alive at a later one.
- Backward compatibility is retained for the single-index case (`K = 1`).

The output consists of prevalence estimates and uncertainty summaries organised by index date and estimation horizon.

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

Additional thesis materials are available in `notebooks/`. These include a script for generating synthetic registry-style datasets and three **Jupyter notebooks** used in the thesis analyses.

- generation and export of synthetic registry-style datasets, used as controlled test data for the validation tasks below,
- consistency and accuracy checks for the multi-index implementation,
- runtime benchmarking against the reference single-index workflow.


The thesis analyses are provided as **Jupyter notebooks**, i.e. interactive documents that combine code, text, and output in a single file. They can be opened either in **VS Code** with the Jupyter extension installed or in **JupyterLab / Jupyter Notebook** in a web browser. To run the R code in them, Jupyter must have access to an R kernel such as **IRkernel**. Setup instructions are available in the official [Jupyter installation guide](https://jupyter.org/install) and the [IRkernel installation documentation](https://irkernel.github.io/installation/).

## Package Installation

Note: the repository name is `rprev-ext`, but the package name remains `rprev`.

Install from GitHub with `devtools`:

```r
# install.packages("devtools")
devtools::install_github("lrado1/rprev-ext", ref = "master")
```

If you cloned the repository locally, you can either load it directly for development or install it from the local checkout:

```r
# install.packages("devtools")
devtools::load_all(".")   # load directly from the local repository
# or
devtools::install(".")    # install the package from the local repository
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
