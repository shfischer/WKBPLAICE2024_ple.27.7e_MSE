ICES WKBPLAICE 2024 – MSE for ple.27.7e
================

## Introduction

This repository contains the data code for a stock-specific MSE for the
ICES category 3 data-limited stock

- plaice (*Pleuronectes platessa*) in Division 7.e (western English
  Channel)

as part of the ICES WKBPLAICE 2024 benchmark.

The MSE framework and code is based on this publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2023). Risk equivalence in data-limited and data-rich fisheries
> management: an example based on the ICES advice framework. Fish and
> Fisheries 24(2): 231-247. <https://doi.org/10.1111/faf.12722>

which included ple.27.7e as a case study. The original code for the
above publication is available from
[`shfischer/MSE_risk_comparison`](https://github.com/shfischer/MSE_risk_comparison).

The operating models (OMs) are created using the SAM
[`stockassessment`](https://github.com/fishfollower/SAM/) R package and
follow the principles developed during the ICES Workshop on North Sea
stocks management strategy evaluation
([WKNSMSE](https://doi.org/10.17895/ices.pub.5090)).

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and its
[`mse`](https://github.com/flr/mse) package.

For exact reproducibility, R packages versions are recorded with
[renv](https://rstudio.github.io/renv/) in a
[renv.lock](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/renv.lock)
file.

## Repository structure

- `funs_*`: Function libraries, defining the functions used in the other
  scripts

  - [`funs.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/funs.R):
    generic function library, including definition of data-limited
    management procedures (MPs)

  - [`funs_analysis.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/funs_analysis.R):
    for analysis of results

  - [`funs_GA.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/funs_GA.R):
    functions used in the optimisation with the genetic algorithm (GA)

  - [`funs_OM.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/funs_OM.R):
    functions for generating the operating models

  - [`funs_WKNSMSE.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/funs_WKNSMSE.R):
    functions required for the ICES MSY rule

- `OM_*`: Scripts for operating models (OMs, including alternative OMs)

  - [`OM_ple.27.7e.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/OM_ple.27.7e.R)
    for plaice

  - [`OM_MSY.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/OM_MSY.R):
    script for estimating MSY

  - [`OM_MSY.pbs`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/OM_MSY.pbs):
    job script, calling `OM_MSY.R`

- `MP_*`: Script for running and analysing the MSE

  - [`MP_analysis.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/MP_analysis.R):
    script for analysing MSE results (summarising, exporting,
    visualisation, …)

  - [`MP_run.R`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/MP_run.R):
    script for running any MP in the MSE and optimising MPs

  - `MP_*.pbs`: job submission scripts, used for running MP_run.R on a
    high-performance computing cluster,
    e.g. [`MP_run_rfb_mult.pbs`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/MP_run_rfb_mult.pbs)
    for optimising the rfb rule with a multiplier

[`input/`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/tree/master/input):
This directory contains all files required for generating the OMs for
the three stocks (`OM_*.R`)

- [`input/ple.27.7e/preparation/`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/tree/master/input/ple.27.7e/preparation):
  for plaice

- [`input/OM_refpts.csv`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/input/OM_refpts.csv):
  summarised OM reference points

[`output/`](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/tree/master/output):
This directory contains some summarised results

## R, R packages and version info

The MSE simulations were run with R:

``` r
> sessionInfo()
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)
```

The package versions and their dependencies are recorded with the R
package [renv](https://rstudio.github.io/renv/) and stored in the file
[renv.lock](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/renv.lock).
The exact package version can be restored by cloning this repository,
navigating into this folder in R (or setting up a project), installing
the renv package

``` r
install.packages("renv")
```

and calling

``` r
renv::restore()
```

See [renv](https://rstudio.github.io/renv/) and the package
documentation for details.

The framework is based on the Fisheries Library in R (FLR) framework and
uses the [FLR packages](https://flr-project.org/)
[`FLCore`,](https://github.com/flr/FLCore)
[`FLasher`](https://github.com/flr/FLasher),
[`FLBRP`](https://github.com/flr/FLBRP),
[`FLAssess`](https://github.com/flr/FLAssess),
[`FLXSA`](https://github.com/flr/FLXSA),
[`ggplotFL`](https://github.com/flr/ggplotFL),
[`mse`](https://github.com/flr/mse), and
[`FLfse`](https://github.com/shfischer/FLfse). See
[renv.lock](https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/master/renv.lock)
for version details and sources.

Also, the R package
[`stockassessment`](https://github.com/fishfollower/SAM)is used.

For running the optimisations on a high-performance computing cluster, a
suitable MPI back-end and the R package
[`Rmpi`](https://cran.r-project.org/web/packages/Rmpi/index.html) are
required.
