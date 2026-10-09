# BCDint

Scripts for block-coordinate descent in linear causal models with interventional
data. This is a script collection, not an R package.

## Dependencies

Use R >= 4.1. Core fitting, simulation, data generation, tests, and `plot_f.R`
need only base R. `sachs.R` loads **tidyverse** (including purrr, dplyr and
ggplot2), **MASS**, and **SEMgraph**. SEMgraph supplies the Sachs dataset and has
Bioconductor dependencies. Install them in R with:

```r
install.packages(c("tidyverse", "MASS", "BiocManager"))
BiocManager::install("SEMgraph", update=FALSE, ask=FALSE)
```

See the [SEMgraph project](https://fernandopalluzzi.github.io/SEMgraph/) and
[BiocManager installation documentation](https://bioconductor.github.io/BiocManager/reference/install.html).
The optional `plotGraph()` helper in `generateModel.R` needs `ggm` and `igraph`:

```r
BiocManager::install(c("ggm", "igraph"), update=FALSE, ask=FALSE)
```

Symbolic computations require Mathematica or Wolfram Engine/Cloud. R scripts do
not install packages automatically.

## Run from the repository root

Use fresh sessions so that global variables from previous analyses are absent:

```sh
Rscript --vanilla tests/regression.R
Rscript --vanilla tests/simulation-smoke.R
Rscript --vanilla simulation.R > simulation-results.txt
Rscript --vanilla sachs.R > sachs-results.txt
Rscript --vanilla plot_f.R
```

`ricf_dg.R`, `ricf_int.R`, and `generateModel.R` define functions and can be
sourced independently. `simulation.R` and `sachs.R` source sibling files using
relative paths. `plot_f.R` creates `Rplots.pdf` in a noninteractive session;
subsequent plotting scripts may replace the default-device output.

`tests/simulation-smoke.R` runs the real simulation driver over all 24
configurations with two replicates each. It checks successful execution,
convergence, finite estimates/metrics, and all saved archives. Scripts and outputs
are isolated in a temporary directory, which is removed on success and retained
for debugging on failure. It does not overwrite the repository's `data/` files
or verify the paper's statistical results.

## Simulation

`simulation.R` runs all 24 configurations, with seed 1 and 1000 replicates per
configuration, and automatically creates `data/` for its `.rda` files. Sourcing
this script also launches the full study. `BCD_REPL` can override the replicate
count; the smoke test sets it only inside its isolated run. The saved `res` matrix has rows
AggBCD/BCD-i and columns convergence count, mean **MSE**, mean observational
population likelihood gap, and total elapsed seconds. Average elapsed milliseconds
per replicate are `1000 * res[,4] / repl`; this is not CPU time.

Each BCD-i replicate uses one default initialization through `ricf_int_()`;
AggBCD similarly fits each environment once. The 1000 replicates are different
graphs/data sets, not repeated starts on one data set. The simulator retains the
legacy `cov()` convention, coefficient magnitudes .3--1, and its original
intervention sample-count sampler. Different hardware/software need not reproduce timings.

## Data generation and fitting

`generateData()` obtains its dimension from B. Supply `target.length` for explicit
environment counts, or `n` for a common count. Without a sample count it stops
with an explanatory error rather than reading a global n:

```r
source("generateModel.R")
B <- matrix(0, 3, 3)
B[2,1] <- B[3,2] <- 1
set.seed(1)
observational <- generateData(B, diag(3), n=100)
interventional <- generateData(B, diag(3),
    targets=list(numeric(0), 2), target.length=c(100, 100))

source("ricf_int.R")
fit <- ricf_int(L=t(B), data=t(interventional$Y),
    targets=interventional$targets,
    target.length=interventional$target.length)
```

`ricf_int_()` performs one fit. `ricf_int()` compares `restarts` random starts,
a zero/identity start, and a supplied/default start using all environments. It
forwards solver controls. `fit$restart_scores` and `fit$selected_restart` expose
the comparison. The single-variable case estimates variance from rows where the
variable was not intervened upon.

## Sachs likelihood convention

`sachs.R` loads the data explicitly from SEMgraph, sets seed 1, log-transforms
and centers each selected environment separately, and requests `covariance="ml"`.
This uses `crossprod(Y)/n` for the zero-mean Gaussian likelihood. Default fitting
uses `covariance="unbiased"` to preserve the simulation results. For ML mode,
inputs must already be centered or have known zero mean.

The shared `llh_int()` score modifies the coefficient matrix in **both** precision
factors for each intervention. Its `covariance` argument must match the fit.
With `covariance="ml"`, it returns twice the Gaussian log-likelihood, with
additive constants omitted. Therefore ordinary log-likelihood differences are
half the score differences. The corrected main Sachs scores are approximately
-19025.09 (DAG) and -19003.54 (cyclic); these replace the original erroneous scores.
They favor the cyclic candidate under this working model, not proof of a
biological feedback loop. Smaller graphs at the end are supplementary comparisons.

## Symbolic ML degree and helper plot

Copy `ML-deg.txt` into a fresh Mathematica/Wolfram Cloud session, or run:

```sh
wolframscript -file ML-deg.txt
```

The covariance log-determinant sign has been corrected. The script uses random
sufficient statistics; the supplied Cloud examples produced elimination degrees
6 and 9. Inspecting `GG[[1]]` alone is not a complete generic ML-degree certificate.

The helper plot in `plot_f.R` illustrates the four cases in Lemma A.1/Figure 6.
