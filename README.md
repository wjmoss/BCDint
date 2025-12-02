# BCDint

A small collection of scripts for block-coordinate descent (BCD) algorithms in linear causal models for interventional data.

## Algorithms
- **`ricf_dg.R`**: BCD algorithm (Drton et al., 2019), simplied version for directed graphs (without bidirected edges)  
- **`ricf_int.R`**: BCD algorithm for interventional data

## Simulations
- **`generateModel.R`**: Data generation process used in simulations  
- **`simulation.R`**: Runs simulation and evaluation for one configuration $(v, n, l, d) = (p, m/p, k, d)$

## Maximum Likelihood Degree Computation
- **`ML-deg.txt`**: Mathematica code for computing maximum likelihood degrees:
  - One for **observational data** only
  - One for **observational + interventional data**

## Notes
- Not a package; just standalone scripts  
- No installation required  
- Core dependencies: **base R**
- Optional (for visualization function in `generateModel.R`): **ggm**, **BiocManager**, **graph**, **igraph**
- For symbolic computations, copy-paste the ML_deg.txt file into a **Wolfram Engine** or **Mathematica** session
