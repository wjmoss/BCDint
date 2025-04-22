# BCDint

A small collection of scripts for block-coordinate descent (BCD) algorithms in linear causal models.

## Algorithms
- **`ricf_dg`**: BCD algorithm (Drton et al., 2019) for directed graphs (without bidirected edges)  
- **`ricf_int`**: Early version of the BCD algorithm for interventional data  
- **`ricf_int2`**: Updated BCD algorithm for interventional data

## Simulations
- **`generateModel`**: Data generation process used in simulations  
- **`simulationcyc`**: Runs simulation and evaluation for one configuration $(v, n, l, d) = (p, m, k, d)$

## Maximum Likelihood Degree Computation
- **`ML_deg`**: Mathematica code for computing maximum likelihood degrees:
  - One for **observational data** only
  - One for **observational + interventional data**

## Notes
- Not a package; just standalone scripts  
- No installation required  
- Dependencies: base R and Mathematica only
