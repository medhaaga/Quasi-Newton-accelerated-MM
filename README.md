##  Quasi-Newton Acceleration of EM and MM Algorithms via Broyden’s Method

### Dependencies
The code requires R pacakages `quasiNewtonMM`, `turboEM`, `SQUAREM`, and `daarem` to run the MM acceleration methods considered in the paper Agarwal and Xu (2021).

### Reproducing Experimental Results

This repository contains the code for all examples in Agarwal and Xu (2021). The four examples are:
1. Landweber’s method for quadratic minimization 
2. Truncated beta binomial 
3. Generalized eigenvalues
4. Multivariate t-distribution
5. Sparse Logistic Regression

To run the optimization routine using all methods, run `code.R` files in individual examples folder. This may take some time. The R objects are stored `Out` folder for each example.

To output all tables and figures from the manuscript, run `OneRun.R`.
