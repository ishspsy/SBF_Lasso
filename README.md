# SBF_Lasso: Efficient Functional Lasso Kernel Smoothing for High-dimensional Additive Regression


## Overview

*SBF_Lasso* is a novel functional Lasso kernel smoothing method, which combines the idea of functional Lasso with the smooth backfitting technique.


## Main functions

- proj_kj: projection operator of the k th component function to the j th component space.

- NW: Nadaraya-Watson estimator

- SBF_update_j: first step of the algorithm based on solving system of equations.

- shrink_update_j: second step of the algoirthm based on thresholding operator.

- sparse_SBF_fit: main SBF algorithm.

- sim: simulation analysis.



## Main analysis

Please follow the links to reproduce the clustering results of TCGA data sets

- [Example_BRCA.m](https://github.com/ishspsy/MKerW-A/blob/master/example_BRCA.m)
: This file includes a brief analysis step included in the paper when patients of Breast Invasive Carcinoma (BRCA) cancer are used in our analysis.  
User may follow the similar analysis step for the other cancer types.

- [Example_BRCA.pdf](https://github.com/ishspsy/MKerW-A/blob/master/example_BRCA.pdf)
: This pdf file includes a brief analysis step as well as a survival curve using Breast Invasive Carcinoma (BRCA) cancer.

-  [Clustering analysis for each cancer type ("main_real.m")](https://github.com/ishspsy/MKerW-A/blob/master/main_real.m)
: Generate clustering results as well as basic survival analysis results for each of the 22 cancer types. All the analysis in the paper are based on
the results from this file.

-  [Clustering analysis for identifying 22 cancer types ("main_simul.m")](https://github.com/ishspsy/MKerW-A/blob/master/main_simul.m)
: Generate clustering results related to identifying 22 cancer types.

-  [Analysis on target cluster number ("main_simul_cls.m")](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_cls.m)
: Choose target cluster numbers.

-  [Sensitivity analysis ("main_simul_robust.m")](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_robust.m)
: Sensitivity analysis with respect to changes of regularization parameters.

-  [Sensitivity test w.r.t. additive noise ("main_simul_pert.m")](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_pert.m)
: Sensitivity analysis with respect to additive noise.


**Note** Most of the simulations and CCLE data applications were implemented on an iMac (3.6 GHz, 10 core, Intel Core i9, 64GB 2667 MHz DDR4) using R. 





## CCLE Data
CCLE dataset, which includes IC50 values and gene expressions for multiple cancer cell lines. This dataset originates from here. 




## Authors

* Eun Ryung Lee (erlee@skku.edu), Seyoung Park (ishspsy@yonsei.ac.kr), Enno Mammen (mammen@math.uni-heidelberg.de), and Byeong U. Park (bupark2000@gmail.com)
  
 Department of Statistics, Sungkyunkwan University, Department of Statistics, Yonsei University, Institute for Applied Mathematics, Heidelberg University, Department of Statistics, Seoul National University.


## License

This project is licensed under the MIT License.




