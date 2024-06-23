# SBF_Lasso: Efficient Functional Lasso Kernel Smoothing for High-dimensional Additive Regression


## Overview

*SBF_Lasso* is a novel functional Lasso kernel smoothing method, which combines the idea of functional Lasso with the smooth backfitting technique.


## Main functions

- proj_kj: projection operator of the k th component function to the j th component space.

- NW: Nadaraya-Watson estimator.

- SBF_update_j: first step of the algorithm based on solving system of equations.

- shrink_update_j: second step of the algoirthm based on thresholding operator.

- sparse_SBF_fit: main SBF algorithm.

- sim: simulation analysis.


## Authors

* Eun Ryung Lee (erlee@skku.edu), Seyoung Park (ishspsy@yonsei.ac.kr), Enno Mammen (mammen@math.uni-heidelberg.de), and Byeong U. Park (bupark2000@gmail.com)
  
 Department of Statistics, Sungkyunkwan University, Department of Statistics, Yonsei University, Institute for Applied Mathematics, Heidelberg University, Department of Statistics, Seoul National University.


## License

This project is licensed under the MIT License.




