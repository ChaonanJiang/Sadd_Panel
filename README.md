# Overview

This repository contains the working paper "Saddlepoint approximations for spatial panel data models", 2021, Chaonan Jiang, 
Davide La Vecchia, Elvezio Ronchetti, Olivier Scaillet. The repository contains R files which reproduce Figure 1-7 and Table 
1-3 in the manuscript.
# Data 

* Investment rate
* Saving rate
* inverse distance weight matrix
* 7 nearest neighbor weight matrix


# Instructions for reproducibility

1- [WeightMatrix.R](https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC/WeightMatrix.R) displays the geometry of three different spatial weight matrices shown in Figure 1. 

2- [SARAR.R](https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SARAR_MC/SARAR.R) generates a MC simulation based on SARAR(1,1) model to compare the distribution of MLE to the Gaussian asymptotic distribution via QQ-plots in Figure 2.

3- In the [SAR_MC folder](https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC), n24_WnRook.R, n24_WnQueen.R, n24_WnQueen_Torus.R, n100_WnRook.R, n100_WnQueen.R, n100_WnQueen_Torus.R make use of the algorithm in Section 5 of the paper to obtain the saddlepoint density approximation of the MLE for two sample sizes n = 24 or n = 100 and three weight matrices: Rook, Queen, Queen with torus. The codes generate Figure 3-6 as available in the paper. The size and weight matrix are indicated by the names of R files.

4- [OECD.R](https://github.com/ChaonanJiang/Sadd_Panel/blob/master/FH_Puzzle/OECD.R) shows London network for inverse distance and 7 nearest neighbours weight matrices, as shown in Figure 7.

5- In the [FH_Puzzle folder](https://github.com/ChaonanJiang/Sadd_Panel/blob/master/FH_Puzzle), FH_lambda.R obtains the MLEs for three parameters and 2 weight matrices in three sub-periods 1960-1970, 1971-1985 and 1986-2000 shown in Table 1. Following the similar Step 2-5 of the algorithm in Section 5, FH_lambda.R, FH_rho.R and FH_LR.R perform the saddlepoint and Wald tests for the composite hypotheses <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda=0" title="\lambda=0" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\rho=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho=0" title="\rho=0" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda=\rho=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda=\rho=0" title="\lambda=\rho=0" /></a> respectively. Then we can compute the p-values available in Table 2 under 2 weight matrices in three sub-periods. 
# Additional info.
* Author of the R codes: Chaonan Jiang.

* Creation date: June 2020.

* Last update: June 2021.

* R version: 3.6.1 (>=).
