# Overview
This repository contains the working paper "Saddlepoint approximations for spatial panel data models", 2019, Chaonan Jiang, Davide La Vecchia, Elvezio Ronchetti, Olivier Scaillet. The repository contains R files which reproduce Figure 1-8 and Table 1-2 in the manuscript.
# Content
1- WeightMatrix.R (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC/WeightMatrix.R) displays the geometry of three different spatial weight matrices shown in Figure 1. 

2- SAR_mle.R (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC/SAR_mle.R) generates a MC simulation based on SAR(1) model to compare the distribution of MLE to the Gaussian asymptotic distribution via QQ-plots in Figure 2.

3- In the SAR_MC folder (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC), n24_WnRook.R, n24_WnQueen.R, n24_WnQueen_Torus.R, n100_WnRook.R, n100_WnQueen.R, n100_WnQueen_Torus.R make use of the algorithm in Section 5 of the paper to obtain the saddlepoint density approximation of the MLE for two sample sizes n=24 or 100 and three weight matrices: Rook, Queen, Queen with torus. The codes generate Figure 3-7 available in the paper. The size and weight matrix are indicated by the names of R files.

4- OECD.R (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/FH_Puzzle/OECD.R) shows London network for inverse distance and 7 nearest neighbours shown in Figure 8.

5- In the FH_Puzzle folder (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/FH_Puzzle), Lambda_60-70.R, Lambda_71-85.R, Lambda_86-20.R get the MLEs for three parameters and 2 weight matrices in three sub-periods: inverse distance and 7 nearest neighbours in Table 1. Following the similar Step 2-5 of the algorithm in Section 5, Lambda_60-70.R, Lambda_71-85.R, Lambda_86-20.R, Beta_60-70.R, Beta_71-85.R, Beta_86-20.R, Rho_60-70.R, Rho_71-85.R and Rho_86-20.R generate the saddlepoint density approximation of the MLE of the parameter for 2 weight matrices in three sub-periods. Then we can compute the p-values of parameters available in Table 2. The name of R files consist of parameter and sub-period. 
# Additional info.
*Author of the R codes: Chaonan Jiang.

*Creation data: June 2020.

*Last update: June 2020.

*R version: 3.6.1 (>=).
