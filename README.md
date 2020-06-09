# Overview
This repository contains the working paper "Saddlepoint approximations for spatial panel data models", 2019, Chaonan Jiang, Davide La Vecchia, Elvezio Ronchetti, Olivier Scaillet. The repository contains R files which reproduce Figure 1-8 and Table 1-2 in the manuscript.
# Content
1- WeightMatrix.R (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC/WeightMatrix.R) displays the geometry of three different spatial weight matrices shown in Figure 1. 

2- SAR_mle.R (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC/SAR_mle.R) generates a MC simulation based on SAR(1) model to compare the distribution of MLE to the Gaussian asymptotic distribution via QQ-plots in Figure 2.

3- In the SAR_MC folder (https://github.com/ChaonanJiang/Sadd_Panel/blob/master/SAR_MC), n24_WnRook.R, n24_WnQueen.R, n24_WnQueen_Torus.R, n100_WnRook.R, n100_WnQueen.R, n100_WnQueen_Torus.R make use of the algorithm in Section 5 of the paper to obtain the saddlepoint density approximation of the MLE for two sample sizes n=24 or 100 and three weight matrices: Rook, Queen, Queen with torus. The size and weight matrix are indicated by the names of R files. The codes generate Figure 3-7 available in the paper.

4- OECD.R ()
