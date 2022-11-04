# I.   Initialize.
rm(list = ls())
setwd("/Users/dingyi/Desktop/OneBNB\ /Code")


# II.  Load packages.
source('R Code/Load_Packages.R')
Load_Package()


# III. Simulation
## 1. Generate random graphs with pre-given degree distribution.
### Use R script file Simulation_Generator.R to generate random graphs.
### Store the generated graphs data in 'Data' folder.

## 2. Perform 1-BNB method to analyse simulation data.
## In this part, we use MATLAB to perform the 1-BNB method to obtain bootstrap samples.
## The bootstrap samples data would be output to 'Data' folder.

## 2(1). Read the bootstrap samples generated via 1-BNB method on MATLAB.
## The mat file store all B bootstrap samples for all M generated networks.
## Explanation to the mat file name:
## 'Powlaw' represents power-law distribution;
## 'alpha2' represents parameter taken to be 2;
## '02p'    represents sampling proportion taken to be 20%.
## For geometric distribution with parameter 0.2 and sampling proportion 20%, 
## file name would be:' GeomSimuBootsAdjusted_prob02_02p_1000sparse.mat'.
## Here take power-law degree distribution with parameter 2 as demonstration.
OneBitPowBoots = readMat("Data/GeomSimuBootsAdjusted_prob04_30p_1000sparse.mat",
                          header = FALSE, sparsematrixClass = 'SparseM')
OneBitPowBoots = OneBitPowBoots$SimuBoots
# OneBitPowBoots1 = readMat("Data/PowlawSimuBootsAdjusted_alpha3_05p_m100_1000sparseI1.mat")
# OneBitPowBoots2 = readMat("Data/PowlawSimuBootsAdjusted_alpha3_05p_m100_1000sparseI2.mat")
# OneBitPowBoots1 = OneBitPowBoots1$S1
# OneBitPowBoots2 = OneBitPowBoots2$S2
# OneBitPowBoots  = rbind(OneBitPowBoots1,OneBitPowBoots2)
## 2(2). Obtain the confidence intervals for degree chosen in 'Degree_Sequence' for a single network.
## This is for the visulization of the 1-BNB method in simulation study.
## (a) Choose the first network.
Yhat = OneBitPowBoots[which(OneBitPowBoots[,ncol(OneBitPowBoots)] == 2), ]
## (b) Delete the last column for it only contains the index of networks.
Yhat = Yhat[, -ncol(Yhat)]
## (c) Choose a subset of the whole degree sequence for better display
## Here we choose degree from 0 to 10.
Degree_Sequence = 0:20
## Use 'OneBNB_Degree_CI' function to obtain the confidence intervals.
source('R Code/OneBNB_Degree_CI.R')
Degree_CI = ## Store confidence intervals in Degree_CI.
  OneBNB_Degree_CI(Degree_Sequence = Degree_Sequence,
                   Bootstrap_Samples = Yhat,
                   Sample_Number = 1000,
                   #Sample_Number = 500,
                   Alpha = 0.05)
## !Warning: this process might take a while.
## (d) Draw the confidence intervals along with the theoretical degree distribution.
source('R Code/Degree_CI_Plot.R')
## First create a dataframe for theoretical distribution
source('R Code/Theoretical_Degree_Distribution.R')
## For power-law distribution with parameter 2 or 3, type would be 'powlaw2' or 'powlaw3'.
## For geometric distribution with parameter 0.2 or 0.4, type would be 'geom02' or 'geom04'.
Degree_Distribution = Theoretical_Degree_Distribution(Degree_Sequence = Degree_Sequence,
                                                      type = 'geom', gprob = 0.4)
## Then plot the CI along with theoretical degree distribution.
Degree_CI_Plot_Powerlaw2 =
  Degree_CI_Plot(Degree_CI = Degree_CI, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_Powerlaw2  

## 2(3). 
source('R Code/OneBNB_Evaluation.R')
Eval_OneBNB = 
  OneBNB_Evaluation(OneBitPowBoots = OneBitPowBoots, LoopNum = 50,
                    Sample_Number = 1000, Alpha = 0.05,
                    Degree_Sequence = seq(0,20),
                    type = 'geom',gprob = 0.4)
Eval_OneBNB
## !Warning: this process might take a long while.
# capture.output(Eval_OneBNB,file = "Eval_OneBNB_Pow22_20p.txt")

## 3. For comparison, perform SnowBoot method to analyse simulation data.
## 3(1). Read the generated simulation data.
Simulation_Data = read.csv('Data/PowlawSimu_alpha3.csv')
Simulation_Data = Simulation_Data[,-1]
## 3(2). USE 'SnowBoot_Evaluation' function to analyse the simulation data.
source('R Code/SnowBoot_Evaluation.R')
Start = Sys.time()
Eval_SnowBoot = SnowBoot_Evaluation(Simulation_Data = Simulation_Data,
                                    Degree_Sequence = Degree_Sequence,
                                    LoopNum = 50, Sample_Number = 1000, Alpha = 0.05,
                                    type = 'powerlaw',palpha = 3)
Eval_SnowBoot
End   = Sys.time()
Eval_SnowBoot_pow3 = Eval_SnowBoot
End - Start

Eval_SnowBoot = list()
Eval_SnowBoot[[6]] = Eval_SnowBoot_pow3
capture.output(Eval_SnowBoot,file = "Eval_SnowBoot_powlaw_2.txt")
