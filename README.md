# OneBNB
This file includes simulation and case study for 1-BNB method.

### R Code

- Main.R includes the whole process of simulation and case study including the comparison of five methods, i.e., 1-BNB, SnowBoot, CVE, MRME and CHTE. 
- For simulation study, we take power-law degree distribution with $\gamma = 3$ for demonstration. 
- For case study, the analysis of all three real datasets is included.
- The adjacency matrix completion processs of 1-BNB method is conducted in the Matlab Code file.


### Matlab Code
- Simu_OneBNB_B1000.m conducts the adjacency matrix completion process of 1-BNB method for simuation study. We take power-law degree distribution with $\gamma = 3$ as example.
- OneBNB.m conducts the adjacency matrix completion process of 1-BNB for case study. We take BK social network as example.
- The other .m file are the necessary functions needed for adjacency matrix completion. 

### Data 
Data.zip contains some data needed for the demonstration. Please unzip it before running any code.
