# I.   Initialize.
rm(list = ls())
setwd("/Users/dingyi/Desktop/OneBNB\ /Code/github")


# II. Load functions and packages
source('R Code/Function.R')


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
## 'alpha3' represents parameter taken to be 3;
## '20p'    represents sampling proportion taken to be 20%.
## For geometric distribution with parameter 0.2 and sampling proportion 20%, 
## file name would be:' GeomSimuBootsAdjusted_prob02_20p_1000sparse.mat'.
## Here take power-law degree distribution with parameter 3 as demonstration.
OneBitPowBoots = readMat("Data/PowlawSimuBootsAdjusted_alpha3_20p_1000sparse.mat",
                         header = FALSE, sparsematrixClass = 'SparseM')
OneBitPowBoots = OneBitPowBoots$SimuBoots
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
Degree_CI = ## Store confidence intervals in Degree_CI.
  OneBNB_Degree_CI(Degree_Sequence = Degree_Sequence,
                   Bootstrap_Samples = Yhat,
                   Sample_Number = 1000,
                   Alpha = 0.05)
## !Warning: this process might take a while.
## (d) Draw the confidence intervals along with the theoretical degree distribution.
## First create a dataframe for theoretical distribution
## For power-law distribution, type would be 'powerlaw', use palpha to set parameter.
## For geometric distribution, type would be 'geom', use gprob to set parameter.
Degree_Distribution = Theoretical_Degree_Distribution(Degree_Sequence = Degree_Sequence,
                                                      type = 'powerlaw', palpha = 3)
## Then plot the CI along with theoretical degree distribution.
Degree_CI_Plot =
  Degree_CI_Plot(Degree_CI = Degree_CI, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot  

## 2(3). 
Eval_OneBNB = 
  OneBNB_Evaluation(OneBitPowBoots = OneBitPowBoots, LoopNum = 50,
                    Sample_Number = 1000, Alpha = 0.05,
                    Degree_Sequence = seq(0,20),
                    type = 'powerlaw',palpha = 3)
Eval_OneBNB
## !Warning: this process might take a long while.

## 3. For comparison, perform SnowBoot method to analyse simulation data.
## 3(1). Read the generated simulation data.
Simulation_Data = read.csv('Data/PowlawSimu_alpha3.csv')
Simulation_Data = Simulation_Data[,-1]
## 3(2). USE 'SnowBoot_Evaluation' function to analyse the simulation data.
Eval_SnowBoot = SnowBoot_Evaluation(Simulation_Data = Simulation_Data,
                                    Degree_Sequence = Degree_Sequence,
                                    LoopNum = 50, Sample_Number = 1000, Alpha = 0.05,
                                    type = 'powerlaw',palpha = 3)
Eval_SnowBoot

## 4. For comparison, perform CVE method to analyse simulation data.
CVE_Cover = NULL
CVE_Width = NULL
LoopNum   = 50
palpha    = 3
#gprob    = 0.4
m         = 250
B         = 1000
p         = 0.2 #Bootstrap sample proportion
delta     = 3
lambda    = 0
maxlen    = 11
numCores <- detectCores()
numCores
registerDoParallel(numCores)
for (n in 1:LoopNum) {
  print(paste0("n=", n))
  set.seed(n)
  degseq = rpldis(m, xmin = 1, alpha = palpha) - 1 
  #degseq = rgeom(m, prob = gprob)
  
  # Generate graph.
  set.seed(n)
  degseq[1] = degseq[1] + (sum(degseq) %% 2)
  g = sample_degseq(degseq)
  g = simplify(g)
  
  registerDoParallel(numCores)
  Degree_estimate_bootstrap <-
    foreach(b = 1:B, .combine=cbind) %dopar% {
      id_sample = Node_Sample(g    = g,
                              prob = p,
                              seed = b)
      id_snow = Snow_Sample(g,id_sample,nwave = 0)
      g_sample = subgraph(g,id_snow)
      Y_sample = as_adjacency_matrix(g_sample)
      Y_sample = as.matrix(Y_sample)
      Degree_sample = Count_Degree(Y_sample)
      Nstar = as.numeric(Degree_sample$Freq)
      Nstar = fill_fun(Nstar, maxlen)
      P = diag(p,length(Nstar))
      C = p*(1-p)*diag(Nstar) + delta*diag(1,length(Nstar)) #Ego-centric sampling
      D = cbind(diag(1,maxlen-2),matrix(0,maxlen-2,2)) + 
        cbind(matrix(0,maxlen-2,1),diag(-2,maxlen-2),matrix(0,maxlen-2,1)) +
        cbind(matrix(0,maxlen-2,2),diag(1,maxlen-2))
      estimate = CVX_Estimate(m,Nstar,P,C,D,lambda)
      estimate$Deg
    }
  
  degree_CI = NULL
  for (j in 1:maxlen) {
    f = Degree_estimate_bootstrap[j,]
    f_sort = sort(f)
    alpha = 0.05
    low = f_sort[round(B*alpha/2)]
    up = f_sort[round(B*(1-alpha/2))]
    
    ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
    temp = c(j, low, up)
    degree_CI = rbind(degree_CI, temp)
  }
  colnames(degree_CI) = c("degree", "low", "up")
  rownames(degree_CI) = 1:maxlen
  degree_CI = as.data.frame(degree_CI)
  Degree_Sequence = 0:(maxlen-1)
  temp = DegreeCI_Evaluation(degree_CI,Degree_Sequence,
                             type = 'powerlaw',palpha)
  # temp = DegreeCI_Evaluation(degree_CI,Degree_Sequence,
  #                            type = 'geom',gprob)
  CVE_Cover = rbind(CVE_Cover,temp$k)
  CVE_Width = rbind(CVE_Width,temp$w)
}
ave = function(x){
  mean(x, na.rm = TRUE)
}
CP = apply(CVE_Cover, 2, ave) ## CP: Coverage Probability.
AW = apply(CVE_Width ,2, ave)      ## AW: Average Width.
Eval_CVE = cbind(Degree_Sequence, CP,AW)
colnames(Eval_CVE) = c('Degree','CP','AW')
Eval_CVE

## 5. For comparison, perform MRME method to analyse simulation data.
CVE_Cover = NULL
CVE_Width = NULL
LoopNum   = 50
palpha    = 3
#gprob    = 0.4
m         = 250
B         = 1000
p         = 0.1 #Bootstrap sample proportion
maxlen    = 11
numCores <- detectCores()
numCores
registerDoParallel(numCores)
for (n in 1:LoopNum) {
  print(paste0("n=", n))
  set.seed(n)
  degseq = rpldis(m, xmin = 1, alpha = palpha) - 1 
  # Generate graph.
  set.seed(n)
  degseq[1] = degseq[1] + (sum(degseq) %% 2)
  g = sample_degseq(degseq)
  g = simplify(g)
  Degree_estimate_bootstrap <-
    foreach(b = 1:B, .combine=cbind,.packages = c("igraph","expm","poweRlaw","plyr","CVXR")) %dopar% {
      id_sample = Node_Sample(g    = g,
                              prob = p,
                              seed = b)
      id_snow = Snow_Sample(g,id_sample,nwave = 0)
      g_sample = subgraph(g,id_snow)
      Y_sample = as_adjacency_matrix(g_sample)
      Y_sample = as.matrix(Y_sample)
      Dstar    = apply(Y_sample, 2, sum)
      Dmat     = Y_sample%*%Y_sample
      diag(Dmat) = Dstar
      Dhat      = 1/p* Dstar %*%t(Dstar) %*% 
        solve(Dstar %*% t(Dstar) + Dmat + 1*diag(1,nrow(Y_sample))) %*%
        Dstar
      Dhat = floor(Dhat/(sum(floor((Dhat))/m)))
      degree_freq = mutate(data.frame(table(list(Dhat))),
                           proportion = Freq/length(Dhat),
                           logodds = log(proportion/(1-proportion)))
      estimate    = degree_freq$proportion
      estimate    = fill_fun(estimate,maxlen)
    }
  
  degree_CI = NULL
  for (j in 1:maxlen) {
    f = Degree_estimate_bootstrap[j,]
    f_sort = sort(f)
    alpha = 0.05
    low = f_sort[round(B*alpha/2)]
    up = f_sort[round(B*(1-alpha/2))]
    
    ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
    temp = c(j, low, up)
    degree_CI = rbind(degree_CI, temp)
  }
  colnames(degree_CI) = c("degree", "low", "up")
  rownames(degree_CI) = 1:maxlen
  degree_CI = as.data.frame(degree_CI)
  Degree_Sequence = 0:(maxlen-1)
  temp = DegreeCI_Evaluation(degree_CI,Degree_Sequence,
                             type = 'powerlaw',palpha)
  CVE_Cover = rbind(CVE_Cover,temp$k)
  CVE_Width = rbind(CVE_Width,temp$w)
}
ave = function(x){
  mean(x, na.rm = TRUE)
}
CP = apply(CVE_Cover, 2, ave) ## CP: Coverage Probability.
AW = apply(CVE_Width ,2, ave)      ## AW: Average Width.
Eval_MRME = cbind(Degree_Sequence, CP,AW)
colnames(Eval_MRME) = c('Degree','CP','AW')
Eval_MRME

## 6. For comparison, perform CHTE method to analyse simulation data.
CVE_Cover = NULL
CVE_Width = NULL
LoopNum   = 50
palpha    = 3
#gprob    = 0.4
m         = 250
B         = 1000
p         = 0.5
for (n in 1:LoopNum) {
  print(paste0("n=", n))
  set.seed(n)
  degseq = rpldis(m, xmin = 1, alpha = palpha) - 1
  # Generate graph.
  set.seed(n)
  degseq[1] = degseq[1] + (sum(degseq) %% 2)
  g = sample_degseq(degseq)
  g = simplify(g)
  EdgeNum = length(E(g))
  
  Degree_estimate_bootstrap <-
    foreach(b = 1:B, .combine=cbind,.packages = c("igraph","expm","poweRlaw","plyr","CVXR","MASS","sets")) %dopar% {
      print(paste0("n=", n, ",b=", b))
      Y_sample = Pair_Sample(g,p,b)
      EdgeSet = list()
      for (i in 1:nrow(Y_sample)) {
        EdgeSet[[i]] = which(Y_sample[i,] == 1)
      }
      maxDeg = max(sapply(EdgeSet,length))
      EdgeNum_sample = sum(Y_sample[upper.tri(Y_sample)])
      gammaCount = function(vec,Count,vecNum){
        maxCount = length(vec)
        output = choose(maxCount,Count)
        return(output)
      }
      L = list()
      for (c in 1:maxlen) {
        print(c)
        L[[c]] = sapply(EdgeSet,gammaCount,Count = c,vecNum = nrow(Y_sample))
      }
      Lvec = sapply(L,sum)
      Pi   = rep(0,maxlen)
      for (c in 1:maxlen) {
        Pi[c] = choose(EdgeNum,EdgeNum_sample)^(-1) * choose(EdgeNum-c,EdgeNum_sample-c)
      }
      Lhat = Pi * Lvec
      TmatAll = matrix(-1,max(maxDeg,maxlen),max(maxDeg,maxlen))
      for (i in 1:max(maxDeg,maxlen)) {
        for (j in 1:max(maxDeg,maxlen)) {
          TmatAll[i,j] = choose(j,i)
        }
      }
      if (rcond(TmatAll)>=1e-16){
        Tmat = solve(TmatAll)[1:maxlen,1:maxlen]
        Dhat = solve(Tmat) %*% Lhat
      }else{
        Tmat = ginv(TmatAll)[1:maxlen,1:maxlen]
        Dhat = ginv(Tmat) %*% Lhat
      }
      fhat = Dhat/sum(Dhat)
    }
  
  degree_CI = NULL
  for (j in 1:maxlen) {
    f = Degree_estimate_bootstrap[j,]
    f_sort = sort(f)
    alpha = 0.05
    low = f_sort[round(B*alpha/2)]
    up = f_sort[round(B*(1-alpha/2))]
    
    ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
    temp = c(j, low, up)
    degree_CI = rbind(degree_CI, temp)
  }
  colnames(degree_CI) = c("degree", "low", "up")
  rownames(degree_CI) = 1:maxlen
  degree_CI = as.data.frame(degree_CI)
  Degree_Sequence = 0:(maxlen-1)
  temp = DegreeCI_Evaluation(degree_CI,Degree_Sequence,
                             type = 'powerlaw',palpha = palpha)
  CVE_Cover = rbind(CVE_Cover,temp$k)
  CVE_Width = rbind(CVE_Width,temp$w)
}
ave = function(x){
  mean(x, na.rm = TRUE)
}
CP = apply(CVE_Cover, 2, ave) ## CP: Coverage Probability.
AW = apply(CVE_Width ,2, ave)      ## AW: Average Width.
Eval_CHTE = cbind(Degree_Sequence, CP,AW)
colnames(Eval_CHTE) = c('Degree','CP','AW')
Eval_CHTE

# IV.Case Study
## Next are the case studies of three real datasets.

############################################################################
# 1. Collaboration Network
source('R Code/Function.R')
## Read data from edgelist and transform it to graph type data and adjacency matrix.
author_edge = read.csv('Data/StatisticianCollabrationNetwork.csv')
author_edge['larger_author'] = apply(author_edge, 1, function(x){max(x)})
author_edge['smaller_author'] = apply(author_edge, 1, function(x){min(x)})
author_edge = author_edge[,-c(1,2)]
author_edge2 = unique(author_edge) 
g = graph_from_data_frame(author_edge2, directed=FALSE)
W = as_adjacency_matrix(g)
## Obtain the network density.
Network_density = graph.density(g)
## Estimate the scaling parameter.
deg = degree(g)
fit_pl = displ$new(deg)
est = estimate_pars(fit_pl)
Scaling_parameter = est$pars 
## Plot the degree distribution.
deg_dist = degree_distribution(g)
deg_dist = deg_dist[-1]
n = seq(1,length(deg_dist))
Deg_dist = cbind(n,deg_dist)
Deg_dist = as.data.frame(Deg_dist)
Degree_Plot = 
  ggplot(Deg_dist, aes(x = n, y = deg_dist)) +
  geom_line() +
  geom_point(size = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,250,50)) +
  xlab("degree") + ylab("proportion of nodes") +
  xlim(c(0,100)) +
  theme(text = element_text(family="Times New Roman"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
Degree_Plot

## K-Core
## Obtain the 30-core.
W.old = W
W = W.old
converg = FALSE
old.nrow = nrow(W)
while(!converg){
  d = colSums(W)
  to.keep = which(d>=30)
  if(old.nrow==length(to.keep)){
    converg = TRUE
  }
  old.nrow = length(to.keep)
  W = W[to.keep,to.keep]
}
Core = as.matrix(W)
## Choose again a subset of the whole degree sequence of 30-core for better display
## Here we choose degree from 30 to 65.
Degree_Sequence = 31:65
Core_Degree_Distribution = 
  Degree_Distribution_Extract(AdjacencyMatrix =  Core, 
                              Degree_Sequence = Degree_Sequence)
colnames(Core_Degree_Distribution) = c('degree','Freq','proportion','logodds')
## For better visualization, interpolate f(k) = 0 for k = 59,60,61,62.
Core_Degree_Distribution_Adjusted =
  rbind(Core_Degree_Distribution,
        c(59,0,0,0),c(60,0,0,0),c(61,0,0,0),c(62,0,0,0))
Degree_Plot_Core = 
  ggplot(Core_Degree_Distribution_Adjusted, aes(x = degree, y = proportion)) +
  geom_line(color = 'red') +
  geom_point(size = 1,color = 'red') +
  theme_bw() +
  scale_x_continuous(breaks = seq(31,65,1)) +
  xlab("degree") + ylab("proportion of nodes") +
  xlim(31,65) +
  theme(text = element_text(family="Times New Roman"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
Degree_Plot_Core

## 1BNB
Core_Yhat = readMat("Data/CollaBoots_30core_30p_B200.mat")
Core_Yhat = Core_Yhat$BootsResult
## !Warning: this process might take a while.
## Obtain the 1-BNB-estimated CI for degree distribution of 30-core
Sample_Degree_CI_OneBNB =
  OneBNB_Degree_CI(Degree_Sequence = Degree_Sequence,
                   Bootstrap_Samples = Core_Yhat,
                   Sample_Number = 200,
                   Alpha = 0.05)
## Then use 'Degree_CI_Plot' function to plot.
Degree_CI_Plot_1BNB =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_OneBNB,
                 Degree_Distribution = Core_Degree_Distribution_Adjusted,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_1BNB

## CVE
numCores <- detectCores()
numCores
registerDoParallel(numCores)
p = 0.1
B = 200
maxlen = 35
lambda = 0.001
delta  = 1
g_core = graph_from_adjacency_matrix(Core,mode = 'undirected')
m = length(g_core)
g = g_core
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    g_sample = subgraph(g,id_sample)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Degree_sample = Count_Degree(Y_sample)
    Nstar = as.numeric(Degree_sample$Freq)
    Nstar = fill_fun(Nstar, maxlen)
    P = diag(p,length(Nstar))
    C = p*(1-p)*diag(Nstar) + delta*diag(1,length(Nstar)) #Ego-centric sampling
    D = cbind(diag(1,maxlen-2),matrix(0,maxlen-2,2)) + 
      cbind(matrix(0,maxlen-2,1),diag(-2,maxlen-2),matrix(0,maxlen-2,1)) +
      cbind(matrix(0,maxlen-2,2),diag(1,maxlen-2))
    estimate = CVX_Estimate(m,Nstar,P,C,D,lambda)
    estimate$Deg
  }
Sample_Degree_CI_CVE = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CVE = rbind(Sample_Degree_CI_CVE, temp)
}
colnames(Sample_Degree_CI_CVE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CVE) = Degree_Sequence
Sample_Degree_CI_CVE = as.data.frame(Sample_Degree_CI_CVE)
Sample_Degree_CI_CVE$degree = Degree_Sequence
Sample_Degree_CI_CVE$up[which(Sample_Degree_CI_CVE$up<1e-04)] = 0
Sample_Degree_CI_CVE$low[which(Sample_Degree_CI_CVE$low<1e-04)] = 0
Degree_CI_Plot_CVE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CVE, 
                 Degree_Distribution = Core_Degree_Distribution_Adjusted,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_CVE

## MRME
p = 0.1
B = 200
maxlen = 35
g_core = graph_from_adjacency_matrix(Core,mode = 'undirected')
m = length(g_core)
g = g_core
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    g_sample = subgraph(g,id_sample)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Dstar    = apply(Y_sample, 2, sum)
    Dmat     = Y_sample%*%Y_sample
    diag(Dmat) = Dstar
    Dhat      = 1/p* Dstar %*%t(Dstar) %*% 
      solve(Dstar %*% t(Dstar) + Dmat + 1*diag(1,nrow(Y_sample))) %*%
      Dstar
    Dhat = floor(Dhat/(sum(floor((Dhat))/m)))
    degree_freq = mutate(data.frame(table(list(Dhat))),
                         proportion = Freq/length(Dhat),
                         logodds = log(proportion/(1-proportion)))
    estimate    = degree_freq$proportion
    estimate    = fill_fun(estimate,maxlen)
  }
Sample_Degree_CI_MRME = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_MRME = rbind(Sample_Degree_CI_MRME, temp)
}
colnames(Sample_Degree_CI_MRME) = c("degree", "low", "up")
rownames(Sample_Degree_CI_MRME) = 1:maxlen
Sample_Degree_CI_MRME = as.data.frame(Sample_Degree_CI_MRME)
Sample_Degree_CI_MRME$degree = Degree_Sequence
Degree_CI_Plot_MRME =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_MRME, 
                 Degree_Distribution = Core_Degree_Distribution_Adjusted,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_MRME

## CHTE
EdgeNum = length(E(g))
p       = 0.166
B       = 200
maxlen  = 36
g_core = graph_from_adjacency_matrix(Core,mode = 'undirected')
m = length(g_core)
g = g_core
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    g_sample = subgraph(g,id_sample)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    EdgeSet = list()
    for (i in 1:nrow(Y_sample)) {
      EdgeSet[[i]] = which(Y_sample[i,] == 1)
    }
    maxDeg = max(sapply(EdgeSet,length))
    EdgeNum_sample = sum(Y_sample[upper.tri(Y_sample)])
    gammaCount = function(vec,Count,vecNum){
      maxCount = length(vec)
      output = choose(maxCount,Count)
      return(output)
    }
    L = list()
    maxlen = 20
    for (c in 1:maxlen) {
      print(c)
      L[[c]] = sapply(EdgeSet,gammaCount,Count = c,vecNum = nrow(Y_sample))
    }
    Lvec = sapply(L,sum)
    Pi   = rep(0,maxlen)
    n    = length(g_sample)
    for (c in 1:maxlen) {
      Pi[c] = choose(m,n)^(-1) * choose(m-c-1,n-c-1)
    }
    Lhat = Pi * Lvec
    TmatAll = matrix(-1,max(maxDeg,maxlen),max(maxDeg,maxlen))
    for (i in 1:max(maxDeg,maxlen)) {
      for (j in 1:max(maxDeg,maxlen)) {
        TmatAll[i,j] = choose(j,i)
      }
    }
    if (rcond(TmatAll)>=1e-16){
      Tmat = solve(TmatAll)[1:maxlen,1:maxlen]
      Dhat = solve(Tmat) %*% Lhat
    }else{
      Tmat = ginv(TmatAll)[1:maxlen,1:maxlen]
      Dhat = ginv(Tmat) %*% Lhat
    }
    fhat = Dhat/sum(Dhat)
  }

Sample_Degree_CI_CHTE = NULL
for (j in 1:maxlen) {
  if (j <= nrow(Degree_estimate_bootstrap)){
    f = Degree_estimate_bootstrap[j,]
  }else{
    f = Degree_estimate_bootstrap[nrow(Degree_estimate_bootstrap),]
  }
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CHTE = rbind(Sample_Degree_CI_CHTE, temp)
}
colnames(Sample_Degree_CI_CHTE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CHTE) = 1:maxlen
Sample_Degree_CI_CHTE = as.data.frame(Sample_Degree_CI_CHTE)
Sample_Degree_CI_CHTE  = abs(Sample_Degree_CI_CHTE[-1,])
Sample_Degree_CI_CHTE$degree = Degree_Sequence
Degree_CI_Plot_CHTE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CHTE, 
                 Degree_Distribution = Core_Degree_Distribution_Adjusted,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_CHTE

## SnowBoot
Sample_Degree_CI_SnowBoot =
  SnowBoot_Degree_CI(Degree_Sequence = Degree_Sequence,
                     Adjacency_Matrix = Core,
                     Sample_Number = 200, 
                     Alpha = 0.05)
Sample_Degree_CI_SnowBoot = rbind(Sample_Degree_CI_SnowBoot,
                                  Sample_Degree_CI_SnowBoot[nrow(Sample_Degree_CI_SnowBoot),])
Degree_CI_Plot_SnowBoot =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_SnowBoot, 
                 Degree_Distribution = Core_Degree_Distribution_Adjusted,
                 Degree_Sequence = Degree_Sequence)
Degree_CI_Plot_SnowBoot

## Plot adjust
temp = rbind(Sample_Degree_CI_OneBNB,
             Sample_Degree_CI_CVE,
             Sample_Degree_CI_MRME,
             Sample_Degree_CI_CHTE,
             Sample_Degree_CI_SnowBoot)
max(temp$up)

Degree_Plot_Core_adj        = Degree_Plot_Core + ylim(c(0,0.35))
Degree_CI_Plot_1BNB_adj     = Degree_CI_Plot_1BNB + ylim(c(0,0.35))
Degree_CI_Plot_CVE_adj      = Degree_CI_Plot_CVE + ylim(c(0,0.35))
Degree_CI_Plot_MRME_adj     = Degree_CI_Plot_MRME + ylim(c(0,0.35))
Degree_CI_Plot_CHTE_adj     = Degree_CI_Plot_CHTE + ylim(c(0,0.35))
Degree_CI_Plot_SnowBoot_adj = Degree_CI_Plot_SnowBoot + ylim(c(0,0.35))

Collab_Analysis_Plot = list()
Collab_Analysis_Plot[[1]] = Degree_Plot_Core_adj
Collab_Analysis_Plot[[2]] = Degree_CI_Plot_1BNB_adj
Collab_Analysis_Plot[[3]] = Degree_CI_Plot_CVE_adj
Collab_Analysis_Plot[[4]] = Degree_CI_Plot_MRME_adj
Collab_Analysis_Plot[[5]] = Degree_CI_Plot_CHTE_adj
Collab_Analysis_Plot[[6]] = Degree_CI_Plot_SnowBoot_adj
names(Collab_Analysis_Plot) = c("Core","1BNB","CVE","MRME","CHTE","SnowBoot")
#saveRDS(Collab_Analysis_Plot,file = "Data/Collab_Analysis_Plot.rds")
Degree_Distribution = Core_Degree_Distribution[,c(1,3)]
Collab_Eval = Case_Study_Evaluation(Degree_Distribution,
                                    Sample_Degree_CI_OneBNB,
                                    Sample_Degree_CI_CVE,
                                    Sample_Degree_CI_MRME,
                                    Sample_Degree_CI_CHTE,
                                    Sample_Degree_CI_SnowBoot)
Collab_Eval
#capture.output(Collab_Eval, file = "Data/Collab_Eval.txt")

############################################################################
# 2. BK Social Network
## Read
BKEdgeList = read.table("Data/BK_EdgeList.txt",sep = "")
max(BKEdgeList$V1)
g = graph_from_data_frame(BKEdgeList,directed = FALSE)
g = simplify(g)
length(E(g))
graph.density(g)
deg = degree(g)
fit_pl = displ$new(deg)
est = estimate_pars(fit_pl)
Scaling_parameter = est$pars
deg_dist = degree_distribution(g)
deg_dist = deg_dist[-1]
n = seq(1,length(deg_dist))
Deg_dist = cbind(n,deg_dist)
Deg_dist = as.data.frame(Deg_dist)
Degree_Plot = 
  ggplot(Deg_dist, aes(x = n, y = deg_dist)) +
  geom_line(color = 'red') +
  geom_point(size = 1,color = 'red') +
  theme_bw() +
  scale_x_continuous(breaks = seq(1,250,50)) +
  xlab("degree") + ylab("proportion of nodes") +
  xlim(c(1,100)) +
  theme(text = element_text(family="Times New Roman"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
Degree_Plot
Degree_Plot + xlim(c(0,30))

## Sample 
vid_sample = Node_Sample(g,300/length(V(g)),1)
vid_snow   = Snow_Sample(g,vid_sample,0)
g_sample   = subgraph(g,vid_snow)
g_sample   = simplify(g_sample)
length(E(g_sample))
length(V(g_sample))
graph.density(g_sample)
deg_dist_sample = degree_distribution(g_sample)
deg_dist_sample = deg_dist_sample[-1]
n = seq(0,length(deg_dist_sample)-1)
Deg_dist_sample = cbind(n,deg_dist_sample)
Deg_dist_sample = as.data.frame(Deg_dist_sample)
W_sample = as_adjacency_matrix(g_sample)
W_sample = as.matrix(W_sample)
#write.csv(W_sample,file = "Data/BK_sample.csv")
Degree_Sequence = 0:30
Degree_Distribution = Deg_dist[1:length(Degree_Sequence),]
Degree_Distribution$degree     = 1:length(Degree_Sequence)
colnames(Degree_Distribution) = c("degree","proportion")

## 1-BNB
Sample_Yhat = readMat("Data/BK_sample3_15p_B200.mat")
Sample_Yhat = Sample_Yhat$BootsResult
Sample_Degree_CI_OneBNB =
  OneBNB_Degree_CI(Degree_Sequence = Degree_Sequence,
                   Bootstrap_Samples = Sample_Yhat,
                   Sample_Number = 200,
                   Alpha = 0.05)
Sample_Degree_CI_OneBNB$degree = 1:length(Degree_Sequence)
Degree_CI_Plot_1BNB =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_OneBNB, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_1BNB

## CVE
numCores <- detectCores()
numCores
registerDoParallel(numCores)
p = 0.1
B = 200
maxlen = 31
lambda = 0
delta  = 5
m = length(g_sample)
g = g_sample 
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    id_snow = Snow_Sample(g,id_sample,nwave = 0)
    g_sample = subgraph(g,id_snow)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Degree_sample = Count_Degree(Y_sample)
    Nstar = as.numeric(Degree_sample$Freq)
    Nstar = fill_fun(Nstar, maxlen)
    P = diag(p,length(Nstar))
    C = p*(1-p)*diag(Nstar) + delta*diag(1,length(Nstar)) #Ego-centric sampling
    D = cbind(diag(1,maxlen-2),matrix(0,maxlen-2,2)) + 
      cbind(matrix(0,maxlen-2,1),diag(-2,maxlen-2),matrix(0,maxlen-2,1)) +
      cbind(matrix(0,maxlen-2,2),diag(1,maxlen-2))
    estimate = CVX_Estimate(m,Nstar,P,C,D,lambda)
    estimate$Deg
  }
Sample_Degree_CI_CVE = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CVE = rbind(Sample_Degree_CI_CVE, temp)
}
colnames(Sample_Degree_CI_CVE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CVE) = 1:maxlen
Sample_Degree_CI_CVE = as.data.frame(Sample_Degree_CI_CVE)
Sample_Degree_CI_CVE$degree = 1:length(Degree_Sequence)
Sample_Degree_CI_CVE$up[which(Sample_Degree_CI_CVE$up<1e-04)] = 0
Sample_Degree_CI_CVE$low[which(Sample_Degree_CI_CVE$low<1e-04)] = 0
Degree_CI_Plot_CVE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CVE, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_CVE

## MRME
p = 0.08
B = 200
maxlen = 31
m = length(g_sample)
g = g_sample 
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    print(paste0("n=", n, ",b=", b))
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    id_snow = Snow_Sample(g,id_sample,nwave = 0)
    g_sample = subgraph(g,id_snow)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Dstar    = apply(Y_sample, 2, sum)
    Dmat     = Y_sample%*%Y_sample
    diag(Dmat) = Dstar
    Dhat      = 1/p* Dstar %*%t(Dstar) %*% 
      solve(Dstar %*% t(Dstar) + Dmat + 1*diag(1,nrow(Y_sample))) %*%
      Dstar
    Dhat = floor(Dhat/(sum(floor((Dhat))/m)))
    degree_freq = mutate(data.frame(table(list(Dhat))),
                         proportion = Freq/length(Dhat),
                         logodds = log(proportion/(1-proportion)))
    estimate    = degree_freq$proportion
    estimate    = fill_fun(estimate,maxlen)
  }
Sample_Degree_CI_MRME = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_MRME = rbind(Sample_Degree_CI_MRME, temp)
}
colnames(Sample_Degree_CI_MRME) = c("degree", "low", "up")
rownames(Sample_Degree_CI_MRME) = 1:maxlen
Sample_Degree_CI_MRME = as.data.frame(Sample_Degree_CI_MRME)
Sample_Degree_CI_MRME$degree = 1:length(Degree_Sequence)
Degree_CI_Plot_MRME =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_MRME, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_MRME

## CHTE
maxlen  = 31
EdgeNum = length(E(g))
maxlen  = 31
p       = 0.15
B       = 200
m = length(g_sample)
g = g_sample 
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    g_sample = subgraph(g,id_sample)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    EdgeSet = list()
    for (i in 1:nrow(Y_sample)) {
      EdgeSet[[i]] = which(Y_sample[i,] == 1)
    }
    maxDeg = max(sapply(EdgeSet,length))
    EdgeNum_sample = sum(Y_sample[upper.tri(Y_sample)])
    gammaCount = function(vec,Count,vecNum){
      maxCount = length(vec)
      output = choose(maxCount,Count)
      return(output)
    }
    L = list()
    maxlen = 30
    for (c in 1:maxlen) {
      print(c)
      L[[c]] = sapply(EdgeSet,gammaCount,Count = c,vecNum = nrow(Y_sample))
    }
    Lvec = sapply(L,sum)
    Pi   = rep(0,maxlen)
    n    = length(g_sample)
    for (c in 1:maxlen) {
      Pi[c] = choose(m,m-c-1)^(-1) * choose(n,n-c-1)
    }
    Lhat = Pi * Lvec
    TmatAll = matrix(-1,max(maxDeg,maxlen),max(maxDeg,maxlen))
    for (i in 1:max(maxDeg,maxlen)) {
      for (j in 1:max(maxDeg,maxlen)) {
        TmatAll[i,j] = choose(j,i)
      }
    }
    if (rcond(TmatAll)>=1e-16){
      Tmat = solve(TmatAll)[1:maxlen,1:maxlen]
      Dhat = solve(Tmat) %*% Lhat
    }else{
      Tmat = ginv(TmatAll)[1:maxlen,1:maxlen]
      Dhat = ginv(Tmat) %*% Lhat
    }
    fhat = Dhat/sum(Dhat)
  }

Sample_Degree_CI_CHTE = NULL
for (j in 1:maxlen) {
  if (j <= nrow(Degree_estimate_bootstrap)){
    f = Degree_estimate_bootstrap[j,]
  }else{
    f = Degree_estimate_bootstrap[nrow(Degree_estimate_bootstrap),]
  }
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CHTE = rbind(Sample_Degree_CI_CHTE, temp)
}
colnames(Sample_Degree_CI_CHTE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CHTE) = 1:maxlen
Sample_Degree_CI_CHTE = as.data.frame(Sample_Degree_CI_CHTE)
colnames(Sample_Degree_CI_CHTE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CHTE) = 1:maxlen
Sample_Degree_CI_CHTE$up[which(Sample_Degree_CI_CHTE$up<1e-04)] = 0
Sample_Degree_CI_CHTE$low[which(Sample_Degree_CI_CHTE$low<1e-04)] = 0
Sample_Degree_CI_CHTE$degree = 1:length(Degree_Sequence)
Degree_CI_Plot_CHTE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CHTE, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_CHTE

## SnowBoot
Sample_Degree_CI_SnowBoot =
  SnowBoot_Degree_CI(Degree_Sequence = Degree_Sequence,
                     Adjacency_Matrix = W_sample,
                     Sample_Number = 200, n.wave = 1,
                     Alpha = 0.05)
Sample_Degree_CI_SnowBoot = rbind(Sample_Degree_CI_SnowBoot,
                                  Sample_Degree_CI_SnowBoot[nrow(Sample_Degree_CI_SnowBoot),])
Sample_Degree_CI_SnowBoot = Sample_Degree_CI_SnowBoot[-1,]
Sample_Degree_CI_SnowBoot$degree = 1:length(Degree_Sequence)
Degree_CI_Plot_SnowBoot =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_SnowBoot, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_SnowBoot

## Plot
temp = rbind(Sample_Degree_CI_OneBNB,
             Sample_Degree_CI_CVE,
             Sample_Degree_CI_MRME,
             Sample_Degree_CI_CHTE,
             Sample_Degree_CI_SnowBoot)
max(temp$up)

Degree_Plot_adj             = Degree_Plot + ylim(c(0,0.5))
Degree_CI_Plot_1BNB_adj     = Degree_CI_Plot_1BNB + ylim(c(0,0.5))
Degree_CI_Plot_CVE_adj      = Degree_CI_Plot_CVE + ylim(c(0,0.5))
Degree_CI_Plot_MRME_adj     = Degree_CI_Plot_MRME + ylim(c(0,0.5))
Degree_CI_Plot_CHTE_adj     = Degree_CI_Plot_CHTE + ylim(c(0,0.5))
Degree_CI_Plot_SnowBoot_adj = Degree_CI_Plot_SnowBoot + ylim(c(0,0.5))
BK_Analysis_Plot = list()
BK_Analysis_Plot[[1]] = Degree_Plot_adj
BK_Analysis_Plot[[2]] = Degree_CI_Plot_1BNB_adj
BK_Analysis_Plot[[3]] = Degree_CI_Plot_CVE_adj
BK_Analysis_Plot[[4]] = Degree_CI_Plot_MRME_adj
BK_Analysis_Plot[[5]] = Degree_CI_Plot_CHTE_adj
BK_Analysis_Plot[[6]] = Degree_CI_Plot_SnowBoot_adj
names(BK_Analysis_Plot) = c("Original","1BNB","CVE","MRME","CHTE","SnowBoot")
#saveRDS(BK_Analysis_Plot,file = "Data/BK_Analysis_Plot.rds")
BK_Eval = Case_Study_Evaluation(Degree_Distribution,
                                Sample_Degree_CI_OneBNB,
                                Sample_Degree_CI_CVE,
                                Sample_Degree_CI_MRME,
                                Sample_Degree_CI_CHTE,
                                Sample_Degree_CI_SnowBoot)
BK_Eval
#capture.output(BK_Eval, file = "Data/BK_Eval.txt")

############################################################################
# 3. Electrical Grid
## Read
g = read.graph("Data/power.gml", format = "gml")
graph.density(g)
length(V(g))
length(E(g))
graph.density(g)
deg = degree(g)
fit_pl = displ$new(deg)
est = estimate_pars(fit_pl)
Scaling_parameter = est$pars
deg_dist = degree_distribution(g)
deg_dist = deg_dist[-1]
n = seq(1,length(deg_dist))
Deg_dist = cbind(n,deg_dist)
Deg_dist = as.data.frame(Deg_dist)
Deg_dist = Deg_dist[-1,]
Deg_dist$n = 1:nrow(Deg_dist)
Degree_Plot = 
  ggplot(Deg_dist, aes(x = n, y = deg_dist)) +
  geom_line(color = 'red') +
  geom_point(size = 1, color = 'red') +
  theme_bw() +
  scale_x_continuous(breaks = c(1,seq(5,20,5))) +
  xlab("degree") + ylab("proportion of nodes") +
  xlim(c(1,20)) +
  theme(text = element_text(family="Times New Roman"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
Degree_Plot

## Sample 
vid_sample = Node_Sample(g,0.005,1)
vid_snow   = Snow_Sample(g,vid_sample,2)
g_sample   = subgraph(g,vid_snow)
graph.density(g)
deg_dist_sample = degree_distribution(g_sample)
deg_dist_sample = deg_dist_sample[-1]
n = seq(0,length(deg_dist_sample)-1)
Deg_dist_sample = cbind(n,deg_dist_sample)
Deg_dist_sample = as.data.frame(Deg_dist_sample)
W_sample = as_adjacency_matrix(g_sample)
W_sample = as.matrix(W_sample)
#write.csv(W_sample,file = "Data/Power_sample1w2.csv") 
Degree_Sequence = 0:17
Degree_Distribution = Deg_dist[Degree_Sequence + 1,]
colnames(Degree_Distribution) = c("degree","proportion")
Degree_Distribution$degree = Degree_Sequence + 1

## 1-BNB
Sample_Yhat = readMat("Data/Power_sample1w2_15p_B200.mat")
Sample_Yhat = Sample_Yhat$BootsResult
Sample_Degree_CI_OneBNB =
  OneBNB_Degree_CI(Degree_Sequence = Degree_Sequence,
                   Bootstrap_Samples = Sample_Yhat,
                   Sample_Number = 200,
                   Alpha = 0.05)
Sample_Degree_CI_OneBNB$degree = Degree_Sequence + 1
Degree_CI_Plot_1BNB =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_OneBNB, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_1BNB

## CVE
numCores <- detectCores()
numCores
registerDoParallel(numCores)
p = 0.15
B = 200
maxlen = 18
lambda = 0
delta  = 0.01
m = length(g_sample)
g = g_sample 
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    id_snow = Snow_Sample(g,id_sample,nwave = 0)
    g_sample = subgraph(g,id_snow)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Degree_sample = Count_Degree(Y_sample)
    Nstar = as.numeric(Degree_sample$Freq)
    Nstar = fill_fun(Nstar, maxlen)
    P = diag(p,length(Nstar))
    C = p*(1-p)*diag(Nstar) + delta*diag(1,length(Nstar)) #Ego-centric sampling
    D = cbind(diag(1,maxlen-2),matrix(0,maxlen-2,2)) + 
      cbind(matrix(0,maxlen-2,1),diag(-2,maxlen-2),matrix(0,maxlen-2,1)) +
      cbind(matrix(0,maxlen-2,2),diag(1,maxlen-2))
    estimate = CVX_Estimate(m,Nstar,P,C,D,lambda)
    estimate$Deg
  }
Sample_Degree_CI_CVE = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CVE = rbind(Sample_Degree_CI_CVE, temp)
}
colnames(Sample_Degree_CI_CVE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CVE) = 1:maxlen
Sample_Degree_CI_CVE = as.data.frame(Sample_Degree_CI_CVE)
Sample_Degree_CI_CVE$degree = Degree_Sequence + 1
Sample_Degree_CI_CVE$up[which(Sample_Degree_CI_CVE$up<1e-04)] = 0
Sample_Degree_CI_CVE$low[which(Sample_Degree_CI_CVE$low<1e-04)] = 0
Degree_CI_Plot_CVE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CVE, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_CVE

## MRME
p = 0.15
B = 200
maxlen = 18
m = length(g_sample)
g = g_sample
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    print(paste0("n=", n, ",b=", b))
    id_sample = Node_Sample(g,
                            prob = p,
                            seed = b)
    id_snow = Snow_Sample(g,id_sample,nwave = 0)
    g_sample = subgraph(g,id_snow)
    Y_sample = as_adjacency_matrix(g_sample)
    Y_sample = as.matrix(Y_sample)
    Dstar    = apply(Y_sample, 2, sum)
    Dmat     = Y_sample%*%Y_sample
    diag(Dmat) = Dstar
    Dhat      = 1/p* Dstar %*%t(Dstar) %*% 
      solve(Dstar %*% t(Dstar) + Dmat + 1*diag(1,nrow(Y_sample))) %*%
      Dstar
    Dhat = floor(Dhat/(sum(floor((Dhat))/m)))
    degree_freq = mutate(data.frame(table(list(Dhat))),
                         proportion = Freq/length(Dhat),
                         logodds = log(proportion/(1-proportion)))
    estimate    = degree_freq$proportion
    estimate    = fill_fun(estimate,maxlen)
  }
Sample_Degree_CI_MRME = NULL
for (j in 1:maxlen) {
  f = Degree_estimate_bootstrap[j,]
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_MRME = rbind(Sample_Degree_CI_MRME, temp)
}
colnames(Sample_Degree_CI_MRME) = c("degree", "low", "up")
rownames(Sample_Degree_CI_MRME) = 1:maxlen
Sample_Degree_CI_MRME = as.data.frame(Sample_Degree_CI_MRME)
Sample_Degree_CI_MRME$degree = Degree_Sequence + 1
Degree_CI_Plot_MRME =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_MRME, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_MRME

## CHTE
maxlen = 18
EdgeNum = length(E(g))
p       = 0.5
B = 200
m = length(g_sample)
g = g_sample
Degree_estimate_bootstrap <-
  foreach(b = 1:B, .combine=cbind) %dopar% {
    print(paste0("n=", n, ",b=", b))
    Y_sample = Pair_Sample(g,p,b)
    EdgeSet = list()
    for (i in 1:nrow(Y_sample)) {
      EdgeSet[[i]] = which(Y_sample[i,] == 1)
    }
    maxDeg = max(sapply(EdgeSet,length))
    EdgeNum_sample = sum(Y_sample[upper.tri(Y_sample)])
    gammaCount = function(vec,Count,vecNum){
      maxCount = length(vec)
      output = choose(maxCount,Count)
      return(output)
    }
    L = list()
    maxlen = 9
    for (c in 1:maxlen) {
      print(c)
      L[[c]] = sapply(EdgeSet,gammaCount,Count = c,vecNum = nrow(Y_sample))
    }
    Lvec = sapply(L,sum)
    Pi   = rep(0,maxlen)
    for (c in 1:maxlen) {
      Pi[c] = choose(EdgeNum,EdgeNum_sample)^(-1) * choose(EdgeNum-c,EdgeNum_sample-c)
    }
    Lhat = Pi * Lvec
    TmatAll = matrix(-1,max(maxDeg,maxlen),max(maxDeg,maxlen))
    for (i in 1:max(maxDeg,maxlen)) {
      for (j in 1:max(maxDeg,maxlen)) {
        TmatAll[i,j] = choose(j,i)
      }
    }
    if (rcond(TmatAll)>=1e-16){
      Tmat = solve(TmatAll)[1:maxlen,1:maxlen]
      Dhat = solve(Tmat) %*% Lhat
    }else{
      Tmat = ginv(TmatAll)[1:maxlen,1:maxlen]
      Dhat = ginv(Tmat) %*% Lhat
    }
    fhat = Dhat/sum(Dhat)
    fhat = Dhat/sum(Dhat)
  }

Sample_Degree_CI_CHTE = NULL
for (j in 1:maxlen) {
  if (j <= nrow(Degree_estimate_bootstrap)){
    f = Degree_estimate_bootstrap[j,]
  }else{
    f = Degree_estimate_bootstrap[nrow(Degree_estimate_bootstrap),]
  }
  f_sort = sort(f)
  alpha = 0.05
  low = f_sort[round(B*alpha/2)]
  up = f_sort[round(B*(1-alpha/2))]
  
  ## Store CI results in 'degree_CI' for all degree in 'Degree_Sequence'.
  temp = c(j, low, up)
  Sample_Degree_CI_CHTE = rbind(Sample_Degree_CI_CHTE, temp)
}
colnames(Sample_Degree_CI_CHTE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CHTE) = 1:maxlen
Sample_Degree_CI_CHTE = as.data.frame(Sample_Degree_CI_CHTE)
colnames(Sample_Degree_CI_CHTE) = c("degree", "low", "up")
rownames(Sample_Degree_CI_CHTE) = 1:maxlen
Sample_Degree_CI_CHTE = as.data.frame(Sample_Degree_CI_CHTE)
Sample_Degree_CI_CHTE$degree = Degree_Sequence + 1
Degree_CI_Plot_CHTE =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_CHTE, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_CHTE

## SnowBoot
Sample_Degree_CI_SnowBoot =
  SnowBoot_Degree_CI2(Degree_Sequence = Degree_Sequence,
                     Adjacency_Matrix = W_sample,
                     Sample_Number = 1000, n.wave = 2,
                     Alpha = 0.05)
Sample_Degree_CI_SnowBoot = rbind(Sample_Degree_CI_SnowBoot,
                                  Sample_Degree_CI_SnowBoot[nrow(Sample_Degree_CI_SnowBoot),])
Sample_Degree_CI_SnowBoot = Sample_Degree_CI_SnowBoot[-1,]
Sample_Degree_CI_SnowBoot$degree = Degree_Sequence + 1
Degree_CI_Plot_SnowBoot =
  Degree_CI_Plot(Degree_CI = Sample_Degree_CI_SnowBoot, 
                 Degree_Distribution = Degree_Distribution,
                 Degree_Sequence = Degree_Sequence + 1)
Degree_CI_Plot_SnowBoot

## Plot
temp = rbind(Sample_Degree_CI_OneBNB,
             Sample_Degree_CI_CVE,
             Sample_Degree_CI_MRME,
             Sample_Degree_CI_CHTE,
             Sample_Degree_CI_SnowBoot)
max(temp$up)

Degree_Plot_adj             = Degree_Plot + ylim(c(0,0.55))
Degree_CI_Plot_1BNB_adj     = Degree_CI_Plot_1BNB + ylim(c(0,0.55))
Degree_CI_Plot_CVE_adj      = Degree_CI_Plot_CVE + ylim(c(0,0.55))
Degree_CI_Plot_MRME_adj     = Degree_CI_Plot_MRME + ylim(c(0,0.55))
Degree_CI_Plot_CHTE_adj     = Degree_CI_Plot_CHTE + ylim(c(0,0.55))
Degree_CI_Plot_SnowBoot_adj = Degree_CI_Plot_SnowBoot + ylim(c(0,0.55))

Power_Analysis_Plot = list()
Power_Analysis_Plot[[1]] = Degree_Plot_adj
Power_Analysis_Plot[[2]] = Degree_CI_Plot_1BNB_adj
Power_Analysis_Plot[[3]] = Degree_CI_Plot_CVE_adj
Power_Analysis_Plot[[4]] = Degree_CI_Plot_MRME_adj
Power_Analysis_Plot[[5]] = Degree_CI_Plot_CHTE_adj
Power_Analysis_Plot[[6]] = Degree_CI_Plot_SnowBoot_adj
names(Power_Analysis_Plot) = c("Original","1BNB","CVE","MRME","CHTE","SnowBoot")
#saveRDS(Power_Analysis_Plot,file = "Data/Power_Analysis_Plot.rds")

Power_Eval = Case_Study_Evaluation(Degree_Distribution,
                                   Sample_Degree_CI_OneBNB,
                                   Sample_Degree_CI_CVE,
                                   Sample_Degree_CI_MRME,
                                   Sample_Degree_CI_CHTE,
                                   Sample_Degree_CI_SnowBoot)
Power_Eval
#capture.output(Power_Eval,file = 'Data/Power_Eval.txt')

