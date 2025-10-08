
library(mgcv)# rmvn
library(splines) #base spline (bs)
library(Matrix) # KhatriRhao matrix
library(LaplacesDemon) #normal inverse wishart
library(lme4)
library(pracma) #tic toc
Rcpp::sourceCpp("/groups/umcg-lifelines/tmp01/projects/ov21_0374/mamestoy/ProfileREgression/Code/ProfileGibbsCPPV2.cpp")
source("paper_utility.R")
source("/groups/umcg-lifelines/tmp01/projects/ov21_0374/mamestoy/ProfileREgression/Code/utility_functions.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("At least 4 arguments expected: 1. Name of model Correct,Wrong,Flat , 
                                       2. nb of Ind,
                                       3. nb of repeats,
                                       4. nb of experiments", call.=FALSE)
} 



nameSetUp = args[1]
nInd = as.numeric(args[2])#1500
nR = as.numeric(args[3])#3
nbSim = as.numeric(args[4])#40

print(nameSetUp)
print(nInd)
print(nR)
print(nbSim)
set.seed(1)
seeds = sample(40000,nbSim)

paramLat = {}
if (nameSetUp=='Correct')
{
paramLat$name = 'Correct'
paramLat$rho = 0.6
paramLat$spread = 0.21} else if(nameSetUp=='Wrong'){
paramLat$name = 'Wrong'
paramLat$rho = 0.6
paramLat$spread = 1.2} else if(nameSetUp=='Flat'){
  paramLat$name = 'Flat'}


for(simu in 1:nbSim){
  set.seed(seeds[simu])
  tic()
  print(simu)
  sim = sim_lif_v2(nInd,nR,paramLat)
  dataPro = process_Data_outcome(sim$covList, sim$df, intercept = list(FE=T,RE=F,Lat =T))
  
  nC = 30
  priorInit = {}
  priorInit$assign = {}
  priorInit$assign$lambda = 1
  priorInit$assign$mu = rep(0,dataPro$params$qU)
  priorInit$assign$nu = dataPro$params$qU+4
  priorInit$assign$Psi = cov(dataPro$d$U)*(priorInit$assign$nu-dataPro$params$qU-1)
  prior = prior_init(dataPro$params, nC, {})
  startTheta = {}
  f_noLAt = as.formula(paste0('Y~',
                              paste0(sim$covList$FE,collapse = '+'),
                              '+(0+',
                              paste0(sim$covList$RE,collapse = '+'),
                              '|Zind)'
  ))
  lmm_noLAt = lmer(f_noLAt,sim$df,REML =F)
  startTheta$betaFE = fixef(lmm_noLAt)
  theta_Start = theta_init(prior,dataPro$params,nC,startTheta)
  library(pracma)
  
  nSim = 20000
  nBurnIn = 5000
  out = profileLMM_Gibbs(dataPro,nSim+nBurnIn,nBurnIn,nC,prior,theta_Start)
  store = {}
  store$sim = sim
  store$PRout = out
  saveRDS(store,paste0('OutSim/paper_simulation_',paramLat$name,'_',simu,'.RData'))
  toc()
}


