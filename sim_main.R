library(mgcv)# rmvn
library(splines) #base spline (bs)
library(Matrix) # KhatriRhao matrix
library(lme4) #Lmer
library(LaplacesDemon) #normal inverse wishart
library(ClusterR)#GMM
library(mcclust)#arandi
library(mclust) #dmvnorm
# library(dbscan)
# library(fpc)#pamk
library(pracma) # tic/toc zeros
library(Spectrum)#estimate_k/cluster_similarity
library(aricode)#cluster comparison
setwd("C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/")
Rcpp::sourceCpp("ProfileGibbsCPPV2.cpp")
source("paper_utility.R")
source("utility_functions.R")

set.seed(20)
seeds = sample(40000,500)
nInd = 1000
nR = 3
nExp = 25


listSU = list('CorHard','Wrong')

SetUps = {}
rhoClus = 0.0
spreadClus = 0.3
SetUps[['CorHard']] = {}
SetUps[['CorHard']][['paramLat']]={}
SetUps[['CorHard']]$paramLat$name = 'Correct'
SetUps[['CorHard']]$paramLat$rho = 0.0
SetUps[['CorHard']]$paramLat$spread = 0.15
SetUps[['CorHard']]$paramLat$rhoClus = rhoClus
SetUps[['CorHard']]$paramLat$spreadClus = spreadClus
SetUps[['CorHard']]$paramLat$noise = 'low'

SetUps[['Wrong']] = {}
SetUps[['Wrong']][['paramLat']]={}
SetUps[['Wrong']]$paramLat$name = 'Wrong'
SetUps[['Wrong']]$paramLat$rho = 0.5
SetUps[['Wrong']]$paramLat$spread = 0.8
SetUps[['Wrong']]$paramLat$rhoClus = rhoClus
SetUps[['Wrong']]$paramLat$spreadClus = spreadClus
SetUps[['Wrong']]$paramLat$noise = 'low'



RMSE_Fe = array(0,dim =c(3,nExp,length(listSU)),dimnames=list(c('ProLMM','TrueClus','True'),
                                                              c(1:nExp),listSU))

RMSE_FeFull = RMSE_Fe
RMSE_FeLat = RMSE_Fe

RMSE_Re = RMSE_Fe
nC_st= array(0,dim =c(nExp,length(listSU)),dimnames=list(c(1:nExp),listSU))
randI_st= array(0,dim =c(2,nExp,length(listSU)),dimnames=list(c('ProLMM','TrueClus'),
                                                              c(1:nExp),listSU))
ami_st = array(0,dim =c(2,nExp,length(listSU)),dimnames=list(c('ProLMM','TrueClus'),
                                                             c(1:nExp),listSU))
purity = array(0,dim =c(2,nExp,length(listSU)),dimnames=list(c('ProLMM','TrueClus'),
                                                             c(1:nExp),listSU))

REnames = c('RE1','RE2')
FEnames = c('FEgender','FE1','FE2','FE3')
Latnames =  c('FE2')
Assignnames = c('Exp1','Exp2')

theta0 = {} # parameters to estimate
theta0$betaFE = rnorm(length(FEnames)+1) #adding intercept
theta0$WLat = diag(length(Latnames)+1)
theta0$SigRE = diag(length(REnames))
theta0$sigma = 0.5
inter = c(-4,-1,2,-3,0,3,-2,1,4)
for(cidx in 1:9){
  theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)] = rmvn(1,rep(0,length(Latnames)+1),theta0$WLat)/10
  theta0$alphaLat[1+(length(Latnames)+1)*(cidx-1)]=inter[cidx]/10
  print(theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)])
}
bbb = (0:8)*(length(Latnames)+1)
for(cidx in 1:(length(Latnames)+1)){
  theta0$alphaLat[bbb+cidx] = theta0$alphaLat[bbb+cidx]-mean(theta0$alphaLat[bbb+cidx])
}

theta0$alphaMat = zeros((length(Latnames)+1),9)
for(cidx in 1:9){
  theta0$alphaMat[,cidx] = theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)]
}
j = 1


FEnames_RMSE = c('FEgender1','FE1','FE3')
beta0_RMSE = theta0$betaFE[c(2,3,5)]

FEnames_RMSEFull = c('(Intercept)','FEgender1','FE1','FE2','FE3')
beta0_RMSEFull = theta0$betaFE

FEnames_RMSELat = c('(Intercept)','FE2')
beta0_RMSELat = theta0$betaFE[c(1,4)]


for(e in 1:nExp){
  for(SU in listSU){
    print(paste(SU,e))
    j=j+1
    set.seed(seeds[j])
    tic()
    sim =  sim_lif_v2(nInd,nR,SetUps[[SU]]$paramLat,theta0)



    dataPro = process_Data_outcome(sim$covList, sim$df, intercept = list(FE=T,RE=F,Lat =T))

    nC = 30
    # prior
    priorInit = {}
    prior = prior_init(dataPro$params, nC, priorInit)

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

    if (SU == 'CorHard' ){
      nSim = 3000
      nBurnIn = 1000
    }else{
      nSim = 3000
      nBurnIn = 2000
    }

    out = profileLMM_Gibbs(dataPro,nSim+nBurnIn+1,nBurnIn,nC,prior,theta_Start)


    # Optimal Cluster Estimation ------------------------------------------------------------------
    coocc = create_co_occurrence_matrix_r(out$Z)
    tmp = estimate_k(coocc, maxk = 15,showplots=F)
    diffs <- diff(tmp$evals)
    diffs <- diffs[-1]
    nclus <- which.max(abs(diffs[1:15 - 1])) + 1
    optClus = cluster_similarity(coocc, k = nclus, specalg = "Ng")
    # table(optClus+10,sim$theta0$Lat)

    # Post process and comparison
    density_list = list()

    trueCen = list()
    trueCovar = list()
    trueCovar2 = list()
    c=1
    for(c1 in c(-1,0,1)){
      for(c2 in c(-1,0,1)){
        trueCovar=SetUps[[SU]]$paramLat$spreadClus*((1-SetUps[[SU]]$paramLat$rhoClus)*diag(2)+SetUps[[SU]]$paramLat$rhoClus)
        density_list[[c]] <- dmvnorm(as.matrix(sim$df[,c('Exp1','Exp2')]), mean = c(c1,c2), sigma =trueCovar)
        c=c+1
      }
    }
    density_matrix <- do.call(cbind, density_list)
    ZTrueClus <- as.factor(apply(density_matrix, 1, which.max))

    f_lmer = as.formula("Y~1+FE2*Zgmm+FEgender+FE1+FE3+(0+RE1+RE2|Zind)")

    df_true = sim$df
    df_true$Zgmm = sim$theta0$Lat
    lmm_true = lmer(f_lmer,df_true, REML = F)

    df_trueClus = sim$df
    df_trueClus$Zgmm = ZTrueClus
    lmm_trueClus = lmer(f_lmer,df_trueClus, REML = F)

    # Inference ------------------------------------------------------------------
    # RMSE FE
    b0 = sqrt(sum(beta0_RMSE^2))
    RMSE_Fe[1,e,SU] = sqrt(sum((apply(out$beta, 1, median, na.rm=TRUE)[c(5,2,4)]-beta0_RMSE)^2))/b0
    RMSE_Fe[2,e,SU] = sqrt(sum((fixef(lmm_trueClus)[FEnames_RMSE]-beta0_RMSE)^2))/b0
    RMSE_Fe[3,e,SU] = sqrt(sum((fixef(lmm_true)[FEnames_RMSE]-beta0_RMSE)^2))/b0


    b0 = sqrt(sum(beta0_RMSEFull^2))
    RMSE_FeFull[1,e,SU] = sqrt(sum((apply(out$beta, 1, median, na.rm=TRUE)[c(1,5,2,3,4)]-beta0_RMSEFull)^2))/b0
    beta0_RMSE_ = beta0_RMSEFull
    beta0_RMSE_[1] = beta0_RMSEFull[1]+theta0$alphaMat[1,1]
    beta0_RMSE_[4] = beta0_RMSEFull[4]+theta0$alphaMat[2,1]

    RMSE_FeFull[2,e,SU] = sqrt(sum((fixef(lmm_trueClus)[FEnames_RMSEFull]-beta0_RMSE_)^2))/b0
    RMSE_FeFull[3,e,SU] = sqrt(sum((fixef(lmm_true)[FEnames_RMSEFull]-beta0_RMSE_)^2))/b0

    b0 = sqrt(sum(beta0_RMSELat^2))
    RMSE_FeLat[1,e,SU] = sqrt(sum((apply(out$beta, 1, median, na.rm=TRUE)[c(1,3)]-beta0_RMSELat)^2))/b0
    beta0_RMSE_ = beta0_RMSELat
    beta0_RMSE_[1] = beta0_RMSELat[1]+theta0$alphaMat[1,1]
    beta0_RMSE_[2] = beta0_RMSELat[2]+theta0$alphaMat[2,1]

    RMSE_FeLat[2,e,SU] = sqrt(sum((fixef(lmm_trueClus)[FEnames_RMSELat]-beta0_RMSE_)^2))/b0
    RMSE_FeLat[3,e,SU] = sqrt(sum((fixef(lmm_true)[FEnames_RMSELat]-beta0_RMSE_)^2))/b0


    # RMSE Var RE
    v0 = sqrt(sum(diag((sim$theta0$SigRE%*%sim$theta0$SigRE))))
    RMSE_Re[1,e,SU] = sqrt(sum(diag(((apply(out$WRE, c(1,2), median, na.rm=TRUE)-sim$theta0$SigRE)%*%(apply(out$WRE, c(1,2), median, na.rm=TRUE)-sim$theta0$SigRE)))))/v0
    RMSE_Re[2,e,SU] = sqrt(sum(diag(((VarCorr(lmm_trueClus)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_trueClus)$Zind-sim$theta0$SigRE)))))/v0
    RMSE_Re[3,e,SU] = sqrt(sum(diag(((VarCorr(lmm_true)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_true)$Zind-sim$theta0$SigRE)))))/v0

    print(RMSE_Fe[,e,SU])

    randI_st[1,e,SU] = arandi( sim$theta0$Lat,optClus, adjust = TRUE)
    randI_st[2,e,SU] = arandi( sim$theta0$Lat,ZTrueClus, adjust = TRUE)

    purity[1,e,SU] = ClusterPurity(optClus,sim$theta0$Lat)
    purity[2,e,SU] = ClusterPurity(ZTrueClus,sim$theta0$Lat)

    ami_st[1,e,SU] = AMI(optClus,sim$theta0$Lat)
    ami_st[2,e,SU] = AMI(ZTrueClus,sim$theta0$Lat)

    nC_st[e,SU] = nclus

    print(purity[,e,SU])
  }
}
apply(RMSE_Fe, 1, colMeans)
apply(RMSE_FeFull, 1, colMeans)
apply(RMSE_FeLat, 1, colMeans)
apply(RMSE_Re, 1, colMeans)
apply(randI_st, 1, colMeans)
apply(purity, 1, colMeans)
apply(ami_st, 1, colMeans)
colMeans(nC_st)

store = {}
store$RMSE_Fe = RMSE_Fe
store$RMSE_FeFull = RMSE_FeFull
store$RMSE_FeLat = RMSE_FeLat
store$RMSE_Re = RMSE_Re
store$randI_st = randI_st
store$purity = purity
store$ami_st = ami_st
store$nC_st = nC_st
store$theta0 = theta0
store$SetUps = SetUps
saveRDS(store,paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaperFinalHard.RData'))
