library(mgcv)# rmvn
library(splines) #base spline (bs)
library(Matrix) # KhatriRhao matrix
library(lme4) #Lmer
library(LaplacesDemon) #normal inverse wishart
library(ClusterR)#GMM
library(mcclust)#arandi
library(mclust)
library(dbscan)
library(fpc)#pamk
library(pracma) # tic/toc
library(Spectrum)#estimate_k/cluster_similarity
library(aricode)#cluster comparison
setwd("C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/")
Rcpp::sourceCpp("ProfileGibbsCPPV2.cpp")
source("paper_utility.R")
source("utility_functions.R")

set.seed(2)
seeds = sample(40000,500)
nInd = 900
nR = 3
nExp = 10


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



RMSE_Fe = array(0,dim =c(5,nExp,length(listSU)),dimnames=list(c('ProLMM','2step','True','TrueClus','ProLMMStar'),
                                                              c(1:nExp),listSU))

RMSE_Re = RMSE_Fe
RMSE_Lat = array(0,dim =c(5,nExp,length(listSU)),dimnames=list(c('ProLMM','2step','True','TrueClus','ProLMMStar'),
                                                               c(1:nExp),listSU))

nC_st= array(0,dim =c(2,nExp,length(listSU)),dimnames=list(c('ProLMM','2step'),
                                                           c(1:nExp),listSU))
randI_st= array(0,dim =c(3,nExp,length(listSU)),dimnames=list(c('ProLMM','2step','TrueClus'),
                                                              c(1:nExp),listSU))
ami_st = array(0,dim =c(3,nExp,length(listSU)),dimnames=list(c('ProLMM','2step','TrueClus'),
                                                             c(1:nExp),listSU))
purity = array(0,dim =c(3,nExp,length(listSU)),dimnames=list(c('ProLMM','2step','TrueClus'),
                                                             c(1:nExp),listSU))
 
REnames = c('RE1','RE2')
FEnames = c('FEgender','FE1','FE2','FE3')
Latnames =  c('FE2')
Assignnames = c('Exp1','Exp2')

theta0 = {} # parameters to estimate
theta0$betaFE = 2*rnorm(length(FEnames)+1) #adding intercept
theta0$WLat = 4*diag(length(Latnames)+1)
theta0$SigRE = diag(length(REnames))
theta0$sigma = 0.5
# inter = c(-4,-1,2,-3,0,3,-2,1,4)*2
for(cidx in 1:9){
  theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)] = rmvn(1,rep(0,length(Latnames)+1),theta0$WLat)
  # theta0$alphaLat[1+(length(Latnames)+1)*(cidx-1)]=inter[cidx]
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
for(e in 1:nExp){
  for(SU in listSU){#list('CorSimp')
    print(paste(SU,e))
    j=j+1
    set.seed(seeds[j])
    tic()
    sim =  sim_lif_v2(nInd,nR,SetUps[[SU]]$paramLat,theta0)
    

    
    dataPro = process_Data_outcome(sim$covList, sim$df, intercept = list(FE=T,RE=F,Lat =T))
    
    nC = 30
    # prior 
    priorInit = {}
    # priorInit$assign = {}
    # priorInit$assign$lambda = 0.15
    # priorInit$assign$mu = rep(0,dataPro$params$qU)
    # priorInit$assign$nu = dataPro$params$qU+1
    # priorInit$assign$Psi = 0.5*(dataPro$params$qU+priorInit$assign$nu+1)*diag(dataPro$params$qU)
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
    # candidate = (1:14)[(abs(diffs[1:15 - 1])>(0.5*max(abs(diffs[1:15 - 1]))))] + 1
    # nclus <- candidate[length(candidate)]
    nclus <- which.max(abs(diffs[1:15 - 1])) + 1
    optClus = cluster_similarity(coocc, k = nclus, specalg = "Ng")
    # table(optClus+10,sim$theta0$Lat)
    # # table(gmm_result$classification,sim$theta0$Lat)
    # # optClus3 = pam(1-coocc/nSim,k = nclus)
    # 
    
   
    
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
    
    gmm_result <- Mclust(sim$df[sim$covList$Assign], G = 2:15)
    Zgmm = gmm_result$classification
    
    # f_2way = as.formula(paste0('Y~',
    #                            paste0(sim$covList$FE,collapse = '+'),
    #                            '+(0+',
    #                            paste0(sim$covList$RE,collapse = '+'),
    #                            '|Zind)+(',
    #                            paste0(sim$covList$Lat,collapse = '+'),
    #                            '|Zgmm)'))
    f_2way = as.formula("Y~1+Zgmm+FE2*Zgmm+FEgender+FE1+FE2+FE3+(0+RE1+RE2|Zind)")
    df_2way = sim$df
    df_2way$Zgmm = Zgmm
    lmm_2way = lmer(f_2way,df_2way, REML = F)
    # summary(lmm_2way)

    df_true = sim$df
    df_true$Zgmm = sim$theta0$Lat
    lmm_true = lmer(f_2way,df_true, REML = F)
    
    
    df_trueClus = sim$df
    df_trueClus$Zgmm = ZTrueClus
    lmm_trueClus = lmer(f_2way,df_trueClus, REML = F)
    # summary(lmm_true)
    
    df_Pro = sim$df
    df_Pro$Zgmm = optClus
    lmm_Pro = lmer(f_2way,df_Pro, REML = F)
    # Inference ------------------------------------------------------------------
    # RMSE FE
    b0 = sqrt(sum(beta0_RMSE^2))
    RMSE_Fe[1,e,SU] = sqrt(sum((apply(out$beta, 1, median, na.rm=TRUE)[c(5,2,4)]-beta0_RMSE)^2))/b0
    RMSE_Fe[2,e,SU] = sqrt(sum((fixef(lmm_2way)[FEnames_RMSE]-beta0_RMSE)^2))/b0
    RMSE_Fe[3,e,SU] = sqrt(sum((fixef(lmm_true)[FEnames_RMSE]-beta0_RMSE)^2))/b0
    RMSE_Fe[4,e,SU] = sqrt(sum((fixef(lmm_trueClus)[FEnames_RMSE]-beta0_RMSE)^2))/b0
    RMSE_Fe[5,e,SU] = sqrt(sum((fixef(lmm_Pro)[FEnames_RMSE]-beta0_RMSE)^2))/b0
    
    
    # windows()
    # hist(out$beta[1,])
    # RMSE Var RE
    v0 = sqrt(sum(diag((sim$theta0$SigRE%*%sim$theta0$SigRE))))
    RMSE_Re[1,e,SU] = sqrt(sum(diag(((apply(out$WRE, c(1,2), median, na.rm=TRUE)-sim$theta0$SigRE)%*%(apply(out$WRE, c(1,2), median, na.rm=TRUE)-sim$theta0$SigRE)))))/v0
    RMSE_Re[2,e,SU] = sqrt(sum(diag(((VarCorr(lmm_2way)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_2way)$Zind-sim$theta0$SigRE)))))/v0
    RMSE_Re[3,e,SU] = sqrt(sum(diag(((VarCorr(lmm_true)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_true)$Zind-sim$theta0$SigRE)))))/v0
    RMSE_Re[4,e,SU] = sqrt(sum(diag(((VarCorr(lmm_trueClus)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_trueClus)$Zind-sim$theta0$SigRE)))))/v0
    RMSE_Re[5,e,SU] = sqrt(sum(diag(((VarCorr(lmm_Pro)$Zind-sim$theta0$SigRE)%*%(VarCorr(lmm_Pro)$Zind-sim$theta0$SigRE)))))/v0
    
    # v0 = sqrt(sum(diag((sim$theta0$WLat%*%sim$theta0$WLat))))
    # RMSE_Lat[1,e,SU] = sqrt(sum(diag(((apply(out$WL, c(1,2), median, na.rm=TRUE)-sim$theta0$WLat)%*%(apply(out$WL, c(1,2), median, na.rm=TRUE)-sim$theta0$WLat)))))/v0
    # RMSE_Lat[2,e,SU] = sqrt(sum(diag(((VarCorr(lmm_2way)$Zgmm-sim$theta0$WLat)%*%(VarCorr(lmm_2way)$Zgmm-sim$theta0$WLat)))))/v0
    # RMSE_Lat[3,e,SU] = sqrt(sum(diag(((VarCorr(lmm_true)$Zgmm-sim$theta0$WLat)%*%(VarCorr(lmm_true)$Zgmm-sim$theta0$WLat)))))/v0
    # RMSE_Lat[4,e,SU] = sqrt(sum(diag(((VarCorr(lmm_trueClus)$Zgmm-sim$theta0$WLat)%*%(VarCorr(lmm_trueClus)$Zgmm-sim$theta0$WLat)))))/v0
    # RMSE_Lat[5,e,SU] = sqrt(sum(diag(((VarCorr(lmm_Pro)$Zgmm-sim$theta0$WLat)%*%(VarCorr(lmm_Pro)$Zgmm-sim$theta0$WLat)))))/v0
    
  
    print(RMSE_Fe[,e,SU])
    
    randI_st[1,e,SU] = arandi( sim$theta0$Lat,optClus, adjust = TRUE)
    randI_st[2,e,SU] = arandi( sim$theta0$Lat,Zgmm, adjust = TRUE)
    randI_st[3,e,SU] = arandi( sim$theta0$Lat,ZTrueClus, adjust = TRUE)
    
    purity[1,e,SU] = ClusterPurity(optClus,sim$theta0$Lat)
    purity[2,e,SU] = ClusterPurity(Zgmm,sim$theta0$Lat)
    purity[3,e,SU] = ClusterPurity(ZTrueClus,sim$theta0$Lat)
    
    ami_st[1,e,SU] = AMI(optClus,sim$theta0$Lat)
    ami_st[2,e,SU] = AMI(Zgmm,sim$theta0$Lat)
    ami_st[3,e,SU] = AMI(ZTrueClus,sim$theta0$Lat)
    
    nC_st[1,e,SU] = nclus
    nC_st[2,e,SU] = gmm_result$G
    
    print(purity[,e,SU])
    # matchIdx = matrix(apply(table(optClus,sim$theta0$Lat), 2, which.max))
    # gamma = array(0,dim =c(3,nclus,nSim))
    # for (it in 1:nSim){
    #   for (c_ in 1:9){
    #     gamma[,c_,it] = apply(out$gamma[,out$Z[,it][which(optClus == c_)]+1,it],1,median)
    #   }
    # }
    # apply(gamma, 2, rowMeans)[,matchIdx]
    # theta0$alphaMat
    
  }
}
apply(RMSE_Fe, 1, colMeans)
apply(RMSE_Re, 1, colMeans)
apply(RMSE_Lat, 1, colMeans)



apply(randI_st, 1, colMeans)
apply(purity, 1, colMeans)
apply(ami_st, 1, colMeans)


store = {}
store$RMSE_Fe = RMSE_Fe
store$RMSE_Re = RMSE_Re
store$RMSE_Lat = RMSE_Lat
store$randI_st = randI_st
store$purity = purity
store$ami_st = ami_st
store$nC_st = nC_st
store$theta0 = theta0
store$SetUps = SetUps
saveRDS(store,paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaper.RData'))
# saveRDS(store,paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaperExt.RData'))

# store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaper.RData'))
# 
# centroids = array(0,dim =c(2,nclus,nSim))
# coVar = array(0,dim = c(2,2,nclus))
# gamma = array(0,dim =c(3,nclus,nSim))
# for (it in 1:nSim){
#   for (c_ in unique(optClus)){
#     centroids[,c_,it] = apply(out$muClus[,out$Z[,it][which(optClus == c_)]+1,it],1,median)
#     gamma[,c_,it] = apply(out$gamma[,out$Z[,it][which(optClus == c_)]+1,it],1,mean)
#     coVar[,,c_] = coVar[,,c_] +apply(out$PhiClus[[it]][,,out$Z[which(optClus == c_),it]+1],c(1,2),median)/nSim
#   }
# }
# 
# cen = apply(centroids, 2, rowMeans)
# 
# 
# 
# trueCen = list()
# trueCovar = list()
# c=1
# for(c1 in c(-1,0,1)){
#   for(c2 in c(-1,0,1)){
#     trueCen[[c]]=c(c1,c2)
#     trueCovar[[c]]=SetUps[[SU]]$paramLat$spread*((1-SetUps[[SU]]$paramLat$rho)*diag(2)+SetUps[[SU]]$paramLat$rho)
#     c = c+1
#   }
# }
# library(tidyr) 
# 
# windows()
# plot_gaussian_density_zones(lapply(seq_len(ncol(cen)),function(x) cen[ , x]),
#                             lapply(seq_len(dim(coVar)[3]),function(x) coVar[ , ,x]),x_range = c(-2.1,2.1),y_range = c(-2.1,2.1))
# 
# 
# 
# windows()
# plot_gaussian_density_zones(trueCen,
#                             trueCovar,x_range = c(-2.1,2.1),y_range = c(-2.1,2.1))
# 


# windows()
# plot(sim$df[,c('Exp1','Exp2')],col=ZTrueClus)
# 
# windows()
# plot(sim$df[,c('Exp1','Exp2')],col=Zgmm)
# 
# 
# windows()
# plot(sim$df[,c('Exp1','Exp2')],col=sim$theta0$Lat)
# 
# ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
# 
# id = 4
# idbis = c(1,5,2,3,4)[id]
# intFil=ma(out$beta[idbis,],10)
# windows()
# plot(intFil)
# lines(0*intFil+sim$theta0$betaFE[id])


# 
# gmm_result <- Mclust(sim$df[sim$covList$Assign], G = 2:15)
# Zgmm = gmm_result$classification
# 
# windows()
# plot(gmm_result, what = "classification")
# summary(gmm_result)
# windows()
# plot(sim$df[,c('Exp1','Exp2')],col=sim$theta0$Lat)
# 
# windows()
# plot(gmm_result$BIC)
