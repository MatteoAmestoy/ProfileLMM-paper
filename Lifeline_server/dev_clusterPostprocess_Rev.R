
args = commandArgs(trailingOnly=TRUE)
appendix = args[1]

print(paste0('Cluster Characteristics-------',appendix))
out = readRDS(paste0(appendix,'.RData'))
optClusStore = readRDS(paste0('optclus_Rev_',appendix,'.RData'))

nbIt = dim(out$Z)[2]
qU = dim(out$muClus)[1]
qL = dim(out$WL)[1]
names = out$names
# clustering and cluster effects -----------------------------------------------
store = {}
store$BetaFE = out$beta
rownames(store$BetaFE) = out$names$FE

for (case in 1:2){
  if (case==1){
    optClus = optClusStore$optClus
    sampleSim = optClusStore$sampleSim
    sampleObs = optClusStore$sampleObs
    nSim = optClusStore$nSim}else
    {
      optClus = optClusStore$optClus2
      sampleSim = optClusStore$sampleSim2
      sampleObs = optClusStore$sampleObs2
      nSim = optClusStore$nSim2}
  
  nobs = length(optClus)
 
  
  coutClus = table(optClus)
  print(coutClus)
  # clusList = which(coutClus>floor(nobs/15))
  clusList = unique(optClus)
  nclus = length(clusList)
  centroids = array(0,dim =c(qU,nclus,nSim))
  coVar = array(0,dim = c(qU,qU,nclus))
  gamma = array(0,dim =c(qL,nclus,nSim))
  print(table(optClus))
  
  j = 1
  for (it in sampleSim){
    cIdx = 1
    for (c_ in clusList){
      tmpMu = out$muClus[,out$Z[sampleObs,it][which(optClus == c_)]+1,it,drop=F]
      centroids[,cIdx,j] = apply(tmpMu,1,mean)
      coVar[,,cIdx] = coVar[,,cIdx] +
        apply(out$PhiClus[[it]][,,out$Z[sampleObs,it][which(optClus == c_)]+1],c(1,2),mean)/nSim+
        drop(tmpMu)%*%t(drop(tmpMu))/nSim/coutClus[c_]
      gamma[,cIdx,j] = apply(out$gamma[,out$Z[sampleObs,it][which(optClus == c_)]+1,it,drop=F],1,mean)
      cIdx = cIdx+1
    }
    j=j+1
  }
  cen = apply(centroids, 2, rowMeans)
  rownames(cen) = out$names$U
  
  for (c_ in 1:nclus){
    coVar[,,c_] = coVar[,,c_]-cen[,c_]%*%t(cen[,c_])
  }
  gam = apply(gamma, 2, rowMeans)
  rownames(gam) = out$names$Lat
  if(case == 1){
  store$gamma = gamma
  store$centroids = centroids
  store$coVar = coVar
  store$optClus = optClus
  store$gam = gam
  store$cen = cen}else{
    store$gamma2 = gamma
    store$centroids2 = centroids
    store$coVar2 = coVar
    store$optClus2 = optClus
    store$gam2 = gam
    store$cen2 = cen
  }
  
}

saveRDS(store,paste0('profileRev_',appendix,'.RData'))

