library(Spectrum)
source("paper_utility.R")
# setwd('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/LifelinesServer/')

args = commandArgs(trailingOnly=TRUE)
appendix = args[1]
print(paste0('Clustering-------',appendix))
out = readRDS(paste0(appendix,'.RData'))

nbIt = dim(out$Z)[2]
qU = dim(out$muClus)[1]
qL = dim(out$WL)[1]
names = out$names
# clustering and cluster effects -----------------------------------------------


nobs = as.numeric(args[2])#dim(out$Z)[1]
nSim = as.numeric(args[3])

sampleObs = sample(1:dim(out$Z)[1],nobs,replace = F)
sampleSim = (nbIt-nSim+1):nbIt

coocc = create_co_occurrence_matrix_r(out$Z[sampleObs,sampleSim])
tmp = estimate_k(coocc, maxk = 25,showplots=F)
diffs <- diff(tmp$evals)
diffs <- diffs[-1]
# candidate = (1:14)[(abs(diffs[1:15 - 1])>(0.5*max(abs(diffs[1:15 - 1]))))] + 1
# nclus <- candidate[length(candidate)]
nclus <- which.max(abs(diffs[1:25 - 1])) + 1
optClus = cluster_similarity(coocc, k = nclus, specalg = "Ng")

print(table(optClus))

optClusStore = {}
optClusStore$optClus = optClus
optClusStore$nSim = nSim
optClusStore$sampleSim = sampleSim
optClusStore$sampleObs = sampleObs

saveRDS(optClusStore,paste0('optclus',appendix,'.RData'))
