library(mgcv)# rmvn
library(splines) #base spline (bs)
library(Matrix) # KhatriRhao matrix
library(lme4) #Lmer
library(LaplacesDemon) #normal inverse wishart
library(ClusterR)#GMM
library(mcclust)#arandi
library(mclust)
library(dbscan)
library(tidyr)
library(fpc)#pamk
library(Spectrum)#estimate_k/cluster_similarity
library(aricode)#cluster comparison
library(matlib)
setwd("C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/")
Rcpp::sourceCpp("ProfileGibbsCPPV2.cpp")
source("paper_utility.R")
source("utility_functions.R")

set.seed(3)

store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaper.RData'))
# store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/whynot.RData'))
listSU = names(store$SetUps)
# listSU = list('Wrong')
theta0 = store$theta0
SetUps = store$SetUps

out = {}
nInd = 1000
nR = 3
for(SU in listSU){#
  sim =  sim_lif_v2(nInd,nR,SetUps[[SU]]$paramLat,theta0)
  SetUps[[SU]]$data = sim


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

  out[[SU]] = profileLMM_Gibbs(dataPro,nSim+nBurnIn+1,nBurnIn,nC,prior,theta_Start)


  # Optimal Cluster Estimation ------------------------------------------------------------------
  coocc = create_co_occurrence_matrix_r(out[[SU]]$Z)
  tmp = estimate_k(coocc, maxk = 15,showplots=F)
  diffs <- diff(tmp$evals)
  diffs <- diffs[-1]
  # candidate = (1:14)[(abs(diffs[1:15 - 1])>(0.9*max(abs(diffs[1:15 - 1]))))] + 1
  # nclus <- candidate[length(candidate)]
  nclus <- which.max(abs(diffs[1:15 - 1])) + 1
  optClus = cluster_similarity(coocc, k = nclus)
  out[[SU]]$optClus = optClus
  centroids = array(0,dim =c(2,nclus,nSim))
  out[[SU]]$coVar = array(0,dim = c(2,2,nclus))
  gamma = array(0,dim =c(2,nclus,nSim))
  transfert = apply(table(out[[SU]]$optClus,SetUps[[SU]]$data$theta0$Lat),2,which.max)

  j = 1
  for (it in 1:nSim){
    cIdx = 1
    for (c_ in transfert){
      centroids[,cIdx,j] = apply(out[[SU]]$muClus[,out[[SU]]$Z[,it][which(optClus == c_)]+1,it,drop=F],1,mean)
      out[[SU]]$coVar[,,cIdx] = out[[SU]]$coVar[,,cIdx] +apply(out[[SU]]$PhiClus[[it]][,,out[[SU]]$Z[,it][which(optClus == c_)]+1],c(1,2),median)/nSim
      gamma[,cIdx,j] = apply(out[[SU]]$gamma[,out[[SU]]$Z[,it][which(optClus == c_)]+1,it,drop=F],1,mean)
      cIdx = cIdx+1
    }
    j=j+1
  }
  out[[SU]]$cen = apply(centroids, 2, rowMeans)
  out[[SU]]$gam = apply(gamma, 2, rowMeans)
  pct = 0.05
  out[[SU]]$gamqL = apply(gamma, c(1,2), quantile,pct/2)
  out[[SU]]$gamqU = apply(gamma, c(1,2), quantile,1-pct/2)

}

for (SU in listSU){
  out[[SU]]$transfert = apply(table(out[[SU]]$optClus,SetUps[[SU]]$data$theta0$Lat),1,which.max)
  out[[SU]]$invtransfert = apply(table(out[[SU]]$optClus,SetUps[[SU]]$data$theta0$Lat),2,which.max)
}

library(RColorBrewer)
colors_for_labels <- brewer.pal(n = 9, name = "Set1")
windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
j= 1
for (SU in listSU){
  plot(SetUps[[SU]]$data$df$Exp1,
       SetUps[[SU]]$data$df$Exp2,
       col = colors_for_labels[SetUps[[SU]]$data$theta0$Lat],
       ylim = c(-2,2),xlim = c(-2,2),
       xlab = 'Exposure 1',
       ylab = 'Exposure 2',pch = 20,
       title(paste0('Scenario ',j)))
  j=j+1
}
j = 1
windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for (SU in listSU){
  plot(SetUps[[SU]]$data$df$Exp1,
       SetUps[[SU]]$data$df$Exp2,
       col = colors_for_labels[SetUps[[SU]]$data$theta0$Lat],
       ylim = c(-2,2),xlim = c(-2,2),
       xlab = 'Exposure 1',
       ylab = 'Exposure 2',pch = 20,
       title(paste0('Scenario ',j)))
  j=j+1
  for(cl in 1:dim(out[[SU]]$coVar)[3]){
    aaa = ellipse(
      out[[SU]]$cen[ , cl],
      out[[SU]]$coVar[ , ,cl]
      ,0.3,100)
    lines(aaa$xx,aaa$yy,
          col = colors_for_labels[cl],
          lwd = 3)
  }
}

trueCen = list()
trueCovar = list()
trueCovar2 = list()
c=1
for(c1 in c(-1,0,1)){
  for(c2 in c(-1,0,1)){
    trueCen[[c]]=c(c1,c2)
    trueCovar[[c]]=SetUps[[SU]]$paramLat$spreadClus*((1-SetUps[[SU]]$paramLat$rhoClus)*diag(2)+SetUps[[SU]]$paramLat$rhoClus)
    trueCovar2[[c]]= diag(2)
    trueCovar2[[c]][1,1] = 0.5
    c = c+1
  }
}

listMeans = list(trueCen,
                 lapply(1:9,function(x) out[['CorHard']]$cen[ , x]),
                 lapply(1:9,function(x) out[['Wrong']]$cen[ , x])
)

listCovar = list(trueCovar,
                 lapply(1:9,function(x) out[['CorHard']]$coVar[ , ,x]),
                 lapply(1:9,function(x) out[['Wrong']]$coVar[ , ,x]))



windows()
plot_multiple_gaussian_density_zones(listMeans, listCovar,
                                     c('True','Scenario 1','Scenario 2'),
                                     x_range = c(-2.1,2.1),y_range = c(-2.1,2.1),
                                     grid_res = 100, contour_levels = NULL,
                                     list_of_show_points = FALSE,
                                     list_of_data_points = NULL,
                                     component_colors = colors_for_labels,
                                     list_of_point_colors = "black",
                                     list_of_x_axis_labels = "Exposure 1",
                                     list_of_y_axis_labels = "Exposure 2" )



GamTabInter = data.frame(matrix(0,nrow = 3,ncol = 9 ))
GamTabVar = data.frame(matrix(0,nrow = 3,ncol = 9 ))
GamTabInter[1,] = SetUps$Wrong$data$theta0$alphaMat[1,]
GamTabVar[1,] = SetUps$Wrong$data$theta0$alphaMat[2,]
j=1
for (SU in listSU){
  GamTabInter[1+j,] = out[[SU]]$gam[1,]
  GamTabVar[1+j,] = out[[SU]]$gam[2,]
  j=j+1}

colnames(GamTabInter)= SetUps$Wrong$data$theta0$label
colnames(GamTabVar)= SetUps$Wrong$data$theta0$label
rownames(GamTabInter)= c('True','Scenario 1','Scenario 2')
rownames(GamTabVar)= c('True','Scenario 1','Scenario 2')

library(xtable)



xtable_Inter <- xtable(GamTabInter,
                       caption = "Estimated latent cluster specific intercept",
                       label = "tab:intercept",
                       align = rep('c',10), # Alignment for each column (l=left, r=right, c=center)
                       digits = rep(2,10)) # Number of decimal places for each column (0 for non-numeric)

xtable_Var <- xtable(GamTabVar,
                       caption = "Estimated latent cluster specific random slope",
                       label = "tab:slope",
                       align = rep('c',10), # Alignment for each column (l=left, r=right, c=center)
                       digits = rep(2,10)) # Number of decimal places for each column (0 for non-numeric)

# Print the LaTeX code
cat("\n--- LaTeX Table using xtable ---\n")
print(xtable_Inter,
      type = "latex",
      include.rownames = T, # Do not include R row names
      caption.placement = "bottom", # Place caption at the top
      hline.after = c(-1, 0, 3), # Add horizontal lines after header and at end
      floating = TRUE, # Create a floating table environment (table, caption, label)
      latex.environments = c("center") # Center the table
)
print(xtable_Var,
      type = "latex",
      include.rownames = T, # Do not include R row names
      caption.placement = "bottom", # Place caption at the top
      hline.after = c(-1, 0, 3), # Add horizontal lines after header and at end
      floating = TRUE, # Create a floating table environment (table, caption, label)
      latex.environments = c("center") # Center the table
)
