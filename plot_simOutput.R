
setwd("C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/")



store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaper.RData'))
# store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaperExt.RData'))
store = readRDS(paste0('C:/Users/VNOB-0731/Documents/GitHub/profileLMM/Restart/simPaperFinal.RData'))


labelSU = c('Scenario 1', 'Scenario 2')
labelMeth = c('Profile LMM', 'True cluster assignment', 'True cluster parameters')



# Exposure distribution---------------------------------------------------


# param Measures
apply(store$RMSE_Fe, 1, colMeans)
apply(store$RMSE_Re, 1, colMeans)
# apply(store$RMSE_Lat, 1, colMeans)
apply(store$purity, 1, colMeans)
store$nC_st






windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(t(store$RMSE_Fe[,,i]),ylim = c(0,0.7),names= labelMeth,sub =labelSU[i] )
}
mtext("RMSE fixed effects", side = 3, line = -3, outer = TRUE, cex=2)

# windows()
# par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
# for(i in 1:length(labelSU)){
#   boxplot(t(store$RMSE_FeFull[,,i]),ylim = c(0,0.5),names= labelMeth,sub =labelSU[i] )
# }
# mtext("RMSE fixed effects Full", side = 3, line = -3, outer = TRUE, cex=2)
windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(t(store$RMSE_FeLat[,,i]),ylim = c(0,0.7),names= labelMeth,sub =labelSU[i] )
}
mtext("RMSE mean latent interraction", side = 3, line = -3, outer = TRUE, cex=2)



windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(t(store$RMSE_Re[,,i]),ylim = c(0,2),names= labelMeth,sub =labelSU[i] )
}
mtext("RMSE  individual random effect covariance matrix", side = 3, line = -3, outer = TRUE, cex=2)



# cluster Measures
apply(store$randI_st, 1, colMeans)
apply(store$purity, 1, colMeans)
apply(store$ami_st, 1, colMeans)

windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(t(store$randI_st[,,i]),ylim = c(0,1),names= labelMeth[1:2],sub =labelSU[i] )
}
mtext("Adjusted Rand index", side = 3, line = -3, outer = TRUE, cex=2)


windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(t(store$purity[,,i]),ylim = c(0,1),names= labelMeth[1:2],sub =labelSU[i] )
}
mtext("Purity score", side = 3, line = -3, outer = TRUE, cex=2)


windows()
par(mfrow = c(1,2), col.sub="blue", cex.sub=2)
for(i in 1:length(labelSU)){
  boxplot(store$nC_st[,i],ylim = c(0,11),sub =labelSU[i] )
}
mtext("Number of clusters", side = 3, line = -3, outer = TRUE, cex=2)



