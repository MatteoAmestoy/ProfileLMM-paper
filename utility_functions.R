############################################
#                                          #
#  Post processing and preprocessing data  #
#                                          #
############################################

## authors: M. Amestoy
# This script centralises all the utility functions of LMM profile regression.
library(caret) # (dummyVars) one hot encoding of factor data

encodeCat = function(dataframe){
  ' Encodes factor variables dummies LAST is reference
-------Input:
      - dataframe: dataframe to convert
-------Output:
      - out: dataframe converted
  Note:
  '
  typeCol = sapply(dataframe, class)
  out = data.frame(dataframe[,typeCol!='factor'])
  colnames(out) = colnames(dataframe)[typeCol!='factor']
  for(col in colnames(dataframe)[typeCol=='factor']){
    dmy <- dummyVars(paste0(" ~ ",col), data = dataframe)
    dfTmp = data.frame(predict(dmy, newdata = dataframe))
    nfact = dim(dfTmp)[2]
    out[colnames(dfTmp)[2:nfact]] = dfTmp[,colnames(dfTmp)[2:nfact]]
  }
  return(out)
}

process_Data_outcome <- function(covList, dataframe, intercept = list(FE=T,RE=T,Lat =T)) {
  ' Preprocess the data from a LMER formula
-------Input:
      - dataframe: dataframe on which to apply the formula
      - covList: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
                    - REunit statistical unit of the RE colomn name/index
                    - Y outcome (Continuous)
      - intercept (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
      - d dictionary with [XFE,XRE,XLat,U,ZRE] design matrices
      - [[params]] list of the parameters of the data
            - n int nb of obs
            - qFE lint, number of covariates of FE
            - nRE int, number of stat units of RE
            - qRE int, number of covariates of RE
  Note:
  '

  d = {}
  d$names = {}
  d$XFE = encodeCat(dataframe[,covList$FE, drop = FALSE])
  d$names$FE = colnames(d$XFE)
  if (intercept$FE){d$XFE = cbind(1,d$XFE)
  d$names$FE= c('Intercept',d$names$FE) }
  d$XFE = as.matrix(d$XFE)

  d$XRE = encodeCat(dataframe[,covList$RE, drop = FALSE])
  d$names$RE = colnames(d$XRE)
  if (intercept$RE){d$XRE = cbind(1,d$XRE)
  d$names$RE= c('Intercept',d$names$RE) }
  d$XRE = as.matrix(d$XRE)

  d$XLat = encodeCat(dataframe[,covList$Lat, drop = FALSE])
  d$names$Lat = colnames(d$XLat)
  if (intercept$Lat){d$XLat = cbind(1,d$XLat)
  d$names$Lat= c('Intercept',d$names$Lat) }
  d$XLat = as.matrix(d$XLat)
  print(1)
  d$U = dataframe[,covList$Assign, drop = FALSE]
  d$U = as.matrix(d$U)
  d$names$U = colnames(d$U)

  d$Y = drop(dataframe[,covList$Y])

  d$ZRE = as.numeric(factor(dataframe[,covList$REunit]))-1

  params = {}
  params$n = dim(d$XFE)[1]
  params$nRE = max(d$ZRE)+1
  params$qFE = dim(d$XFE)[2]
  params$qRE = dim(d$XRE)[2]
  params$qLat = dim(d$XLat)[2]
  params$qU = dim(d$U)[2]


  return(list(d = d, params = params))

}

profileLMM_Gibbs = function(data,nIt,nBurnIn,nCmax,prior,theta_init){
  ' Preprocess the data from a LMER formula
-------Input:
      - data: output from process_Data_outcome
      - covList: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - REunit statistical unit of the RE colomn name/index
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
      - intercept (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
      - d dictionary with [XFE,XRE,XLat,U,ZRE] design matrices
      - [[params]] list of the parameters of the data
            - n int nb of obs
            - qFE lint, number of covariates of FE
            - nRE int, number of stat units of RE
            - qRE int, number of covariates of RE
  Note:
  '

  gibbs_out = GSLoopCPP(nIt, nBurnIn, nCmax,
                        data$d$Y, data$d$XFE, data$d$XRE, data$d$XLat, data$d$U, data$d$ZRE,
                        theta_init$betaFE,theta_init$sig2, theta_init$SigRE, theta_init$muClus,
                        theta_init$SigmaClus, theta_init$SigLat, theta_init$gammaLat,
                        prior$FE$a, prior$FE$b, prior$FE$lambda, prior$RE$Phi,
                        prior$RE$eta, prior$assign$lambda, prior$assign$mu, prior$assign$nu,
                        prior$assign$Psi, prior$Lat$Phi, prior$Lat$eta, prior$DP$scale, prior$DP$shape)


  return(gibbs_out)
}

prior_init = function(params, nC, priorInit){
  ' initialise the priors of the model
-------Input:
      - params: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - REunit statistical unit of the RE colomn name/index
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
      - priorInit (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
      - prior list of the parameters of the data
            - n int nb of obs
            - qFE lint, number of covariates of FE
            - nRE int, number of stat units of RE
            - qRE int, number of covariates of RE
  Note:
  '
  prior ={}
  if('FE' %in% names(priorInit)){
    prior$FE = priorInit$FE
  }else{
    prior$FE = {}
    prior$FE$lambda = 10**(-6)
    prior$FE$a = 10**(-6)
    prior$FE$b = 10**(-6)
  }
  if('RE' %in% names(priorInit)){
    prior$RE = priorInit$RE
  }else{
    prior$RE = {}
    prior$RE$Phi = diag(params$qRE)
    prior$RE$eta = params$qRE
  }
  if('assign' %in% names(priorInit)){
    prior$assign = priorInit$assign
  }else{
    prior$assign = {}
    prior$assign$lambda = 1
    prior$assign$mu = rep(0,params$qU)
    prior$assign$nu = params$qU+1
    prior$assign$Psi = diag(params$qU)
  }
  if('Lat' %in% names(priorInit)){
    prior$Lat = priorInit$Lat
  }else{
    prior$Lat = {}
    prior$Lat$eta = params$qLat
    prior$Lat$Phi = diag(params$qLat)


  }
  if('DP' %in% names(priorInit)){
    prior$DP = priorInit$DP
  }else{
    prior$DP = {}
    prior$DP$scale = sqrt(nC)
    prior$DP$shape = sqrt(nC)
  }

  return(prior)
}

theta_init = function(prior,params,nC,init){
  ' Preprocess the data from a LMER formula
-------Input:
      - data: output from process_Data_outcome
      - covList: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - REunit statistical unit of the RE colomn name/index
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
      - intercept (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
      - d dictionary with [XFE,XRE,XLat,U,ZRE] design matrices
      - [[params]] list of the parameters of the data
            - n int nb of obs
            - qFE lint, number of covariates of FE
            - nRE int, number of stat units of RE
            - qRE int, number of covariates of RE
  Note:
  '
  theta = {}

  if('sig2' %in% names(init)){
    theta$sig2 = init$sig2
  }else{
    theta$sig2 = rgamma(1,1,1)#rgamma(1,prior$FE$a,prior$FE$b)
  }
  if('betaFE' %in% names(init)){
    theta$betaFE = init$betaFE
  }else{
    theta$betaFE = rnorm(params$qFE,0,theta$sig2/sqrt(prior$FE$lambda))
  }

  if('SigRE' %in% names(init)){
    theta$SigRE = init$SigRE
  }else{
    theta$SigRE = rinvwishart(prior$RE$eta,prior$RE$Phi) #
  }


  if('SigLat' %in% names(init)){
    theta$SigLat = init$SigLat
  }else{
    theta$SigLat = rinvwishart(prior$Lat$eta,prior$Lat$Phi) #
  }
  if('betaLat' %in% names(init)){
    theta$betaLat = init$betaLat
  }else{
    theta$betaLat = rnorm(params$qLat,0,theta$sig2/sqrt(prior$FE$lambda))
  }
  if('gammaLat' %in% names(init)){
    theta$gammaLat = init$gammaLat
  }else{
    theta$gammaLat = matrix(0,nrow = params$qLat,ncol = nC)
    for (c in 1:nC){
      theta$gammaLat[,c] = rmvn(1,rep(0,params$qLat),theta$SigLat)
    }
  }


  if('SigmaClus' %in% names(init)){
    theta$SigmaClus = init$SigmaClus
  }else{
    theta$SigmaClus = array(0, dim = c(params$qU,params$qU,nC))
    for (c in 1:nC){
      theta$SigmaClus[,,c] = rinvwishart(prior$assign$nu,prior$assign$Psi)}
  }
  if('muClus' %in% names(init)){
    theta$muClus = init$muClus
  }else{
    theta$muClus = matrix(0,nrow = params$qU,ncol = nC)
    for (c in 1:nC){
      theta$muClus[,c] = rmvn(1,prior$assign$mu,theta$SigmaClus[,,c]/prior$assign$lambda)*0}
  }


  return(theta)
}
