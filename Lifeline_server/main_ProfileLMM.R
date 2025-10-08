source("utility_functions.R")
library(splines)
library(LaplacesDemon)

Rcpp::sourceCpp("ProfileGibbsCPPV2.cpp")

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("At least 4 arguments expected: 1. session number (1 or 2),
                                       2. nb of individuals,
                                       3. nb max of clusters,
                                       4. nb of iterations,
                                       5. burnin (optional) default 10 percent of iterations,
                                       6. Appendix to defaul name (optional) ", call.=FALSE)
  }


FEnames = c("gender","age","age2","maritalstatus_current_adu_q_1",
            "ethnicity_category_adu_q_1", "birthplace_country_adu_q_1",
            "employment_current_adu_q_1", "degree_highest_adu_q_1",
            "current_smoker_adu_c_2","bmi") # ,"inhouse_partner_adu_q_1", "ex_smoker_adu_c_2""movehouse_followup_adu_q_2",

Latnames = c("gender","age","current_smoker_adu_c_2")#, "ever_smoker_adu_c_2")

Assignnames = c("NO2","OZO","P25" ,"NDV_MD3")#"P10",,"NDV_MD5","NDV_MD1")

#load('/groups/umcg-lifelines/tmp01/projects/ov21_0374/mamestoy/Rcode/df_merge_long_MICE.Rda')
load('df_merge_long_MICE_Dias.Rda')
df_merge_long_MICE$bpavg_diastolic_all_m_1 = as.numeric(levels(df_merge_long_MICE$bpavg_diastolic_all_m_1))[df_merge_long_MICE$bpavg_diastolic_all_m_1]

# filter missing values in outcome or exposures-----------------------
tot = (df_merge_long_MICE[[Assignnames[1]]]!=0)
for (u in Assignnames){
  tot = tot * (df_merge_long_MICE[[u]]!=0)
}
tot = tot * (1-is.na(df_merge_long_MICE[['bpavg_diastolic_all_m_1']]))

df_filt = df_merge_long_MICE[tot==1,]
print(paste0('Nb of observations with full Exposures = ',sum(tot)))
print(paste0('out of  ',dim(df_merge_long_MICE)[1]))


bsession = as.numeric(args[1])
Nobs = as.numeric(args[2])
nC = as.numeric(args[3])
nSim = as.numeric(args[4])
if (length(args)>4){nBurnIn = as.numeric(args[5])
}else{nBurnIn = int(nSim*0.1) }
if (length(args)>5){appendix = args[6]
}else{appendix = "" }
print(args)

set.seed(2)
indiv_ = df_filt$project_pseudo_id[df_filt$waveNum=='t2']
indiv_ = sample(indiv_,Nobs)

if(bsession ==2){
  print('Robustness check simulation')
  indiv_= indiv_[floor(Nobs/2):Nobs]
}else{
  print('Main simulation')
  indiv_= indiv_[1:floor(Nobs/2)]
}

idx = (df_filt$project_pseudo_id %in% indiv_)
df = df_filt[idx,]
print(dim(df))

# ,"wegverkeer" ,
#               "bc2010","UFP_RUN_1km2_Raster","stedelijkheid_adressen_per_km2", "bevolkingsdichtheid_inwoners_per_km2" ,
#               "percentage_personen_0_tot_15_jaar", "percentage_personen_15_tot_25_jaar","percentage_personen_25_tot_45_jaar",
#               "percentage_personen_45_tot_65_jaar","percentage_gehuwd","percentage_gescheid", "percentage_verweduwd" ,
#               "percentage_eenpersoonshuishoudens" ,"percentage_westerse_migratieachtergrond",
#               "dekkingspercentage","gemiddelde_woningwaarde","opleidingsniveau_hoog","opleidingsniveau_middelbaar","percentage_huishoudens_met_laag_inkomen", "percentage_huishoudens_met_hoog_inkomen",
#               "ziekenhuis_excl_buitenpolikliniek_gem_afst_in_km","treinstation_gemiddelde_afstand_in_km" ,
#               "af_healthy","af_nonhealt","af_med","af_educ","af_rec","popdensity",
#               "retailservice" ,"publictransport" )

#,"walkability" extreme values errors
df$age = as.numeric(df$age)
norma = {}
norma$age = {}
norma$age$mean = mean(df$age)
norma$age$std =sqrt(var(df$age))
df$age = (df$age-mean(df$age))/sqrt(var(df$age))
df$age2 = df$age*df$age
norma$bmi = {}
norma$bmi$mean = mean(df$bmi)
norma$bmi$std =sqrt(var(df$bmi))
df$bmi = (df$bmi-mean(df$bmi))/sqrt(var(df$bmi))
for (u in Assignnames ){
  norma[[u]] = {}
  norma[[u]]$mean = mean(df[[u]])
  norma[[u]]$std = sqrt(var(df[[u]]))
  df[[u]] = (df[[u]]-mean(df[[u]]))/sqrt(var(df[[u]]))
}
# df[,Assignnames] = scale(df[,Assignnames])

datetime_objects <- as.POSIXct(paste0(df$date, "-01 00:00:00"),
                               format = "%Y-%m-%d %H:%M:%S",
                               tz = "UTC")

min_datetime <- min(datetime_objects, na.rm = TRUE)
df$dateNum <- as.numeric(difftime(datetime_objects, min_datetime, units = "weeks"))
df$dateNum = (df$dateNum-mean(df$dateNum))/sqrt(var(df$dateNum))
Bt = bs(df$dateNum,df = 3,degree = 3)
attributes(Bt) <- attributes(Bt)["dim"]
df[c('RE1','RE2','RE3')] = Bt
rm(Bt)
REnames = c('RE1','RE2','RE3')

covList={}
covList$FE = FEnames
covList$RE = REnames
covList$Lat = Latnames
covList$Assign = Assignnames
covList$REunit = "project_pseudo_id"
covList$Y = "bpavg_diastolic_all_m_1"

dataPro = process_Data_outcome(covList, df, intercept = list(FE=T,RE=F,Lat =T))

priorInit = {}
# priorInit$DP = {}
# priorInit$DP$scale = 1
# priorInit$DP$shape = 20
# priorInit$assign = {}
# priorInit$assign$lambda = 1
# priorInit$assign$mu = rep(0,dataPro$params$qU)
# priorInit$assign$nu = dataPro$params$qU+1
# priorInit$assign$Psi =3*diag(dataPro$params$qU)

prior = prior_init(dataPro$params, nC, priorInit)
theta_Start = theta_init(prior,dataPro$params,nC,{})


outPro = profileLMM_Gibbs(dataPro,nSim+nBurnIn,nBurnIn,nC,prior,theta_Start)
outPro$names = dataPro$d$names
outPro$norma = norma
saveRDS(outPro,paste0(appendix,'.RData'))
#
# FEnames = c("gender","age","partner_presence_adu_q_1","inhouse_partner_adu_q_1",
#             "maritalstatus_current_adu_q_1", "ethnicity_category_adu_q_1", "birthplace_country_adu_q_1",
#             "employment_current_adu_q_1", "degree_highest_adu_q_1", "movehouse_followup_adu_q_2",
#             "current_smoker_adu_c_2", "ever_smoker_adu_c_2", "ex_smoker_adu_c_2")
#
# Assignnames = c("NO2cum", "NO2","OZOcum","OZO","P10cum","P10","P25cum","P25" , "NDV_MD1cum" ,
#                 "NDV_MD1","NDV_MD5cum","NDV_MD5","NDV_MD3cum","NDV_MD3","wegverkeer" ,
#                 "bc2010","UFP_RUN_1km2_Raster","stedelijkheid_adressen_per_km2", "bevolkingsdichtheid_inwoners_per_km2" ,
#                 "percentage_personen_0_tot_15_jaar", "percentage_personen_15_tot_25_jaar","percentage_personen_25_tot_45_jaar",
#                 "percentage_personen_45_tot_65_jaar","percentage_personen_65_jaar_en_ouder" ,
#                 "percentage_ongehuwd","percentage_gehuwd","percentage_gescheid", "percentage_verweduwd" ,
#                 "percentage_eenpersoonshuishoudens" ,"percentage_westerse_migratieachtergrond","percentage_niet_westerse_migratieachtergrond",
#                 "dekkingspercentage","gemiddelde_woningwaarde","opleidingsniveau_hoog","opleidingsniveau_middelbaar",
#                 "opleidingsniveau_laag","percentage_huishoudens_met_laag_inkomen", "percentage_huishoudens_met_hoog_inkomen",
#                 "ziekenhuis_excl_buitenpolikliniek_gem_afst_in_km","treinstation_gemiddelde_afstand_in_km" ,
#                 "af_healthy","af_nonhealt","af_med","af_educ","af_rec","popdensity","popdensitycum",
#                 "retailservice","retailservicecum","walkability" ,"walkabilitycum","publictransport", "publictransportcum" )
#
#
