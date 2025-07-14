library(Hmisc) 
library(rms)
library(survival)
library(glmnet)
library(ggplot2)
library(pROC)
library(caret)
library(aod)
library(mRMRe)
library(survival)
library(xlsx)
library(flexsurv)
library(compareC)
library(dplyr) 
library(survminer)
library(prodlim)
library(ResourceSelection)
library(corrplot)
library(survivalROC)
library(survIDINRI)
library(PredictABEL)
library(nricens)


# determine the optinal threshold of seperating GGO and solid component
# GGO and solid nodule threshold list
t_ggo_list <- seq(from=-800, to=-500, by=50)
t_soild_list <- seq(from=-500, to=0, by=50)

