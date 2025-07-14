
# train
# obtain cT (9th version) according to solid diameter
train_cT <- train_mannual_soild_d
train_cT[which(train_mannual_soild_d<=10)] <- 1
train_cT[which(train_mannual_soild_d<=20  &  train_mannual_soild_d>10)] <- 2
train_cT[which(train_mannual_soild_d<=30 &  train_mannual_soild_d>20)] <- 3
train_cT[which(train_mannual_soild_d>30)] <- 4

cT <- train_cT
DFS_time <- DFS_time_train
DFS_event <- DFS_event_train
OS_time <- OS_time_train
OS_event <- OS_event_train
info_feas_train_combine <- cbind(info_feas_train, cT, DFS_time, DFS_event, OS_time, OS_event)

# adjust the T stage using the optimal feature by definition of two thresholds of upstage and downstage
# optimal thresholds when maximum the C-index
# downstage_th_list<- seq(from=0.01, to=0.49, by=0.01)
# upstage_th_list<- seq(from=0.50, to=0.99, by=0.01)
Cindex_adjust_results <- c()
group_cutoff_list <- seq(from=0, to=1, by=0.1)
for(group_cutoff in group_cutoff_list){

downstage_th_list<- seq(from=0, to=group_cutoff, by=0.05)
upstage_th_list<- seq(from=group_cutoff, to=1, by=0.05)

for (downstage_th in downstage_th_list){
  for (upstage_th in upstage_th_list){
    
    cT_adjust <- train_cT
    cT_adjust[which(train_cT==4 & solid_mass_percentage_train<downstage_th)] <- 3
    cT_adjust[which(train_cT==3 & solid_mass_percentage_train>upstage_th)] <- 4
    cT_adjust[which(train_cT==3 & solid_mass_percentage_train<downstage_th)] <- 2
    cT_adjust[which(train_cT==2 & solid_mass_percentage_train>upstage_th)] <- 3
    cT_adjust[which(train_cT==2 & solid_mass_percentage_train<downstage_th)] <- 1
    cT_adjust[which(train_cT==1 & solid_mass_percentage_train>upstage_th)] <- 2
    
    Cindex_pathology_adjust_tmp <- roc(as.factor(y_train),as.numeric(cT_adjust), plot=F, 
                          print.thres=F, print.Cindex=F, smooth = F)
    Cindex_pathology_adjust_val_tmp <- round(Cindex_pathology_adjust_tmp$Cindex,4)
    
    Cindex_DFS_adjust_tmp <-  round(CstatisticCI(rcorr.cens(cT_adjust, DFS_y_train))[1],4)
    
    Cindex_OS_adjust_tmp <-  round(CstatisticCI(rcorr.cens(cT_adjust, OS_y_train))[1],4)
    
    Cindex_mean_adjust_tmp <- round((Cindex_pathology_adjust_val_tmp+Cindex_DFS_adjust_tmp+Cindex_OS_adjust_tmp)/3,4)
    
    
    Cindex_adjust_results <- rbind(Cindex_adjust_results, c(downstage_th, upstage_th, 
                                                            Cindex_pathology_adjust_val_tmp, Cindex_DFS_adjust_tmp, Cindex_OS_adjust_tmp, Cindex_mean_adjust_tmp))
    print(paste('downstage_th = ', downstage_th, '; upstage_th = ', upstage_th,
                '; C_index_pathology = ', Cindex_pathology_adjust_val_tmp,
                '; C_index_DFS = ', Cindex_DFS_adjust_tmp,
                '; C_index_OS = ', Cindex_DFS_adjust_tmp,
                '; C_index_mean = ', Cindex_mean_adjust_tmp))
  }
}
}

Cindex_adjust_results <- as.data.frame(Cindex_adjust_results)
colnames(Cindex_adjust_results) <- c('downstage_th', 'upstage_th', 'C_index_pathology', 'C_index_DFS', 'C_index_OS', 'C_index_mean')

# print thresholds according to the max mean C-index
Cindex_adjust_results[which(Cindex_adjust_results$C_index_mean==max(Cindex_adjust_results$C_index_mean)),]


###train set###
# application of the best adjustment thresholds to modify current cT
cT_adjust_train <- train_cT
downstage_th <- 0.45
upstage_th <- 0.75
cT_adjust_train[which(train_cT==4 & solid_mass_percentage_train<downstage_th)] <- 3
cT_adjust_train[which(train_cT==3 & solid_mass_percentage_train>upstage_th)] <- 4
cT_adjust_train[which(train_cT==3 & solid_mass_percentage_train<downstage_th)] <- 2
cT_adjust_train[which(train_cT==2 & solid_mass_percentage_train>upstage_th)] <- 3
cT_adjust_train[which(train_cT==2 & solid_mass_percentage_train<downstage_th)] <- 1
cT_adjust_train[which(train_cT==1 & solid_mass_percentage_train>upstage_th)] <- 2
cT_adjust<- cT_adjust_train
info_feas_train_combine_adjust <- cbind(info_feas_train_combine, cT_adjust)


# compare pathological risk performance
y_pathology_train <- info_feas_train$pathology
Cindex_train <- roc(as.factor(y_pathology_train),as.numeric(train_cT), plot=T, 
                        print.thres=T, print.Cindex=T, smooth = F)
Cindex_adjust_train <- roc(as.factor(y_pathology_train),as.numeric(cT_adjust_train), plot=T, 
                               print.thres=T, print.Cindex=T, smooth = F)
ci.Cindex(Cindex_adjust_train)
roc.test(Cindex_adjust_train, Cindex_train, alternative ='greater')
# IDI # NRI
Cindex_train_data <- cbind(y_pathology_train, train_cT, cT_adjust_train)
mstd_train <- glm(y_pathology_train~.,binomial(logit),data.frame(y_pathology_train, train_cT), x=T) 
mnew_train <- glm(y_pathology_train~.,binomial(logit),data.frame(y_pathology_train, cT_adjust_train), x=T) 
reclassification(data=Cindex_train_data, cOutcome = 1, 
                 predrisk1 = mstd_train$fitted.values, predrisk2 = mnew_train$fitted.values,
                 cutoff=c(0,0.2,0.4,1))


# compare c-index of DFS
CstatisticCI(rcorr.cens(train_cT, DFS_y_train))
CstatisticCI(rcorr.cens(cT_adjust_train, DFS_y_train))
compareC(DFS_time_train, DFS_event_train, cT_adjust_train, train_cT)
# IDI # NRI
DFS_train_data <- cbind(DFS_time_train, DFS_event_train)
DFS_train_data <- as.data.frame(DFS_train_data)
t0 = 365*5
IDI_train <- IDI.INF(DFS_train_data, train_cT, cT_adjust_train, t0, npert = 100 )
IDI.INF.OUT(IDI_train)
IDI.INF.GRAPH(IDI_train)

# KM plot of original cT for DFS
# ggplot
DFS_kmfit_train_cT <- survfit(DFS_y_train~train_cT)
names(DFS_kmfit_train_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_train_cT, data=DFS_y_train,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)

# KM plot of new cT for DFS
DFS_kmfit_train_cT <- survfit(Surv(DFS_time, DFS_event)~cT_adjust, data = info_feas_train_combine_adjust)
names(DFS_kmfit_train_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_train_cT, data=info_feas_train_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# compare c-index of OS
CstatisticCI(rcorr.cens(train_cT, OS_y_train))
CstatisticCI(rcorr.cens(cT_adjust_train, OS_y_train))
compareC(OS_time_train, OS_event_train, cT_adjust_train, train_cT)

# IDI # NRI
OS_train_data <- cbind(OS_time_train, OS_event_train)
OS_train_data <- as.data.frame(OS_train_data)
t0 = 365*5
IDI_train <- IDI.INF(OS_train_data, train_cT, cT_adjust_train, t0, npert = 100 )
IDI.INF.OUT(IDI_train)
IDI.INF.GRAPH(IDI_train)

#cT original
OS_kmfit_train_cT <- survfit(Surv(OS_time, OS_event)~cT, data = info_feas_train_combine_adjust)
names(OS_kmfit_train_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_train_cT, data=info_feas_train_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# cT adjust
OS_kmfit_train_cT <- survfit(Surv(OS_time, OS_event)~cT_adjust, data = info_feas_train_combine_adjust)
names(OS_kmfit_train_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_train_cT, data=info_feas_train_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)



###test1 set###
# obtain cT (9th version) according to solid diameter
test1_cT <- test1_mannual_soild_d
test1_cT[which(test1_mannual_soild_d<=10)] <- 1
test1_cT[which(test1_mannual_soild_d<=20  &  test1_mannual_soild_d>10)] <- 2
test1_cT[which(test1_mannual_soild_d<=30 &  test1_mannual_soild_d>20)] <- 3
test1_cT[which(test1_mannual_soild_d>30)] <- 4

cT <- test1_cT
DFS_time <- DFS_time_test1
DFS_event <- DFS_event_test1
OS_time <- OS_time_test1
OS_event <- OS_event_test1
info_feas_test1_combine <- cbind(info_feas_test1, cT, DFS_time, DFS_event, OS_time, OS_event)

# application of the best adjustment thresholds to modify current cT
cT_adjust_test1 <- test1_cT
downstage_th <- 0.45
upstage_th <- 0.75
cT_adjust_test1[which(test1_cT==4 & solid_mass_percentage_test1<downstage_th)] <- 3
cT_adjust_test1[which(test1_cT==3 & solid_mass_percentage_test1>upstage_th)] <- 4
cT_adjust_test1[which(test1_cT==3 & solid_mass_percentage_test1<downstage_th)] <- 2
cT_adjust_test1[which(test1_cT==2 & solid_mass_percentage_test1>upstage_th)] <- 3
cT_adjust_test1[which(test1_cT==2 & solid_mass_percentage_test1<downstage_th)] <- 1
cT_adjust_test1[which(test1_cT==1 & solid_mass_percentage_test1>upstage_th)] <- 2
cT_adjust<- cT_adjust_test1
info_feas_test1_combine_adjust <- cbind(info_feas_test1_combine, cT_adjust)


# compare pathological risk performance
y_pathology_test1 <- info_feas_test1$pathology
Cindex_test1 <- roc(as.factor(y_pathology_test1),as.numeric(test1_cT), plot=T, 
                    print.thres=T, print.Cindex=T, smooth = F)
Cindex_adjust_test1 <- roc(as.factor(y_pathology_test1),as.numeric(cT_adjust_test1), plot=T, 
                           print.thres=T, print.Cindex=T, smooth = F)
ci.Cindex(Cindex_adjust_test1)
roc.test(Cindex_adjust_test1, Cindex_test1, alternative ='greater')
# IDI # NRI
Cindex_test1_data <- cbind(y_pathology_test1, test1_cT, cT_adjust_test1)
mstd_test1 <- glm(y_pathology_test1~.,binomial(logit),data.frame(y_pathology_test1, test1_cT), x=T) 
mnew_test1 <- glm(y_pathology_test1~.,binomial(logit),data.frame(y_pathology_test1, cT_adjust_test1), x=T) 
reclassification(data=Cindex_test1_data, cOutcome = 1, 
                 predrisk1 = mstd_test1$fitted.values, predrisk2 = mnew_test1$fitted.values,
                 cutoff=c(0,0.2,0.4,1))


# compare c-index of DFS
CstatisticCI(rcorr.cens(test1_cT, DFS_y_test1))
CstatisticCI(rcorr.cens(cT_adjust_test1, DFS_y_test1))
compareC(DFS_time_test1, DFS_event_test1, cT_adjust_test1, test1_cT)
# IDI # NRI
DFS_test1_data <- cbind(DFS_time_test1, DFS_event_test1)
DFS_test1_data <- as.data.frame(DFS_test1_data)
t0 = 365*5
IDI_test1 <- IDI.INF(DFS_test1_data, test1_cT, cT_adjust_test1, t0, npert = 100 )
IDI.INF.OUT(IDI_test1)
IDI.INF.GRAPH(IDI_test1)

# KM plot of original cT for DFS
# ggplot
DFS_kmfit_test1_cT <- survfit(DFS_y_test1~test1_cT)
names(DFS_kmfit_test1_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_test1_cT, data=DFS_y_test1,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)

# KM plot of new cT for DFS
DFS_kmfit_test1_cT <- survfit(Surv(DFS_time, DFS_event)~cT_adjust, data = info_feas_test1_combine_adjust)
names(DFS_kmfit_test1_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_test1_cT, data=info_feas_test1_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# compare c-index of OS
CstatisticCI(rcorr.cens(test1_cT, OS_y_test1))
CstatisticCI(rcorr.cens(cT_adjust_test1, OS_y_test1))
compareC(OS_time_test1, OS_event_test1, cT_adjust_test1, test1_cT)

# IDI # NRI
OS_test1_data <- cbind(OS_time_test1, OS_event_test1)
OS_test1_data <- as.data.frame(OS_test1_data)
t0 = 365*5
IDI_test1 <- IDI.INF(OS_test1_data, test1_cT, cT_adjust_test1, t0, npert = 100 )
IDI.INF.OUT(IDI_test1)
IDI.INF.GRAPH(IDI_test1)

#cT original
OS_kmfit_test1_cT <- survfit(Surv(OS_time, OS_event)~cT, data = info_feas_test1_combine_adjust)
names(OS_kmfit_test1_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_test1_cT, data=info_feas_test1_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# cT adjust
OS_kmfit_test1_cT <- survfit(Surv(OS_time, OS_event)~cT_adjust, data = info_feas_test1_combine_adjust)
names(OS_kmfit_test1_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_test1_cT, data=info_feas_test1_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


###test2 set###
# obtain cT (9th version) according to solid diameter
test2_cT <- test2_mannual_soild_d
test2_cT[which(test2_mannual_soild_d<=10)] <- 1
test2_cT[which(test2_mannual_soild_d<=20  &  test2_mannual_soild_d>10)] <- 2
test2_cT[which(test2_mannual_soild_d<=30 &  test2_mannual_soild_d>20)] <- 3
test2_cT[which(test2_mannual_soild_d>30)] <- 4

cT <- test2_cT
DFS_time <- DFS_time_test2
DFS_event <- DFS_event_test2
OS_time <- OS_time_test2
OS_event <- OS_event_test2
info_feas_test2_combine <- cbind(info_feas_test2, cT, DFS_time, DFS_event, OS_time, OS_event)

# application of the best adjustment thresholds to modify current cT
cT_adjust_test2 <- test2_cT
downstage_th <- 0.45
upstage_th <- 0.75
cT_adjust_test2[which(test2_cT==4 & solid_mass_percentage_test2<downstage_th)] <- 3
cT_adjust_test2[which(test2_cT==3 & solid_mass_percentage_test2>upstage_th)] <- 4
cT_adjust_test2[which(test2_cT==3 & solid_mass_percentage_test2<downstage_th)] <- 2
cT_adjust_test2[which(test2_cT==2 & solid_mass_percentage_test2>upstage_th)] <- 3
cT_adjust_test2[which(test2_cT==2 & solid_mass_percentage_test2<downstage_th)] <- 1
cT_adjust_test2[which(test2_cT==1 & solid_mass_percentage_test2>upstage_th)] <- 2
cT_adjust<- cT_adjust_test2
info_feas_test2_combine_adjust <- cbind(info_feas_test2_combine, cT_adjust)


# compare pathological risk performance
y_pathology_test2 <- info_feas_test2$pathology
Cindex_test2 <- roc(as.factor(y_pathology_test2),as.numeric(test2_cT), plot=T, 
                    print.thres=T, print.Cindex=T, smooth = F)
Cindex_adjust_test2 <- roc(as.factor(y_pathology_test2),as.numeric(cT_adjust_test2), plot=T, 
                           print.thres=T, print.Cindex=T, smooth = F)
ci.Cindex(Cindex_adjust_test2)
roc.test(Cindex_adjust_test2, Cindex_test2, alternative ='greater')
# IDI # NRI
Cindex_test2_data <- cbind(y_pathology_test2, test2_cT, cT_adjust_test2)
mstd_test2 <- glm(y_pathology_test2~.,binomial(logit),data.frame(y_pathology_test2, test2_cT), x=T) 
mnew_test2 <- glm(y_pathology_test2~.,binomial(logit),data.frame(y_pathology_test2, cT_adjust_test2), x=T) 
reclassification(data=Cindex_test2_data, cOutcome = 1, 
                 predrisk1 = mstd_test2$fitted.values, predrisk2 = mnew_test2$fitted.values,
                 cutoff=c(0,0.2,0.4,1))


# compare c-index of DFS
CstatisticCI(rcorr.cens(test2_cT, DFS_y_test2))
CstatisticCI(rcorr.cens(cT_adjust_test2, DFS_y_test2))
compareC(DFS_time_test2, DFS_event_test2, cT_adjust_test2, test2_cT)
# IDI # NRI
DFS_test2_data <- cbind(DFS_time_test2, DFS_event_test2)
DFS_test2_data <- as.data.frame(DFS_test2_data)
t0 = 365*5
IDI_test2 <- IDI.INF(DFS_test2_data, test2_cT, cT_adjust_test2, t0, npert = 100 )
IDI.INF.OUT(IDI_test2)
IDI.INF.GRAPH(IDI_test2)

# KM plot of original cT for DFS
# ggplot
DFS_kmfit_test2_cT <- survfit(DFS_y_test2~test2_cT)
names(DFS_kmfit_test2_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_test2_cT, data=DFS_y_test2,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)

# KM plot of new cT for DFS
DFS_kmfit_test2_cT <- survfit(Surv(DFS_time, DFS_event)~cT_adjust, data = info_feas_test2_combine_adjust)
names(DFS_kmfit_test2_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_test2_cT, data=info_feas_test2_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# compare c-index of OS
CstatisticCI(rcorr.cens(test2_cT, OS_y_test2))
CstatisticCI(rcorr.cens(cT_adjust_test2, OS_y_test2))
compareC(OS_time_test2, OS_event_test2, cT_adjust_test2, test2_cT)

# IDI # NRI
OS_test2_data <- cbind(OS_time_test2, OS_event_test2)
OS_test2_data <- as.data.frame(OS_test2_data)
t0 = 365*5
IDI_test2 <- IDI.INF(OS_test2_data, test2_cT, cT_adjust_test2, t0, npert = 100 )
IDI.INF.OUT(IDI_test2)
IDI.INF.GRAPH(IDI_test2)

#cT original
OS_kmfit_test2_cT <- survfit(Surv(OS_time, OS_event)~cT, data = info_feas_test2_combine_adjust)
names(OS_kmfit_test2_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_test2_cT, data=info_feas_test2_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# cT adjust
OS_kmfit_test2_cT <- survfit(Surv(OS_time, OS_event)~cT_adjust, data = info_feas_test2_combine_adjust)
names(OS_kmfit_test2_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_test2_cT, data=info_feas_test2_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)








###TCIA set###
# obtain cT (9th version) according to solid diameter
TCIA_cT <- TCIA_mannual_soild_d
TCIA_cT[which(TCIA_mannual_soild_d<=10)] <- 1
TCIA_cT[which(TCIA_mannual_soild_d<=20  &  TCIA_mannual_soild_d>10)] <- 2
TCIA_cT[which(TCIA_mannual_soild_d<=30 &  TCIA_mannual_soild_d>20)] <- 3
TCIA_cT[which(TCIA_mannual_soild_d>30)] <- 4

cT <- TCIA_cT
DFS_time <- DFS_time_TCIA
DFS_event <- DFS_event_TCIA
OS_time <- OS_time_TCIA
OS_event <- OS_event_TCIA
info_feas_TCIA_combine <- cbind(info_feas_TCIA, cT, DFS_time, DFS_event, OS_time, OS_event)

# application of the best adjustment thresholds to modify current cT
cT_adjust_TCIA <- TCIA_cT
downstage_th <- 0.45
upstage_th <- 0.75
cT_adjust_TCIA[which(TCIA_cT==4 & solid_mass_percentage_TCIA<downstage_th)] <- 3
cT_adjust_TCIA[which(TCIA_cT==3 & solid_mass_percentage_TCIA>upstage_th)] <- 4
cT_adjust_TCIA[which(TCIA_cT==3 & solid_mass_percentage_TCIA<downstage_th)] <- 2
cT_adjust_TCIA[which(TCIA_cT==2 & solid_mass_percentage_TCIA>upstage_th)] <- 3
cT_adjust_TCIA[which(TCIA_cT==2 & solid_mass_percentage_TCIA<downstage_th)] <- 1
cT_adjust_TCIA[which(TCIA_cT==1 & solid_mass_percentage_TCIA>upstage_th)] <- 2
cT_adjust<- cT_adjust_TCIA
info_feas_TCIA_combine_adjust <- cbind(info_feas_TCIA_combine, cT_adjust)


# compare pathological risk performance
y_pathology_TCIA <- info_feas_TCIA$pathology
Cindex_TCIA <- roc(as.factor(y_pathology_TCIA),as.numeric(TCIA_cT), plot=T, 
                   print.thres=T, print.Cindex=T, smooth = F)
Cindex_adjust_TCIA <- roc(as.factor(y_pathology_TCIA),as.numeric(cT_adjust_TCIA), plot=T, 
                          print.thres=T, print.Cindex=T, smooth = F)
ci.Cindex(Cindex_adjust_TCIA)
roc.test(Cindex_adjust_TCIA, Cindex_TCIA, alternative ='greater')
# IDI # NRI
Cindex_TCIA_data <- cbind(y_pathology_TCIA, TCIA_cT, cT_adjust_TCIA)
mstd_TCIA <- glm(y_pathology_TCIA~.,binomial(logit),data.frame(y_pathology_TCIA, TCIA_cT), x=T) 
mnew_TCIA <- glm(y_pathology_TCIA~.,binomial(logit),data.frame(y_pathology_TCIA, cT_adjust_TCIA), x=T) 
reclassification(data=Cindex_TCIA_data, cOutcome = 1, 
                 predrisk1 = mstd_TCIA$fitted.values, predrisk2 = mnew_TCIA$fitted.values,
                 cutoff=c(0,0.2,0.4,1))


# compare c-index of DFS
CstatisticCI(rcorr.cens(TCIA_cT, DFS_y_TCIA))
CstatisticCI(rcorr.cens(cT_adjust_TCIA, DFS_y_TCIA))
compareC(DFS_time_TCIA, DFS_event_TCIA, cT_adjust_TCIA, TCIA_cT)
# IDI # NRI
DFS_TCIA_data <- cbind(DFS_time_TCIA, DFS_event_TCIA)
DFS_TCIA_data <- as.data.frame(DFS_TCIA_data)
t0 = 365*5
IDI_TCIA <- IDI.INF(DFS_TCIA_data, TCIA_cT, cT_adjust_TCIA, t0, npert = 100 )
IDI.INF.OUT(IDI_TCIA)
IDI.INF.GRAPH(IDI_TCIA)

# KM plot of original cT for DFS
# ggplot
DFS_kmfit_TCIA_cT <- survfit(DFS_y_TCIA~TCIA_cT)
names(DFS_kmfit_TCIA_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_TCIA_cT, data=DFS_y_TCIA,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)

# KM plot of new cT for DFS
DFS_kmfit_TCIA_cT <- survfit(Surv(DFS_time, DFS_event)~cT_adjust, data = info_feas_TCIA_combine_adjust)
names(DFS_kmfit_TCIA_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(DFS_kmfit_TCIA_cT, data=info_feas_TCIA_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Recurrence-free survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# compare c-index of OS
CstatisticCI(rcorr.cens(TCIA_cT, OS_y_TCIA))
CstatisticCI(rcorr.cens(cT_adjust_TCIA, OS_y_TCIA))
compareC(OS_time_TCIA, OS_event_TCIA, cT_adjust_TCIA, TCIA_cT)

# IDI # NRI
OS_TCIA_data <- cbind(OS_time_TCIA, OS_event_TCIA)
OS_TCIA_data <- as.data.frame(OS_TCIA_data)
t0 = 365*5
IDI_TCIA <- IDI.INF(OS_TCIA_data, TCIA_cT, cT_adjust_TCIA, t0, npert = 100 )
IDI.INF.OUT(IDI_TCIA)
IDI.INF.GRAPH(IDI_TCIA)

#cT original
OS_kmfit_TCIA_cT <- survfit(Surv(OS_time, OS_event)~cT, data = info_feas_TCIA_combine_adjust)
names(OS_kmfit_TCIA_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_TCIA_cT, data=info_feas_TCIA_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)


# cT adjust
OS_kmfit_TCIA_cT <- survfit(Surv(OS_time, OS_event)~cT_adjust, data = info_feas_TCIA_combine_adjust)
names(OS_kmfit_TCIA_cT$strata) <- c('cT1a', 'cT1b', 'cT1c', 'cT2a')
ggsurvplot(OS_kmfit_TCIA_cT, data=info_feas_TCIA_combine_adjust,
           risk.table = TRUE, # Add risk table
           risk.table.y.text.col = F,
           fontsize =  5,
           tables.theme = theme_cleantable(),
           pval = TRUE,
           #conf.int = TRUE,
           font.legend = c(16),
           # font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(16),
           xlab = "Time (days)",
           ylab = "Overall survival rate",
           # risk.table.col = "strata", # Change risk table color by group
           linetype = "strata", # Change line type by groups
)
