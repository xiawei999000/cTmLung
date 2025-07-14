# evaluate the performance of each feature in risk level prediction, DFS and OS prediction
x_train <- info_feas_train[, 12:30]

# risk level prediction
fea_pathology_cindex_list <- c()
fea_names <- colnames(x_train)
for(fea_id in seq(from=1, to=ncol(x_train), by=1)){
  fea_temp <- x_train[,fea_id]
  fea_temp_cindex <- roc(as.factor(y_train),as.numeric(fea_temp), plot=F , print.thres=TRUE, print.auc=TRUE, smooth = F)
  fea_temp_cindex_val <- round(fea_temp_auc$auc,4)
  fea_pathology_cindex_list <- rbind(fea_pathology_cindex_list, fea_temp_cindex_val)
}
rownames(fea_pathology_cindex_list) <- colnames(x_train)

# DFS prediction
DFS_time_train <- as.numeric(as.character(info_feas_SYSUCC_train$DFS))
DFS_event_train <- info_feas_SYSUCC_train$DFS_status
DFS_y_train <- survival::Surv(DFS_time_train, DFS_event_train)
DFS_kmfit_train <- survfit(DFS_y_train~1)
# c-index
fea_DFS_cindex_list <- c()
for(fea_id in seq(from=1, to=ncol(x_train), by=1)){
  fea_temp <- x_train[,fea_id]
  cindex_tmp <- CstatisticCI(rcorr.cens(fea_temp, DFS_y_train))[1]
  
  if(cindex_tmp<0.5){
    cindex_tmp <- 1 - cindex_tmp
  }
  
  cindex_tmp <- round(cindex_tmp,4)
  
  fea_DFS_cindex_list <- rbind(fea_DFS_cindex_list, cindex_tmp)
}

# OS prediction
OS_time_train <- as.numeric(as.character(info_feas_SYSUCC_train$OS))
OS_event_train <- info_feas_SYSUCC_train$os_status
OS_y_train <- survival::Surv(OS_time_train, OS_event_train)
OS_kmfit_train <- survfit(OS_y_train~1)
# c-index
fea_OS_cindex_list <- c()
for(fea_id in seq(from=1, to=ncol(x_train), by=1)){
  fea_temp <- x_train[,fea_id]
  cindex_tmp <- CstatisticCI(rcorr.cens(fea_temp, OS_y_train))[1]
  
  if(cindex_tmp<0.5){
    cindex_tmp <- 1 - cindex_tmp
  }
  cindex_tmp <- round(cindex_tmp,4)
  fea_OS_cindex_list <- rbind(fea_OS_cindex_list, cindex_tmp)
}

# sum of c-index
mean_cindex_list <- round((fea_pathology_cindex_list+fea_DFS_cindex_list+fea_OS_cindex_list)/3, 4)
train_fea_cindex_list <- cbind(fea_AUC_list, fea_DFS_cindex_list, fea_OS_cindex_list, mean_cindex_list)
train_fea_cindex_list <- as.data.frame(train_fea_cindex_list)
colnames(train_fea_cindex_list) <- c('pathology_cindex', 'DFS_cindex', 'OS_cindex', 'mean_cindex')

# the best feature
rownames(train_fea_cindex_list)[which(train_fea_cindex_list$mean_cindex==max(train_fea_cindex_list$mean_cindex))]
