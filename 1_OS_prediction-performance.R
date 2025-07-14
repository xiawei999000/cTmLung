# evaluating the OS prediction performance of model under each combination of GGO and solid thresholds
results_cindex_os <- c()

fix_seed <- 666

for (t_ggo in t_ggo_list) {
  for (t_soild in t_soild_list) {
    feas_file_name <- paste(t_ggo, t_soild, sep = "")

    feas_feas_file_dir <- paste('D:/train_features/', feas_file_name, '.xlsx', sep = "")
    info_feas_SYSUCC_train <- read.xlsx(feas_feas_file_dir, 1, header=TRUE) 

    # exclude the patient withsolid part diameter of node > 40 mm
    info_feas_SYSUCC_train <- info_feas_SYSUCC_train[-which(info_feas_SYSUCC_train$Imaging_ID %in% exclude_id_SYSUCC),]
    
    x_train <- info_feas_SYSUCC_train[, 13:31]

    # LASSO cox modelling 
    OS_time_train <- as.numeric(as.character(info_feas_SYSUCC_train$OS))
    OS_event_train <- info_feas_SYSUCC_train$OS_status
    OS_y_train <- survival::Surv(OS_time_train, OS_event_train)
    set.seed(fix_seed)
    cvfit <- cv.glmnet(as.matrix(x_train),OS_y_train,
                       nfolds=5,family = "cox",type.measure="C")
    s_val <- cvfit$lambda.min
    Cindex_cv_temp <- cvfit$cvm[which(cvfit$lambda==s_val)]

    # print
    print(paste(feas_file_name, '5 fold mean C-index = ', round(Cindex_cv_temp,4)))

    # combine all results
    results_cindex_os <- rbind(results_cindex_os, c(t_ggo, t_soild, Cindex_cv_temp))
  }
}

colnames(results_cindex_os) <- c('t_ggo', 't_solid', 'C_index_5_fold_mean')

# show the cindexs corresponding to thresholds
cindexs <- results_cindex_os[,3]
cindexs_matrix <-matrix(cindexs, ncol = 7)
rownames(cindexs_matrix) <- as.character(t_soild_list)
colnames(cindexs_matrix) <- as.character(t_ggo_list)
corrplot(cindexs_matrix, method='number',number.digits =4, is.corr=F, col = COL1('YlOrRd'), title = "C_index_5_fold_mean_OS")





