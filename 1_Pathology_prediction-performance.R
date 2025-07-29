# evaluating the pathology grade prediction performance of model under each combination of GGO and solid thresholds
# cindex (AUC) for each combination of thresholds
# cindex and AUC is the same value under the binary classification situation
results_cindex_pathology <- c()

fix_seed <- 666

for (t_ggn in t_ggn_list) {
  for (t_soild in t_soild_list) {
    feas_file_name <- paste( t_ggn, t_soild, sep = "")
    feas_file_dir <- paste('D:/train_features/', feas_file_name, '.xlsx', sep = "")
    info_feas_train <- read.xlsx(feas_file_dir, 1, header=TRUE)    

    x_train <- info_feas_train[, 13:31]
    y_train <- info_feas_train$pathology
  
    # LASSO modelling
    set.seed(fix_seed)
    cvfit <- cv.glmnet(as.matrix(x_train),y_train,
                       nfolds=5,family = "binomial",type.measure="cindex")
    s_val <- cvfit$lambda.1se
    cindex_cv_temp <- cvfit$cvm[which(cvfit$lambda==s_val)]
    
    # print
    print(paste(feas_file_name, '5 fold mean C-index = ', round(Cindex_cv_temp,4)))

    # combine all results
    results_cindex_pathology <- rbind(results_cindex_pathology, c(t_ggn, t_soild, cindex_cv_temp))
    
  }
}

colnames(results_cindex_pathology) <- c('t_ggn', 't_solid', 'C_index_5_fold_mean')

# show the cindexs corresponding to thresholds
cindexs <- results_cindex_pathology[,3]
cindex_matrix <-matrix(results_cindex_pathology[,3], ncol = 7)
rownames(cindex_matrix) <- as.character(t_soild_list)
colnames(cindex_matrix) <- as.character(t_ggn_list)
corrplot(cindex_matrix, method='number',number.digits =4, is.corr=F, col = COL1('YlOrRd'), title = "C_index_5_fold_mean_pathology")



