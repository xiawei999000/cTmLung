# sum and calculate the mean cindex to find the optimal threshold
sum_cindex_train_list <- results_cindex_os$C_index_5_fold_mean + results_cindex_DFS$C_index_5_fold_mean + results_cindex_pathology$C_index_5_fold_mean

# show the cindex corresponding to thresholds
mean_matrix <-matrix(sum_cindex_train_list/3, ncol = 7)
rownames(mean_matrix) <- as.character(t_soild_list)
colnames(mean_matrix) <- as.character(t_ggn_list)
corrplot(mean_matrix, method='number',number.digits =4, is.corr=F, col = COL1('YlOrRd'))

# optimal threshold with the max mean cindex
optimal_results <- results_bind[which(mean_matrix==max(mean_matrix)),]

t_ggn <- optimal_results[1]
t_soild <- optimal_results[2]

feas_optimal_threshold <- paste( t_ggn, t_soild, sep = "")
# optimal threshold
feas_optimal_threshold



# load the features under the optimal threshold
# train
feas_file_dir <- paste('D:/train_features/', feas_optimal_threshold, '.xlsx', sep = "")
info_feas_train <- read.xlsx(feas_file_dir, 1, header=TRUE)

# test1
feas_file_dir_test1 <- paste('D:/test1_features/', feas_optimal_threshold, '.xlsx', sep = "")
info_feas_test1 <- read.xlsx(feas_file_dir_test1, 1, header=TRUE)

# test2
feas_file_dir_test2 <- paste('D:/test2_features/', feas_optimal_threshold, '.xlsx', sep = "")
info_feas_test2 <- read.xlsx(feas_file_dir_test2, 1, header=TRUE)

# TCIA
feas_file_dir_TCIA <- paste('D:/TCIA_features/', feas_optimal_threshold, '.xlsx', sep = "")
info_feas_TCIA <- read.xlsx(feas_file_dir_TCIA, 1, header=TRUE)
