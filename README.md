# AI-driven discovery of a simple and generalizable computer tomography imaging biomarker for improved clinical T staging in early-stage lung adenocarcinoma


# The lung nodule automatic segmentation model:
The weights of the lung nodule automatic segmentation model based on nnUNet: https://zenodo.org/records/16608410


# code structure:
> [quantitative features calculation.py/](https://github.com/xiawei999000/cTmLung/blob/main/quantitative%20features%20calculation.py): functions of quantitative features calculation based on different thresholds (tggo and ts).
>
> 
> Determination of the optimal component thresholds combination and the optimal quantitative feature to obtain the modified clinical T stage (cTm):
> >[0_import_lib_data.R/](https://github.com/xiawei999000/cTmLung/blob/main/0_import_lib_data.R): load required libs and define the thresholds ranges of extracting GGO and solid components.
> >
> >[1_DFS_prediction-performance.R/](https://github.com/xiawei999000/cTmLung/blob/main/1_DFS_prediction-performance.R): LASSO regression for DFS prediction based on 5 fold cross-validation.
> >
> >[1_OS_prediction-performance.R/](https://github.com/xiawei999000/cTmLung/blob/main/1_OS_prediction-performance.R): LASSO regression for OS prediction based on 5 fold cross-validation.  
> >
> >[1_Pathology_prediction-performance.R/](https://github.com/xiawei999000/cTmLung/blob/main/1_Pathology_prediction-performance.R): LASSO regression for pathology prediction based on 5 fold cross-validation.
> >
> >[2_Optimal_threshold_determination.R/](https://github.com/xiawei999000/cTmLung/blob/main/2_Optimal_threshold_determination.R): Determination of optimal thresholds based on the highest mean C-index of pathology, DFS, and OS prediction.
> >
> >[3_Optimal_feature_determination.R](https://github.com/xiawei999000/cTmLung/blob/main/3_Optimal_feature_determination.R): Determination of optimal quantitative feature based on the highest mean C-index of pathology, DFS, and OS prediction.
> >
> >[4_cTm.R/](https://github.com/xiawei999000/cTmLung/blob/main/4_cTm.R): The modification of current clinical T stage according to the value of SM%.
> 
> 
> [cTm_application.py/](https://github.com/xiawei999000/cTmLung/blob/main/cTm_application.py): The code for building an executable program for cTm demonstration.


# An executable program for the demonstration of cTm: 
https://zenodo.org/records/16608410


