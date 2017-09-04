# R codes

## Classification modeling
This folder contain R codes used for the classification modelling performed in this study. 

Four different implementations are provided: 

* oversampling for multiple classes (*oversampling_all_classes*)
* oversampling for binary class (*oversampling_dual_vs_mono_epitopes*)
* undersampling for multiple classes (*undersampling_all_classes*)
* undersampling for binary class (*undersampling_dual_vs_mono_epitopes*)

The four implementatons are provided because differences in algorithms for different classification instances were hardcoded without user selection interface as to simplify repeated code tests.

To run, only the R script needs to be started where nothing else needs to be done.

The four implementations provided accounts for all possible classification settings used in this study.

The input files to the R scripts contains less descriptors than those mentioned in Descriptor Legends. This is due to the removal of descriptors with near zero variance done as a separate step from the main classification modeling.

## Bar plots on Feature importance

The folder *Figure_3_bar_plot_feature_importance* contains the R codes, input files and output files for re-creating Figure 3 which shows the Gini index of significant descriptors derived from the random forest model.
