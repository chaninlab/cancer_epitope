This folder contain R codes used for the classification modelling performed in this study. 

Four different implementations are provided: 

* oversampling for multiple classes 
* oversampling for binary class
* undersampling for multiple classes 
* undersampling for binary class

The four implementatons are provided because differences in algorithms for different classification instances were hardcoded without user selection interface as to simplify repeated code tests.

To run, only the R script needs to be started where nothing else needs to be done.

The four implementations provided accounts for all possible classification settings used in this study.

The input files to the R scripts contains less descriptors than those mentioned in Descriptor Legends. This is due to the removal of descriptors with near zero variance done as a separate step from the main classification modeling.
