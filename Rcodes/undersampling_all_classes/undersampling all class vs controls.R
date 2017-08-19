library(readxl)
library(caret)
require(sampling)
library(ranger)
#require(data.table)

setwd("d:/R_test/undersampling_all_classes")

#Change this to the correct input file name
df <- read.csv("full_epitopes_rm_zero_var.csv")

df<-na.omit(df)
Class <- df$Class
Class <- as.factor(Class)
descriptors <- df[, 2:ncol(df)]
whole_data_set <- cbind(Class,descriptors)



#a vector to store random seed for each of the repeat runs of a 10 fold CV
seed_vec<-c()
#IMPORTANT! the number of seeds in there also determines how many repeats are made for a modelling run
seed_vec<-c(seed_vec,650)  #,495,123,789,754)

#storage for the results, NOTE that each of the test type needs its own storage frame or overwritting will occur
instance_results_recall<-data.frame()
kappa_results_recall<-data.frame()
class_results_recall<-data.frame()

instance_results_CV<-data.frame()
kappa_results_CV<-data.frame()
class_results_CV<-data.frame()

instance_results_test<-data.frame()
kappa_results_test<-data.frame()
class_results_test<-data.frame()


#lists used to store variable importance calculation from each repeat this section of data frames is used to store results for grand average claculation
variable_importance_list <- list(length(seed_vec))



#outer LOOP for repeated runs 
for(REPEATS in 1:length(seed_vec))
{
 set.seed(seed_vec[REPEATS])
 #IMPORTANT!! the data set needs to be reassigned for each repeat, otherwise the create datapartitioning will create smaller
 #and smaller set with every repeat, first repeat its 0.8 of the initial dataset, 2nd repeat its 0.8 of the first repeat
 df<-whole_data_set
  
  
  
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing START *-*-*-*-*-*-*-*-*-*-  
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing START *-*-*-*-*-*-*-*-*-*-   
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing START *-*-*-*-*-*-*-*-*-*-  
class_level<-levels(Class)
class_list<-list()

#before K-means clustering can commence, one needs to determine the number of clusters that needs to be created
#based on the smallest cluster ther is. For ease of implementation, do the class counting twice
max_class=-1
min_class=1e5
for(i in 1:length(class_level))
{
  indiv_class<-class_level[i]
  temp_class<-subset(df,Class==indiv_class)
  class_length<-nrow(temp_class)
  if(class_length<min_class)
    min_class<-class_length
  if(class_length>max_class)
    max_class<-class_length
}
#print("min and max class")
#print(min_class)
#print(max_class)  
  
#this is to store all compressed classes in a single frame, like 
balanced_class<-data.frame()

for(i in 1:length(class_level))
{
  indiv_class<-class_level[i]
  temp_class<-subset(df,Class==indiv_class)

  temp_class_descriptors<-temp_class[,2:ncol(temp_class)]
  #print("the descriptors of this class")
  #print(temp_class_descriptors)
  
  if(nrow(temp_class)==min_class)
  {
    balanced_class<-rbind(balanced_class,temp_class)
    next
  }
  
  cl<-kmeans(temp_class_descriptors,min_class)
  centroids<-cl$centers
  
 
  #print("the centroids of this cluster")
  #print(centroids)
  
  centroids<-data.frame(centroids)
  #append the class name
  temp_append_class<-subset(Class,Class==indiv_class)
  temp_append_class<-data.frame(temp_append_class)
  temp_append_class<-temp_append_class[1:nrow(centroids),]

  centroids<-cbind(temp_append_class,centroids)
  names(centroids)[1]<-"Class"
 
  #print("the formatted compressed class")
  #print(centroids)
  
  #line<-readline()
  balanced_class<-rbind(balanced_class,centroids)
}


#print("the finalized, balanced classes")
#print(balanced_class)  
  
df<-balanced_class  
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing END *-*-*-*-*-*-*-*-*-*-  
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing END *-*-*-*-*-*-*-*-*-*-   
#*-*-*-*-*-*-*-*-*-*-*-*-*- code for undersampling class balancing END *-*-*-*-*-*-*-*-*-*-  

  


#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  

  #IMPORTANT!!! this model training on the training set is to be used for the test set as well, do not overwrite it
  model <- ranger::ranger(Class ~., data = df, write.forest  = TRUE, save.memory = TRUE, importance = "impurity",num.trees=500)
  actual <- df$Class
 # prediction <- ranger::predictions(model, df)
  prediction <- predict(model, df)
  prediction <- prediction$predictions
  
  #print("recall test actual")
  #print(actual)
  #line<-readline()
  
  #print("predicted")
  #print(prediction)
  #line<-readline()
  
  results <- caret::confusionMatrix(prediction, actual)
  #print("recall result")
  #print(results)
  #line<-readline()
  
  #code for saving the results_instances -*-*-*-*-*-*-*-*-*-*-*-*-*
  instance<-data.frame(as.matrix(results))
  instance_results_recall<-rbind(instance_results_recall,instance)
  #add a blank row behind the results of each REPAT as spacer
  #first create a one row matrix of the same length as the existing data frame, instance_results_recall in this case
  temprow<-matrix(c(rep.int(NA,length(instance_results_recall))),nrow=1,ncol=length(instance_results_recall))
  #make this one row matrix a dataframe and give its cols the same names as instance_results_recall
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(instance_results_recall)
  #bind the empty row to the existing dataframe
  instance_results_recall<-rbind(instance_results_recall,newrow)
  
  #code for saving the kappa results -*-*-*-*-*-*-* start -*-*-*-*-*-*
  kappa_res<-data.frame(as.matrix(results,what="overall"))
  kappa_results_recall<-rbind(kappa_results_recall,kappa_res)
  temprow<-matrix(c(rep.int(NA,length(kappa_results_recall))),nrow=1,ncol=length(kappa_results_recall))
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(kappa_results_recall)
  kappa_results_recall<-rbind(kappa_results_recall,newrow)
  
  #code for saving the class specificity results -*-*-*-*-*-*-* start -*-*-*-*-*-*
  class_res<-data.frame(as.matrix(results,what="classes"))
  class_results_recall<-rbind(class_results_recall,class_res)
  temprow<-matrix(c(rep.int(NA,length(class_results_recall))),nrow=1,ncol=length(class_results_recall))
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(class_results_recall)
  class_results_recall<-rbind(class_results_recall,newrow)
  
  write.csv(as.matrix(instance_results_recall), file=paste0("Recall_results_instances.csv"))
  write.csv(as.matrix(kappa_results_recall), file=paste0("Recall_results_kappa.csv"))	
  write.csv(as.matrix(class_results_recall), file=paste0("Recall_results_class.csv"))
  
  
  #*-*-*--**-*- code for storing variable importance 
  variable_importance<-model$variable.importance
  variable_importance_list[[REPEATS]]<-variable_importance
  
  
#-*-*-*-*-*-*-*-*-* code for model recall performance END -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  	
#-*-*-*-*-*-*-*-*-* code for model recall performance END -*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
#-*-*-*-*-*-*-*-*-* code for model recall performance END -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
  

#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  folds<-createFolds(df$Class,min_class,list=TRUE,returnTrain=FALSE)
	totalCV_result<-data.frame()
	
	print("split into how many folds?")
	print(min_class)
	for (i in 1:length(folds))
	{
	  print("fold no.")
	  print(i)
	  avector<-data.frame(folds[i])
	  avector<-avector[,1]     
	  test<-df[avector,]
	  train<-df[-avector,]
	  
	  model <- ranger::ranger(Class ~., data = train, write.forest  = TRUE, save.memory = TRUE, importance = "impurity", num.trees = 500)
	  # this uses the model build by the above line to predict the test set  
	  prediction <- predict(model, test)
	  #this line extracts the predicted result for each sample of the test set
	  prediction <- prediction$predictions
	  
	  
	  # extracting the actual target variable value for each sample of the test set to compare against predicted value
	  actual <- test$Class
	  # binding a 2 column list of predicted vs actual Y-value   
	  data <- data.frame(prediction, actual)

	  # storing the predicte vs actual results for each of the folds in a very long df, appended in row fashion
	  totalCV_result <- rbind(totalCV_result,data)
	}
	
	prediction <- totalCV_result$prediction
	actual <- totalCV_result$actual
	results <- caret::confusionMatrix(prediction, actual)
	#WARNING!!!! removal of this object must NOT happen before all applications of it is completed
	rm(totalCV_result)
	
	#code for saving the results_instances -*-*-*-*-*-*-* start -*-*-*-*-*-*
	instance<-data.frame(as.matrix(results))
	instance_results_CV<-rbind(instance_results_CV,instance)
	#add a blank row behind the results of each REPAT as spacer
	#first create a one row matrix of the same length as the existing data frame, instance_results_CV in this case
	temprow<-matrix(c(rep.int(NA,length(instance_results_CV))),nrow=1,ncol=length(instance_results_CV))
	#make this one row matrix a dataframe and give its cols the same names as instance_results_CV
	newrow<-data.frame(temprow)
	colnames(newrow)<-colnames(instance_results_CV)
	#bind the empty row to the existing dataframe
	instance_results_CV<-rbind(instance_results_CV,newrow)
	#print(instance_results_CV)
	#code for saving the results_instances -*-*-*-*-*-*-* end -*-*-*-*-*-*
	
	
	#code for saving the kappa results -*-*-*-*-*-*-* start -*-*-*-*-*-*
	kappa_res<-data.frame(as.matrix(results,what="overall"))
	kappa_results_CV<-rbind(kappa_results_CV,kappa_res)
	#print(kappa_results_cv)
	
	temprow<-matrix(c(rep.int(NA,length(kappa_results_CV))),nrow=1,ncol=length(kappa_results_CV))
	newrow<-data.frame(temprow)
	colnames(newrow)<-colnames(kappa_results_CV)
	kappa_results_CV<-rbind(kappa_results_CV,newrow)
	#code for saving the kappa results -*-*-*-*-*-*-* end -*-*-*-*-*-*
	
	
	#code for saving the class specificity results -*-*-*-*-*-*-* start -*-*-*-*-*-*
	class_res<-data.frame(as.matrix(results,what="classes"))
	class_results_CV<-rbind(class_results_CV,class_res)
	temprow<-matrix(c(rep.int(NA,length(class_results_CV))),nrow=1,ncol=length(class_results_CV))
	newrow<-data.frame(temprow)
	colnames(newrow)<-colnames(class_results_CV)
	class_results_CV<-rbind(class_results_CV,newrow)
	#code for saving the class specificity results -*-*-*-*-*-*-* end -*-*-*-*-*-*
	
	#line<-readline()
	#write the table of class instances prediction
	write.csv(as.matrix(instance_results_CV), file=paste0("CV_results_instances.csv"))
	write.csv(as.matrix(kappa_results_CV), file=paste0("CV_results_kappa.csv"))	
	write.csv(as.matrix(class_results_CV), file=paste0("CV_results_class.csv"))
#-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	
#-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	
#-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	
#-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	
#-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	
	
	
	print(paste0("repeat no: ",REPEATS," complete"))
#lower bracket for the REPEAT runs loop	
}

#*-*-*-*-*--* calculating and saving feature importance grand average START *-*-*--*-*
variable_importance_list <- data.frame(variable_importance_list)
mean_Importance <- apply(variable_importance_list, MARGIN = 1, FUN = mean)
sd_Importance <- apply(variable_importance_list, MARGIN = 1, FUN = sd)

mean_Importance <- data.frame(round(mean_Importance, digits = 3))
sd_Importance <- data.frame(round(sd_Importance, digits = 3))
feature_importance <- cbind(mean_Importance, sd_Importance)
feature_importance <- data.frame(feature_importance)

write.csv(feature_importance, file = paste0("Feature Importance.csv"))
#*-*-*-*-*--* calculating and saving feature importance grand average END *-*-*--*-*



print("programe complete")