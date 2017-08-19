library(readxl)
library(caret)
require(sampling)
library(ranger)
#require(data.table)


#Change this to the working directory
setwd("d:/R_test/all_classes")

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
seed_vec<-c(seed_vec,867)


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
  
 #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing start *-*-*-*-*-*-*-*-*-*-*-
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing start *-*-*-*-*-*-*-*-*-*-*-
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing start *-*-*-*-*-*-*-*-*-*-*-
class_level<-levels(Class)

max_class=-1
min_class=1e5
class_list<-list()
for(i in 1:length(class_level))
{
  indiv_class<-class_level[i]
  #print(indiv_class)
  temp_class<-subset(df,Class==indiv_class)
  class_list<-append(class_list,list(temp_class))
  
  class_length<-nrow(temp_class)
 # print(temp_class)
#  print(nrow(temp_class))
  if(class_length<min_class)
    min_class<-class_length
  if(class_length>max_class)
    max_class<-class_length
}

#print("class_list")
#print(class_list)
#print("min and max class")
#print(min_class)
#print(max_class)

#oversampling creates mulitple training data sets, each training data set should be the size of the smallest dataset?
training_cluster_count<-(max_class/min_class)
rounded_training_cluster_count<-round(training_cluster_count)

#print("training_cluster_count")
#print(training_cluster_count)
#print("rounded_training_cluster_count")
#print(rounded_training_cluster_count)

if(rounded_training_cluster_count<training_cluster_count){
  rounded_training_cluster_count<-rounded_training_cluster_count+1
}

#print("readjusted cluster count")
#print(rounded_training_cluster_count)

total_sample_required_per_class<-min_class*rounded_training_cluster_count

#print("total_sample_required_per_class")
#print(total_sample_required_per_class)
#stop()

#finalized_balanced_classes<-list()
finalized_balanced_classes<-data.frame()
  
for(i in 1:length(class_level))
{
  temp_class<-data.frame(class_list[i])
  
  balanced_class_temp<-data.frame()
  #in all case, a first filling of this temporal balanced frame with one scrambled copie of the current class is needed
  append_samples<-sample(nrow(temp_class),replace = FALSE)
  balancing_samples<-temp_class[append_samples,]
  balanced_class_temp<-rbind(balanced_class_temp,balancing_samples)

  missing_sample_count<-total_sample_required_per_class-(nrow(balanced_class_temp))
  while(nrow(balanced_class_temp)<total_sample_required_per_class)
  {
    #remember, its the current class count that needs to be compared to the number of missing samples in the BALANCED class
    #NOT the current appended balanced class sample count
    if(nrow(temp_class)<missing_sample_count)
    {
      #the missing sample number is greater than the number of samples in this class, scramble the entire current class and append it
      append_samples<-sample(nrow(temp_class),replace = FALSE)
      balancing_samples<-temp_class[append_samples,]
      balanced_class_temp<-rbind(balanced_class_temp,balancing_samples)
    } else
    {
      #before sampling a lower number of samples from a class, the whole class index needs to be scrampled
      append_samples<-sample(nrow(temp_class),replace = FALSE)
      balancing_samples<-temp_class[append_samples,]
      
      #now pick only the missing number of index from the scrambled index

      
      #missing_sample_count<-total_sample_required_per_class-(nrow(balanced_class_temp))
      balancing_samples<-balancing_samples[1:missing_sample_count,]
      balanced_class_temp<-rbind(balanced_class_temp,balancing_samples)
    }
    
    #the missing sample count changes at each loop and needs to be recalculated
    missing_sample_count<-total_sample_required_per_class-(nrow(balanced_class_temp))
  }
 
  #print("the balanced class")
  #print(balanced_class_temp)
  
  #at this point each iteration has oversampled a class to a desireable number samples store them in a list for further processing
  #afterwards, the memory needs to be cleared
  
  #this finalized_balanced_classes is a list that contains the balanced classes as elements, to
  #use them, one needs to first extract each element and then convert it into a data frame
  
 # finalized_balanced_classes<-append(finalized_balanced_classes,list(balanced_class_temp))
  finalized_balanced_classes<-rbind(finalized_balanced_classes,balanced_class_temp)
  

} #lower bracket for iterating through class_list

#for(i in 1:length(finalized_balanced_classes))
#{
 # print("balanced class no.")
 # print(i)
 # print("*-*-*--**-")
 # print(finalized_balanced_classes[i])
#}


#cleaning the memory
rm(class_level)
rm(max_class)
rm(min_class)
rm(class_list)
rm(indiv_class)
rm(temp_class)
rm()
rm(training_cluster_count)
rm(rounded_training_cluster_count)
rm(balanced_class_temp)
rm(append_samples)
rm(balancing_samples)


#print(finalized_balanced_classes)
#line<-readline()

df<-finalized_balanced_classes
rm(finalized_balanced_classes)

#write.csv(df, file=paste0("finalized_balanced_class.csv"))
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing end *-*-*-*-*-*-*-*-*-*-*-
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing end *-*-*-*-*-*-*-*-*-*-*-
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- code for oversample class balancing end *-*-*-*-*-*-*-*-*-*-*- 
  
#print("got to here")

#print(df)
#line<-readline()
#-*-*-*-*-*-*-*-*-* code for splitting the complete dataset into test and train set START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  training_set_index<-caret::createDataPartition(df$Class,p=0.80,list=FALSE)
  training_set<-df[training_set_index,]
  test_set<-df[-training_set_index,]
  
  df<-training_set
#  print("the df frame reassigned as training set")
 # print(df)
#-*-*-*-*-*-*-*-*-* code for splitting the complete dataset into test and train set END -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  

#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
#-*-*-*-*-*-*-*-*-* code for model recall performance START -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  

  #IMPORTANT!!! this model training on the training set is to be used for the test set as well, do not overwrite it
  model <- ranger::ranger(Class ~., data = df, write.forest  = TRUE, save.memory = TRUE, importance = "impurity", num.trees = 500)
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
 
  
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*		
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  actual<- test_set$Class
  
  prediction <- predict(model, test_set)
  prediction <- prediction$predictions
  
  results <- caret::confusionMatrix(prediction, actual)
  
  #print("predictions")
  #print(prediction)
  #print("actual")
  #print(actual)
  #line<-readline()
  #code for saving the results_instances -*-*-*-*-*-*-*-*-*-*-*-*-*
  instance<-data.frame(as.matrix(results))
  instance_results_test<-rbind(instance_results_test,instance)
  temprow<-matrix(c(rep.int(NA,length(instance_results_test))),nrow=1,ncol=length(instance_results_test))
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(instance_results_test)
  instance_results_test<-rbind(instance_results_test,newrow)
  
  #code for saving the kappa results -*-*-*-*-*-*-* start -*-*-*-*-*-*
  kappa_res<-data.frame(as.matrix(results,what="overall"))
  kappa_results_test<-rbind(kappa_results_test,kappa_res)
  temprow<-matrix(c(rep.int(NA,length(kappa_results_test))),nrow=1,ncol=length(kappa_results_test))
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(kappa_results_test)
  kappa_results_test<-rbind(kappa_results_test,newrow)
  
  #code for saving the class specificity results -*-*-*-*-*-*-* start -*-*-*-*-*-*
  class_res<-data.frame(as.matrix(results,what="classes"))
  class_results_test<-rbind(class_results_test,class_res)
  temprow<-matrix(c(rep.int(NA,length(class_results_test))),nrow=1,ncol=length(class_results_test))
  newrow<-data.frame(temprow)
  colnames(newrow)<-colnames(class_results_test)
  class_results_test<-rbind(class_results_test,newrow)
  
  write.csv(as.matrix(instance_results_test), file=paste0("Test_results_instances.csv"))
  write.csv(as.matrix(kappa_results_test), file=paste0("Test_results_kappa.csv"))	
  write.csv(as.matrix(class_results_test), file=paste0("Test_results_class.csv"))
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*		
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#-*-*-*-*-*-*-* Code to run test set -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* end -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
  
  
  
  
  

#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
#-*-*-*-*-*-*-*-*-* Code to run 10CV for ONE time -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* start -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  folds<-createFolds(df$Class,10,list=TRUE,returnTrain=FALSE)
	totalCV_result<-data.frame()
	
	
	for (i in 1:length(folds))
	{
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
