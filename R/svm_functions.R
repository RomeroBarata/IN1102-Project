library(kernlab)

svmEnsemble <- function(data_sets_list, C = 1){
    # Transform the last column of each data set to be a factor
    data_sets_list <- lapply(data_sets_list, function(data){
        data$Class <- factor(data$Class)
        data
    })
    
    # Train the ensemble of SVMs
    svm_ens <- lapply(data_sets_list, function(training){
        ksvm(Class ~ ., data = training, kernel = "rbfdot", C = C)
    })
    
    class(svm_ens) <- "svm_ensemble"
    svm_ens
}

predict.svm_ensemble <- function(object, newdata){
    # object: A list of svms
    # newdata: A list of data sets to predict
    majorityVote(object, newdata)
}