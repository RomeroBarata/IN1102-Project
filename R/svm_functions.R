library(kernlab)

svmEnsemble <- function(data_sets_list, 
                        params_list = list(list(C = 1))){
    # Transform the last column of each data set on a factor
    data_sets_list <- lapply(data_sets_list, function(data){
        data$Class <- factor(data$Class)
        data
    })
    
    # Training function
    train <- function(training, params_list){
        C <- if (!is.null(params_list$C)) params_list$C else 1
        ksvm(Class ~ ., data = training, kernel = "rbfdot", C = C)
    }
    svm_ens <- mapply(train, data_sets_list, params_list, SIMPLIFY = FALSE)
    
    class(svm_ens) <- "svm_ensemble"
    svm_ens
}

predict.svm_ensemble <- function(object, newdata){
    # object: A list of svms
    # newdata: A list of data sets to predict
    majorityVote(object, newdata)
}