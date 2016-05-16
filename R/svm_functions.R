library(kernlab)

svmEnsemble <- function(data_sets_list, 
                        params_list = list(list(C = 1)), 
                        pre_process = c("center", "scale")){
    # Training function
    train <- function(training, params_list){
        C <- if (!is.null(params_list$C)) params_list$C else 1
        ksvm(Class ~ ., data = training, 
             scaled = if (!is.null(pre_process)) TRUE else FALSE, 
             kernel = "rbfdot", C = C)
    }
    svm_ens <- mapply(train, data_sets_list, params_list, SIMPLIFY = FALSE)
    
    class(svm_ens) <- "svm_ensemble"
    svm_ens
}

predict.svm_ensemble <- function(object, newdata, ...){
    # object: A list of svms
    # newdata: A list of data sets to predict
    majorityVote(object, newdata)
}