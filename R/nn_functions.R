library(nnet)

# Function to train the ensemble of neural networks.
nnEnsemble <- function(data_sets_list, 
                       params_list = list(list(size = 3, decay = 0, maxit = 100)), 
                       pre_process = c("center", "scale")){
    # Pre-process the training sets
    pre_processed_data_sets_list <- lapply(data_sets_list,
                                          preProcessTraining, pre_process = pre_process)
    data_sets_list <- extractTrainingSets(pre_processed_data_sets_list)
    
    # Training function
    train <- function(training, params_list){
        size <- if (!is.null(params_list$size)) params_list$size else 3
        decay <- if (!is.null(params_list$decay)) params_list$decay else 0
        maxit <- if (!is.null(params_list$maxit)) params_list$maxit else 100
        nnet(Class ~ ., data = training, size = size,
             decay = decay, maxit = maxit, MaxNWts = 1500)
    }
    nn_ens <- mapply(train, data_sets_list, params_list, SIMPLIFY = FALSE)
    
    pre_process_args <- extractPreProcessArgs(pre_processed_data_sets_list)
    nn_ens <- c(nn_ens, pre_process_args)
    
    class(nn_ens) <- "nn_ensemble"
    nn_ens
}

# Generic predict function for the ensemble of neural networks.
predict.nn_ensemble <- function(object, newdata, ...){
    # object: A list of nns
    # newdata: A list of data sets to predict
    num_classifiers <- length(object) / 2
    newdata <- mapply(preProcessTesting, 
                      newdata, object[(num_classifiers + 1):length(object)], 
                      MoreArgs = list(...), SIMPLIFY = FALSE)
    object <- object[1:num_classifiers]
    
    majorityVote(object, newdata, type = "class")
}