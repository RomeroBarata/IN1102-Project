library(nnet)

nnEnsemble <- function(data_sets_list, 
                       params_list = list(list(size = 3, decay = 0, maxit = 100))){
    # Transform the last column of each data set to be a factor
    data_sets_list <- lapply(data_sets_list, function(data){
        data$Class <- factor(data$Class)
        data
    })
    
    # Training function
    train <- function(training, params_list){
        size <- if (!is.null(params_list$size)) params_list$size else 3
        decay <- if (!is.null(params_list$decay)) params_list$decay else 0
        maxit <- if (!is.null(params_list$maxit)) params_list$maxit else 100
        nnet(Class ~ ., data = training, size = size,
             decay = decay, maxit = maxit, MaxNWts = 1500)
    }
    nn_ens <- mapply(train, data_sets_list, params_list, SIMPLIFY = FALSE)
    
    class(nn_ens) <- "nn_ensemble"
    nn_ens
}

predict.nn_ensemble <- function(object, newdata){
    # object: A list of nns
    # newdata: A list of data sets to predict
    majorityVoteNN(object, newdata)
}

majorityVoteNN <- function(objects_list, newdata_list){
    predictions <- mapply(predict, objects_list, newdata_list, SIMPLIFY = FALSE)
    predictions <- sapply(predictions, function(mat){
        apply(mat, 1, which.max) - 1
    })
    predictions <- apply(predictions, 1, function(row){
        idx <- which.max(table(row))
        as.integer(names(table(row))[idx])
    })
}