library(nnet)

nnEnsemble <- function(data_sets_list, size = 3, decay = 0, maxit = 100){
    # Transform the last column of each data set to be a factor
    data_sets_list <- lapply(data_sets_list, function(data){
        data$Class <- factor(data$Class)
        data
    })
    
    # Train the ensemble of Neural Networks
    nn_ens <- lapply(data_sets_list, function(training){
        nnet(Class ~ ., data = training, size = size, 
             decay = decay, maxit = maxit, MaxNWts = 1500)
    })
    
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