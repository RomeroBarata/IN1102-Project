library(plyr)

# Function to create the stratified folds given the classes
# of the examples.
createStratifiedFolds <- function(y, folds = 10, repeats = 5){
    sps <- list()
    classes_dist <- daply(y, .(Class), nrow)
    for (i in 1:repeats){
        partition <- vector(mode = "numeric", length = nrow(y))
        for(j in 1:length(classes_dist)){
            max_sample <- ceiling(classes_dist[j] / folds) * folds
            class_partition <- sample(rep(1:folds, length.out = max_sample))[1:classes_dist[j]]
            class_id <- as.integer(names(classes_dist)[j])
            partition[as.vector(y == class_id)] <- class_partition
        }
        sps[[paste("rep_", i, sep = "")]] <- partition
    }
    sps
}

# Function that implements the repeated cross-validated
# training process.
repeatedCVTrain <- function(method, data_list, 
                            method_args = list(list()), 
                            seed = NULL, pre_process = NULL, ...){
    y <- data_list[[1]][ncol(data_list[[1]])]
    if (!is.null(seed)) set.seed(seed)
    sps <- createStratifiedFolds(y, ...)
    results <- as.vector(laply(sps, cvTrain, data_list, method, method_args, pre_process))
    
    final_results <- list(Accuracy = results, 
                          Accuracy_MEAN = mean(results), Accuracy_SD = sd(results))
}

# Auxiliary function for the repeatedCVTrain function.
cvTrain <- function(spartition, data_list, 
                    method, method_args = list(list()), pre_process = NULL){
    nfolds <- length(unique(spartition))
    method_args <- if (sum(sapply(method_args, length)) > 0) method_args else NULL
    accuracy <- vector(mode = "numeric", length = nfolds)
    for(i in 1:nfolds){
        training_list <- splitCVTrain(data_list, spartition, i)
        # Train the model
        if (!is.null(method_args))
            model <- do.call(method, list(training_list, method_args, pre_process))
        else
            model <- do.call(method, list(training_list, pre_process))
        
        testing_list <- splitCVTest(data_list, spartition, i)
        # Make the predictions
        predictions <- predict(model, testing_list, pre_process = pre_process)
        
        accuracy[i] <- mean(predictions == data_list[[1]][spartition == i, ]$Class)
    }
    accuracy
}

# Auxiliary function for the cvTrain function to extract the training
# set from a data set given the stratified partition.
splitCVTrain <- function(data_sets_list, spartition, i){
    lapply(data_sets_list, function(data, spartition, i){
        data[!(spartition == i), ]
    }, spartition = spartition, i = i)
}

# Auxiliary function for the cvTrain function to extract the testing
# set from a data set given the stratified partition.
splitCVTest <- function(data_sets_list, spartition, i){
    lapply(data_sets_list, function(data, spartition, i){
        data[spartition == i, -ncol(data)]
    }, spartition = spartition, i = i)
}