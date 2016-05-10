library(plyr)

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

repeatedCVTrain <- function(method, data_list, 
                            method_args = list(), seed = NULL, pre_process = NULL, ...){
    y <- data_list[[1]][ncol(data_list[[1]])]
    if (!is.null(seed)) set.seed(seed)
    sps <- createStratifiedFolds(y, ...)
    results <- as.vector(laply(sps, cvTrain, data_list, method, method_args, pre_process))
    
    final_results <- list(Accuracy = results, 
                          Accuracy_MEAN = mean(results), Accuracy_SD = sd(results))
}

cvTrain <- function(spartition, data_list, 
                    method, method_args = list(), pre_process = NULL){
    nfolds <- length(unique(spartition))
    accuracy <- vector(mode = "numeric", length = nfolds)
    for(i in 1:nfolds){
        training_list <- splitCVTrain(data_list, spartition, i)
        # Pre-process the training sets
        pre_processed_training_list <- lapply(training_list, 
                                              preProcessTraining, pre_process = pre_process)
        training_list <- extractTrainingSets(pre_processed_training_list)
        # Train the model
        model <- do.call(method, c(list(training_list), method_args))
        
        testing_list <- splitCVTest(data_list, spartition, i)
        # Pre-process the training sets
        testing_list <- mapply(preProcessTesting, 
                               testing_list, pre_processed_training_list, 
                               MoreArgs = list(pre_process = pre_process), SIMPLIFY = FALSE)
        # Make the predictions
        predictions <- predict(model, testing_list)
        
        accuracy[i] <- mean(predictions == data_list[[1]][spartition == i, ]$Class)
    }
    accuracy
}

splitCVTrain <- function(data_sets_list, spartition, i){
    lapply(data_sets_list, function(data, spartition, i){
        data[!(spartition == i), ]
    }, spartition = spartition, i = i)
}

splitCVTest <- function(data_sets_list, spartition, i){
    lapply(data_sets_list, function(data, spartition, i){
        data[spartition == i, -ncol(data)]
    }, spartition = spartition, i = i)
}