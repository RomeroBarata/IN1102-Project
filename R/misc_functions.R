# Function to compute the majority vote of a set of classifiers
# on a set of data sets.
majorityVote <- function(objects_list, newdata_list, ...){
    predictions <- mapply(predict, objects_list, newdata_list, MoreArgs = list(...))
    predictions <- apply(predictions, 1, function(row){
        idx <- which.max(table(row))
        as.integer(names(table(row))[idx])
    })
}

# Function to perform a brute-force optimization of the parameters
# of a SVM.
bruteForceOptimSVM <- function(..., params){
    mean_acc <- vector(mode = "numeric", length = nrow(params))
    for (i in 1:nrow(params)){
        svms_params <- list(list(C = params[i, 1]), 
                            list(C = params[i, 2]), 
                            list(C = params[i, 3]))
        mean_acc[i] <- repeatedCVTrain(..., method_args = svms_params)$Accuracy_MEAN
    }
    cbind(params, mean_acc)
}

# Function to perform a brute-force optimization of the parameters
# of a NN.
bruteForceOptimNN <- function(..., params){
    mean_acc <- vector(mode = "numeric", length = nrow(params))
    for (i in 1:nrow(params)){
        nns_params <- list(list(size = params[i, 1], decay = params[i, 2], maxit = 1500), 
                           list(size = params[i, 3], decay = params[i, 4], maxit = 1500), 
                           list(size = params[i, 5], decay = params[i, 6], maxit = 1500))
        mean_acc[i] <- repeatedCVTrain(..., method_args = nns_params)$Accuracy_MEAN
    }
    cbind(params, mean_acc)
}

# Function to create a confidence interval for a set of real values
# using a t-distribution (thus, not assuming normality of the values).
createConfidenceInterval <- function(values, confidence = 0.99){
    x <- mean(values)
    s <- sd(values)
    n <- length(values)
    
    error <- qt(confidence, df = n - 1) * s / sqrt(n)
    c(error, x - error, x + error)
}