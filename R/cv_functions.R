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

repeatedCVTrainBayes <- function(data_list, ...){
    y <- data_list[[1]][ncol(data_list[[1]])]
    sps <- createStratifiedFolds(y, ...)
    results <- as.vector(laply(sps, cvTrainBayes, data_list))
    
    final_results <- list(Accuracy = results, 
                          Accuracy_MEAN = mean(results), Accuracy_SD = sd(results))
}

cvTrainBayes <- function(spartition, data_list){
    nfolds <- length(unique(spartition))
    accuracy <- vector(mode = "numeric", length = nfolds)
    for(i in 1:nfolds){
        bayes_model1 <- bayes(data_list[[1]][!(spartition == i), ])
        bayes_model2 <- bayes(data_list[[2]][!(spartition == i), ])
        bayes_model3 <- bayes(data_list[[3]][!(spartition == i), ])
        
        predictions <- majorityVote(list(bayes_model1, bayes_model2, bayes_model3), 
                                    list(data_list[[1]][spartition == i, ], 
                                         data_list[[2]][spartition == i, ], 
                                         data_list[[3]][spartition == i, ]))
        
        accuracy[i] <- mean(predictions == data_list[[1]][spartition == i, ]$Class)
    }
    accuracy
}