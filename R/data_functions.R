readDataSets <- function(data_path, files_names){
    lapply(files_names, function(file_name, data_path){
        read.table(file.path(data_path, file_name), header = FALSE)
    }, data_path = data_path)
}

addClasses <- function(data_sets_list){
    lapply(data_sets_list, transform, Class = rep(0:9, each = 200))
}

cleanDataSets <- function(data_sets_list){
    lapply(data_sets_list, function(data) removeInconsistencies(unique(data)))
}

removeInconsistencies <- function(data){
    idx <- rep(FALSE, nrow(data))
    data_atts <- data[, -ncol(data)]
    hashs <- as.character(apply(data_atts, 1, digest::digest))
    duplicated_hashs <- hashs[duplicated(hashs)]
    for (dhash in duplicated_hashs){
        idx <- idx | (dhash == hashs) 
    }
    data[!idx, ]
}

preProcessTraining <- function(training, pre_process = NULL){
    if (is.null(pre_process)) return(list(training = training))
    
    training_class <- training[, ncol(training), drop = FALSE]
    training <- training[, -ncol(training), drop = FALSE]
    
    center <- scale_factor <- NULL
    if ("center" %in% pre_process){
        training <- scale(training, center = TRUE, scale = FALSE)
        center <- attributes(training)$`scaled:center`
        attributes(training)$`scaled:center` <- NULL
    }
    if ("scale" %in% pre_process){
        training <- scale(training, center = FALSE, 
                          scale = apply(training, 2, function(x) sqrt(sum(x^2))))
        scale_factor <- attributes(training)$`scaled:scale`
        attributes(training)$`scaled:scale` <- NULL
    }
    
    result <- list()
    result$training <- cbind(training, training_class)
    result$center <- center
    result$scale_factor <- scale_factor
    result
}

extractTrainingSets <- function(pre_processed_training_list){
    lapply(pre_processed_training_list, function(x) x$training)
}

preProcessTesting <- function(testing, pre_processed_training, pre_process = NULL){
    if (is.null(pre_process)) return(testing)
    
    if ("center" %in% pre_process){
        center <- pre_processed_training$center
        testing <- t(apply(testing, 1, 
                           function(x, center) x - center, 
                           center = center))
    }
    if ("scale" %in% pre_process){
        scale_factor <- pre_processed_training$scale_factor
        testing <- t(apply(testing, 1, 
                           function(x, scale_factor) x / scale_factor, 
                           scale_factor = scale_factor))
    }
    
    testing
}