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