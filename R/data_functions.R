readDataSets <- function(data_path, files_names){
    lapply(files_names, function(file_name, data_path){
        read.table(file.path(data_path, file_name), header = FALSE)
    }, data_path = data_path)
}

addClasses <- function(data_sets_list){
    lapply(data_sets_list, transform, Class = rep(0:9, each = 200))
}