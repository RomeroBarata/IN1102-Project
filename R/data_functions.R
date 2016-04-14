readDataSets <- function(data_path, files_names){
    lapply(files_names, function(file_name, data_path){
        read.table(file.path(data_path, file_name), header = FALSE)
    }, data_path = data_path)
}