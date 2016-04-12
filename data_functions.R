readDataSets <- function(data_path, files_names){
  lapply(files_names, function(file_name, data_path){
    read.table(file = paste(data_path, file_name, sep = ""), header = FALSE)
  }, data_path = data_path)
}