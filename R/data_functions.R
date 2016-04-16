readDataSets <- function(data_path, files_names){
    lapply(files_names, function(file_name, data_path){
        read.table(file.path(data_path, file_name), header = FALSE)
    }, data_path = data_path)
}

# This function loads the dissimilarity matrices. If any matrix was
# not generated, it will generate all of them and save into a file.
# 'data_files' and 'diss_data_files' must have the same length.
loadDissMtcs <- function(data_path, data_files, diss_data_files) {
#    input_files <- c("mfeat-fac", "mfeat-fou", "mfeat-kar")

#    output_files <- c(file.path(data_path, "mfeat-fac-diss"),
#                      file.path(data_path, "mfeat-fou-diss"),
#                      file.path(data_path, "mfeat-kar-diss"))

    input_files <- data_files

    output_files <- diss_data_files

    # Gets only those whose dissimilarity matrix does not exist
    if(sum(!sapply(output_files, file.exists))) {
        # Read the data sets
        data_sets_list <- readDataSets(data_path, input_files)

        # Compute the dissimilarity matrix for each data set
        diss_matrix_list <- lapply(data_sets_list, function(data_set){
            as.matrix(cluster::daisy(data_set, metric = "euclidean"))
        })

        apply(rbind(output_files, diss_matrix_list), 2,
              function(tuple, data_path){
                  filename <- file.path(data_path, tuple[1])
                  if(!file.exists(filename)) {
                      print(filename)
                      write.table(tuple[2], file = filename,
                                  quote = FALSE, row.names = FALSE,
                                  col.names = FALSE)
                  }
              }, data_path = data_path)
        return(diss_matrix_list)
    } else {
        return(readDataSets(".", output_files))
    }
}
