R_PATH <- "R"
DATA_PATH <- "data"
FILES_NAMES <- c("mfeat-fou", "mfeat-kar", "mfeat-zer")

source(file.path(R_PATH, "data_functions.R"))

diss_data_files <- sapply(FILES_NAMES, paste, "-diss", sep = "")

genDissMtcs(DATA_PATH, FILES_NAMES, diss_data_files)
