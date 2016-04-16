# Define constants
R_PATH <- "R"
DATA_PATH <- "data"
FILES_NAMES <- c("mfeat-fac", "mfeat-fou", "mfeat-kar")
DISS_FILES_NAMES <- c("mfeat-fac-diss", "mfeat-fou-diss",
                      "mfeat-kar-diss")

# Creates the dissimilarity matrices if they were not created
source(file.path(R_PATH, "data_functions.R"))

print("Loading dissimilarity matrices.", quote = FALSE)
diss_matrix_list <- loadDissMtcs(DATA_PATH, FILES_NAMES,
                                 DISS_FILES_NAMES)
