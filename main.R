# Define constants
R_PATH <- "R"
DATA_PATH <- "data"
FILES_NAMES <- c("mfeat-fou", "mfeat-kar", "mfeat-zer")

# Source files
source(file.path(R_PATH, "bayes_functions.R"))
source(file.path(R_PATH, "cv_functions.R"))
source(file.path(R_PATH, "data_functions.R"))
source(file.path(R_PATH, "nn_functions.R"))
source(file.path(R_PATH, "svm_functions.R"))

# Read the data sets
data_sets_list <- readDataSets(DATA_PATH, FILES_NAMES)

# Compute the dissimilarity matrix for each data set
# diss_matrix_list <- lapply(data_sets_list, function(data_set){
#     as.matrix(cluster::daisy(data_set, metric = "euclidean"))
# })

# Add the Class variable to each data set
data_sets_list <- addClasses(data_sets_list)
# Remove repeated and inconsistent examples
data_sets_list <- cleanDataSets(data_sets_list)

# Question 2(a)
# Train an ensemble of bayes classifiers
bayes_ens <- repeatedCVTrain(method = "bayesEnsemble", 
                             data_list = data_sets_list, 
                             seed = 1235, folds = 5, repeats = 3)

# Question 2(b)
svms_params <- list(list(C = 0.5), 
                    list(C = 2), 
                    list(C = 0.25))
svm_ens <- repeatedCVTrain(method = "svmEnsemble", 
                           data_list = data_sets_list, 
                           method_args = svms_params, 
                           seed = 1235, folds = 5, repeats = 3)

nn_params <- list(list(size = 9, decay = 5e-2, maxit = 1500), 
                  list(size = 9, decay = 5e-2, maxit = 1500), 
                  list(size = 9, decay = 5e-2, maxit = 1500))
nn_ens <- repeatedCVTrain(method = "nnEnsemble", 
                          data_list = data_sets_list, 
                          method_args = nn_params, 
                          seed = 1235, pre_process = c("center", "scale"), 
                          folds = 5, repeats = 3)

# Print the results
print(bayes_ens)
print(svm_ens)
print(nn_ens)
