library(PMCMR)
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
seed = 1235
folds <- 40
repeats <- 1
# Train an ensemble of bayes classifiers
bayes_ens <- repeatedCVTrain(method = "bayesEnsemble", 
                             data_list = data_sets_list, 
                             seed = seed, folds = folds,
                             repeats = repeats)

# Question 2(b)
svm_ens <- repeatedCVTrain(method = "svmEnsemble", 
                           data_list = data_sets_list, 
                           method_args = list(C = 1), 
                           seed = seed, folds = folds,
                           repeats = repeats)

nn_ens <- repeatedCVTrain(method = "nnEnsemble", 
                          data_list = data_sets_list, 
                          method_args = list(size = 9, decay = 5e-2, maxit = 1500), 
                          seed = seed, pre_process = c("center", "scale"), 
                          folds = folds, repeats = repeats)

# Print the results
print(bayes_ens)
print(svm_ens)
print(nn_ens)

# Table with the accuracies for each fold of each model
accuracy_matrix <- matrix(c(bayes_ens$Accuracy, svm_ens$Accuracy,
                           nn_ens$Accuracy), nrow = folds, ncol = 3,
                         dimnames = list(NULL, c("bayes", "svm", "nn")
                                         ))

print(accuracy_matrix)

# Performing statistical test
test <- friedman.test(accuracy_matrix);
print(test);

# Performing statistical post-test
post_test <- posthoc.friedman.nemenyi.test(accuracy_matrix)
print(post_test)
