if (!require(needs)) install.packages('needs')
# Import the necessary packages
library(needs)
needs(PMCMR, plyr, nnet, ggplot2, kernlab)

# Define constants
R_PATH <- "R"
DATA_PATH <- "data"
FILES_NAMES <- c("mfeat-fou", "mfeat-kar", "mfeat-zer")
FOLDS <- 40
REPEATS <- 1
SEED <- 1235

# Source files
source(file.path(R_PATH, "bayes_functions.R"))
source(file.path(R_PATH, "cv_functions.R"))
source(file.path(R_PATH, "data_functions.R"))
source(file.path(R_PATH, "misc_functions.R"))
source(file.path(R_PATH, "mixture_functions.R"))
source(file.path(R_PATH, "nn_functions.R"))
source(file.path(R_PATH, "plot_functions.R"))
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
                             seed = SEED, folds = FOLDS, repeats = REPEATS)

# Question 2(b)
svms_params <- list(list(C = 32), 
                    list(C = 4), 
                    list(C = 8))
svm_ens <- repeatedCVTrain(method = "svmEnsemble", 
                           data_list = data_sets_list, 
                           method_args = svms_params, 
                           seed = SEED, folds = FOLDS, repeats = REPEATS)

nns_params <- list(list(size = 12, decay = 5e-2, maxit = 1500), 
                   list(size = 12, decay = 5e-2, maxit = 1500), 
                   list(size = 12, decay = 5e-2, maxit = 1500))
nn_ens <- repeatedCVTrain(method = "nnEnsemble", 
                          data_list = data_sets_list, 
                          method_args = nns_params, 
                          seed = SEED, pre_process = c("center", "scale"), 
                          folds = FOLDS, repeats = REPEATS)

# Question 2(c)
mix_params <- list(svm_ensemble = svms_params, 
                   nn_ensemble = nns_params)
mix_ens <- repeatedCVTrain(method = "mixture", 
                           data_list = data_sets_list, 
                           method_args = mix_params, 
                           seed = SEED, pre_process = c("center", "scale"), 
                           folds = FOLDS, repeats = REPEATS)

# Print the results
print(bayes_ens)
print(svm_ens)
print(nn_ens)
print(mix_ens)

# Table with the accuracies for each fold of each model
accuracy_matrix <- matrix(c(bayes_ens$Accuracy, svm_ens$Accuracy,
                            nn_ens$Accuracy, mix_ens$Accuracy), nrow = FOLDS, ncol = 4,
                          dimnames = list(NULL, c("bayes", "svm", "nn", "mix")
                          ))

print(accuracy_matrix)

# Performing statistical test
test <- friedman.test(accuracy_matrix)
print(test)

# Performing statistical post-test
post_test <- posthoc.friedman.nemenyi.test(accuracy_matrix)
print(post_test)

# Build a data frame for the plotting function
results_df <- buildPlotDataFrame(bayes_ens, svm_ens, nn_ens, mix_ens)
# Plot the results
createDotPlot(results_df)
