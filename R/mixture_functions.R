# Function to train the mixture of ensembles.
mixture <- function(data_sets_list, params_list = NULL, 
                    pre_process = c("center", "scale")){
    bayes_ens <- bayesEnsemble(data_sets_list)
    svm_ens <- svmEnsemble(data_sets_list, 
                           params_list = params_list[["svm_ensemble"]], 
                           pre_process = pre_process)
    nn_ens <- nnEnsemble(data_sets_list, 
                         params_list = params_list[["nn_ensemble"]], 
                         pre_process = pre_process)
    mix_ens <- list(bayes_ens = bayes_ens, svm_ens = svm_ens, nn_ens = nn_ens)
    class(mix_ens) <- "mix_ens"
    mix_ens
}

# Generic predict function for the mix of ensembles.
predict.mix_ens <- function(object, newdata, ...){
    predictions <- sapply(object, predict, newdata = newdata, ...)
    apply(predictions, 1, function(row){
        idx <- which.max(table(row))
        as.integer(names(table(row))[idx])
    })
}