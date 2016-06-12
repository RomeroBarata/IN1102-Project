library(ggplot2)

# Function to create an appropriate data frame for 
# the ggplot function of the ggplot2 package.
buildPlotDataFrame <- function(bayes, svm, nn, mix){
    accuracy_matrix <- cbind(bayes$Accuracy, svm$Accuracy, 
                             nn$Accuracy, mix$Accuracy)
    data.frame(Classifier = c("Bayes", "SVM", "NN", "Mix"), 
               Acc = c(bayes$Accuracy_MEAN, svm$Accuracy_MEAN, 
                       nn$Accuracy_MEAN, mix$Accuracy_MEAN), 
               Err = apply(accuracy_matrix, 2, createConfidenceInterval)[1, ])
}

# Function to create dotplot of the classifiers' results.
createDotPlot <- function(results_df){
    p <- ggplot(results_df, aes(x = Classifier, y = Acc)) + 
        scale_x_discrete(limits = c("Bayes", "NN", "Mix", "SVM")) + theme_minimal()
    p + geom_point(size = 2.5) + 
        geom_errorbar(aes(ymin = Acc - Err, ymax = Acc + Err, width = 0.1)) + 
        xlab("Classificador") + 
        ylab("Taxa de Acerto") + 
        coord_flip()
}