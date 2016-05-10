bayes <- function(training){
    mus <- ddply(training, .(Class), colwise(mean))
    cov_mats <- ddply(training, .(Class), function(x){
        x$Class <- NULL
        x <- cov(x)
    })
    a_priori_probs <- ddply(training, .(Class), nrow)
    a_priori_probs[, 2] <- a_priori_probs[, 2] / sum(a_priori_probs[, 2])
    
    bayes_model <- buildModel(mus, cov_mats, a_priori_probs)
    class(bayes_model) <- "bayes"
    bayes_model
}

bayesEnsemble <- function(data_sets_list){
    bayes_ens <- lapply(data_sets_list, function(training) bayes(training))
    class(bayes_ens) <- "bayes_ensemble"
    bayes_ens
}

buildModel <- function(mus, cov_mats, a_priori_probs){
    bayes_model <- list()
    for(i in 1:nrow(mus)){
        current_name <- paste("dis_fun_", i - 1, sep = "")
        current_mu <- subset(mus, Class == (i - 1))
        current_mu$Class <- NULL
        current_cov <- subset(cov_mats, Class == (i - 1))
        current_cov$Class <- NULL
        current_priori <- subset(a_priori_probs, Class == (i - 1))
        current_priori$Class <- NULL
        bayes_model[[current_name]]$mu <- as.numeric(current_mu)
        bayes_model[[current_name]]$cov_mat <- data.matrix(current_cov)
        bayes_model[[current_name]]$a_priori <- as.numeric(current_priori)
    }
    bayes_model
}

predict.bayes <- function(object, newdata){
    predictions <- sapply(object, function(dis_fun, newdata){
        apply(newdata, 1, dmvnorm, 
              mu = dis_fun$mu, Sigma = dis_fun$cov_mat) + log(dis_fun$a_priori)
    }, newdata = newdata)
    predictions <- apply(predictions, 1, which.max) - 1
}

predict.bayes_ensemble <- function(object, newdata){
    # object: A list of bayes classifiers
    # newdata: A list of data sets to predict
    majorityVote(object, newdata)
}

dmvnorm <- function(x, mu, Sigma, log = TRUE){
    if(log){
        p1 <- -0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu)
        p2 <- -0.5 * as.numeric(unlist(determinant(Sigma))["modulus"])
        return (p1 + p2)
    } else{
        d <- length(mu)
        p1 <- 1 / ((2 * pi)^(d / 2) * sqrt(abs(det(Sigma))))
        p2 <- exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
        return (p1 * p2)
    }
}

majorityVote <- function(objects_list, newdata_list){
    predictions <- mapply(predict, objects_list, newdata_list)
    predictions <- apply(predictions, 1, function(row){
        idx <- which.max(table(row))
        as.integer(names(table(row))[idx])
    })
}