majorityVote <- function(objects_list, newdata_list, ...){
    predictions <- mapply(predict, objects_list, newdata_list, MoreArgs = list(...))
    predictions <- apply(predictions, 1, function(row){
        idx <- which.max(table(row))
        as.integer(names(table(row))[idx])
    })
}