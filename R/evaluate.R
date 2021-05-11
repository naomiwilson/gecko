#' Make Confusion Matrix of Results
#' @param p1v2 second matrix in output of perform_NBC aka "sample_by_sample"
#' @export
make_confusion_matrix = function(p1v2) {
  final = as.factor(p1v2["result",])
  truth = as.factor(p1v2["truth",])
  output = caret::confusionMatrix(final, truth)
}

#' Get list of features included in model
#' @param model output from train_NBC
#' @export
get_modeled_features = function(model) {
  modeled_features = colnames(model[[1]])
  return(modeled_features)
}
