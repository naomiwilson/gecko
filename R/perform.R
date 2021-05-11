#' Use trained model to predict on data
#'
#' @param model (object) direct output from train_NBC
#' @param asv_table (double) rows= ASV names, cols=sample names
#' @param sample_data (character matrix) column 1 = sample names; column 2 = class labels (e.g. "asthmatic" and "healthy")
#' @param prob_class1 (integer) number samples of test class 1 (class 1 is alphabetically first, class 2 is next) / total number samples see "sample_classes" from train_NBC to check order
#' @export
perform_NBC = function(model, asv_table, sample_data, prob_class1) {
  savedParams_overall = model[[3]]
  savedParams_class1 = model[[1]]
  savedParams_class2 = model[[2]]
  subset_features = intersect(colnames(savedParams_overall), row.names(asv_table))
  #filter asv_table
  if(is.matrix(sample_data)){
    asv_table = asv_table[subset_features, sample_data[,1]]
  } else {
    asv_table = asv_table[subset_features, sample_data[1]]}
  #acquire class names
  sample_classes = model[[4]]
  ##Calculate P(X|Y = class1)
  p_given_class1 = p_giv_class(features = subset_features, asv_table = asv_table, param_table = savedParams_class1, sample_data = sample_data) # matrix cols=samples, rows=ASVs 90x95
  p_given_class2 = p_giv_class(features = subset_features, asv_table = asv_table, param_table = savedParams_class2, sample_data = sample_data) # matrix cols=samples, rows=ASVs 90x95
  p_x_independent_of_y = p_giv_class(features = subset_features, asv_table = asv_table, param_table = savedParams_overall, sample_data = sample_data) # matrix cols=samples, rows=ASVs 90x95
  ##Take the log of probabilities
  p_given_class1 = log(p_given_class1)
  p_given_class2 = log(p_given_class2)
  p_x_independent_of_y = log(p_x_independent_of_y)
  ##Perform log(P(X|Y)) - log(P(X))
  p_given_class1 = p_given_class1 - p_x_independent_of_y    # = log( P(class1 | X )
  p_given_class2 = p_given_class2 - p_x_independent_of_y    # = log( P(class2 | X )
  ##Get ASV asthma scores = log( P(class1 | X ) - log( P(class2 | X )
  #NOTE: class 1 is assumed as the case, class 2 is assumed as the control (we can make this user defined)
  feature_scores = p_given_class1 - p_given_class2
  ##Sum log probabilities for each sample to acquire log(P(Y|X)) - log(P(Y))
  p_1v2 = matrix(nrow = 4, ncol = dim(as.matrix(asv_table))[2])
  if(is.matrix(sample_data)){
    colnames(p_1v2) = colnames(asv_table)
  } else {
    colnames(p_1v2) = sample_data[1]
  }
  rownames(p_1v2) = c(sample_classes[1], sample_classes[2], "result", "truth")
  p_1v2[sample_classes[1],] = colSums(p_given_class1)
  p_1v2[sample_classes[2],] = colSums(p_given_class2)
  ##Finally, I just need to add the log(P(Y)) to each value in this table
  prob_class2 = 1 - prob_class1
  prob_class1 = log(prob_class1)
  prob_class2 = log(prob_class2)
  p_1v2[sample_classes[1],] = p_1v2[sample_classes[1],]+prob_class1
  p_1v2[sample_classes[2],] = p_1v2[sample_classes[2],]+prob_class2
  ##Solve argmax.
  p_1v2 = data.frame(p_1v2)
  if(is.matrix(sample_data)){
    p_1v2["truth",] = sample_data[,2]
  } else {
    p_1v2["truth",] = sample_data[2]
  }
  if(is.matrix(sample_data)){
    for(p in sample_data[,1]) {
      q = as.numeric(p_1v2[sample_classes, p])
      p_1v2["result",p] = sample_classes[which(q == max(q))]
    }
  } else {
    for(p in sample_data[1]) {
      q = as.numeric(p_1v2[sample_classes, p])
      p_1v2["result",p] = sample_classes[which(q == max(q))]
    }
  }

  classification_rate = sum(p_1v2["result",]==p_1v2["truth",])/ncol(p_1v2)
  output = list(classification_rate = as.numeric(classification_rate),
                sample_by_sample = p_1v2,
                model = model,
                score_per_feature_per_sample = feature_scores)
  return(output)
}
