#' Make Confusion Matrix of Results
#' @param p1v2 second matrix in output of perform_NBC aka "sample_by_sample"
#' @export
make_confusion_matrix = function(p1v2) {
  final = as.factor(p1v2["result",])
  truth = as.factor(p1v2["truth",])
  output = caret::confusionMatrix(final, truth)
  return(output)
}

#' Get list of features included in model
#' @param train_NBC_output output from train_NBC
#' @export
get_modeled_features = function(train_NBC_output) {
  modeled_features = colnames(train_NBC_output[[1]])
  return(modeled_features)
}

#' Get feature AUCs using pROC
#'
#' @param train_NBC_output train_NBC output
#' @param perform_NBC_output perform_NBC output
#' @param my_sample_data character table col1=sample names, col2=class names
#' @return table rownames=feature names, col "feature_aucs" containing numeric values from 0.5 to 1
#' @export
getFeatureAUCs <- function(train_NBC_output, perform_NBC_output, my_sample_data){
  feature_aucs <- c()
  feature_names_list <- c()
  for (asv in row.names(perform_NBC_output$score_per_feature_per_sample)){
    scores_by_feature <- data.frame(t(perform_NBC_output$score_per_feature_per_sample))
    scores_by_feature$class = my_sample_data[,2][my_sample_data[,1]==row.names(scores_by_feature)]
    scores_by_feature$binaryClass = ifelse(scores_by_feature$class==train_NBC_output[[4]][1], 1, 0) # make binary: class1 = 1 (case), class2 = 0 (control)
    roc_obj <- pROC::roc(scores_by_feature$binaryClass, eval(parse(text=paste("scores_by_feature$", asv, sep = ""))))
    feature_aucs <- c(feature_aucs, pROC::auc(roc_obj))
    feature_names_list <- c(feature_names_list, asv)
  } # don't be alarmed, there will be a lot of print-outs...
  feature_aucs <- data.frame(feature_aucs, row.names = feature_names_list)
  return(feature_aucs)
}

#' Histograms of sample scores by class (makes 3 plots)
#' @param scores.df dataframe with column named "score" containing numeric score values and "class" column containing class names per sample by row
#' @param donor0 control class donor sample name (e.g. healthyPatient32)
#' @param donor1 case class donor sample name (e.g. asthmaPatient19)
#' @param title0 title for control plot
#' @param title1 title for case plot
#' @param class0.color hex value, or R color name
#' @param class1.color hex value, or R color name
#' @param both.color hex value, or R color name
#' @param x.label x-axis label
#' @param numbins breaks for hist
#' @return 3 plots inline
#' @export
make_sample_score_histograms <- function(scores.df,
                                         donor0 = "MARS0022.F.16s", donor1 = "MARS0043.F.16s",
                                         title0 = "Healthy (59) Sample Scores",
                                         title1 = "Asthmatic (36) Sample Scores",
                                         class1.color = "#DE4968FF", class0.color = "#51127CFF", both.color = "lightblue",
                                         x.label = "Sample Score",
                                         numbins=10)
{
  hist(scores.df$score, main="ALL Sample Scores", xlab = x.label, col = "lightblue", breaks=numbins)
  abline(v = scores.df[donor1,]$score, col = class1.color, lwd = 2, lty="dashed")
  abline(v = scores.df[donor0,]$score, col = class0.color, lwd=2, lty="dashed")

  hist(scores.df$score[scores.df$class == 1],
       main=title1, xlab = x.label,
       breaks =15, col = class1.color)
  abline(v = scores.df[donor1,]$score, col = "black", lwd = 2, lty="dashed")
  median(scores.df$score[scores.df$class == 1])

  hist(scores.df$score[scores.df$class == 0],
       main=title0, xlab = x.label,
       breaks = 20, col = class0.color)
  abline(v = scores.df[donor0,]$score, col = "black", lwd=2, lty="dashed")
}

#' Make Histogram of sample scores with both classes overlapping
#'
#' @param scores.df dataframe with column named "score" containing numeric score values and "class" column containing class names per sample by row
#' @param donor0 control class donor sample name (e.g. healthyPatient32)
#' @param donor1 case class donor sample name (e.g. asthmaPatient19)
#' @param class0.color hex value, or R color name
#' @param class1.color hex value, or R color name
#' @param x.label x-axis label
#' @param y.label y-axis label
#' @param plot.title for ggtitle
#' @param numbins breaks for hist
#' @return ggplot
#' @export
make_overlapping_score_histogram <- function(scores.df,
                                             donor0 = "MARS0022.F.16s", donor1 = "MARS0043.F.16s",
                                             class1.color = "#DE4968FF", class0.color = "#51127CFF",
                                             x.label = "Sample Score", y.label = "Frequency", plot.title="Score Histogram",
                                             numbins = 30){
    ggplot2::ggplot(scores.df, aes(x = score)) +
    geom_histogram(aes(color = class, fill = class),
                   position = "identity", bins = numbins, alpha = 0.7) +
    scale_color_manual(values = c(class1.color, class0.color)) +
    scale_fill_manual(values = c(class1.color, class0.color)) +
    theme_classic(
    ) +
    ggtitle(plot.title)+
    scale_x_continuous(name=x.label) +
    scale_y_continuous(name=y.label) +
    geom_vline(xintercept = scores.df[donor1,]$score,
               linetype="dashed",
               color = class1.color, size=1.2) +
    geom_vline(xintercept = scores.df[donor0,]$score,
               linetype="dashed",
               color = class0.color, size=1.2)
}

#' Log Loss as defined on kaggle forum
#'
#' @param act actual value
#' @param pred predicted value
#' @export
LogLoss<-function(act, pred)
{
  eps = 1e-15;
  nr = length(pred)
  pred = matrix(sapply( pred, function(x) max(eps,x)), nrow = nr)
  pred = matrix(sapply( pred, function(x) min(1-eps,x)), nrow = nr)
  ll = sum(act*log(pred) + (1-act)*log(1-pred))
  ll = ll * -1/(length(act))
  return(ll);
}

#' Grab best thresholds per feature from ROC curve
#'
#' @param sample_data character matrix: column 1 = sample names; column 2 = class labels (e.g. "asthmatic" and "healthy")
#' @param perform_NBC_output entire output from perform_NBC
#' @param train_NBC_output entire output from train_NBC
#' @param featureAUCsRanked vector of AUCs in rank order
#' @param idx1 starting rank for resulting list
#' @param idx2 ending rank for resulting list
#' @export
getBestThresholds <- function(sample_data, perform_NBC_output, train_NBC_output, featureAUCsRanked, idx1, idx2) {
  # asvThresholds = getBestThresholds(idx1 = 1, idx2 = 392,
  #                     featureAUCsRanked = featureAUCsForHist,
  #                     sample_data = my_sample_data,
  #                     perform_NBC_output = my_output,
  #                     train_NBC_output = my_model)
  # transform scores matrix:
  scores_by_feature <- data.frame(t(perform_NBC_output$score_per_feature_per_sample))
  scores_by_feature$class = sample_data[,2][sample_data[,1]==row.names(scores_by_feature)]
  scores_by_feature$binaryClass = ifelse(scores_by_feature$class==train_NBC_output[[4]][1], 1, 0) # make binary: class1 = 1 (case), class2 = 0 (control)
  thresholds = seq(min(scores_by_feature[,1:(length(names(scores_by_feature))-2)])-1,
                   max(scores_by_feature[,1:(length(names(scores_by_feature))-2)])+1,
                   by = 0.1)
  bestThresholds = data.frame(row.names = row.names(featureAUCsRanked[idx1:idx2,]))
  for (asv in row.names(featureAUCsRanked[idx1:idx2,])){
    roc_obj <- roc(scores_by_feature$binaryClass, eval(parse(text=paste("scores_by_feature$", asv, sep = ""))))
    AUC = auc(roc_obj)
    # TPR=rev(roc_obj$sensitivities)
    # FPR=rev(1 - roc_obj$specificities)
    # labels=roc_obj$response
    # scores=roc_obj$predictor
    thresholdBest = coords(roc_obj, "best", ret = "threshold")
    bestThresholds[asv,1] = thresholdBest
    bestThresholds[asv,2] = AUC
  }
  return(bestThresholds)
}

#' Get Concordance using ROC NB Score Thresholds
#'
#' @param perform_NBC_output NBC output
#' @param asvThresholds best thresholds table with column "threshold", rownames feature names
#' @param ctrl_class_name default "healthy"
#' @param case_class_name default "asthmatic"
#' @return binary matrix features vs samples 0=not concordant, 1=concordant, meaning healthy sample scored well for healthy for just that ASV value
#' @export
getCountMatrixByROCThreshold <- function(perform_NBC_output,
                                         asvThresholds,
                                         ctrl_class_name = "healthy",
                                         case_class_name = "asthmatic"){
  featureScoreMatrix = perform_NBC_output$score_per_feature_per_sample
  sample_table = data.frame(t(perform_NBC_output$sample_by_sample))
  countMatrixByROCThreshold <- matrix(nrow = dim(featureScoreMatrix)[1], ncol = dim(featureScoreMatrix)[2])
  count = 0
  for (i in 1:dim(featureScoreMatrix)[2]) {
    if (as.character(sample_table$truth[colnames(featureScoreMatrix)[i]]) == case_class_name) {
      # print("this is asthmatic")
      print(colnames(featureScoreMatrix)[i])
      countMatrixByROCThreshold[,i] <-
        as.numeric(featureScoreMatrix[,i] > asvThresholds[row.names(featureScoreMatrix),"threshold"])
      count = count +1
    }
    if (as.character(sample_table$truth[colnames(featureScoreMatrix)[i]]) == ctrl_class_name) {
      # print("this is healthy")
      print(colnames(featureScoreMatrix)[i])
      countMatrixByROCThreshold[,i] <-
        as.numeric(featureScoreMatrix[,i] < asvThresholds[row.names(featureScoreMatrix),"threshold"])
      count = count +1
    }
    if (as.character(sample_table$truth[colnames(featureScoreMatrix)[i]]) != ctrl_class_name && as.character(sample_table$truth[colnames(featureScoreMatrix)[i]]) != case_class_name) {break}
  }
  if(count != dim(featureScoreMatrix)[2]){break}
  colnames(countMatrixByROCThreshold) = colnames(featureScoreMatrix)
  row.names(countMatrixByROCThreshold) = row.names(featureScoreMatrix)
  return(countMatrixByROCThreshold)
}

#' Get proportion of pairwise concordant taxa per dyad
#'
#' @param mat1_i integer column i for matrix 1
#' @param mat2_j integer column j for matrix 2
#' @param asv_table features x samples: relative abundance
#' @param feature_list list of features to include
#' @param mat1 class 1 binary sample concordance matrix (based on ROC thresholds from getCountMatrixByROCThreshold)
#' @param mat2 class 2 binary sample concordance matrix (based on ROC thresholds from getCountMatrixByROCThreshold)
#' @param limit_of_detection relative abundance lower cutoff
#' @param type character option specifying proportion "proportion" or absolute "absolute" counts
#' @export
count_concordant_pairs <- function(mat1_i, mat2_j, asv_table, feature_list, mat1, mat2, limit_of_detection, type="proportion") {
  # don't count any where the ASV is ND in both:
  EITHERPRESENT = as.numeric(asv_table[feature_list,colnames(mat1)[mat1_i]] > limit_of_detection |
                               asv_table[feature_list,colnames(mat2)[mat2_j]] > limit_of_detection)
  BOTHCONCORDANT = as.numeric((mat1[feature_list,mat1_i] + mat2[feature_list,mat2_j]) == 2)
  if(type=="proportion"){
    return(sum(EITHERPRESENT*BOTHCONCORDANT)/sum(EITHERPRESENT))
  }
  if(type=="absolute"){
    return(sum(EITHERPRESENT*BOTHCONCORDANT))
  }
}

#' Get full matrix of pairwise concordant features between 2 samples
#'
#' @param mat1_i integer column i for matrix 1
#' @param mat2_j integer column j for matrix 2
#' @param asv_table features x samples: relative abundance
#' @param feature_list list of features to include
#' @param mat1 class 1 binary sample concordance matrix (based on ROC thresholds from getCountMatrixByROCThreshold)
#' @param mat2 class 2 binary sample concordance matrix (based on ROC thresholds from getCountMatrixByROCThreshold)
#' @param limit_of_detection relative abundance lower cutoff
#' @export
get_pairwise_concordant_features <- function(mat1_i, mat2_j, asv_table, feature_list, mat1, mat2, limit_of_detection) {
  # don't count any where the ASV is ND in both:
  EITHERPRESENT = as.numeric(asv_table[feature_list,colnames(mat1)[mat1_i]] > limit_of_detection | asv_table[feature_list,colnames(mat2)[mat2_j]] > limit_of_detection)
  BOTHCONCORDANT = as.numeric((mat1[feature_list,mat1_i] + mat2[feature_list,mat2_j]) == 2)
  concordant_list = EITHERPRESENT*BOTHCONCORDANT
  return(concordant_list)
}

#' Get Pairwise Concordance
#'
#' @param perform_NBC_output entire output from perform_NBC
#' @param sampleConcordanceMatrix binary matrix of all samples per taxon notifying sample concordance
#' @param featureAUCsForConcordanceList list of feature names: subset of ASVs we want to use to count pairwise concordance (top X, or AUC>0.6)
#' @param allHumanASVs ASV table for all samples
#' @param limit_of_detection relative abundance cutoff
#' @param PWCtype type of pairwise concordance count (proportion or absolute)
#' @return list of matrices 1) proportion of pairwise concordant taxa per pair and 2) multidimensional binary matrix- one per case sample
#' @export
get_concordance_matrices <- function(perform_NBC_output,
                                     sampleConcordanceMatrix,
                                     featureAUCsForConcordanceList,
                                     allHumanASVs,
                                     limit_of_detection,
                                     case_class_name = "asthmatic",
                                     ctrl_class_name = "healthy",
                                     PWCtype="proportion"){
  sample_table=data.frame(t(perform_NBC_output$sample_by_sample))
  if(length(featureAUCsForConcordanceList) == 1){
    countMatrixByROCThresholdAsthmatics <-
      t(as.matrix(sampleConcordanceMatrix[featureAUCsForConcordanceList,
                                          row.names(sample_table[sample_table$truth == case_class_name,])]))
    row.names(countMatrixByROCThresholdAsthmatics) <- featureAUCsForConcordanceList
    countMatrixByROCThresholdHealthy <-
      t(as.matrix(sampleConcordanceMatrix[featureAUCsForConcordanceList,
                                          row.names(sample_table[sample_table$truth == ctrl_class_name,])]))
    row.names(countMatrixByROCThresholdHealthy) <- featureAUCsForConcordanceList
  } else{
    countMatrixByROCThresholdAsthmatics <-
      sampleConcordanceMatrix[featureAUCsForConcordanceList,
                              row.names(sample_table[sample_table$truth == case_class_name,])]
    countMatrixByROCThresholdHealthy <-
      sampleConcordanceMatrix[featureAUCsForConcordanceList,
                              row.names(sample_table[sample_table$truth == ctrl_class_name,])]

  }
  countMatrixByROCThresholdPairs <- matrix(nrow = dim(countMatrixByROCThresholdAsthmatics)[2],
                                           ncol = dim(countMatrixByROCThresholdHealthy)[2])
  row.names(countMatrixByROCThresholdPairs) = colnames(countMatrixByROCThresholdAsthmatics)
  colnames(countMatrixByROCThresholdPairs) = colnames(countMatrixByROCThresholdHealthy)
  countMatrixByROCThresholdPairsVerbose <- list()
  for (i in 1:dim(countMatrixByROCThresholdAsthmatics)[2]) {
    for (j in 1:dim(countMatrixByROCThresholdHealthy)[2]){
      countMatrixByROCThresholdPairs[i, j] <- count_concordant_pairs(mat1_i=i,
                                                                     mat2_j=j,
                                                                     asv_table = allHumanASVs,
                                                                     feature_list = featureAUCsForConcordanceList,
                                                                     mat1 = countMatrixByROCThresholdAsthmatics,
                                                                     mat2 = countMatrixByROCThresholdHealthy,
                                                                     limit_of_detection = limit_of_detection,
                                                                     type = PWCtype)
      # save which asvs these are:
      if (j == 1){
        countMatrixByROCThresholdPairsVerbose[[colnames(countMatrixByROCThresholdAsthmatics)[i]]] <-
          matrix(nrow = dim(countMatrixByROCThresholdAsthmatics)[1],
                 ncol = dim(countMatrixByROCThresholdHealthy)[2])
        row.names(countMatrixByROCThresholdPairsVerbose[[colnames(countMatrixByROCThresholdAsthmatics)[i]]]) =
          row.names(countMatrixByROCThresholdAsthmatics)
        colnames(countMatrixByROCThresholdPairsVerbose[[colnames(countMatrixByROCThresholdAsthmatics)[i]]]) =
          colnames(countMatrixByROCThresholdHealthy)
      }
      countMatrixByROCThresholdPairsVerbose[[colnames(countMatrixByROCThresholdAsthmatics)[i]]][,colnames(countMatrixByROCThresholdHealthy)[j]] <-
        get_pairwise_concordant_features(mat1_i=i,
                                  mat2_j=j,
                                  asv_table = allHumanASVs,
                                  feature_list = featureAUCsForConcordanceList,
                                  mat1 = countMatrixByROCThresholdAsthmatics,
                                  mat2 = countMatrixByROCThresholdHealthy,
                                  limit_of_detection = limit_of_detection)
    }
  }
  return(list(countMatrixByROCThresholdPairs, countMatrixByROCThresholdPairsVerbose))
}

