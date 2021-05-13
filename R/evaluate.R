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
  }
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

#' Plot a feature ROC Curve
#'
#' @param sample_data character table col1=sample names, col2=class names
#' @param perform_NBC_output perform_NBC output
#' @param train_NBC_output train_NBC output
#' @param featureAUCsRanked vector of AUCs in rank order
#' @param idx1 integer rank of first feature to plot
#' @param idx2 integer rank of last feature to plot
#' @param title plot title
#' @param want_label Boolean value - TRUE will print feature name and AUC on the plot (can only fit up to ~6)
#' @export
PlotAllROCByFeature <- function(sample_data, perform_NBC_output, train_NBC_output, featureAUCsRanked, idx1, idx2, title, want_label=TRUE) {
  # transform scores matrix:
  scores_by_feature <- data.frame(t(perform_NBC_output$score_per_feature_per_sample))
  scores_by_feature$class = sample_data[,2][sample_data[,1]==row.names(scores_by_feature)]
  scores_by_feature$binaryClass = ifelse(scores_by_feature$class==train_NBC_output[[4]][1], 1, 0) # make binary: class1 = 1 (case), class2 = 0 (control)
  # get lotsa unique colors
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # Plot em
  plot(1, type="n", xlab="FPR", ylab="TPR", xlim=c(0, 1), ylim=c(0, 1), main=title)
  thresholds = seq(min(scores_by_feature[,1:(length(names(scores_by_feature))-2)])-1, max(scores_by_feature[,1:(length(names(scores_by_feature))-2)])+1, by = 0.1)
  # for (asv in names(scores_by_feature[,1:(length(names(scores_by_feature))-2)])){
  x_coord = 0.9
  y_coord = 0.75
  for (asv in row.names(featureAUCsRanked[idx1:idx2,])){
    roc_obj <- pROC::roc(scores_by_feature$binaryClass, eval(parse(text=paste("scores_by_feature$", asv, sep = ""))))
    AUC = pROC::auc(roc_obj)
    TPR=rev(roc_obj$sensitivities)
    FPR=rev(1 - roc_obj$specificities)
    labels=roc_obj$response
    scores=roc_obj$predictor
    label_color = sample(col_vector,1)
    thresholdBest = pROC::coords(roc_obj, "best", ret = "threshold", transpose = FALSE)
    lines(FPR, TPR, type="l", col=label_color)
    if (want_label == TRUE){
      if (sum(row.names(featureAUCsRanked[idx1:idx2,]) %in% asv)==1) {
        text(x_coord, y_coord, asv, col = label_color)
        text(1-x_coord, 1-y_coord+0.4, round(AUC, digits = 2), col = label_color)
        # text(1-x_coord+0.2, 1-y_coord+0.4, round(thresholdBest, digits = 2), col = label_color)
        x_coord = x_coord
        y_coord = y_coord - 0.075
      }
    }
  }
  abline(a = 0, b = 1, col="black", lw = 3, lty='dashed') # add AUC=0.5 line
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
