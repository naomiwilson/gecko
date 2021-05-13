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
      # print(colnames(featureScoreMatrix)[i])
      countMatrixByROCThreshold[,i] <-
        as.numeric(featureScoreMatrix[,i] > asvThresholds[row.names(featureScoreMatrix),"threshold"])
      count = count +1
    }
    if (as.character(sample_table$truth[colnames(featureScoreMatrix)[i]]) == ctrl_class_name) {
      # print("this is healthy")
      # print(colnames(featureScoreMatrix)[i])
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

#' Get proportion or count of pairwise concordant taxa per dyad
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

#' Get Dyads needed to cover given list of features
#'
#' @param ASVs_to_cover list of features for which you'd like to get coverage in a gnotobiotic exptt
#' @param perform_NBC_output perform_NBC_output
#' @param sampleConcordanceMatrix output of getCountMatrixByROCThreshold
#' @param allHumanASVs filtered ASV Table used for model training
#' @param limit_of_detection limit_of_detection (default=0)
#' @param case_class_name case_class_name (default = "asthmatic")
#' @param ctrl_class_name ctrl_class_name (default = "healthy")
#' @return matrix of dyads best for full coverage of input feature list. Columns: Case (disease cohort) and Ctrl (healthy/control cohort)
#' @export
dyads_for_coverage <- function(ASVs_to_cover,  # list of names
                               perform_NBC_output,
                               sampleConcordanceMatrix,
                               allHumanASVs, # filtered asv table from NBC models
                               limit_of_detection=0,
                               case_class_name = "asthmatic",
                               ctrl_class_name = "healthy"){
  # Initiate a data frame:
  bestPairsForCoverage <- data.frame(Case=character(),
                                     Ctrl=character(),
                                     stringsAsFactors=FALSE)
  countDet <- 0  # set counter

  # drop asvs that are below LOD
  rare_and_lowRA = row.names(allHumanASVs)[rowSums(allHumanASVs < my_min_rel_abund) == dim(allHumanASVs)[2]]
  ASVs_to_cover = ASVs_to_cover[!(ASVs_to_cover %in% rare_and_lowRA)]

  while(length(ASVs_to_cover)>0){
    conc_matrix_list = get_concordance_matrices(perform_NBC_output = perform_NBC_output,
                                                sampleConcordanceMatrix=sampleConcordanceMatrix,
                                                featureAUCsForConcordanceList=ASVs_to_cover,
                                                allHumanASVs = allHumanASVs,
                                                limit_of_detection=limit_of_detection,
                                                case_class_name =case_class_name,
                                                ctrl_class_name = ctrl_class_name,
                                                PWCtype="proportion")
    best_pairs_tmp <- reshape2::melt(conc_matrix_list[[1]])
    best_a_tmp <- as.character(best_pairs_tmp[which.max(best_pairs_tmp$value),]$Var1) #asthmatic
    best_h_tmp <- as.character(best_pairs_tmp[which.max(best_pairs_tmp$value),]$Var2) #healthy

    winning_pair <- data.frame(conc_matrix_list[[2]][best_a_tmp])
    winning_pair <- winning_pair[paste0(best_a_tmp, ".", best_h_tmp)]

    countDet <- countDet + 1

    bestPairsForCoverage[countDet, "Case"] <- best_a_tmp
    bestPairsForCoverage[countDet, "Ctrl"] <- best_h_tmp

    asvs_covered <- row.names(winning_pair)[winning_pair > 0]
    ASVs_to_cover <- ASVs_to_cover[!(ASVs_to_cover %in% asvs_covered)]

    print(paste0("pair: ", best_a_tmp, ", ", best_h_tmp))
    print(ASVs_to_cover)

  }
  return(bestPairsForCoverage) # include all for testing
}
