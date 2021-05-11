#' Estimate Binomial Parameter
#'
#' This function estimates probability of non-detection per taxon:
#' Number of samples below the limit of detection divided by
#' the total number of samples.
#' A pseudocount of 1/2 is also included.
#'
#' @param taxa
#' @param limit_of_detection numeric relative abundance value
#' @return numeric weight value
#' @export
estimate_binomial_parameter = function(taxa, limit_of_detection) {
  return((sum((taxa<=limit_of_detection))+0.5)/(length(taxa)+1))
}

#' Calculate P(B)
#'
#' @param x numeric : relative abundance of the feature of interest
#' @param feature string : feature of interest with representation in the provided parameter table
#' @param pararms_table data frame / matrix : parameter table containing P(not detected) and shape parameters for beta distribution:
#' @return numeric value
#' @export
calculate_probability = function(x, feature, params_table) {
  #x: [numeric] the relative abundance of the feature of interest
  #feature: [string] feature of interest with representation in the provided parameter table
  #pararms_table: [parameter x feature data.frame/matrix] parameter table containing P(not detected) and shape parameters for beta distribution.
  if(x == 0) {
    #Get P(not detected)
    return(params_table["pNotDetected",feature])
  }  else {
    #Prob detected from beta div
    return(c(1- params_table["pNotDetected",feature])*stats::dbeta(x = x, shape1 = params_table["shape1",feature], shape2 = params_table["shape2",feature]))
  }
}

#' @export
prior = function(alpha, mu, sig) {
  # Gaussian distribution
  res = (mu-alpha)/sig^2
  return(res)
}

#' Maximum a posteriori estimation
#'
#' @param taxa vector of relative abundances for a given taxon
#' @param sig_a Gaussian prior parameter for alpha (default=30)
#' @param sig_b Gaussian prior parameter for beta (default=100)
#' @export
getMapEstimates = function(taxa, sig_a = 30, sig_b = 100) {
  taxa = taxa[taxa>0]
  n = length(taxa)
  if (n == 0) {return(c(1,1))}
  mu_a = 1 # default starting place looks like presence vs absence - can make user-defined if they have a good idea of the presence/absence distribution
  mu_b = 1 # default starting place looks like presence vs absence
  i = 0
  alpha_old = stats::runif(1, min = 0, max = 2)
  beta_old = stats::runif(1, min = 0, max = 100)
  while(TRUE) {
    if (i>1000) {break}
    g_1 = digamma(alpha_old+beta_old) - digamma(alpha_old) + mean(log(taxa)) + prior(alpha_old, mu = mu_a, sig = sig_a)/n
    g_2 = digamma(beta_old+alpha_old) - digamma(beta_old) + mean(log(1-taxa)) + prior(beta_old, mu = mu_b, sig = sig_b)/n
    #hessian functions
    g_11 = trigamma(alpha_old+beta_old) - trigamma(alpha_old) - 1/sig_a^2/n
    g_12 = trigamma(alpha_old+beta_old)
    g_21 = trigamma(alpha_old+beta_old)
    g_22 = trigamma(alpha_old+beta_old) - trigamma(beta_old) - 1/sig_b^2/n
    #Solve
    matrix_a =  matrix(data = c(g_11, g_21, g_12, g_22), 2, 2)
    #Ensure Solution is Not Singular
    matrix_a = try(solve(matrix_a), silent = TRUE)
    if (class(matrix_a) == "try-error") {
      alpha_old = stats::runif(1, min = 0, max = 2)
      beta_old = stats::runif(1, min = 0, max = 100)
      i = i+1
      next
    }
    vector_b = c(g_1, g_2)
    results = c(alpha_old, beta_old) - matrix_a%*%c(g_1, g_2)
    alpha_new = results[1]
    beta_new = results[2]
    if (alpha_new<0 | beta_new<0) {
      alpha_old = stats::runif(1, min = 0, max = 2)  # max values were chosen by Ariel to keep alpha from exploding (another thing that could be user-defined)
      beta_old = stats::runif(1, min = 0, max = 100)
      i = i+1
      next
    }
    convergence = all(abs(results - c(alpha_old, beta_old))/results<0.001) # a chosen error value that could be subject to optimization
    if (convergence) {
      output = c(alpha_new, beta_new)
      names(output) = c("shape1", "shape2")
      return(output)
    } else {
      alpha_old = alpha_new
      beta_old = beta_new
      i = i+1
      next
    }
  }
  return(c("Failure", alpha_new, beta_new, alpha_old, beta_old))
}


#' #General Acquire Parameters Function
#' @export
acquire_parameters = function(asv_table, limit_of_detection) {
  pNotDetected = apply(X = asv_table, MARGIN = 1, estimate_binomial_parameter, limit_of_detection)
  betaParams = apply(X = asv_table, MARGIN = 1, getMapEstimates)
  savedParams = rbind(pNotDetected, betaParams)
  row.names(savedParams) = c("pNotDetected", "shape1", "shape2")
  return(savedParams)
}

#' General Calculate P(X|Y) function)
#' @export
p_giv_class = function(features, asv_table, param_table, sample_data) {
  p_given_class = matrix(nrow = length(features), ncol = dim(as.matrix(asv_table))[2])
  if (dim(as.matrix(asv_table))[2]>1){
    colnames(p_given_class) = colnames(asv_table)
  } else {
    colnames(p_given_class) = sample_data[1]
  }
  rownames(p_given_class) = features
  for (feature in features)  {
    if (length(sample_data) > 2) {
      for(sample in colnames(asv_table)) {
        p_given_class[feature,sample] = calculate_probability(x = asv_table[feature,sample],
                                                              feature = feature,
                                                              params_table = param_table)
      }
    } else {
      p_given_class[feature,sample_data[1]] = calculate_probability(x = asv_table[feature],
                                                                    feature = feature,
                                                                    params_table = param_table)
    }
  }
  return(p_given_class)
}

#' @export
make_confusion_matrix = function(p1v2) {
  final = as.factor(p1v2["result",])
  truth = as.factor(p1v2["truth",])
  output = caret::confusionMatrix(final, truth)
}

#' @export
get_modeled_features = function(model) {
  modeled_features = colnames(model[[1]])
  return(modeled_features)
}

#' NBC Training Function
#'   # "asv_table" = double
#'         rows= ASV names, cols=sample names
#' "sample_data" = character matrix
#'         column 1 = sample names; column 2 = class labels (e.g. "asthmatic" and "healthy")
#' "minimum_detection" = integer
#'         is how many samples in which the ASV has to be present in to be included in the training
#' "min_rel_abund" = integer
#'         sequencing detection limit - the relative abundance value at which we no longer trust the ASV is truly present
#'
#'
#' @export
train_NBC = function(asv_table, sample_data, minimum_detection, min_rel_abund) {
  #filter asv_table
  asv_table = asv_table[,sample_data[,1]]
  # print(dim(asv_table))
  asv_table = asv_table[rowSums(asv_table>min_rel_abund)>=minimum_detection,]
  #acquire class names
  sample_classes = sort(unique(sample_data[,2])) # class 1 is always alphabetically first
  #acquire overall parameters
  savedParams_overall = acquire_parameters(asv_table = asv_table, limit_of_detection = min_rel_abund)
  #acquire parameters for class 1
  sample_data_1 = sample_data[sample_data[,2]==sample_classes[1],]
  asv_table_1 = asv_table[,sample_data_1[,1]]
  # asv_table_1 = asv_table_1[rowSums(asv_table_1>0)>=minimum_detection,] # getting rid of this for md=7
  savedParams_class1 = acquire_parameters(asv_table = asv_table_1, limit_of_detection = min_rel_abund)
  #acquire parameters for class 2
  sample_data_2 = sample_data[sample_data[,2]==sample_classes[2],]
  asv_table_2 = asv_table[,sample_data_2[,1]]
  # asv_table_2 = asv_table_2[rowSums(asv_table_2>0)>=minimum_detection,] # getting rid of this for md=7
  savedParams_class2 = acquire_parameters(asv_table = asv_table_2, limit_of_detection = min_rel_abund)
  #Identify shared parameters
  class_1_features = colnames(savedParams_class1)
  class_2_features = colnames(savedParams_class2)
  overall_features = colnames(savedParams_overall)
  subset_features = intersect(class_1_features, class_2_features)
  #Filter down to shared parameters
  savedParams_overall = savedParams_overall[,subset_features]
  savedParams_class1 = savedParams_class1[,subset_features]
  savedParams_class2 = savedParams_class2[,subset_features]
  prob_class1 = sum(sample_data[,2]==sample_classes[1])/nrow(sample_data)
  model = list(savedParams_class1, savedParams_class2, savedParams_overall, sample_classes, prob_class1)
  return(model)
}
