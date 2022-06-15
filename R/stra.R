#' Estimate Population Mean
#'
#' This function is used to estimate the population mean in stratified simple random sampling.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the sample mean of each strata.
#' @param s2h a vector, the sample variance of each strata.
#' @param alpha a number, the confidence level.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
stra_srs_mean <- function(Nh = NULL,Wh = NULL,N = NULL,nh,ybarh,s2h,
                               alpha = 0.05,stra_est = FALSE){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      Wh <- Nh/sum(Nh)
    }
  }

  # the estimated results of each strata
  fh <- nh/Nh
  ybarh_var_est <- (1-fh)/nh*s2h
  if(stra_est){
    ybarh_se_est <- sqrt(ybarh_var_est)
    ybarh_ci_est <- conf_interval(ybarh,ybarh_se_est,alpha = alpha)
    ybarh_ci_left <- ybarh_ci_est$ci_left_bound
    ybarh_ci_right <- ybarh_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,ybarh,ybarh_var_est,
                                       ybarh_se_est,ybarh_ci_left,
                                       ybarh_ci_right))
  }
  # the overall estimated results
  ybar <- sum(Wh*ybarh)
  ybar_var_est <- sum(Wh^2*ybarh_var_est)
  ybar_se_est <- sqrt(ybar_var_est)
  ybar_ci_est <- conf_interval(ybar,ybar_se_est,alpha = alpha)
  ybar_ci_left <- ybar_ci_est$ci_left_bound
  ybar_ci_right <- ybar_ci_est$ci_right_bound
  mean_result <- matrix(c(ybar,ybar_var_est,ybar_se_est,ybar_ci_left,
                          ybar_ci_right),nrow = 1)
  colnames(mean_result) <- c("ybar","ybar_var_est","ybar_se_est",
                          "ybar_ci_left","ybar_ci_right")
  rownames(mean_result) <- "Stra"

  if(stra_est){
    return(list(stra_result = stra_result,
                mean_result = as.data.frame(mean_result)))
  }else{
    return(list(mean_result = as.data.frame(mean_result)))
  }
}

#' Sample Information Extraction
#'
#' This function is an auxiliary function of \code{stra_srs_mean} and
#' \code{stra_srs_prop}, which is used to extract some information from the sample.
#'
#' @param sample_data a list, each element denotes a strata sample.
#' @param type string, "mean" for \code{stra_srs_mean}，"prop" for \code{stra_srs_prop}.
stra_sample_extract <- function(sample_data,type = "mean"){
  nh <- sapply(sample_data, length)
  if(type == "mean"){
    ybarh <- sapply(sample_data, mean)
    s2h <- sapply(sample_data, var)
    names(ybarh) <- NULL
    names(s2h) <- NULL
    names(nh) <- NULL
    return(list(nh = nh,ybarh = ybarh,s2h = s2h))
  }
  if(type == "prop"){
    ph <- sapply(sample_data, mean)
    names(ph) <- NULL
    names(nh) <- NULL
    return(list(nh = nh,ph = ph))
  }
}


#' Sampling simulation
#'
#' A function for simulation calculations that draws samples from the population.
#'
#' @param pop_data a list, each element denotes a strata population.
#' @param n a number, the sample size one wishes.
#' @param w a vector, the allocation proportion of each stratum sample. It will be normalized.
#' If it is not given, then the function will use proportional allocation.
#' @note Since there is rounding error, the final sample size may not be \code{n}.
stra_sample_generate <- function(pop_data,n,wh = NULL){
  stra_num <- length(pop_data)
  Nh <- map_int(pop_data,nrow)
  if(is.null(wh)){
    wh <- Nh/sum(Nh)
  }else{
    wh <- wh/sum(wh)
  }
  size <- round(n*wh)
  re <- list(NULL)
  for (i in 1:stra_num) {
    index <- sample(1:Nh[i],size = size[i])
    re[[i]] <- pop_data[[i]][index,]
  }
  return(re)
}

#' Estimate Population Proportion
#'
#' This function is used to estimate The proportion of units in the population
#' that have a certain characteristic in stratified simple random sampling.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ph a vector, the sample proportion of each strata.
#' @param alpha a number, the confidence level.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
stra_srs_prop <- function(Nh = NULL,Wh = NULL,N = NULL,nh,ph,
                               alpha = 0.05,stra_est = FALSE){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      Wh <- Nh/sum(Nh)
    }
  }

  sing_index <- which(ph == 0)
  s2h <- nh/(nh-1)*ph*(1-ph)
  s2h[sing_index] = 0 # in case some strata doesn't have such unit
  fh <- nh/Nh
  # the estimation for each strata
  ph_var_est <- (1-fh)/nh*s2h
  if(stra_est){
    ph_se_est <- sqrt(ph_var_est)
    ph_ci_est <- conf_interval(ph,ph_se_est,alpha = alpha)
    ph_ci_left <- ph_ci_est$ci_left_bound
    ph_ci_right <- ph_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,ph,ph_var_est,
                                       ph_se_est,ph_ci_left,
                                       ph_ci_right))
  }
  # the overall estimation
  p <- sum(Wh*ph)
  p_var_est <- sum(Wh^2*ph_var_est)
  p_se_est <- sqrt(p_var_est)
  p_ci_est <- conf_interval(p,p_se_est,alpha = alpha)
  p_ci_left <- p_ci_est$ci_left_bound
  p_ci_right <- p_ci_est$ci_right_bound
  p_result <- matrix(c(p,p_var_est,p_se_est,p_ci_left,
                          p_ci_right),nrow = 1)
  colnames(p_result) <- c("p","p_var_est","p_se_est",
                             "p_ci_left","p_ci_right")
  rownames(p_result) <- "Stra"

  if(stra_est){
    return(list(stra_result = stra_result,
                p_result = as.data.frame(p_result)))
  }else{
    return(list(p_result = as.data.frame(p_result)))
  }
}


#' Sample Size Allocation
#'
#' This function is used to calculate the sample size of each strata under the given sample allocation method.
#' @param Nh a vector, the size of each population strata.
#' @param S2h a vector, the variance of each population strata, which can be
#' estimated by historical data. It can be ignored when \code{type} is "prop".
#' @param Ph a vector, the proportion of each population strata, which can be
#' estimated by historical data. It can be ignored when \code{type} is "mean".
#' @param ch a vector, the unit cost of each strata. It can be just the cost weights instead of specific values.
#' @param n a number, the sample size one wishes.
#' @param allocation string, the allocation method. "Prop" denotes proportional
#' allocation, "Opt" denotes optimal allocation, "Neyman" demotes Neyman allocation.
#' @param type string, the type of the parameter to be estimated,
#' "mean" for population mean and "prop" for population proportion.
stra_allocation <- function(Nh,S2h = NULL,Ph = NULL,ch = NULL,n,
                            allocation,type = "mean"){
  # proportional allocation requires only Nh
  if(allocation == "Prop"){
    nh <- n*Nh/sum(Nh)
  }else{
    if(type == "prop"){
      if(is.null(Ph)){
        stop("Please give Ph when type is 'prop'!")
      }
      S2h <- sqrt(Nh/(Nh-1)*Ph*(1-Ph))
    }else{
      if(is.null(S2h)){
        stop("Please give Ph when type is 'mean'!")
      }
    }
    if(allocation == "Opt"){
      if(is.null(ch)){
        stop("Please give ch when allocation is optimal allocation!")
      }
      nh <- n*Nh*sqrt(S2h)/sqrt(ch)/sum(Nh*sqrt(S2h)/sqrt(ch))
    }
    if(allocation == "Neyman"){
      nh <- n*Nh*sqrt(S2h)/sum(Nh*sqrt(S2h))
    }
  }
  wh <- nh/n
  nh <- round(nh)
  return(list(n = n,allocation = allocation,nh = nh,wh = wh,
              size = sum(nh)))
}


#' Estimate Sample size
#'
#' #' This function estimates the sample size required for a given mean/proportion estimator precision.
#' @param Nh a vector, the size of each population strata.
#' @param Wh a vector, the weight of each population strata. It can be ignored
#' since it is for the case N is infinity where \code{Nh} cannot be specified.
#' @param varh_his a vector, historical estimates of each strata population variance.
#' @param mean_his a vector, historical estimates of each strata population mean.
#' @param Ph_his a vector, historical estimates of each strata population proportion.
#' @param P_his a number, a vector, historical estimates of population prportion.
#' @param ch a vector, the unit cost of each strata. It can be just the cost weights instead of specific values.
#' @param method string, the form of precision requirement.
#' For upper limit of variance, use "V", for absolute error limit, use "d",
#' for relative error limit, use "r", for upper limit of variation coefficient,
#' use "CV".
#' @param bound a number, the upper limit value of the precision requirement.
#' @param alpha a number, the confidence level.
#' @return a list contain the form of precision requirement, uncorrected sample
#' size and corrected sample size.
#' @param allocation string, the allocation method. "Prop" denotes proportional
#' allocation, "Opt" denotes optimal allocation, "Neyman" demotes Neyman allocation.
#' @param type string, the type of the parameter to be estimated,
#' "mean" for population mean and "prop" for population proportion.
stra_size <- function(Nh,Wh = NULL,varh_his = NULL,mean_his = NULL,
                      Ph_his,ch = NULL,P_his = NULL,method,bound,
                      alpha = 0.05,allocation,type = "mean"){
  quan <- qnorm(1-alpha/2)
  N <- sum(Nh)
  if(is.null(Wh)){
    Wh <- Nh/N
  }

  if(type == "prop"){
    varh_his <- Nh/(Nh-1)*Ph_his*(1-Ph_his)
    if(method == "r" | method == "CV"){
      if(is.null(P_his)){
        stop("The current accuracy is specified as relative error limit or
             upper limit of variation coefficient, please give the historical
             population proportion！")
      }
      mean_his <- P_his
    }
  }

  # allocation proportion for proportional allocation
  if(allocation == "Prop"){
    wh <- Wh
  }
  # allocation proportion for optimal allocation
  if(allocation == "Opt"){
    if(is.null(varh_his)){
      stop("Please give historical strata variancewhen allocation is optimal
           allocation!")
    }
    if(is.null(ch)){
      stop("Please give strata unit cost when allocation is optimal
           allocation!")
    }
    wh <- Wh*sqrt(varh_his)/sqrt(ch)/sum(Wh*sqrt(varh_his)/sqrt(ch))
  }
  # allocation proportion for Neyman allocation
  if(allocation == "Neyman"){
    if(is.null(varh_his)){
      stop("Please give historical strata variancewhen allocation is Neyman
           allocation!")
    }
    wh <- Wh*sqrt(varh_his)/sum(Wh*sqrt(varh_his))
  }

  # precision requirement transform
  if(method == "d"){
    bound <- (bound/quan)^2
  }
  if(method == "r"){
    if(is.null(mean_his)){
      stop("The current accuracy is specified as relative error limit,
           please give the historical strata mean!")
    }
    bound <- (bound*mean_his/quan)^2
  }
  if(method == "CV"){
    if(is.null(mean_his)){
      stop("The current accuracy is specified as upper limit of variation coefficient,
           please give the historical strata mean!")
    }
    bound <- (bound*mean_his)^2
  }
  # calculate sample size and strata sample size
  n <- sum(Wh^2*varh_his/wh)/(bound+sum(Wh*varh_his)/N)
  nh <- ceiling(n*wh)
  return(list(n = n,nh = nh,wh = wh,size = sum(nh)))
}
