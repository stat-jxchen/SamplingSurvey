#' Confidence Interval
#'
#' This function is used to estimate the confidence interval for the parameter.
#'
#' @param stat_obs a number, the the observed value of the estimator.
#' @param se_est a number, the standard error estimate of the estimator.
#' @param alpha a number, the confidence level.
#' @return A list containing:\tabular{ll}{
#'    \code{absolute_error} \tab absolute error limit. \cr
#'    \tab \cr
#'    \code{relative_error} \tab relative error limit. \cr
#'    \tab \cr
#'    \code{ci_left_bound} \tab lower confidence bound. \cr
#'    \tab \cr
#'    \code{ci_right_bound} \tab upper confidence bound. \cr
#' }
#' @examples
#' x <- rnorm(100)
#' xbar <- mean(x)
#' se <- sqrt(var(x)/100)
#' conf_interval(xbar,se,alpha = 0.05)
conf_interval <- function(stat_obs,se_est,alpha = 0.05){
  quan <- qnorm(1-alpha/2)
  # calculate estimates of the absolute error limit (d) and the relative error limit (r)
  d <- quan*se_est
  r <- d/stat_obs
  # two endpoints of the confidence interval
  left_bound <- stat_obs-d
  right_bound <- stat_obs+d
  # return
  result <- list(absolute_error = d,relative_error = r,
                 ci_left_bound = left_bound,ci_right_bound = right_bound)
  return(result)
}



#' Estimate Population Mean
#'
#' This function is used to estimate the population mean in simple random sampling.
#'
#' @param sample_data a vector, the sample data.
#' @param N a number, the size of population, If it is unknown or very large, it can be omitted.
#' By default, infinity is used for calculation.
#' @param alpha a number, the confidence level.
#' @return A list containing:\tabular{ll}{
#'    \code{ybar_result} \tab a list containing estimation of the population mean, estimation of the variance of the estimator and the estimation of the standard error of the estimator. \cr
#'    \tab \cr
#'    \code{ci_est} \tab a list contain absolute/relative error limit and the estimation of lower/upper confidence bound of the estimator.\cr
#'    \tab \cr
#' }
srs_mean <- function(sample_data,N = Inf,alpha = 0.05){
  n <- length(sample_data) # sample size
  f <- n/N # sample ratio
  if(f>1|f<0){
    stop("The sampling ratio is greater than 1 or less than 0,
         please check that the population size N is entered correctly!")
  }
  ybar <- mean(sample_data)
  s2 <- var(sample_data)
  ybar_var_est <- (1-f)*s2/n
  ybar_se_est <- sqrt(ybar_var_est) # estimated se of ybar
  # CI
  ci_est <- conf_interval(ybar,ybar_se_est,alpha)
  ybar_result <- list(ybar = ybar,ybar_var_est = ybar_var_est,
                      ybar_se_est = ybar_se_est)
  result <- c(ybar_result,ci_est)
  return(result)
}

#' Estimate Population Total
#'
#' This function is used to estimate the population total in simple random sampling.
#'
#' @param sample_data a vector, the sample data.
#' @param N a number, the size of population.
#' @param alpha a number, the confidence level.
#' @return A list containing:\tabular{ll}{
#'    \code{ytot_result} \tab a list containing estimation of the population total, estimation of the variance of the estimator and the estimation of the standard error of the estimator. \cr
#'    \tab \cr
#'    \code{ci_est} \tab a list contain absolute/relative error limit and the estimation of lower/upper confidence bound of the estimator.\cr
#'    \tab \cr
#' }
srs_total <- function(sample_data,N,alpha = 0.05){
  n <- length(sample_data)
  f <- n/N
  if(f>1|f<0){
    stop("The sampling ratio is greater than 1 or less than 0,
         please check that the population size N is entered correctly!")
  }
  ybar <- mean(sample_data)
  s2 <- var(sample_data)
  ytot <- N*ybar
  ytot_var_est <- N^2*(1-f)*s2/n
  ytot_se_est <- sqrt(ytot_var_est) # estimated se of ytot
  # CI
  ci_est <- conf_interval(ytot,ytot_se_est,alpha)
  ytot_result <- list(ytot = ytot,ytot_var_est = ytot_var_est,
                      ytot_se_est = ytot_se_est)
  result <- c(ytot_result,ci_est)
  return(result)
}

#' Estimate Population Proportion
#'
#' This function is used to estimate The proportion of units in the population that have a certain characteristic.
#'
#' @param event_num a number, the number of units in the sample with a certain characteristic.
#' @param n a number, the sample size
#' @param N a number, the size of population, If it is unknown or very large, it can be omitted.
#' By default, infinity is used for calculation.
#' @param alpha a number, the confidence level.
#' @return A list containing:\tabular{ll}{
#'    \code{p_result} \tab a list containing estimation of the population proportion, estimation of the variance of the estimator and the estimation of the standard error of the estimator. \cr
#'    \tab \cr
#'    \code{ci_est} \tab a list contain absolute/relative error limit and the estimation of lower/upper confidence bound of the estimator.\cr
#'    \tab \cr
#' }
srs_prop <- function(event_num,n,N = Inf,alpha = 0.05){
  f <- n/N
  if(f>1|f<0){
    stop("The sampling ratio is greater than 1 or less than 0,
         please check that the population size N is entered correctly!")
  }
  if(n<event_num){
    stop("event_num can not exceed n!")
  }
  p <-  event_num/n # the estimation of the proportion
  p_var_est <- (1-f)/(n-1)*p*(1-p)
  p_se_est <- sqrt(p_var_est) # estimated se of p
  # CI
  ci_est <- conf_interval(p,p_se_est,alpha)
  p_result <- list(P_est = p,p_var_est = p_var_est,p_se_est = p_se_est)
  result <- c(p_result,ci_est)
  return(result)
}

#' Estimate the Number of Units with a Certain Type
#'
#' This function is used to estimate The number of units in the population that have a certain characteristic.
#'
#' @param event_num a number, the number of cells in the sample with a certain characteristic.
#' @param n a number, the sample size
#' @param N a number, the size of population.
#' @param alpha a number, the confidence level.
#' @return A list containing:\tabular{ll}{
#'    \code{a_result} \tab a list containing estimation of the number of units in the population that have a certain characteristic, estimation of the variance of the estimator and the estimation of the standard error of the estimator. \cr
#'    \tab \cr
#'    \code{ci_est} \tab a list contain absolute/relative error limit and the estimation of lower/upper confidence bound of the estimator.\cr
#'    \tab \cr
#' }
srs_num <- function(event_num,n,N,alpha = 0.05){
  f <- n/N
  if(f>1|f<0){
    stop("The sampling ratio is greater than 1 or less than 0,
         please check that the population size N is entered correctly!")
  }
  if(n<event_num){
    stop("event_num can not exceed n!")
  }
  p <-  event_num/n
  a <- N*p # the estimation of the number
  a_var_est <- N^2*(1-f)/(n-1)*p*(1-p)
  a_se_est <- sqrt(a_var_est) # estimated se of a
  # CI
  ci_est <- conf_interval(a,a_se_est,alpha)
  a_result <- list(A_est = a,a_var_est = a_var_est,a_se_est = a_se_est)
  result <- c(a_result,ci_est)
  return(result)
}

#' Sample Size for Mean
#'
#' This function estimates the sample size required for a given mean estimator precision.
#'
#' @param N a number, the size of population, If it is unknown or very large, it can be omitted.
#' By default, infinity is used for calculation.
#' @param mean_his a number, historical estimates of population mean. It need be
#' specified if the precision requirement is a relative error limit or an upper
#' coefficient of variation limit.
#' @param var_his a number, historical estimates of population variation.
#' @param method string, the form of precision requirement.
#' For upper limit of variance, use "V", for absolute error limit, use "d",
#' for relative error limit, use "r", for upper limit of variation coefficient,
#' use "CV".
#' @param bound a number, the upper limit value of the precision requirement.
#' @param alpha a number, the confidence level.
#' @return a list contain the form of precision requirement, uncorrected sample
#' size and corrected sample size.
#' @note If the total estimator precision is given, convert it to the mean estimator
#' precision then use this function. This function will convert all forms of upper precision to
#' upper limit of variance in calculation.
srs_size_for_mean <- function(var_his,mean_his = NULL,N = Inf,
                              method,bound,alpha = 0.05){

  quan <- qnorm(1-alpha/2)

  if(method == "d"){
    bound <- (bound/quan)^2
  }
  if(method == "r"){
    if(is.null(mean_his)){
      stop("The current accuracy is specified as relative error limit,
           please give the historical mean!")
    }
    bound <- (bound*mean_his/quan)^2
  }
  if(method == "CV"){
    if(is.null(mean_his)){
      stop("The current accuracy is specified as upper limit of variation coefficient,
           please give the historical mean!")
    }
    bound <- (bound*mean_his)^2
  }
  n0 <- var_his/bound
  size <- ifelse(is.infinite(N),n0,n0/(1+n0/N))
  return(list(method = method,n0 = ceiling(n0),size = ceiling(size)))
}


#' Sample Size for Proportion
#'
#' This function estimates the sample size required for a given proportion estimator precision.
#'
#' @param N a number, the size of population, If it is unknown or very large, it can be omitted.
#' By default, infinity is used for calculation.
#' @param p_his a number, historical estimates of population proportion.
#' @param method string, the form of precision requirement.
#' For upper limit of variance, use "V", for absolute error limit, use "d",
#' for relative error limit, use "r", for upper limit of variation coefficient,
#' use "CV".
#' @param bound a number, the upper limit value of the precision requirement.
#' @param alpha a number, the confidence level.
#' @return a list contain the form of precision requirement, uncorrected sample
#' size and corrected sample size.
#' @note If the number estimator precision is given, convert it to the proportion estimator
#' precision then use this function. This function will convert all forms of upper precision to
#' upper limit of variance in calculation.
srs_size_for_prop <- function(p_his,N = Inf,method,bound,alpha = 0.05){

  quan <- qnorm(1-alpha/2)

  if(method == "d"){
    bound <- (bound/quan)^2
  }
  if(method == "r"){
    bound <- (bound*p_his/quan)^2
  }
  if(method == "CV"){
    bound <- (bound*p_his)^2
  }
  n0 <- p_his*(1-p_his)/bound
  size <- ifelse(is.infinite(N),n0,n0/(1+(n0-1)/N))
  return(list(method = method,n0 = ceiling(n0),size = ceiling(size)))
}

