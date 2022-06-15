#' Estimate Population Mean Under Double Sampling with Stratification
#'
#' This function estimates the population mean when double sampling that the first
#' sampling method is simple random sample sample and the second sampling method
#' is stratified sample is used.
#'
#' @param nh1 a vector, the strata sample size of the first phase sample.
#' @param nh2 a vector, the strata sample size of the second phase sample.
#' @param ybarh a vector, the strata sample mean of the second phase sample.
#' @param s2h a vector, the strata sample variance of the second phase sample.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
twophase_stra_mean <- function(nh1, nh2, ybarh, s2h, N = Inf, alpha = 0.05){
  n1 <- sum(nh1)
  wh1 <- nh1/n1
  ybar_stD <- sum(wh1*ybarh)
  var_est_depos1 <- sum((1/nh2-1/nh1)*wh1^2*s2h)
  var_est_depos2 <- (1/n1-1/N)*sum(wh1*(ybarh-ybar_stD)^2)
  ybar_stD_var_est <- var_est_depos1+var_est_depos2
  ybar_stD_se_est <- sqrt(ybar_stD_var_est)
  ybar_stD_ci_est <- conf_interval(ybar_stD,ybar_stD_se_est,alpha = alpha)
  ybar_stD_ci_left <- ybar_stD_ci_est$ci_left_bound
  ybar_stD_ci_right <- ybar_stD_ci_est$ci_right_bound
  result <- list(ybar_stD = ybar_stD,
                 ybar_stD_var_est = ybar_stD_var_est,
                 ybar_stD_se_est = ybar_stD_se_est,
                 ybar_stD_ci_left = ybar_stD_ci_left,
                 ybar_stD_ci_right = ybar_stD_ci_right)
  return(result)
}

#' Estimate Population Total Under Double Sampling with Stratification
#'
#' This function estimates the population total when double sampling that the first
#' sampling method is simple random sample sample and the second sampling method
#' is stratified sample is used.
#'
#' @param nh1 a vector, the strata sample size of the first phase sample.
#' @param nh2 a vector, the strata sample size of the second phase sample.
#' @param ybarh a vector, the strata sample mean of the second phase sample.
#' @param s2h a vector, the strata sample variance of the second phase sample.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
twophase_stra_total <- function(nh1, nh2, ybarh, s2h, N, alpha = 0.05){
  n1 <- sum(nh1)
  wh1 <- nh1/n1
  ybar_stD <- sum(wh1*ybarh)
  y_stD <- ybar_stD*N
  var_est_depos1 <- sum((1/nh2-1/nh1)*wh1^2*s2h)*N^2
  var_est_depos2 <- (1/n1-1/N)*sum(wh1*(ybarh-ybar_stD)^2)*N^2
  y_stD_var_est <- var_est_depos1+var_est_depos2
  y_stD_se_est <- sqrt(y_stD_var_est)
  y_stD_ci_est <- conf_interval(y_stD,y_stD_se_est,alpha = alpha)
  y_stD_ci_left <- y_stD_ci_est$ci_left_bound
  y_stD_ci_right <- y_stD_ci_est$ci_right_bound
  result <- list(y_stD = y_stD,
                 y_stD_var_est = y_stD_var_est,
                 y_stD_se_est = y_stD_se_est,
                 y_stD_ci_left = y_stD_ci_left,
                 y_stD_ci_right = y_stD_ci_right)
  return(result)
}

#' Estimate Population Mean Under Double Sampling with Ratio Estimation
#'
#' This function estimates the population mean when double sampling that the first
#' sample is only used for the auxiliary variable x and the second sample is used
#' for the main variable y.
#'
#' @param n1 a vector, the sample size of the first phase sample.
#' @param n2 a vector, the sample size of the second phase sample.
#' @param xbar1 a number, the sample mean of the auxiliary variable x in the
#' first phase sample.
#' @param ybar a number, the sample mean of the main variable y in the second
#' phase sample.
#' @param xbar a number, the sample mean of the auxiliary variable x in the
#' second phase sample
#' @param sx2 a number, the sample variance of the main variable y in the second
#' phase sample.
#' @param sy2 a number, the sample variance of the auxiliary variable x in the
#' second phase sample.
#' @param syx a number, the sample covariance of y and x in the second phase sample.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
twophase_ratio_mean <- function(n1, n2, xbar1, ybar, xbar, sx2, sy2, syx,
                                N = Inf, alpha = 0.05){
  ratio_est <- ybar/xbar
  ybar_RD <- ratio_est*xbar1
  var_est_depos1 <- (1/n1-1/N)*sy2
  var_est_depos2 <- (1/n2-1/n1)*(sy2+ratio_est^2*sx2-2*ratio_est*syx)
  ybar_RD_var_est <- var_est_depos1+var_est_depos2
  ybar_RD_se_est <- sqrt(ybar_RD_var_est)
  ybar_RD_ci_est <- conf_interval(ybar_RD,ybar_RD_se_est,alpha = alpha)
  ybar_RD_ci_left <- ybar_RD_ci_est$ci_left_bound
  ybar_RD_ci_right <- ybar_RD_ci_est$ci_right_bound
  result <- list(ratio_est = ratio_est,
                 ybar_RD = ybar_RD,
                 ybar_RD_var_est = ybar_RD_var_est,
                 ybar_RD_se_est = ybar_RD_se_est,
                 ybar_RD_ci_left = ybar_RD_ci_left,
                 ybar_RD_ci_right = ybar_RD_ci_right)
  return(result)
}

#' Estimate Population Total Under Double Sampling with Ratio Estimation
#'
#' This function estimates the population total when double sampling that the first
#' sample is only used for the auxiliary variable x and the second sample is used
#' for the main variable y.
#'
#' @param n1 a vector, the sample size of the first phase sample.
#' @param n2 a vector, the sample size of the second phase sample.
#' @param xbar1 a number, the sample mean of the auxiliary variable x in the
#' first phase sample.
#' @param ybar a number, the sample mean of the main variable y in the second
#' phase sample.
#' @param xbar a number, the sample mean of the auxiliary variable x in the
#' second phase sample
#' @param sx2 a number, the sample variance of the main variable y in the second
#' phase sample.
#' @param sy2 a number, the sample variance of the auxiliary variable x in the
#' second phase sample.
#' @param syx a number, the sample covariance of y and x in the second phase sample.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
twophase_ratio_total <- function(n1, n2, xbar1, ybar, xbar, sx2, sy2, syx,
                                N, alpha = 0.05){
  ratio_est <- ybar/xbar
  y_RD <- ratio_est*xbar1*N
  var_est_depos1 <- (1/n1-1/N)*sy2*N^2
  var_est_depos2 <- (1/n2-1/n1)*(sy2+ratio_est^2*sx2-2*ratio_est*syx)*N^2
  y_RD_var_est <- var_est_depos1+var_est_depos2
  y_RD_se_est <- sqrt(y_RD_var_est)
  y_RD_ci_est <- conf_interval(y_RD,y_RD_se_est,alpha = alpha)
  y_RD_ci_left <- y_RD_ci_est$ci_left_bound
  y_RD_ci_right <- y_RD_ci_est$ci_right_bound
  result <- list(ratio_est = ratio_est,
                 y_RD = y_RD,
                 y_RD_var_est = y_RD_var_est,
                 y_RD_se_est = y_RD_se_est,
                 y_RD_ci_left = y_RD_ci_left,
                 y_RD_ci_right = y_RD_ci_right)
  return(result)
}

#' Sample Information Extraction
#'
#' This function is an auxiliary function which can extract the sample mean, sample variance and sample covariance of the
#' main variable y and the auxiliary variable x.
#'
#' @param yxsamplea matrix or dataframe, which have two columns. One for y,
#' another for x.
#' @param yloc a number, the location of y column in the dataframe.
extract_srs_yxsample <- function(yxsample,yloc){
  ysample <- yxsample[[yloc]]
  xsample <- yxsample[,-yloc][[1]]
  ybar <- mean(ysample)
  xbar <- mean(xsample)
  sy2 <- var(ysample)
  sx2 <- var(xsample)
  syx <- cov(ysample,xsample)
  return(list(ybar = ybar,
              xbar = xbar,
              sy2 = sy2,
              sx2 = sx2,
              syx = syx))
}

