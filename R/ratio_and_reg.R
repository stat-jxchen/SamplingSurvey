#' Estimate Ratio
#'
#' This function estimates the ratio between the main variable y and the auxiliary variable x.
#'
#' @param ybar a number, the sample mean of the main variable y.
#' @param xbar a number, the sample mean of the auxiliary variable x.
#' @param sy2 a number, the sample variance of the main variable y.
#' @param sx2 a number, the sample variance of the auxiliary variable x.
#' @param syx a number, the sample covariance of y and x.
#' @param n a number, the sample size.
#' @param N a number, the population size.
#' @param var_method string, the method to obtain a estimation of the variance
#' of the estimator. If choose "xbar",then will use \code{xbar} directly, else if
#' choose "Xbar", then \code{Xbar} needs be given.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @export
ratio <- function(ybar, xbar, sy2, sx2, syx, n, N = Inf,
                       var_method = "xbar", ci_method = "ordinary",
                       Xbar = NULL, alpha = 0.05){
  ratio_est <- ybar/xbar
  f <- n/N
  if(var_method == "xbar"){
    ratio_var_est <- (1-f)/(n*xbar^2)*(sy2+ratio_est^2*sx2-2*ratio_est*syx)
  }
  if(var_method == "Xbar"){
    ratio_var_est <- (1-f)/(n*Xbar^2)*(sy2+ratio_est^2*sx2-2*ratio_est*syx)
  }
  ratio_se_est <- sqrt(ratio_var_est)
  if(ci_method == "ordinary"){
    ratio_ci_est <- conf_interval(ratio_est,ratio_se_est,alpha = alpha)
    ratio_ci_left <- ratio_ci_est$ci_left_bound
    ratio_ci_right <- ratio_ci_est$ci_right_bound
  }
  if(ci_method == "equation"){
    quan = qnorm(1-alpha/2)
    nf <- (1-f)/n
    c2_ybar <- nf*sy2/ybar^2
    c2_xbar <- nf*sx2/xbar^2
    c_yxbar <- nf*syx/(xbar*ybar)
    mediate1 <- 1-quan^2*c_yxbar
    mediate2 <- quan*sqrt((c2_ybar+c2_xbar-2*c_yxbar)-
                            quan^2*(c2_ybar*c2_xbar-c_yxbar^2))
    mediate3 <- 1-quan^2*c2_xbar
    ratio_ci_left <- ratio_est*(mediate1-mediate2)/mediate3
    ratio_ci_right <- ratio_est*(mediate1+mediate2)/mediate3
  }
  result <- list(ratio_est = ratio_est,
                 ratio_var_est = ratio_var_est,
                 ratio_se_est = ratio_se_est,
                 ratio_ci_left = ratio_ci_left,
                 ratio_ci_right = ratio_ci_right)
  return(result)

}


#' Estimate Population Mean with Ratio Estimation
#'
#' This function uses ratio estimation to estimate the population mean.
#'
#' @param ybar a number, the sample mean of the main variable y.
#' @param xbar a number, the sample mean of the auxiliary variable x.
#' @param sy2 a number, the sample variance of the main variable y.
#' @param sx2 a number, the sample variance of the auxiliary variable x.
#' @param syx a number, the sample covariance of y and x.
#' @param n a number, the sample size.
#' @param N a number, the population size.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @export
ratio_mean <- function(ybar, xbar, sy2, sx2, syx, n, Xbar, N = Inf,
                            ci_method = "ordinary", alpha = 0.05){
  ratio_est <- ybar/xbar
  ybar_ratio <- ratio_est*Xbar
  f <- n/N
  ybar_var_est <- (1-f)/n*(sy2+ratio_est^2*sx2-2*ratio_est*syx)
  ybar_se_est <- sqrt(ybar_var_est)
  if(ci_method == "ordinary"){
    ybar_ci_est <- conf_interval(ybar_ratio,ybar_se_est,alpha = alpha)
    ybar_ci_left <- ybar_ci_est$ci_left_bound
    ybar_ci_right <- ybar_ci_est$ci_right_bound
  }
  if(ci_method == "equation"){
    quan = qnorm(1-alpha/2)
    nf <- (1-f)/n
    c2_ybar <- nf*sy2/ybar^2
    c2_xbar <- nf*sx2/xbar^2
    c_yxbar <- nf*syx/(xbar*ybar)
    mediate1 <- 1-quan^2*c_yxbar
    mediate2 <- quan*sqrt((c2_ybar+c2_xbar-2*c_yxbar)-
                            quan^2*(c2_ybar*c2_xbar-c_yxbar^2))
    mediate3 <- 1-quan^2*c2_xbar
    ybar_ci_left <- ybar_ratio*(mediate1-mediate2)/mediate3
    ybar_ci_right <- ybar_ratio*(mediate1+mediate2)/mediate3
  }
  result <- list(ybar_ratio = ybar_ratio,
                 ybar_var_est = ybar_var_est,
                 ybar_se_est = ybar_se_est,
                 ybar_ci_left = ybar_ci_left,
                 ybar_ci_right = ybar_ci_right)
  return(result)

}


#' Estimate Population Total with Ratio Estimation
#'
#' This function uses ratio estimation to estimate the population total.
#'
#' @param ybar a number, the sample mean of the main variable y.
#' @param xbar a number, the sample mean of the auxiliary variable x.
#' @param sy2 a number, the sample variance of the main variable y.
#' @param sx2 a number, the sample variance of the auxiliary variable x.
#' @param syx a number, the sample covariance of y and x.
#' @param n a number, the sample size.
#' @param N a number, the population size.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @export
ratio_total <- function(ybar, xbar, sy2, sx2, syx, n, Xbar, N,
                            ci_method = "ordinary", alpha = 0.05){

  mean_result <- ratio_mean(ybar = ybar, xbar = xbar, sy2 = sy2, sx2 = sx2,
                            syx = syx, n = n, Xbar = Xbar, N = N,
                            ci_method = ci_method, alpha = alpha)
  y_ratio <- mean_result$ybar_ratio*N
  y_var_est <- mean_result$ybar_var_est*N^2
  y_se_est <- mean_result$ybar_se_est*N
  y_ci_left <- mean_result$ybar_ci_left*N
  y_ci_right <- mean_result$ybar_ci_right*N
  result <- list(y_ratio = y_ratio,
                 y_var_est = y_var_est,
                 y_se_est = y_se_est,
                 y_ci_left = y_ci_left,
                 y_ci_right = y_ci_right)
  return(result)

}


#' Sample Information Extraction
#'
#' This function is an auxiliary function for ratio estimation and regression estimation
#' and can extract the sample mean, sample variance and sample covariance of the
#' main variable y and the auxiliary variable x.
#'
#' @param y_sample a vector, the sample of the main variable y.
#' @param x_sample a vector, the sample of the auxiliary variable x.
#' @export
ratio_sample_extract <- function(y_sample,x_sample){
  xbar <- mean(x)
  ybar <- mean(y)
  n <- length(y)
  sx2 <- var(x)
  sy2 <- var(y)
  syx <- cov(x,y)
  result <- list(ybar = ybar,
                 xbar = xbar,
                 sy2 = sy2,
                 sx2 = sx2,
                 syx = syx,
                 n = n)
  return(result)
}


#' Estimate Population Mean with Regression Estimation
#'
#' This function uses regression estimation to estimate the population mean.
#'
#' @param ybar a number, the sample mean of the main variable y.
#' @param xbar a number, the sample mean of the auxiliary variable x.
#' @param sy2 a number, the sample variance of the main variable y.
#' @param sx2 a number, the sample variance of the auxiliary variable x.
#' @param syx a number, the sample covariance of y and x.
#' @param n a number, the sample size.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param N a number, the population size.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param beta0 a number.
#' @param alpha a number, the confidence level.
#' @export
reg_mean <- function(ybar, xbar, sy2, sx2, syx, n, Xbar, N = Inf,
                          beta_method = "reg_coef", beta0 = NULL,
                          alpha = 0.05){
  f <- n/N
  if(beta_method == "reg_coef"){
    beta <- syx/sx2
    ybar_var_est <- (1-f)/n*(n-1)/(n-2)*(sy2-syx^2/sx2)
  }

  if(beta_method == "const"){
    beta = beta0
    ybar_var_est <- (1-f)/n*(sy2+beta^2*sx2-2*beta*syx)
  }
  ybar_reg <- ybar+beta*(Xbar-xbar)
  ybar_se_est <- sqrt(ybar_var_est)
  ybar_ci_est <- conf_interval(ybar_reg,ybar_se_est,alpha = alpha)
  ybar_ci_left <- ybar_ci_est$ci_left_bound
  ybar_ci_right <- ybar_ci_est$ci_right_bound
  result <- list(ybar_reg = ybar_reg,
                 ybar_var_est = ybar_var_est,
                 ybar_se_est = ybar_se_est,
                 ybar_ci_left = ybar_ci_left,
                 ybar_ci_right = ybar_ci_right)
  return(result)
}


#' Estimate Population Total with Regression Estimation
#'
#' This function uses regression estimation to estimate the population total.
#'
#' @param ybar a number, the sample mean of the main variable y.
#' @param xbar a number, the sample mean of the auxiliary variable x.
#' @param sy2 a number, the sample variance of the main variable y.
#' @param sx2 a number, the sample variance of the auxiliary variable x.
#' @param syx a number, the sample covariance of y and x.
#' @param n a number, the sample size.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param N a number, the population size.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param beta0 a number.
#' @param alpha a number, the confidence level.
#' @export
reg_total <- function(ybar, xbar, sy2, sx2, syx, n, Xbar, N,
                          beta_method = "reg_coef", beta0 = NULL,
                          alpha = 0.05){
  mean_result <- reg_mean(ybar = ybar, xbar = xbar, sy2 = sy2, sx2 = sx2,
                               syx = syx, n = n, Xbar = Xbar, N = N,
                               beta_method = beta_method, beta0 = beta0,
                               alpha = alpha)
  y_reg <- mean_result$ybar_reg*N
  y_var_est = mean_result$ybar_var_est*N^2
  y_se_est = mean_result$ybar_se_est*N
  y_ci_left = mean_result$ybar_ci_left*N
  y_ci_right = mean_result$ybar_ci_right*N
  result <- list(y_reg = y_reg,
                 y_var_est = y_var_est,
                 y_se_est = y_se_est,
                 y_ci_left = y_ci_left,
                 y_ci_right = y_ci_right)
  return(result)
}


#' Separate Ratio Estimator for Population Mean
#'
#' This function estimates the population mean using the separate ratio estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbarh a vector, the historical strata mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
#' @export
sep_ratio_mean <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                           sy2h, sx2h, syxh, Xbarh, alpha = 0.05,
                           stra_est = FALSE){
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
  Rh_est <- ybarh/xbarh
  if(stra_est){
    ybarh_ratio <- Rh_est*Xbarh
    ybarh_ratio_var_est <- (1-fh)/nh*(sy2h+Rh_est^2*sx2h-2*Rh_est*syxh)
    ybarh_ratio_se_est <- sqrt(ybarh_ratio_var_est)
    ybarh_ratio_ci_est <- conf_interval(ybarh_ratio,ybarh_ratio_se_est,
                                        alpha = alpha)
    ybarh_ratio_ci_left <- ybarh_ratio_ci_est$ci_left_bound
    ybarh_ratio_ci_right <- ybarh_ratio_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,ybarh_ratio,
                                       ybarh_ratio_var_est,
                                       ybarh_ratio_se_est,
                                       ybarh_ratio_ci_left,
                                       ybarh_ratio_ci_right))
  }

  # the overall estimated results
  ybar_RS <- sum(Wh*Rh_est*Xbarh)
  ybar_RS_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h+Rh_est^2*sx2h-2*Rh_est*syxh))
  ybar_RS_se_est <- sqrt(ybar_RS_var_est)
  ybar_RS_ci_est <- conf_interval(ybar_RS,ybar_RS_se_est,alpha = alpha)
  ybar_RS_ci_left <- ybar_RS_ci_est$ci_left_bound
  ybar_RS_ci_right <- ybar_RS_ci_est$ci_right_bound
  mean_result <- matrix(c(ybar_RS,ybar_RS_var_est,ybar_RS_se_est,
                          ybar_RS_ci_left,ybar_RS_ci_right),nrow = 1)
  colnames(mean_result) <- c("ybar","ybar_var_est","ybar_se_est",
                             "ybar_ci_left","ybar_ci_right")
  rownames(mean_result) <- "Sep_Ratio"

  if(stra_est){
    return(list(stra_result = stra_result,
                mean_result = as.data.frame(mean_result)))
  }else{
    return(list(mean_result = as.data.frame(mean_result)))
  }
}

#' Separate Ratio Estimator for Population Total
#'
#' This function estimates the population total using the separate ratio estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbarh a vector, the historical strata mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
#' @export
sep_ratio_total <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                            sy2h, sx2h, syxh, Xbarh, alpha = 0.05,
                            stra_est = FALSE){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }

  # the estimated results of each strata
  fh <- nh/Nh
  Rh_est <- ybarh/xbarh
  if(stra_est){
    yh_ratio <- Rh_est*Xbarh*Nh
    yh_ratio_var_est <- (1-fh)/nh*(sy2h+Rh_est^2*sx2h-2*Rh_est*syxh)*Nh^2
    yh_ratio_se_est <- sqrt(yh_ratio_var_est)
    yh_ratio_ci_est <- conf_interval(yh_ratio,yh_ratio_se_est,
                                     alpha = alpha)
    yh_ratio_ci_left <- yh_ratio_ci_est$ci_left_bound
    yh_ratio_ci_right <- yh_ratio_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,yh_ratio,
                                       yh_ratio_var_est,
                                       yh_ratio_se_est,
                                       yh_ratio_ci_left,
                                       yh_ratio_ci_right))
  }

  # the overall estimated results
  y_RS <- sum(Wh*Rh_est*Xbarh)*N
  y_RS_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h+Rh_est^2*sx2h-2*Rh_est*syxh))*N^2
  y_RS_se_est <- sqrt(y_RS_var_est)
  y_RS_ci_est <- conf_interval(y_RS,y_RS_se_est,alpha = alpha)
  y_RS_ci_left <- y_RS_ci_est$ci_left_bound
  y_RS_ci_right <- y_RS_ci_est$ci_right_bound
  total_result <- matrix(c(y_RS,y_RS_var_est,y_RS_se_est,
                           y_RS_ci_left,y_RS_ci_right),nrow = 1)
  colnames(total_result) <- c("y","y_var_est","y_se_est",
                              "y_ci_left","y_ci_right")
  rownames(total_result) <- "Sep_Ratio"

  if(stra_est){
    return(list(stra_result = stra_result,
                total_result = as.data.frame(total_result)))
  }else{
    return(list(total_result = as.data.frame(total_result)))
  }
}



#' Combined Ratio Estimator for Population Mean
#'
#' This function estimates the population mean using the combined ratio estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @export
com_ratio_mean <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                           sy2h, sx2h, syxh, Xbar, alpha = 0.05){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }
  fh <- nh/Nh
  ybar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = ybarh,
                            s2h = sy2h, alpha = alpha)$mean_result$ybar
  xbar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = xbarh,
                            s2h = sx2h, alpha = alpha)$mean_result$ybar
  R_est <- ybar_est/xbar_est
  ybar_RC <- R_est*Xbar
  ybar_RC_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h+R_est^2*sx2h-2*R_est*syxh))
  ybar_RC_se_est <- sqrt(ybar_RC_var_est)
  ybar_RC_ci_est <- conf_interval(ybar_RC,ybar_RC_se_est,alpha = alpha)
  ybar_RC_ci_left <- ybar_RC_ci_est$ci_left_bound
  ybar_RC_ci_right <- ybar_RC_ci_est$ci_right_bound
  result <- list(R_est = R_est,
                 ybar_RC = ybar_RC,
                 ybar_RC_var_est = ybar_RC_var_est,
                 ybar_RC_se_est = ybar_RC_se_est,
                 ybar_RC_ci_left = ybar_RC_ci_left,
                 ybar_RC_ci_right = ybar_RC_ci_right)
  return(result)
}

#' Combined Ratio Estimator for Population Total
#'
#' This function estimates the population total using the combined ratio estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @export
com_ratio_total <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                            sy2h, sx2h, syxh, Xbar, alpha = 0.05){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }
  fh <- nh/Nh
  ybar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = ybarh,
                            s2h = sy2h, alpha = alpha)$mean_result$ybar
  xbar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = xbarh,
                            s2h = sx2h, alpha = alpha)$mean_result$ybar
  R_est <- ybar_est/xbar_est
  y_RC <- R_est*Xbar*N
  y_RC_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h+R_est^2*sx2h-2*R_est*syxh))*N^2
  y_RC_se_est <- sqrt(y_RC_var_est)
  y_RC_ci_est <- conf_interval(y_RC,y_RC_se_est,alpha = alpha)
  y_RC_ci_left <- y_RC_ci_est$ci_left_bound
  y_RC_ci_right <- y_RC_ci_est$ci_right_bound
  result <- list(R_est = R_est,
                 y_RC = y_RC,
                 y_RC_var_est = y_RC_var_est,
                 y_RC_se_est = y_RC_se_est,
                 y_RC_ci_left = y_RC_ci_left,
                 y_RC_ci_right = y_RC_ci_right)
  return(result)
}



#' Separate Regression Estimator for Population Mean
#'
#' This function estimates the population mean using the separate regression estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the strata sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param Xbarh a vector, the historical strata mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @param betah0 a number or a vector. If the coefficient for each stata is the
#' same, then can just input a number.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
#' @export
sep_reg_mean <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                         sy2h, sx2h, syxh, Xbarh, beta_method = "reg_coef",
                         alpha = 0.05, betah0 = NULL, stra_est = FALSE){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }

  # the estimated results of each strata
  fh <- nh/Nh
  if(beta_method == "reg_coef"){
    betah <- syxh/sx2h
    ybar_regs_var_est <- sum(Wh^2*(1-fh)/nh*(nh-1)/(nh-2)*(sy2h-syxh^2/sx2h))
  }
  if(beta_method == "const"){
    betah <- betah0
    ybar_regs_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h-2*betah*syxh+betah^2*sx2h))
  }
  if(stra_est){
    ybarh_regs <- ybarh+betah*(Xbarh-xbarh)
    ybarh_regs_var_est <- sum((1-fh)/nh*(nh-1)/(nh-2)*(sy2h-syxh^2/sx2h))
    ybarh_regs_se_est <- sqrt(ybarh_regs_var_est)
    ybarh_regs_ci_est <- conf_interval(ybarh_regs,ybarh_regs_se_est,
                                        alpha = alpha)
    ybarh_regs_ci_left <- ybarh_regs_ci_est$ci_left_bound
    ybarh_regs_ci_right <- ybarh_regs_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,betah,
                                       ybarh_regs,
                                       ybarh_regs_var_est,
                                       ybarh_regs_se_est,
                                       ybarh_regs_ci_left,
                                       ybarh_regs_ci_right))
  }
  # the overall estimated results
  ybar_regs <- sum(Wh*(ybarh+betah*(Xbarh-xbarh)))
  ybar_regs_se_est <- sqrt(ybar_regs_var_est)
  ybar_regs_ci_est <- conf_interval(ybar_regs,ybar_regs_se_est,alpha = alpha)
  ybar_regs_ci_left <- ybar_regs_ci_est$ci_left_bound
  ybar_regs_ci_right <- ybar_regs_ci_est$ci_right_bound
  mean_result <- matrix(c(ybar_regs,ybar_regs_var_est,ybar_regs_se_est,
                          ybar_regs_ci_left,ybar_regs_ci_right),nrow = 1)
  colnames(mean_result) <- c("ybar","ybar_var_est","ybar_se_est",
                             "ybar_ci_left","ybar_ci_right")
  rownames(mean_result) <- "Sep_Reg"

  if(stra_est){
    return(list(stra_result = stra_result,
                mean_result = as.data.frame(mean_result)))
  }else{
    return(list(mean_result = as.data.frame(mean_result)))
  }
}

#' Separate Regression Estimator for Population Total
#'
#' This function estimates the population total using the separate regression estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the strata sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param Xbarh a vector, the historical strata mean of the auxiliary variable x.
#' @param alpha a number, the confidence level.
#' @param betah0 a number or a vector. If the coefficient for each stata is the
#' same, then can just input a number.
#' @param stra_est logical, whether to calculate and output the estimated results
#' of each strata.
#' @export
sep_reg_total <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                          sy2h, sx2h, syxh, Xbarh, beta_method = "reg_coef",
                          alpha = 0.05, betah0 = NULL, stra_est = FALSE){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }

  # the estimated results of each strata
  fh <- nh/Nh
  if(beta_method == "reg_coef"){
    betah <- syxh/sx2h
    y_regs_var_est <- sum(Wh^2*(1-fh)/nh*(nh-1)/(nh-2)*(sy2h-syxh^2/sx2h))*N^2
  }
  if(beta_method == "const"){
    betah <- betah0
    y_regs_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h-2*betah*syxh+betah^2*sx2h))*N^2
  }
  if(stra_est){
    yh_regs <- (ybarh+betah*(Xbarh-xbarh))*Nh
    yh_regs_var_est <- sum((1-fh)/nh*(nh-1)/(nh-2)*(sy2h-syxh^2/sx2h))*Nh^2
    yh_regs_se_est <- sqrt(yh_regs_var_est)
    yh_regs_ci_est <- conf_interval(yh_regs,yh_regs_se_est,
                                    alpha = alpha)
    yh_regs_ci_left <- yh_regs_ci_est$ci_left_bound
    yh_regs_ci_right <- yh_regs_ci_est$ci_right_bound
    stra_result <- as.data.frame(cbind(Nh,nh,Wh,betah,
                                       yh_regs,
                                       yh_regs_var_est,
                                       yh_regs_se_est,
                                       yh_regs_ci_left,
                                       yh_regs_ci_right))
  }
  # the overall estimated results
  y_regs <- sum(Wh*(ybarh+betah*(Xbarh-xbarh)))*N
  y_regs_se_est <- sqrt(y_regs_var_est)
  y_regs_ci_est <- conf_interval(y_regs,y_regs_se_est,alpha = alpha)
  y_regs_ci_left <- y_regs_ci_est$ci_left_bound
  y_regs_ci_right <- y_regs_ci_est$ci_right_bound
  total_result <- matrix(c(y_regs,y_regs_var_est,y_regs_se_est,
                           y_regs_ci_left,y_regs_ci_right),nrow = 1)
  colnames(total_result) <- c("y","y_var_est","y_se_est",
                              "y_ci_left","y_ci_right")
  rownames(total_result) <- "Sep_Reg"

  if(stra_est){
    return(list(stra_result = stra_result,
                total_result = as.data.frame(total_result)))
  }else{
    return(list(total_result = as.data.frame(total_result)))
  }
}



#' Combined Regression Estimator for Population Mean
#'
#' This function estimates the population mean using the combined regression estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the strata sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param alpha a number, the confidence level.
#' @param beta0 a number.
#' @export
com_reg_mean <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                         sy2h, sx2h, syxh, Xbar, beta_method = "reg_coef",
                         alpha = 0.05, beta0 = NULL){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }
  fh <- nh/Nh
  ybar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = ybarh,
                            s2h = sy2h, alpha = alpha)$mean_result$ybar
  xbar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = xbarh,
                            s2h = sx2h, alpha = alpha)$mean_result$ybar
  if(beta_method == "reg_coef"){
    beta <- sum(Wh^2*(1-fh)*syxh/nh)/sum(Wh^2*(1-fh)*sx2h/nh)
  }
  if(beta_method == "const"){
    beta = beta0
  }
  ybar_regc <- ybar_est+beta*(Xbar-xbar_est)
  ybar_regc_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h-2*beta*syxh+beta^2*sx2h))
  ybar_regc_se_est <- sqrt(ybar_regc_var_est)
  ybar_regc_ci_est <- conf_interval(ybar_regc,ybar_regc_se_est,alpha = alpha)
  ybar_regc_ci_left <- ybar_regc_ci_est$ci_left_bound
  ybar_regc_ci_right <- ybar_regc_ci_est$ci_right_bound
  result <- list(beta = beta,
                 ybar_regc = ybar_regc,
                 ybar_regc_var_est = ybar_regc_var_est,
                 ybar_regc_se_est = ybar_regc_se_est,
                 ybar_regc_ci_left = ybar_regc_ci_left,
                 ybar_regc_ci_right = ybar_regc_ci_right)
  return(result)
}

#' Combined Regression Estimator for Population Total
#'
#' This function estimates the population total using the combined regression estimator.
#'
#' @param Nh a vector, the population size of each strata.
#' If it is be provided, then \code{Wh} and \code{N} can be ignored.
#' @param Wh a vector, the weight of each strata.
#' @param N a number, the population size.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the strata sample mean of the main variable y.
#' @param xbarh a vector, the strata sample mean of the auxiliary variable x.
#' @param sy2h a vector, the strata sample variance of the main variable y.
#' @param sx2h a vector, the strata sample variance of the auxiliary variable x.
#' @param syxh a vector, the strata sample covariance of y and x.
#' @param Xbar a number, the historical mean of the auxiliary variable x.
#' @param beta_method string, coefficient of auxiliary variables x in regression estimation.
#' If choose "reg_coef", then the strata sample regression coefficient will be used.
#' If choose "const", then \code{beta0} need be given and be used.
#' @param alpha a number, the confidence level.
#' @param beta0 a number.
#' @export
com_reg_total <- function(Nh = NULL, Wh = NULL, N = NULL, nh, ybarh, xbarh,
                          sy2h, sx2h, syxh, Xbar, beta_method = "reg_coef",
                          alpha = 0.05, beta0 = NULL){
  condition <- c(is.null(Nh),is.null(Wh),is.null(N))
  if(condition[1]){
    if(sum(condition)<2){
      stop("Nh is not given, then Wh and N must be given together!")
    }else{
      Nh <- Wh*N
    }
  }else{
    if(condition[2]){
      N <- sum(Nh)
      Wh <- Nh/N
    }
  }
  fh <- nh/Nh
  ybar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = ybarh,
                            s2h = sy2h, alpha = alpha)$mean_result$ybar
  xbar_est <- stra_srs_mean(Nh = Nh,Wh = Wh,N = N,nh = nh, ybarh = xbarh,
                            s2h = sx2h, alpha = alpha)$mean_result$ybar
  if(beta_method == "reg_coef"){
    beta <- sum(Wh^2*(1-fh)*syxh/nh)/sum(Wh^2*(1-fh)*sx2h/nh)
  }
  if(beta_method == "const"){
    beta = beta0
  }
  y_regc <- (ybar_est+beta*(Xbar-xbar_est))*N
  y_regc_var_est <- sum(Wh^2*(1-fh)/nh*(sy2h-2*beta*syxh+beta^2*sx2h))*N^2
  y_regc_se_est <- sqrt(y_regc_var_est)
  y_regc_ci_est <- conf_interval(y_regc,y_regc_se_est,alpha = alpha)
  y_regc_ci_left <- y_regc_ci_est$ci_left_bound
  y_regc_ci_right <- y_regc_ci_est$ci_right_bound
  result <- list(beta = beta,
                 y_regc = y_regc,
                 y_regc_var_est = y_regc_var_est,
                 y_regc_se_est = y_regc_se_est,
                 y_regc_ci_left = y_regc_ci_left,
                 y_regc_ci_right = y_regc_ci_right)
  return(result)
}

