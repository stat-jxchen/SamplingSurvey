#' Estimate Population Mean Under Cluster Sampling
#'
#' This function estimates the population mean when cluster sampling is used and
#' the cluster sizes are equal.
#'
#' @param M a number, the size of each cluster.
#' @param s2b a number, the between-sample variance.
#' @param N a number, the number of cluster in population.
#' @param ybari a vector, the sample mean of each cluster in sample.
#' @param ybar a number, the mean of the sample sum of each cluster.
#' @param n a number, the number of cluster in sample.
#' @param alpha a number, the confidence level.
cluster_eqM_mean <- function(M, s2b, N = Inf, ybari = NULL, ybar = NULL,
                            n = NULL, alpha = 0.05){
  if(!is.null(ybari)){
    n <- length(ybari)
    f <- n/N
    ybar_cluster_est <- mean(ybari)
  }else if(!is.null(ybar)){
    if(is.null(n)){
      stop("n must be given when use ybar!")
    }else{
      f <- n/N
      ybar_cluster_est <- ybar/M
    }
  }else{
    stop("One of ybari and ybar must be given!")
  }
  ybar_cluster_var_est <- (1-f)/(n*M)*s2b
  ybar_cluster_se_est <- sqrt(ybar_cluster_var_est)
  ybar_cluster_ci_est <- conf_interval(ybar_cluster_est, ybar_cluster_se_est,
                                       alpha = alpha)
  ybar_cluster_ci_left <- ybar_cluster_ci_est$ci_left_bound
  ybar_cluster_ci_right <- ybar_cluster_ci_est$ci_right_bound
  result <- list(ybar_cluster_est = ybar_cluster_est,
                 ybar_cluster_var_est = ybar_cluster_var_est,
                 ybar_cluster_se_est = ybar_cluster_se_est,
                 ybar_cluster_ci_left = ybar_cluster_ci_left,
                 ybar_cluster_ci_right = ybar_cluster_ci_right)
  return(result)

}

#' Estimate Population Total Under Cluster Sampling
#'
#' This function estimates the population total when cluster sampling is used and
#' the cluster sizes are equal.
#'
#' @param M a number, the size of each cluster.
#' @param s2b a number, the between-sample variance.
#' @param N a number, the number of cluster in population.
#' @param ybari a vector, the sample mean of each cluster in sample.
#' @param ybar a number, the mean of the sample sum of each cluster.
#' @param n a number, the number of cluster in sample.
#' @param alpha a number, the confidence level.
cluster_eqM_total <- function(M, s2b, N, ybari = NULL, ybar = NULL,
                             n = NULL, alpha = 0.05){
  if(!is.null(ybari)){
    n <- length(ybari)
    f <- n/N
    y_cluster_est <- mean(ybari)*N*M
  }else if(!is.null(ybar)){
    if(is.null(n)){
      stop("n must be given when use ybar!")
    }else{
      f <- n/N
      y_cluster_est <- ybar*N
    }
  }else{
    stop("One of ybari and ybar must be given!")
  }

  y_cluster_var_est <- (1-f)/n*s2b*N^2*M
  y_cluster_se_est <- sqrt(y_cluster_var_est)
  y_cluster_ci_est <- conf_interval(y_cluster_est, y_cluster_se_est,
                                    alpha = alpha)
  y_cluster_ci_left <- y_cluster_ci_est$ci_left_bound
  y_cluster_ci_right <- y_cluster_ci_est$ci_right_bound
  result <- list(y_cluster_est = y_cluster_est,
                 y_cluster_var_est = y_cluster_var_est,
                 y_cluster_se_est = y_cluster_se_est,
                 y_cluster_ci_left = y_cluster_ci_left,
                 y_cluster_ci_right = y_cluster_ci_right)
  return(result)

}

#' Estimate Population Proportion Under Cluster Sampling
#'
#' This function estimates the population proportion when cluster sampling is used.
#'
#' @param m a number or a vector, the size of each cluster, when the size
#' is the same, then \code{m} can be simplified as a number.
#' @param p a vector, the sample proportion of each cluster in sample.
#' @param N a number, the number of cluster in population.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param alpha a number, the confidence level.
cluster_prop <- function(m, p, N = Inf, ci_method = "ordinary", alpha = 0.05){
  n <- length(p)
  f <- n/N
  m_type <- length(unique(m))
  if(m_type == 1){
    p_cluster <- mean(p)
    p_var_est <- (1-f)/n*var(p)
    p_se_est <- sqrt(p_var_est)
    p_ci_est <- conf_interval(p_cluster,p_se_est,alpha = alpha)
    p_ci_left <- p_ci_est$ci_left_bound
    p_ci_right <- p_ci_est$ci_right_bound
  }else{
    a <- m*p
    abar <- mean(a)
    mbar <- mean(m)
    sa2 <- var(a)
    sm2 <- var(m)
    sam <- cov(a,m)
    ratio_re <- ratio(ybar = abar, xbar = mbar, sy2 = sa2, sx2 = sm2,
                      syx = sam, n = n, N = N, var_method = "xbar",
                      ci_method = ci_method, alpha = 0.05)
    p_cluster <- ratio_re$ratio_est
    p_var_est <- ratio_re$ratio_var_est
    p_se_est <- ratio_re$ratio_se_est
    p_ci_left <- ratio_re$ratio_ci_left
    p_ci_right <- ratio_re$ratio_ci_right
  }
  result <- list(p_cluster = p_cluster,
                 p_var_est = p_var_est,
                 p_se_est = p_se_est,
                 p_ci_left = p_ci_left,
                 p_ci_right = p_ci_right)
  return(result)

}

#' Estimate Population Total Under Cluster Sampling
#'
#' This function estimates the population total when cluster sampling is used and
#' the cluster sizes are unequal.
#'
#' @param m a vector, the size of each cluster in sample.
#' @param yi a vector, the sample sum of each cluster.
#' @param N a number, the number of cluster in population.
#' @param M0 a number, the number of secondary units in population.
#' @param method string, the estimate method. "srs" means simple estimation,
#' "ratio" means ratio estimation.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param alpha a number, the confidence level.
cluster_uneqM_total <- function(m, yi, N, M0 = NULL, method = "ratio",
                                ci_method = "ordinary", alpha = 0.05){
  n <- length(yi)
  f <- n/N
  if(method == "srs"){
    y_cluster <- N*mean(yi)
    y_cluster_var_est <- N^2*(1-f)/n*var(yi)
    y_cluster_se_est <- sqrt(y_cluster_var_est)
    y_cluster_ci_est <- conf_interval(y_cluster, y_cluster_se_est,
                                      alpha = alpha)
    y_cluster_ci_left <- y_cluster_ci_est$ci_left_bound
    y_cluster_ci_right <- y_cluster_ci_est$ci_right_bound
  }

  if(method == "ratio"){
    if(is.null(M0)){
      stop("M0 must be given when use ratio estimation!")
    }
    ybar <- mean(yi)
    mbar <- mean(m)
    sy2 <- var(yi)
    sm2 <- var(m)
    sym <- cov(yi,m)
    ratio_re <- ratio_total(ybar = ybar, xbar = mbar, sy2 = sy2, sx2 = sm2,
                            syx = sym, n = n, Xbar = M0/N, N = N,
                            ci_method = ci_method, alpha = alpha)
    y_cluster <- ratio_re$y_ratio
    y_cluster_var_est <- ratio_re$y_var_est
    y_cluster_se_est <- ratio_re$y_se_est
    y_cluster_ci_left <- ratio_re$y_ci_left
    y_cluster_ci_right <- ratio_re$y_ci_right

  }
  result <- list(y_cluster = y_cluster,
                 y_cluster_var_est = y_cluster_var_est,
                 y_cluster_se_est = y_cluster_se_est,
                 y_cluster_ci_left = y_cluster_ci_left,
                 y_cluster_ci_right = y_cluster_ci_right)
  return(result)
}


#' Estimate Population Mean Under Cluster Sampling
#'
#' This function estimates the population mean when cluster sampling is used and
#' the cluster sizes are unequal.
#'
#' @param m a vector, the size of each cluster in sample.
#' @param yi a vector, the sample sum of each cluster.
#' @param N a number, the number of cluster in population.
#' @param M0 a number, the number of primary units in population.
#' @param method string, the estimate method. "srs" means simple estimation,
#' "ratio" means ratio estimation.
#' @param var_method string, the method to obtain a estimation of the variance
#' of the ratio estimator. If choose "mbar",then will use the mean of \code{m}
#' directly, else if choose "Mbar", then will use \code{M0}/\code{N}.
#' @param ci_method string, the method to obtain a estimation of the confidence
#' interval of the estimator, "ordinary" for normal approximation, and "equation"
#' for equation method.
#' @param alpha a number, the confidence level.
cluster_uneqM_mean <- function(m, yi, N, M0 = NULL, method = "ratio",
                               var_method = "mbar", ci_method = "ordinary",
                               alpha = 0.05){
  n <- length(yi)
  f <- n/N
  if(method == "srs"){
    if(is.null(M0)){
      stop("M0 must be given when use simple estimation!")
    }
    ybar_cluster <- N*mean(yi)/M0
    ybar_cluster_var_est <- N^2*(1-f)/n*var(yi)/M0^2
    ybar_cluster_se_est <- sqrt(ybar_cluster_var_est)
    ybar_cluster_ci_est <- conf_interval(ybar_cluster, ybar_cluster_se_est,
                                         alpha = alpha)
    ybar_cluster_ci_left <- ybar_cluster_ci_est$ci_left_bound
    ybar_cluster_ci_right <- ybar_cluster_ci_est$ci_right_bound
  }

  if(method == "ratio"){
    ybar <- mean(yi)
    mbar <- mean(m)
    sy2 <- var(yi)
    sm2 <- var(m)
    sym <- cov(yi,m)
    if(var_method == "mbar"){
      ratio_re <- ratio(ybar = ybar, xbar = mbar, sy2 = sy2, sx2 = sm2,
                        syx = sym, n = n, N = N, ci_method = ci_method,
                        alpha = alpha)
    }
    if(var_method == "Mbar"){
      if(is.null(M0)){
        stop("M0 must be given when use ratio estimation and use Mbar to estimate
             the variance of estimator!")
      }else{
        ratio_re <- ratio(ybar = ybar, xbar = mbar, sy2 = sy2, sx2 = sm2,
                          syx = sym, n = n, N = N, Xbar = M0/N,
                          ci_method = ci_method, var_method = "Xbar",
                          alpha = alpha)
      }
    }
    ybar_cluster <- ratio_re$ratio_est
    ybar_cluster_var_est <- ratio_re$ratio_var_est
    ybar_cluster_se_est <- ratio_re$ratio_se_est
    ybar_cluster_ci_left <- ratio_re$ratio_ci_left
    ybar_cluster_ci_right <- ratio_re$ratio_ci_right

  }
  result <- list(ybar_cluster = ybar_cluster,
                 ybar_cluster_var_est = ybar_cluster_var_est,
                 ybar_cluster_se_est = ybar_cluster_se_est,
                 ybar_cluster_ci_left = ybar_cluster_ci_left,
                 ybar_cluster_ci_right = ybar_cluster_ci_right)
  return(result)
}
