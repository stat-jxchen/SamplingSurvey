#' Estimating the Subpopulation Mean
#'
#' This function is used to estimate the subpopulation mean.
#'
#' @param y_sub a vector, a sample of a subpopulation.
#' @param n a number, the sample size of the entire population.
#' @param N a number, the size of the entire population.
#' @param alpha a number, the confidence level.
#' @note This function is for the data that have origin data.
#' @export
srs_subpop_mean1 <- function(y_sub, n, N = Inf, alpha = 0.05){
  f <- n/N
  n_sub <- length(y_sub)
  ybar_sub <- mean(y_sub)
  ybar_sub_var_est <- (1-f)/n_sub*var(y_sub)
  ybar_sub_se_est <- sqrt(ybar_sub_var_est)
  ybar_sub_ci_est <- conf_interval(ybar_sub,ybar_sub_se_est,alpha = alpha)
  ybar_sub_ci_left <- ybar_sub_ci_est$ci_left_bound
  ybar_sub_ci_right <- ybar_sub_ci_est$ci_right_bound
  result <- list(ybar_sub = ybar_sub,
                 ybar_sub_var_est = ybar_sub_var_est,
                 ybar_sub_se_est = ybar_sub_se_est,
                 ybar_sub_ci_left = ybar_sub_ci_left,
                 ybar_sub_ci_right = ybar_sub_ci_right)
  return(result)
}

#' Estimating the Subpopulation Total
#'
#' This function is used to estimate the subpopulation total.
#'
#' @param y_sub a vector, a sample of a subpopulation.
#' @param n a number, the sample size of the entire population.
#' @param N a number, the size of the entire population.
#' @param alpha a number, the confidence level.
#' @note This function is for the data that have origin data.
#' @export
srs_subpop_total1 <- function(y_sub, n, N, alpha = 0.05){
  f <- n/N
  n_sub <- length(y_sub)
  p_sub <- n_sub/n
  q_sub <- 1-p_sub
  ybar_sub <- mean(y_sub)
  y_sub <- sum(y_sub)*N/n
  y_sub_var_est <- (1-f)*N^2/(n*(n-1))*
    ((n_sub-1)*var(y_sub)+n*p_sub*q_sub*ybar_sub^2)
  y_sub_se_est <- sqrt(y_sub_var_est)
  y_sub_ci_est <- conf_interval(y_sub,y_sub_se_est,alpha = alpha)
  y_sub_ci_left <- y_sub_ci_est$ci_left_bound
  y_sub_ci_right <- y_sub_ci_est$ci_right_bound
  result <- list(y_sub = y_sub,
                 y_sub_var_est = y_sub_var_est,
                 y_sub_se_est = y_sub_se_est,
                 y_sub_ci_left = y_sub_ci_left,
                 y_sub_ci_right = y_sub_ci_right)
  return(result)
}

#' Estimating the Subpopulation Mean
#'
#' This function is used to estimate the subpopulation mean.
#'
#' @param ybar_sub a number, the sample mean of a subpopulation.
#' @param s2_sub a number, the sample variance of a subpopulation.
#' @param n_sub a number, the sample size of a subpopulation.
#' @param n a number, the sample size of the entire population.
#' @param N a number, the size of the entire population.
#' @param alpha a number, the confidence level.
#' @note This function is for the data that only have some information from the
#' origin data.
#' @export
srs_subpop_mean2 <- function(ybar_sub, s2_sub, n_sub, n, N = Inf,
                             alpha = 0.05){
  f <- n/N
  ybar_sub_var_est <- (1-f)/n_sub*s2_sub
  ybar_sub_se_est <- sqrt(ybar_sub_var_est)
  ybar_sub_ci_est <- conf_interval(ybar_sub,ybar_sub_se_est,alpha = alpha)
  ybar_sub_ci_left <- ybar_sub_ci_est$ci_left_bound
  ybar_sub_ci_right <- ybar_sub_ci_est$ci_right_bound
  result <- list(ybar_sub = ybar_sub,
                 ybar_sub_var_est = ybar_sub_var_est,
                 ybar_sub_se_est = ybar_sub_se_est,
                 ybar_sub_ci_left = ybar_sub_ci_left,
                 ybar_sub_ci_right = ybar_sub_ci_right)
  return(result)
}

#' Estimating the Subpopulation Total
#'
#' This function is used to estimate the subpopulation total.
#'
#' @param ybar_sub a number, the sample mean of a subpopulation.
#' @param s2_sub a number, the sample variance of a subpopulation.
#' @param n_sub a number, the sample size of a subpopulation.
#' @param n a number, the sample size of the entire population.
#' @param N a number, the size of the entire population.
#' @param alpha a number, the confidence level.
#' @note This function is for the data that only have some information from the
#' origin data.
#' @export
srs_subpop_total2 <- function(ybar_sub, s2_sub, n_sub, n, N, alpha = 0.05){
  f <- n/N
  p_sub <- n_sub/n
  q_sub <- 1-p_sub
  y_sub <- n_sub*ybar_sub*N/n
  y_sub_var_est <- (1-f)*N^2/(n*(n-1))*
    ((n_sub-1)*s2_sub+n*p_sub*q_sub*ybar_sub^2)
  y_sub_se_est <- sqrt(y_sub_var_est)
  y_sub_ci_est <- conf_interval(y_sub,y_sub_se_est,alpha = alpha)
  y_sub_ci_left <- y_sub_ci_est$ci_left_bound
  y_sub_ci_right <- y_sub_ci_est$ci_right_bound
  result <- list(y_sub = y_sub,
                 y_sub_var_est = y_sub_var_est,
                 y_sub_se_est = y_sub_se_est,
                 y_sub_ci_left = y_sub_ci_left,
                 y_sub_ci_right = y_sub_ci_right)
  return(result)
}

#' Estimate Population Mean with Poststratification
#'
#' This function is used to estimate the population mean in simple random sampling with poststratification method.
#'
#' @param Wh a vector, the weight of each strata.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the sample mean of each strata.
#' @param s2h a vector, the sample variance of each strata.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
#' @export
poststra_mean <- function(Wh, nh, ybarh, s2h, N = Inf, alpha = 0.05){
  n <- sum(nh)
  f <- n/N
  ybar_pst <- sum(Wh*ybarh)
  ybar_pst_var_est <- (1-f)/n*sum(Wh*s2h)+1/n^2*sum((1-Wh)*s2h)
  ybar_pst_se_est <- sqrt(ybar_pst_var_est)
  ybar_pst_ci_est <- conf_interval(ybar_pst,ybar_pst_se_est,alpha = alpha)
  ybar_pst_ci_left <- ybar_pst_ci_est$ci_left_bound
  ybar_pst_ci_right <- ybar_pst_ci_est$ci_right_bound
  result <- list(ybar_pst = ybar_pst,
                 ybar_pst_var_est = ybar_pst_var_est,
                 ybar_pst_se_est = ybar_pst_se_est,
                 ybar_pst_ci_left = ybar_pst_ci_left,
                 ybar_pst_ci_right = ybar_pst_ci_right)
  return(result)
}

#' Estimate Population Total with Poststratification
#'
#' This function is used to estimate the population total in simple random sampling with poststratification method.
#'
#' @param Wh a vector, the weight of each strata.
#' @param nh a vector, the sample size of each strata.
#' @param ybarh a vector, the sample mean of each strata.
#' @param s2h a vector, the sample variance of each strata.
#' @param N a number, the population size.
#' @param alpha a number, the confidence level.
#' @export
poststra_total <- function(Wh, nh, ybarh, s2h, N, alpha = 0.05){
  n <- sum(nh)
  f <- n/N
  y_pst <- sum(Wh*ybarh)*N
  y_pst_var_est <- ((1-f)/n*sum(Wh*s2h)+1/n^2*sum((1-Wh)*s2h))*N^2
  y_pst_se_est <- sqrt(y_pst_var_est)
  y_pst_ci_est <- conf_interval(y_pst,y_pst_se_est,alpha = alpha)
  y_pst_ci_left <- y_pst_ci_est$ci_left_bound
  y_pst_ci_right <- y_pst_ci_est$ci_right_bound
  result <- list(y_pst = y_pst,
                 y_pst_var_est = y_pst_var_est,
                 y_pst_se_est = y_pst_se_est,
                 y_pst_ci_left = y_pst_ci_left,
                 y_pst_ci_right = y_pst_ci_right)
  return(result)
}
