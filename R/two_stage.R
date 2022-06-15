#' Estimate Population Mean Under Two-Stage Sampling
#'
#' This function estimates the population mean when two-stage sampling is used and
#' the cluster sizes are equal.
#'
#' @param M a number, the number of secondary units in each primary unit in population.
#' @param m a number, the number of secondary units in each primary unit in sample.
#' @param ybari a vector, the sample mean of each primary unit in sample.
#' @param s2i a vector, the sample variance of each primary unit in sample.
#' @param N a number, the number of primary units in population.
#' @param alpha a number, the confidence level.
two_stage_eqM_mean <- function(M, m, ybari, s2i, N = Inf, alpha = 0.05){
  n <- length(ybari)
  f1 <- n/N
  f2 <- m/M
  s12 <- var(ybari)
  s22 <- mean(s2i)
  ybar_ts_est <- mean(ybari)
  ybar_ts_var_est <- (1-f1)/n*s12+f1*(1-f2)/(n*m)*s22
  ybar_ts_se_est <- sqrt(ybar_ts_var_est)
  ybar_ts_ci_est <- conf_interval(ybar_ts_est,ybar_ts_se_est,alpha = alpha)
  ybar_ts_ci_left <- ybar_ts_ci_est$ci_left_bound
  ybar_ts_ci_right <- ybar_ts_ci_est$ci_right_bound
  result <- list(ybar_ts_est = ybar_ts_est,
                 ybar_ts_var_est = ybar_ts_var_est,
                 ybar_ts_se_est = ybar_ts_se_est,
                 ybar_ts_ci_left = ybar_ts_ci_left,
                 ybar_ts_ci_right = ybar_ts_ci_right)
  return(result)
}

#' Estimate Population Total Under Two-Stage Sampling
#'
#' This function estimates the population total when two-stage sampling is used and
#' the cluster sizes are equal.
#'
#' @param M a number, the number of secondary units in each primary unit in population.
#' @param m a number, the number of secondary units in each primary unit in sample.
#' @param ybari a vector, the sample mean of each primary unit in sample.
#' @param s2i a vector, the sample variance of each primary unit in sample.
#' @param N a number, the number of primary units in population.
#' @param alpha a number, the confidence level.
two_stage_eqM_total <- function(M, m, ybari, s2i, N, alpha = 0.05){
  n <- length(ybari)
  f1 <- n/N
  f2 <- m/M
  s12 <- var(ybari)
  s22 <- mean(s2i)
  y_ts_est <- mean(ybari)*N
  y_ts_var_est <- ((1-f1)/n*s12+f1*(1-f2)/(n*m)*s22)*N^2
  y_ts_se_est <- sqrt(y_ts_var_est)
  y_ts_ci_est <- conf_interval(y_ts_est,y_ts_se_est,alpha = alpha)
  y_ts_ci_left <- y_ts_ci_est$ci_left_bound
  y_ts_ci_right <- y_ts_ci_est$ci_right_bound
  result <- list(y_ts_est = y_ts_est,
                 y_ts_var_est = y_ts_var_est,
                 y_ts_se_est = y_ts_se_est,
                 y_ts_ci_left = y_ts_ci_left,
                 y_ts_ci_right = y_ts_ci_right)
  return(result)
}

#' Estimate Population Proportion Under Two-Stage Sampling
#'
#' This function estimates the population proportion when two-stage sampling is used and
#' the cluster sizes are equal.
#'
#' @param M a number, the number of secondary units in each primary unit in population.
#' @param m a number, the number of secondary units in each primary unit in sample.
#' @param s2i a vector, the sample proportion of each primary unit in sample.
#' @param N a number, the number of primary units in population.
#' @param alpha a number, the confidence level.
two_stage_eqM_prop <- function(M, m, pj, N = Inf, alpha = 0.05){
  n <- length(pj)
  f1 <- n/N
  f2 <- m/M
  qj <- 1-pj
  p <- mean(pj)
  p_var_est <- (1-f1)/n*var(pj)+f1*(1-f2)/(n^2*(m-1))*sum(pj*qj)
  p_se_est <- sqrt(p_var_est)
  p_ci_est <- conf_interval(p,p_se_est,alpha = alpha)
  p_ci_left <- p_ci_est$ci_left_bound
  p_ci_right <- p_ci_est$ci_right_bound
  result <- list(p = p,
                 p_var_est = p_var_est,
                 p_se_est = p_se_est,
                 p_ci_left = p_ci_left,
                 p_ci_right = p_ci_right)
  return(result)
}

#' Estimate Population Total Under Two-Stage Sampling
#'
#' This function estimates the population total when two-stage sampling is used and
#' the cluster sizes are unequal.
#'
#' @param M0 a number, the number of secondary units in population.
#' @param M2i a vector, the number of secondary units in each primary unit that is drawn into sample.
#' @param mi a vector, the number of secondary units in each primary unit in sample.
#' @param ybari a vector, the sample mean of each primary unit in sample.
#' @param s2i a vector, the sample variance of each primary unit in sample.
#' @param N a number, the number of primary units in population.
#' @param alpha a number, the confidence level.
two_stage_uneqM_total <- function(M0, M2i, mi, ybari, s2i, N, alpha = 0.05){
  n <- length(ybari)
  f1 <- n/N
  f2i <- mi/M2i
  ybar_ts_est <- sum(M2i*ybari)/sum(M2i)
  y_ts_est <- M0*ybar_ts_est
  var_est_dep1 <- N^2*(1-f1)/n*sum(M2i^2*(ybari-ybar_ts_est)^2)/(n-1)
  var_est_dep2 <- N/n*sum(M2i^2*(1-f2i)*s2i/mi)
  y_ts_var_est <- var_est_dep1+var_est_dep2
  y_ts_se_est <- sqrt(y_ts_var_est)
  y_ts_ci_est <- conf_interval(y_ts_est,y_ts_se_est,alpha = alpha)
  y_ts_ci_left <- y_ts_ci_est$ci_left_bound
  y_ts_ci_right <- y_ts_ci_est$ci_right_bound
  result <- list(y_ts_est = y_ts_est,
                 ybar_ts_est = ybar_ts_est,
                 y_ts_var_est = y_ts_var_est,
                 y_ts_se_est = y_ts_se_est,
                 y_ts_ci_left = y_ts_ci_left,
                 y_ts_ci_right = y_ts_ci_right)
  return(result)
}

#' Estimate Population Mean Under Two-Stage Sampling
#'
#' This function estimates the population mean when two-stage sampling is used and
#' the cluster sizes are unequal.
#'
#' @param M0 a number, the number of secondary units in population.
#' @param M2i a vector, the number of secondary units in each primary unit that is drawn into sample.
#' @param mi a vector, the number of secondary units in each primary unit in sample.
#' @param ybari a vector, the sample mean of each primary unit in sample.
#' @param s2i a vector, the sample variance of each primary unit in sample.
#' @param N a number, the number of primary units in population.
#' @param alpha a number, the confidence level.
two_stage_uneqM_mean <- function(M0, M2i, mi, n, ybari, s2i, N = Inf,
                                  alpha = 0.05){
  n <- length(ybari)
  f1 <- n/N
  f2i <- mi/M2i
  ybar_ts_est <- sum(M2i*ybari)/sum(M2i)
  var_est_dep1 <- N^2*(1-f1)/n*sum(M2i^2*(ybari-ybar_ts_est)^2)/(n-1)
  var_est_dep2 <- N/n*sum(M2i^2*(1-f2i)*s2i/mi)
  ybar_ts_var_est <- (var_est_dep1+var_est_dep2)/M0^2
  ybar_ts_se_est <- sqrt(ybar_ts_var_est)
  ybar_ts_ci_est <- conf_interval(ybar_ts_est,ybar_ts_se_est,alpha = alpha)
  ybar_ts_ci_left <- ybar_ts_ci_est$ci_left_bound
  ybar_ts_ci_right <- ybar_ts_ci_est$ci_right_bound
  result <- list(ybar_ts_est = ybar_ts_est,
                 ybar_ts_var_est = ybar_ts_var_est,
                 ybar_ts_se_est = ybar_ts_se_est,
                 ybar_ts_ci_left = ybar_ts_ci_left,
                 ybar_ts_ci_right = ybar_ts_ci_right)
  return(result)
}

# two_stage_var <- function(Mi, mi, Ybari, S2i, n){
#   N <- length(Ybari)
#   f1 <- n/N
#   f2i <- mi/Mi
#   Ybarbar <- sum(Mi*Ybari)/sum(Mi)
#   var_dep1 <- N^2*(1-f1)/(n*(N-1))*sum(Mi^2*(Ybari-Ybarbar)^2)
#   var_dep2 <- N/n*sum(Mi^2*(1-f2i)*S2i/mi)
#   Var <- var_dep1+var_dep2
#   return(Var)
# }
