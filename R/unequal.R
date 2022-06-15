#' Hansen-Hurwitz Estimator for Population Total
#'
#' This function estimates the population total with Hansen-Hurwitz estimator
#' when PPS sampling is used.
#'
#' @param y a vector, the sample data.
#' @param z a vector, the probability that the corresponding sample is selected
#' in a single sampling.
#' @param alpha a number, the confidence level.
#' @export
PPS_HH_total <- function(y, z, alpha = 0.05){
  n <- length(y)
  y_HH <- mean(y/z)
  y_HH_var_est <- var(y/z)/(n-1)
  y_HH_se_est <- sqrt(y_HH_var_est)
  y_HH_ci_est <- conf_interval(y_HH,y_HH_se_est,alpha = alpha)
  y_HH_ci_left <- y_HH_ci_est$ci_left_bound
  y_HH_ci_right <- y_HH_ci_est$ci_right_bound
  result <- list(y_HH = y_HH,
                 y_HH_var_est = y_HH_var_est,
                 y_HH_se_est = y_HH_se_est,
                 y_HH_ci_left = y_HH_ci_left,
                 y_HH_ci_right = y_HH_ci_right)
  return(result)
}

#' Hansen-Hurwitz Estimator for Population Mean
#'
#' This function estimates the population mean with Hansen-Hurwitz estimator
#' when PPS sampling is used.
#'
#' @param y a vector, the sample data.
#' @param z a vector, the probability that the corresponding sample is selected
#' in a single sampling.
#' @param N a number, the size of population.
#' @param alpha a number, the confidence level.
#' @export
PPS_HH_mean <- function(y,z,N,alpha = 0.05){
  n <- length(y)
  ybar_HH <- mean(y/z)/N
  ybar_HH_var_est <- var(y/z)/(N^2*(n-1))
  ybar_HH_se_est <- sqrt(ybar_HH_var_est)
  ybar_HH_ci_est <- conf_interval(ybar_HH,ybar_HH_se_est,alpha = alpha)
  ybar_HH_ci_left <- ybar_HH$ci_left_bound
  ybar_HH_ci_right <- ybar_HH$ci_right_bound
  result <- list(ybar_HH = ybar_HH,
                 ybar_HH_var_est = ybar_HH_var_est,
                 ybar_HH_se_est = ybar_HH_se_est,
                 ybar_HH_ci_left = ybar_HH_ci_left,
                 ybar_HH_ci_right = ybar_HH_ci_right)
  return(result)
}

#' Horvitz-Thompson Estimator for Population Total
#'
#' This function estimates the population total with Horvitz-Thompson estimator
#' when PPS sampling is used.
#'
#' @param y a vector, the sample data.
#' @param PI a vector, the inclusion probability probability.
#' @param PIM a matrix, the joint inclusion probability.
#' @param alpha a number, the confidence level.
#' @export
PiPS_HT_total <- function(y,PI,PIM,alpha = 0.05){
  n <- length(y)
  y_HT <- sum(y/PI)
  # this is only the first term of the estimated variance
  y_HT_var_est <- sum((1-PI)/PI^2*y^2)
  # the second term
  mediate <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      mediate <- mediate+
        2*(PIM[i,j]-PI[i]*PI[j])/(PI[i]*PI[j]*PIM[i,j])*y[i]*y[j]
    }
  }
  y_HT_var_est <- y_HT_var_est+mediate
  y_HT_se_est <- sqrt(y_HT_var_est)
  y_HT_ci_est <- conf_interval(y_HT,y_HT_se_est,alpha = alpha)
  y_HT_ci_left <- y_HT_ci_est$ci_left_bound
  y_HT_ci_right <- y_HT_ci_est$ci_right_bound
  result <- list(y_HT = y_HT,
                 y_HT_var_est = y_HT_var_est,
                 y_HT_se_est = y_HT_se_est,
                 y_HT_ci_left = y_HT_ci_left,
                 y_HT_ci_right = y_HT_ci_right)
  return(result)
}

#' Horvitz-Thompson Estimator for Population Total
#'
#' This function estimates the population total with Horvitz-Thompson estimator
#' when PPS sampling is used.
#'
#' @param y a vector, the sample data.
#' @param PI a vector, the inclusion probability.
#' @param PIM a matrix, the joint inclusion probability.
#' @param N a number, the size of population.
#' @param alpha a number, the confidence level.
#' @export
PiPS_HT_mean <- function(y,PI,PIM,N,alpha = 0.05){
  n <- length(y)
  ybar_HT <- sum(y/PI)/N
  # this is only the first term of the estimated variance
  ybar_HT_var_est <- sum((1-PI)/PI^2*y^2)
  # the second term
  mediate <- 0
  for (i in 1:n) {
    for (j in (i+1):n) {
      mediate <- mediate+
        2*(PIM[i,j]-PI[i]*PI[j])/(PI[i]*PI[j]*PIM[i,j])*y[i]*y[j]
    }
  }
  ybar_HT_var_est <- (ybar_HT_var_est+mediate)/N^2
  ybar_HT_se_est <- sqrt(ybar_HT_var_est)
  ybar_HT_ci_est <- conf_interval(ybar_HT,ybar_HT_se_est,alpha = alpha)
  ybar_HT_ci_left <- ybar_HT_ci_est$ci_left_bound
  ybar_HT_ci_right <- ybar_HT_ci_est$ci_right_bound
  result <- list(ybar_HT = ybar_HT,
                 ybar_HT_var_est = ybar_HT_var_est,
                 ybar_HT_se_est = ybar_HT_se_est,
                 ybar_HT_ci_left = ybar_HT_ci_left,
                 ybar_HT_ci_right = ybar_HT_ci_right)
  return(result)
}


#' Estimate Population Total Under Rao-Hartley-Cochran sampling
#'
#' This function estimates the population total when RHC sampling is used.
#'
#' @param ysample a vector, the sample data.
#' @param yprob a vector, the probability that the corresponding sample is selected
#' into sample.
#' @param group_prob_sum a vector, the sum of probability in a group where the
#' sample is selected.
#' @param group_N a vector, the size of the group where the sample is selected.
#' @param alpha a number, the confidence level.
#' @export
RHC_total <- function(ysample, yprob, group_prob_sum, group_N, alpha = 0.05){
  N <- sum(group_N)
  y_RHC <- sum(ysample/yprob*group_prob_sum)
  var_coef <- (sum(group_N^2)-N)/(N^2-sum(group_N^2))
  y_RHC_var_est <- var_coef*(sum(group_prob_sum*(ysample/yprob-y_RHC)^2))
  y_RHC_se_est <- sqrt(y_RHC_var_est)
  y_RHC_ci_est <- conf_interval(y_RHC, y_RHC_se_est, alpha = alpha)
  y_RHC_ci_left <- y_RHC_ci_est$ci_left_bound
  y_RHC_ci_right <- y_RHC_ci_est$ci_right_bound
  result <- list(y_RHC = y_RHC,
                 y_RHC_var_est = y_RHC_var_est,
                 y_RHC_se_est = y_RHC_se_est,
                 y_RHC_ci_left = y_RHC_ci_left,
                 y_RHC_ci_right = y_RHC_ci_right)
  return(result)
}

#' Estimate Population Mean Under Rao-Hartley-Cochran sampling
#'
#' This function estimates the population mean when RHC sampling is used.
#'
#' @param ysample a vector, the sample data.
#' @param yprob a vector, the probability that the corresponding sample is selected
#' into sample.
#' @param group_prob_sum a vector, the sum of probability in a group where the
#' sample is selected.
#' @param group_N a vector, the size of the group where the sample is selected.
#' @param alpha a number, the confidence level.
#' @export
RHC_mean <- function(ysample, yprob, group_prob_sum, group_N, alpha = 0.05){
  N <- sum(group_N)
  y_RHC <- sum(ysample/yprob*group_prob_sum)
  ybar_RHC <- y_RHC/N
  var_coef <- (sum(group_N^2)-N)/(N^2*(N^2-sum(group_N^2)))
  ybar_RHC_var_est <- var_coef*(sum(group_prob_sum*(ysample/yprob-y_RHC)^2))
  ybar_RHC_se_est <- sqrt(ybar_RHC_var_est)
  ybar_RHC_ci_est <- conf_interval(ybar_RHC, ybar_RHC_se_est, alpha = alpha)
  ybar_RHC_ci_left <- ybar_RHC_ci_est$ci_left_bound
  ybar_RHC_ci_right <- ybar_RHC_ci_est$ci_right_bound
  result <- list(ybar_RHC = ybar_RHC,
                 ybar_RHC_var_est = ybar_RHC_var_est,
                 ybar_RHC_se_est = ybar_RHC_se_est,
                 ybar_RHC_ci_left = ybar_RHC_ci_left,
                 ybar_RHC_ci_right = ybar_RHC_ci_right)
  return(result)
}

# # Example test
# ysample <- c(312.43,531.78,572.25)
# yprob <- c(40,60,93)/795
# group_prob_sum <- c(208,259,328)/795
# group_N <- c(5,5,5)
# RHC_mean(ysample = ysample, yprob = yprob,group_prob_sum = group_prob_sum,
#          group_N = group_N)
# RHC_total(ysample = ysample, yprob = yprob,group_prob_sum = group_prob_sum,
#          group_N = group_N)



