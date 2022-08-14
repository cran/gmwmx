## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval=F,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#  # libraries
#  library(simts)
#  library(gmwmx)
#  library(dplyr)

## -----------------------------------------------------------------------------
#  check_if_theta_in_ci <- function(theta_vec, ci_mat) {
#    vec_emp_coverage <- vector(mode = "logical", length = length(theta_vec))
#    for (i in seq(length(theta_vec))) {
#      if (dplyr::between(theta_vec[i], ci_mat[i, 1], ci_mat[i, 2])) {
#        vec_emp_coverage[i] <- T
#      } else {
#        vec_emp_coverage[i] <- F
#      }
#    }
#    return(as.numeric(vec_emp_coverage))
#  }

## -----------------------------------------------------------------------------
#  # we consider example the model considered in model 3
#  phase <- 0.45
#  amplitude <- 2.5
#  sigma2_wn <- 15
#  sigma2_powerlaw <- 10
#  d <- 0.4
#  bias <- 0
#  trend <- 5 / 365.25
#  cosU <- amplitude * cos(phase)
#  sinU <- amplitude * sin(phase)

## -----------------------------------------------------------------------------
#  # consider n years of daily observations
#  year <- 20
#  n <- year * 365

## ---- eval=F------------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("SMAC-Group/simts")

## -----------------------------------------------------------------------------
#  # define model for generating gaussian white noise + PLP
#  model_gaussian_wn_plp <- WN(sigma2 = sigma2_wn) + PLP(sigma2 = sigma2_powerlaw, d = d)

## -----------------------------------------------------------------------------
#  # generate data
#  # define time at which there are jumps
#  jump_vec <- c(200, 300, 500)
#  jump_height <- c(10, 15, 20)

## -----------------------------------------------------------------------------
#  # define seed
#  myseed <- 123
#  
#  # add trend, gaps and sin
#  nbr_sin <- 1 # define number of sinusoidal process

## -----------------------------------------------------------------------------
#  # define A matrix
#  A <- create_A_matrix(1:n, jump_vec, n_seasonal = nbr_sin)

## -----------------------------------------------------------------------------
#  # define beta
#  x_0 <- c(bias, trend, cosU, sinU, jump_height)

## -----------------------------------------------------------------------------
#  # define number of Monte Carlo simulation
#  n_simu <- 1000
#  
#  # define number of parameter estimated (depends on the model)
#  # bias + trend + height * nbr of jump vec + 2* nbr of sin process + wn + powerlaw parameters
#  nbr_param_check_coverage <- 4
#  nbr_param <- 2 + length(jump_vec) + 2 * nbr_sin + 3
#  dim_mat_results <- (nbr_param + nbr_param_check_coverage) * 3 + 1 * 3
#  
#  # define matrix of results
#  mat_results <- matrix(NA, ncol = dim_mat_results, nrow = n_simu)

## -----------------------------------------------------------------------------
#  for (simu_b in seq(n_simu)) {
#  
#    # fix seed for reproducibility
#    set.seed(myseed + simu_b)
#  
#    # generate residuals from a Gaussian White noise + Power law process
#    eps <- simts::gen_gts(model = model_gaussian_wn_plp, n = n)
#  
#    # create time series
#    yy <- A %*% x_0 + eps
#  
#    # create gnssts
#    gnssts_obj <- create.gnssts(t = 1:length(yy), y = yy, jumps = jump_vec)
#  
#    # MLE
#    fit_simu_b_mle <- estimate_hector(
#      x = gnssts_obj,
#      model = "wn+powerlaw",
#      n_seasonal = 1
#    )
#  
#    # gmwmx 1 step
#    fit_simu_b_gmwm_1_step <- estimate_gmwmx(
#      x = gnssts_obj,
#      model = "wn+powerlaw",
#      theta_0 = c(0.1, 0.1, 0.1),
#      n_seasonal = 1,
#      ci = T,
#      k_iter = 1
#    )
#  
#    # gmwmx 2 steps
#    fit_simu_b_gmwm_2_step <- estimate_gmwmx(
#      x = gnssts_obj,
#      model = "wn+powerlaw",
#      theta_0 = c(0.1, 0.1, 0.1),
#      n_seasonal = 1,
#      ci = T,
#      k_iter = 2
#    )
#  
#  
#    # define vector of estimators
#    fit_mle <- c(fit_simu_b_mle$beta_hat, fit_simu_b_mle$theta_hat)
#    fit_gmwm_1_step <- c(fit_simu_b_gmwm_1_step$beta_hat, fit_simu_b_gmwm_1_step$theta_hat)
#    fit_gmwm_2_step <- c(fit_simu_b_gmwm_2_step$beta_hat, fit_simu_b_gmwm_2_step$theta_hat)
#  
#    # compute coverage for each method
#    alpha <- .05
#    z_val <- qnorm(1 - alpha / 2)
#  
#    mat_ci_mle_beta <- matrix(c(
#      fit_simu_b_mle$beta_hat - z_val * fit_simu_b_mle$beta_std,
#      fit_simu_b_mle$beta_hat + z_val * fit_simu_b_mle$beta_std
#    ),
#    byrow = F, ncol = 2
#    )
#  
#    mat_ci_gmwm_1_step_beta <- matrix(c(
#      fit_simu_b_gmwm_1_step$beta_hat - z_val * fit_simu_b_gmwm_1_step$beta_std,
#      fit_simu_b_gmwm_1_step$beta_hat + z_val * fit_simu_b_gmwm_1_step$beta_std
#    ),
#    byrow = F, ncol = 2
#    )
#  
#    mat_ci_gmwm_2_step_beta <- matrix(c(
#      fit_simu_b_gmwm_2_step$beta_hat - z_val * fit_simu_b_gmwm_2_step$beta_std,
#      fit_simu_b_gmwm_2_step$beta_hat + z_val * fit_simu_b_gmwm_2_step$beta_std
#    ),
#    byrow = F, ncol = 2
#    )
#  
#    # save empirical coverage
#    inside_ci_mle <- check_if_theta_in_ci(x_0, mat_ci_mle_beta)[1:4]
#    inside_ci_gmwm_1 <- check_if_theta_in_ci(x_0, mat_ci_gmwm_1_step_beta)[1:4]
#    inside_ci_gmwm_2 <- check_if_theta_in_ci(x_0, mat_ci_gmwm_2_step_beta)[1:4]
#  
#    # define vector of estimated parameters as object
#    res <- c(
#      fit_mle, fit_gmwm_1_step,
#      fit_gmwm_2_step,
#      inside_ci_mle,
#      inside_ci_gmwm_1, inside_ci_gmwm_2,
#      fit_simu_b_mle$estimation_time[3], fit_simu_b_gmwm_1_step$estimation_time[3], fit_simu_b_gmwm_2_step$estimation_time[3]
#    )
#  
#    # save in matrix of results
#    mat_results[simu_b, ] <- res
#  
#    # print status
#    cat(paste("Completed simulation", simu_b, "\n",sep = " "))
#  }

## -----------------------------------------------------------------------------
#  # define names of mat results
#  name_param_functionnal <- c("bias", "trend", "cosU", "sinU")
#  colnames(mat_results) <- c(
#      paste("mle", names(fit_mle), sep = "_"),
#      paste("gmwm_1_step", names(fit_gmwm_1_step), sep = "_"),
#      paste("gmwm_2_step", names(fit_gmwm_2_step), sep = "_"),
#      paste0("mle_inside_ci_", name_param_functionnal),
#      paste0("gmwm_1_inside_ci_", name_param_functionnal),
#      paste0("gmwm_2_inside_ci_", name_param_functionnal),
#      c("time_mle", "time_gmwm_1", "time_gmwm_2")
#    )

