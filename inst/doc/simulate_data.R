## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## -----------------------------------------------------------------------------
# libraries
library(simts)
library(gmwmx)

## ---- eval=T, echo=T, message=F-----------------------------------------------
phase =     0.45
amplitude = 2.5
sigma2_wn =       15
sigma2_powerlaw = 10
d =               0.4
bias =            0
trend =           5/365.25
cosU =            amplitude*cos(phase)
sinU =            amplitude*sin(phase)

## ---- eval=T, echo=T, message=F-----------------------------------------------
n = 5*365

## ---- eval=F------------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("SMAC-Group/simts")

## ---- eval=F, echo=T, message=F-----------------------------------------------
#  model_i = WN(sigma2 = sigma2_wn) + PLP(sigma2 = sigma2_powerlaw, d = d)

## ---- eval=T, echo=T, message=F-----------------------------------------------
# define time at which there are jumps
jump_vec =  c(600, 1200)
jump_height = c(20, 30)

# define myseed
myseed=123

## ---- eval=F, echo=T, message=F-----------------------------------------------
#  # generate residuals
#  eps = simts::gen_gts(model = model_i, n= n)

## ---- evaL=T, echo=F----------------------------------------------------------
file_path_eps_simulate_data = system.file("extdata", "eps.rda", package = "gmwmx", mustWork = T)
load(file_path_eps_simulate_data)

## ---- eval=T, echo=T, message=F-----------------------------------------------
# add trend and sin
A = gmwmx::create_A_matrix(t_nogap = 1:length(eps), n_seasonal =  1, jumps = NULL)

# define beta
x_0 = c(bias, trend,  cosU,  sinU)
  
# create time series
deterministic_signal = A %*% x_0

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
plot(deterministic_signal, type="l")
plot(eps)

## -----------------------------------------------------------------------------
# add trend, gaps and sin
A = gmwmx::create_A_matrix(t_nogap = 1:length(eps), n_seasonal =  1, jumps = jump_vec)

# define beta
x_0 = c(bias, trend,  cosU,  sinU, jump_height)
  
# create time series
deterministic_signal = A %*% x_0

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
plot(deterministic_signal, type="l")
plot(eps)

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
yy = deterministic_signal + eps
plot(yy)

## ---- eval=T------------------------------------------------------------------
# save signal in temp
gnssts_obj = create.gnssts(t = 1:length(yy), y = yy, jumps = jump_vec)

## ---- eval=T------------------------------------------------------------------
class(gnssts_obj)

## ---- eval=F, echo=T----------------------------------------------------------
#  write.gnssts(gnssts_obj, filename = "simulated_data.mom")

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
fit_gmwmx = gmwmx::estimate_gmwmx(x = gnssts_obj,
                                  model_string = "wn+powerlaw",
                                  n_seasonal = 1,
                                  theta_0 = c(0.1, 0.1, 0.1),
                                  k_iter = 1)

fit_gmwmx
plot(fit_gmwmx)

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
fit_gmwmx_2 = gmwmx::estimate_gmwmx(x = gnssts_obj,
                                    model_string = "wn+powerlaw",
                                    n_seasonal = 1,
                                    theta_0 = c(0.1, 0.1, 0.1),
                                    k_iter = 2)

fit_gmwmx_2
plot(fit_gmwmx_2)

## ---- eval=F, echo=T----------------------------------------------------------
#  fit_mle_hector = gmwmx::estimate_hector(x = gnssts_obj,
#                                      model_string = "wn+powerlaw",
#                                      n_seasonal = 1
#                                      )
#  
#  

## ---- eval=F, echo=F----------------------------------------------------------
#  save(fit_mle_hector,file="fit_mle_hector.rda")

## ---- eval=T, echo=F----------------------------------------------------------
# fit_mle_hector is save in inst/extdata
file_path = system.file("extdata", "fit_mle_hector.rda", package = "gmwmx", mustWork = T)
load(file = file_path)

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
fit_mle_hector
plot(fit_mle_hector)

