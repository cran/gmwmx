## ---- warning=F, message=F----------------------------------------------------
library(gmwmx)

## ---- echo = F, eval=F--------------------------------------------------------
#  dobs = gmwmx::PBO_get_station("DOBS", column = "dN")
#  write.gnssts(dobs, filename = "data_dobs.mom")

## ---- echo=F------------------------------------------------------------------
file_path = system.file("extdata", "data_dobs.mom", package = "gmwmx", mustWork = T)
data_dobs = read.gnssts(filename = file_path)

## -----------------------------------------------------------------------------
data_dobs = read.gnssts(filename = file_path)

## -----------------------------------------------------------------------------
class(data_dobs)

## -----------------------------------------------------------------------------
str(data_dobs)

## ---- fig.height=5, fig.align='center', fig.width=6---------------------------
plot(data_dobs$t, data_dobs$y, type="l")

## -----------------------------------------------------------------------------
fit_dobs_wn_plp_gmwmx = estimate_gmwmx(x = data_dobs, theta_0 = c(0.1, 0.1, 0.1), 
                                       model_string = "wn+powerlaw", 
                                       n_seasonal = 1, ci = T)

## -----------------------------------------------------------------------------
class(fit_dobs_wn_plp_gmwmx)

## -----------------------------------------------------------------------------
print(fit_dobs_wn_plp_gmwmx)
fit_dobs_wn_plp_gmwmx$beta_hat
fit_dobs_wn_plp_gmwmx$theta_hat

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
plot(fit_dobs_wn_plp_gmwmx)

## ---- eval=F, echo=T----------------------------------------------------------
#  fit_dobs_wn_plp_gmwmx_2 = estimate_gmwmx(x = data_dobs, theta_0 = c(0.1, 0.1, 0.1), model_string = "wn+powerlaw", n_seasonal = 1, k_iter = 2)

## ---- eval=F, echo=T----------------------------------------------------------
#  fit_dobs_wn_plp_mle = estimate_hector(x = data_dobs,
#                                        model_string = "wn+powerlaw",
#                                        n_seasonal = 1)
#  
#  
#  
#  

## ---- eval=T, echo=F----------------------------------------------------------
file_path_mle = system.file("extdata", "fit_dobs_wn_plp_mle.rda", package = "gmwmx", mustWork = T)
load(file_path_mle)

## ---- fig.height=8, fig.width=6, fig.align='center'---------------------------
plot(fit_dobs_wn_plp_mle)
fit_dobs_wn_plp_mle$beta_hat
fit_dobs_wn_plp_mle$theta_hat

## ---- eval=F, echo=T----------------------------------------------------------
#  cola = PBO_get_station("COLA", column = "dE")

## ---- eval=F, echo=F----------------------------------------------------------
#  save(cola, file="cola.rda")

## ---- eval=T, echo=F----------------------------------------------------------
cola_path = system.file("extdata", "cola.rda", package = "gmwmx", mustWork = T)
load(cola_path)

## ---- fig.height=8, fig.width=6, fig.align='center', eval=T-------------------
fit_cola_wn_plp = estimate_gmwmx(cola, model_string = "wn+powerlaw", 
                                 theta_0 = c(0.1,0.1,0.1),
                                 n_seasonal = 1, 
                                 ci = T)
plot(fit_cola_wn_plp)

fit_cola_wn_fgn = estimate_gmwmx(cola, model_string = "wn+fgn", theta_0 = c(0.1,0.1,0.2),
                                                          n_seasonal = 1, 
                                                          ci = T)
plot(fit_cola_wn_fgn)

fit_cola_wn_matern = estimate_gmwmx(cola, model_string = "wn+matern", 
                                    theta_0 = c(0.1,0.1,0.1,0.1),
                                                          n_seasonal = 1, 
                                                          ci = T)
plot(fit_cola_wn_matern)


## ---- fig.height=8, fig.width=8, fig.align='center', eval=T-------------------
compare_fits(fit_cola_wn_plp, fit_cola_wn_matern)

