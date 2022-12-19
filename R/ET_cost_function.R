#' Creates a cost function for different simulation setups based on RMSE
#' 
#' Creates a cost function for parameter calibration, keeping non-calibrated
#' parameter values fixed and calibrating the parameters corresponding to setups
#' \code{BRC} and \code{FULL} from Stocker et al., 2020 GMD. The cost function
#' computes the root mean squared error (RMSE) on the calibrated parameters.
#' 
#' @param params_modl A list of model parameter values, including \code{'kphio',
#' 'soilm_par_a', 'soilm_par_b', 'tau_acclim_tempstress' }and \code{'par_shape_tempstress'}
#' in that order.
#' @param setup A character string (\code{'BRC'} or \code{'FULL'}) indicating which
#' parameters are calibrated. For \code{setup = 'BRC'} only the quantum yield
#' efficiency \code{kphio} is calibrated; for \code{setup = 'FULL'} it also includes
#' the soil moisture stress parameters \code{soilm_par_a} and \code{soilm_par_b}
#' for calibration.
#' @param method A character string indicating the optimization method that will
#' be used, either \code{'BayesianTools'} or \code{'GenSA'}.
#' 
#' 
#' @importFrom magrittr '%>%'
#' 
#' @return A cost function which computes the RMSE of the simulated GPP by the P-model 
#' versus the observed GPP. This cost function has as arguments a list of calibratable
#' model parameters \code{par}, a data frame of observations \code{obs}, and a
#' data frame of driver data \code{drivers}.
#' 
#' @details The resulting cost function performs a P-model run for the value of
#' \code{par} given as argument and the remaining non-calibratable parameters
#' are held constant (specified via \code{params_modl}).
#' 
#' Since the calibration routine in \code{BayesianTools} is based on maximizing 
#' a cost function and we want to minimize the RMSE, the opposite value, 
#' \code{(-1)*RMSE}, is returned if \code{method = 'BayesianTools'}. \code{GenSA}
#' minimizes the given objective function, so the plain RMSE is returned when
#' \code{method = 'GenSA'}.
#' 
#' @export
#' 
#' @examples \dontrun{
#' # Set model parameters
#' pars <- list(
#'   kphio          = 0.04,
#'   soilm_par_a    = 2.8,
#'   soilm_par_b    = 1.7,
#'   tau_acclim_tempstress  = 7.3,
#'   par_shape_tempstress   = 0.1
#'   )
#' 
#' # Write cost function
#' cost_rmse_kphio <- create_cost_rmse_pmodel(
#'   params_modl = pars,
#'   setup = 'BRC',
#'   method = 'BayesianTools'
#'   )
#' }



cost_mae_ET <- function(
  par,
  obs,
  drivers,
  inverse = FALSE
){
  
  # predefine variables for CRAN check compliance
  sitename <- data <- NULL
  
  ## execute model for this parameter set tttwwwww
  ## For calibrating quantum yield efficiency only
  params_modl <- list(
    # kphio           = par[[1]],
    whc           = par[[1]],
    kalb_sw       = par[[2]],
    kw            = par[[3]],
    kCw           = par[[4]],
    kphio         = 0.04998,
    soilm_par_a     = 0.33349283,
    soilm_par_b     = 1.45602286,
    tau_acclim_tempstress = 10,
    par_shape_tempstress  = 0.0
    # soilm_par_a     = par[[2]],
    # soilm_par_b     = par[[3]],
    # tau_acclim_tempstress = par[[4]],
    # par_shape_tempstress  = par[[5]],
  )
  
  # run the model (with the previous model parameter list)
  df <- runread_pmodel_f(
    drivers, 
    par = params_modl,
    makecheck = TRUE,
    parallel = FALSE
  )
  
  # cleanup
  df <- df %>%
    dplyr::select(sitename, data) %>% 
    tidyr::unnest(data) %>%
    dplyr::rename(
      'LE_mod' = 'latenth'
    )
  
  obs <- obs %>%
    dplyr::select(sitename, data) %>% 
    tidyr::unnest(data) %>%
    dplyr::rename('LE_obs' = 'LE_F_MDS') # rename for later
  
  # left join with observations
  df <- dplyr::left_join(df, obs, by = c("sitename", "date"))
  
  # Calculate cost (RMSE)
  cost <- sqrt( mean( (df$LE_obs - df$LE_mod )^2, na.rm = TRUE ) )
  
  if (inverse) cost <- 1.0 / cost
  
  return(cost)
}
  