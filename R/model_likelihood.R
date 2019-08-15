#' Outpout model likelihood
#'
#' @description Function to output the model likelihhood given a dataset, model, and vector of parameters. Note this does not include the prior.
#'
#' @param dataList dataList created by \code{create_dataList}
#' @param paramsVec a parameter vector of length 4 or 5 created by the first column output by \code{plotmod}. If model == 1 or 3, the order will be: DN, K600, sigma2, logpost. If model == 2 or 4 will be: K600, DN, Nfix, sigma2, logpost.
#' @param model Model number for fitting algorithm (see details bleow).
#'
#' @details Model determines which model to estimate. 0 is the oxygen model (Eq. 2 from Nifong et al), 1 is the single station model without N consumption (DN base model Eq. 3 from Nifong et al.), 2 is the single station model with N consumption (DN + Nconsume Eq. 4 from Nifong et al.), the  3 being the two-station model without N consumption (DN base Eq. 5 from Nifong et al.), and 4 being the two station model with N consumption (DN N consume Eq. 6 from Nifong et al.).
#'
#' @export
model_likelihood <- function(dataList, paramsVec, model, incl_DN = TRUE){

  #################################################################
  # FUNCTIONS
  #################################################################
  # Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  KcorN2 <- function( temp,  K600) {
    K600/(600/(1615-(temp*92.15)+(2.349*temp^2)-(0.0240*temp^3)))^-0.5
  }


  # Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  KcorO2 <- function( temp,  K600) {
    K600/(600/(1568-(temp*86.04)+(2.142*temp^2)-(0.0216*temp^3)))^-0.5
  }


  #################################################################
  # DATA
  #################################################################
  # Model specifications
  mod = model
  nobs = dataList$nobs
  ncol = dataList$ncol
  incl_DN = incl_DN

  # Data
  data_obj = dataList$data_obj

  # For models 0 data_obj is the following
  # -- Column 1 = tempup (celcius)
  # -- Column 2 = tempdown (celcius)
  # -- Column 3 = o2equilup (g m^-3)
  # -- Column 4 = o2equildown (g m^-3)
  # -- Column 5 = o2up (g m^-3)
  # -- Column 6 = O2down (g m^-3)

  # For models 1 and 2 data_obj is the following
  # -- Column 1 = temp (celcius)
  # -- Column 2 = n2equil (g m^-3)
  # -- Column 3 = n2 (g m^-3)
  # -- Column 4 = delta t (days^-1)

  # For models 3 and 4 data_obj is the following
  # -- Column 1 = tempup (celcius)
  # -- Column 2 = tempdown (celcius)
  # -- Column 3 = n2equilup (g m^-3)
  # -- Column 4 = n2equildown (g m^-3)
  # -- Column 5 = n2up (g m^-3)
  # -- Column 6 = n2down (g m^-3)

  # Fixed values
  z = dataList$z # Water depth (m)
  tt = dataList$tt # Travel time between stations
  lag = dataList$lag # Time lag (if mod < 3 it is 0)

  PPFD = dataList$PPFD; # Photosynthetic photon flux density (PPFD/light) (mol m-2 s-1)
  PPFDtotal = dataList$PPFDtotal # Daily total of PPFD (mol m-2 s-1 d-1)

  # Prior specifications
  Kmean = dataList$Kmean
  Ksd = dataList$Ksd


  #################################################################
  # Parameters
  #################################################################
  # -- For 02 models
  # params[1] = K = Gas exchange rates (d-1)
  # params[2] = GPP = Gross primary production  (g m-2 d-1)
  # params[3] = ER = Ecosystem respiration  (g m-2 d-1)

  # -- For N2 models
  # params[1] = K600 = Gas exchange rates (d-1)
  # params[2] = DN = Denitrification rate (g m-2 d-1)
  # params[3] = Nfix = N2 fixation rate (g m-2 d-1)

  if(mod %in% c(1,3)){
    params <- paramsVec[1:2]
    sigma <- paramsVec[3]
  }

  if(mod %in% c(0,2,4)){
    params <- paramsVec[1:3]
    sigma <- paramsVec[4]
  }

  #################################################################
  # Models
  #################################################################
  estVal <- c()

  # -- Two-station oxygen models (Models 0)
  if(mod == 0){
    # Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5] + params[3] * tt / z
      estVal[i] = estVal[i] + params[2] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal
      estVal[i] = estVal[i] + KcorO2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2
      estVal[i] = estVal[i] / (1 + KcorO2(data_obj[i, 1], params[1]) * tt / 2)
    }

    # Likelihood
    data_obj[,6] ~ normal( estVal, sigma)
  }


  # -- One-station models (Models 1 and 2)
  if(mod != 0){
    if(mod < 3){
      # Initialize
      estVal[1] = data_obj[1, 3]

      # Model loop
      for(i in 2:nobs){
        estVal[i] = estVal[i - 1]

        # Remove DN
        if(incl_DN == 1){
          estVal[i] = estVal[i] + params[2] * data_obj[i, 4] / z
        }

        estVal[i] = estVal[i] + KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] * (data_obj[i, 2] + data_obj[i - 1, 2] - estVal[i - 1])/2

        # Add N consumption (DN + Nconsume) for model 2
        if(mod == 2){
          estVal[i] = estVal[i] + params[2+incl_DN] / z * PPFD[i] / PPFDtotal
        }

        estVal[i] = estVal[i] / (1 + KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] / 2)
      }

      # Likelihood
      nll = -sum(dnorm(data_obj[2:nrow(data_obj),3], estVal[2:length(estVal)], sigma, TRUE))
    }
  }


  # -- Two-station models (Models 3 and 4)
  if(mod > 2){
    # Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5]

      # Remove DN
      if(incl_DN == 1){
        estVal[i] = estVal[i] +params[2] * tt / z
      }

      estVal[i] = estVal[i] + KcorN2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2

      # Add N consumption (DN + Nconsume) for model 4
      if(mod == 4){
        estVal[i] = estVal[i] +params[2+incl_DN] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal
      }

      estVal[i] = estVal[i] / (1 + KcorN2(data_obj[i, 1], params[1]) * tt / 2)
    }

    # Likelihood
    nll = -sum(dnorm(data_obj[,6], estVal, sigma, TRUE))
  }

  return(nll)
}
