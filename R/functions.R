#' Water density
#'
#' @description Function to calculate the water density as a function of temperature.
#'
#' @param temp Water temperature (Celcius)
#'
#' @export
#'
#' @examples
#'
#' watdens(29)
#'
watdens<-function(temp){

  t<-temp
  A <- 7.0132e-5
  B <- 7.926295e-3
  C <-  -7.575477e-5
  D<- 7.314701e-7
  E <-  -3.596363e-9
  to<- 3.9818

  dens<- (999.97358- (A*(t-to) + B*(t-to)^2 +C*(t-to)^3 + D*(t-to)^4+E*(t-to)^5) ) -4.873e-3 + 1.708e-4*t - 3.108e-6 * t^2
  dens/1000
}


#' Oxygen saturation
#'
#' @description From Garcia and Gordon 1992 L&O. This calculation is a slight approximation because we do not account for water density, but the difference is much smaller (0.01 mg/L at 20 deg) than anyone's ability to calibrate an oxygen sonde.
#'
#' @param temp Water temperature (Celcius)
#' @param bp Barymetric pressure (elevation-corrected mm Hg)
#' @export
#'
#' @examples
#'
#' osat(26, 0)
osat <- function(temp,bp){
  sato<-(exp(2.00907 + 3.22014 * (log((298.15-temp) / (273.15 + temp))) + 4.0501 * (log((298.15 - temp) / (273.15 + temp))) ^ 2 + 4.94457 * (log((298.15 - temp) / (273.15 + temp))) ^ 3 - 0.256847 * (log((298.15 - temp) / (273.15 + temp))) ^ 4 + 3.88767 * (log((298.15 - temp) / (273.15 + temp))) ^ 5)) * 1.4276 * bp / 760

  return(sato)
}


#' N2 saturation
#'
#' @description Function to calculate the N2 saturation.
#'
#' @param temp Temperature (Celcius)
#' @param salinity Salinity
#' @param bp Barymetric pressure (elevation-corrected mm Hg)
#'
#' @export
#'
#' @examples
#'
#' n2satfc(15, 0, 760)
#'
n2satfc <- function(temp,salinity,bp){
  Ts=log((298.15-temp)/(273.15+temp))

  A0=6.42931
  A1=2.92704
  A2=4.32531
  A3=4.69149
  B0=-0.00744129
  B1=-0.00802566
  B2=-0.0146775

  ln.C = 6.42931 + 2.92704 * Ts + 4.32531 * Ts^2 + 4.69149 * Ts^3 + salinity*(-0.00744129 + -0.00802566 * Ts + -0.0146775 * Ts^2)
  u<-10^(8.10765-(1750.286/(235+temp)))
  C<-exp(ln.C)*((bp-u)/(760-u))
  converted<-C*watdens(temp)*(28.014/1000) # converts N2 from umolN2/kg to gN/m3
  return(converted)
}



#' Ar Saturation
#'
#' @description Function to calculate the ar saturation.
#'
#' @param temp Temperature (Celcius)
#' @param salinity Salinity
#' @param bp Barymetric pressure (elevation-corrected mm Hg)
#'
#' @export
#'
#' @examples
#'
#' arsatfc(15, 0, 760)
#'
arsatfc<-function(temp, salinity, bp){
  Ts=log((298.15-temp)/(273.15+temp))

  A0=2.79150
  A1=3.17609
  A2=4.13116
  A3=4.90379
  B0=-0.00696233
  B1=-0.00766670
  B2=-0.0116888

  ln.C = A0 + A1 * Ts + A2 * Ts^2 + A3 * Ts^3 + salinity*(B0 + B1 * Ts + B2 * Ts^2)
  u<-10^(8.10765-(1750.286/(235+temp)))
  C<-exp(ln.C)*((bp-u)/(760-u))
  converted<-C*watdens(temp)*(39.948/1000)# converts Ar from umol/kg to g/m3
  return(converted)
}


#' Data list
#'
#' @description Creates a list of data objects to be used by stan
#'
#' @param data data.frame object with columnds for station, temperature, n2.ar, light, bp, z, and dtime, arsat, and n2sat, or oxy in the case of model == 0.
#' @param model Model number for fitting algorithm (see details bleow).
#' @param Kmean Mean of normal prior distribution for K600
#' @param Ksd Sd of normal prior distribution for K600
#' @param up Name indicating the up-river station name
#' @param down Name indicating the down-river station name
#' @param tt Time between stations for two station model
#' @param lag_divisor The value for tt to be divided by to estimate lag time between stations for two-station model
#' @param depth Depth (m)
#' @param PPFDstart Start time (hours) for calculating total of photosynthetic photon flux density (PPFD). Calculated by summing across column light for the next number of hours specified from \code{nhrs} from the first PPFDstart time.
#' @param nhrs The number of hours to calculate the total of photosynthetic photon flux density (PPFD) from \code{PPFDstart}. For example, if \code{PPFDstart = 10} and \code{nhrs = 4}, \code{PPFDtotal} will be the sum of PPFD from 10 am to 1:59 pm.
#'
#' @details Model determines which model to estimate. 0 is the oxygen model (Eq. 2 from Nifong et al), 1 is the single station model without N consumption (DN base model; Eq. 3 from Nifong et al.), 2 is the single station model with N consumption (DN + Nconsume; Eq. 4 from Nifong et al.), the  3 being the two-station model without N consumption (DN base; Eq. 5 from Nifong et al.), and 4 being the two station model with N consumption (DN N consume; Eq. 6 from Nifong et al.).
#'
#' @export
#'
create_dataList <- function(data, model = 1, Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, lag_divisor = 0.0636573611, depth = 0.5588, PPFDstart = 14, nhrs = 24){

  # Convert data
  if(model > 0){
    data$arsat <- arsatfc(temp=data$temp, salinity=0, bp=data$bp)
    data$n2sat <- n2satfc(temp=data$temp, salinity=0, bp=data$bp)
    data$n2convert <- data$n2.ar * (28.014 / 1000) / (39.948 / 1000) * data$arsat
  }

  data$dtime<-chron::chron(dates=as.character(data$date), times=as.character(data$time))


  # Create list of data objects
  data_list <- list()

  ############################################
  # Oxygen models
  if(model == 0){

    # Subset data
    updata <- data[data$station == up,]
    downdata <- data[data$station == down,]
    nup  <- nrow(updata)
    ndown <- nrow(downdata)

    if(nrow(updata) != nrow(downdata)){
      stop("Number of up and down stations do not match")
    }

    # Data objects and size
    data_list$lag <- round(tt / lag_divisor) # 0.010947721
    up_sub <- 1:(nup - data_list$lag)
    down_sub <- (1 + data_list$lag):ndown

    data_list$nobs  <- length(up_sub)
    data_obj <- matrix(NA, ncol = 6, nrow = data_list$nobs)

    # Temperature
    data_obj[,1] <- updata$temp[up_sub]
    data_obj[,2] <- downdata$temp[down_sub]

    # O2 equilibrium concentration
    data_obj[,3] <- osat(data_obj[,1], updata$bp[up_sub] )
    data_obj[,4] <- osat(data_obj[,2], downdata$bp[down_sub])

    # O2
    data_obj[,5] <- updata$oxy[up_sub]
    data_obj[,6] <- downdata$oxy[down_sub]

    # Light
    data_list$PPFD <- downdata$light # Photosynthetic photon flux density (PPFD/light)


    # Calculate daily total of PPFD
    PPFDstart <- which(chron::hours(downdata$dtime) >= PPFDstart)[1] # Get the first observation to begin calculating daily total of PPFD
    end_ppfd <- which(downdata$dtime < (downdata$dtime[PPFDstart] + nhrs/24)) # Get all observations within 24 hours of start
    end_ppfd <- end_ppfd[length(end_ppfd)] # Whats the last observation

    ppfdtime <- downdata$dtime[end_ppfd] - downdata$dtime[PPFDstart] + (downdata$dtime[2] - downdata$dtime[1])
    print(paste0("Assuming equal sampling intervals, PPFDtotal is calculated over:" ))
    print(ppfdtime)

    data_list$PPFDtotal <- sum(downdata$light[PPFDstart:end_ppfd])

    # Parameters
    data_list$Kmean = Kmean
    data_list$Ksd = Ksd

    # Depth and time
    data_list$z = depth
    data_list$tt = tt
    data_list$time <- downdata$dtime[down_sub]

    # Data size
    data_list$data_obj <- data_obj
    data_list$ncol  <- ncol(data_obj)

    data_list$observed <- data_obj[,6]
  }


  ############################################
  # 1 Station models
  if(model %in% c(1,2)){

    # Data objects and size
    data_list$nobs  <- nrow(data)

    data_obj <- matrix(NA, ncol = 4, nrow = data_list$nobs)

    # Assign data to obj
    data_obj[,1] <- data$temp     # Temperature
    data_obj[,2] <- n2satfc(data_obj[,1], 0, data$bp )    # N2 equilibrium concentration
    data_obj[,3] <- data$n2convert     # N2
    data_obj[,4] <- c(-999, data$dtime[2:length(data$dtime)] - data$dtime[1:(length(data$dtime) - 1)])    # delta t

    # Parameters
    data_list$Kmean = Kmean
    data_list$Ksd = Ksd
    data_list$lag <- 0 # Not used

    # Depth and time
    data_list$z = depth
    data_list$tt = 0 # Not used
    data_list$time <- data$dtime

    # Light
    data_list$PPFD <- data$light # Photosynthetic photon flux density (PPFD/light)

    # Calculate daily total of PPFD
    PPFDstart <- which(chron::hours(data$dtime) >= PPFDstart)[1] # Get the first observation to begin calculating daily total of PPFD
    end_ppfd <- which(data$dtime < (data$dtime[PPFDstart] + nhrs/24)) # Get all observations within 24 hours of start
    end_ppfd <- end_ppfd[length(end_ppfd)] # Whats the last observation

    ppfdtime <- data$dtime[end_ppfd] - data$dtime[PPFDstart] + (data$dtime[2] - data$dtime[1])
    print(paste0("Assuming equal sampling intervals, PPFDtotal is calculated over:" ))
    print(ppfdtime)


    data_list$PPFDtotal <- sum(data$light[PPFDstart:end_ppfd])


    # Data size
    data_list$data_obj <- data_obj
    data_list$ncol  <- ncol(data_obj)
    data_list$observed <- data_obj[,3]

  }

  ############################################
  # 2 Station models
  if(model %in% c(3,4)){

    # Subset data
    updata <- data[data$station == up,]
    downdata <- data[data$station == down,]
    nup  <- nrow(updata)
    ndown <- nrow(downdata)

    if(nrow(updata) != nrow(downdata)){
      stop("Number of up and down stations do not match")
    }

    # Data objects and size
    data_list$lag <- round(tt / lag_divisor) # 3 for measured Ditch 2 HRT  ~4.58333/24= 0.1909720833/3 =0.063657
    up_sub <- 1:(nup - data_list$lag)
    down_sub <- (1 + data_list$lag):ndown

    data_list$nobs  <- length(up_sub)
    data_obj <- matrix(NA, ncol = 6, nrow = data_list$nobs)

    # Temperature
    data_obj[,1] <- updata$temp[up_sub]
    data_obj[,2] <- downdata$temp[down_sub]

    # N2 equilibrium concentration
    data_obj[,3] <- n2satfc(data_obj[,1], 0, updata$bp[up_sub] )
    data_obj[,4] <- n2satfc(data_obj[,2], 0, downdata$bp[down_sub])

    # N2
    data_obj[,5] <- updata$n2convert[up_sub]
    data_obj[,6] <- downdata$n2convert[down_sub]

    # Light
    data_list$PPFD <- downdata$light # Photosynthetic photon flux density (PPFD/light)


    # Calculate daily total of PPFD
    PPFDstart <- which(chron::hours(downdata$dtime) >= PPFDstart)[1] # Get the first observation to begin calculating daily total of PPFD
    end_ppfd <- which(downdata$dtime < (downdata$dtime[PPFDstart] + nhrs/24)) # Get all observations within 24 hours of start
    end_ppfd <- end_ppfd[length(end_ppfd)] # Whats the last observation

    ppfdtime <- downdata$dtime[end_ppfd] - downdata$dtime[PPFDstart] + (downdata$dtime[2] - downdata$dtime[1])
    print(paste0("Assuming equal sampling intervals, PPFDtotal is calculated over:" ))
    print(ppfdtime)


    data_list$PPFDtotal <- sum(downdata$light[PPFDstart:end_ppfd])

    # Parameters
    data_list$Kmean = Kmean
    data_list$Ksd = Ksd

    # Depth and time
    data_list$z = depth
    data_list$tt = tt
    data_list$time <- downdata$dtime[down_sub]

    # Data size
    data_list$data_obj <- data_obj
    data_list$ncol  <- ncol(data_obj)

    data_list$observed <- data_obj[,6]
  }

  return(data_list)
}



#' Estimate model
#'
#' @description fits a stan model and returns a stan model object
#'
#' @param dataList dataList created by \code{create_dataList}
#' @param model Model number for fitting algorithm (see details bleow).
#' @param nChains Number of chains to run
#' @param niter Number of iterations to run for each chain (including burnin)
#' @param burnin A positive integer specifying the number of warmup (aka burnin) iterations per chain.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console, which might be helpful for model debugging.
#'
#' @details Model determines which model to estimate. 0 is the oxygen model (Eq. 2 from Nifong et al), 1 is the single station model without N consumption (DN base model; Eq. 3 from Nifong et al.), 2 is the single station model with N consumption (DN + Nconsume; Eq. 4 from Nifong et al.), the  3 being the two-station model without N consumption (DN base; Eq. 5 from Nifong et al.), and 4 being the two station model with N consumption (DN N consume; Eq. 6 from Nifong et al.).
#'
#' @examples
#'
#' # Run model Eq. 5
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, depth = 0.5588)
#' mod <- fitmod(dataList, model = 3)
#' plotmod(mod, dataList = dataList, model = 3, file = NULL)
#'
#' # Run model Eq. 6
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.19097290833, depth = 0.5588)
#' mod2 <- fitmod(dataList, model = 4, verbose = FALSE)
#' plotmod(StanFit = mod2, dataList = dataList, model = 4)
#'
#'@export
fitmod <- function(dataList, model = 3, nChains = 2, niter = 5000, burnin = 1000, verbose = FALSE, incl_DN = TRUE){

  # Model set up
  dataList$mod = model

  # 02 Model
  if(model == 0 ){
    dataList$nparam = 3
  }

  # N2 Model
  if(model != 0){
    if(model %in% c(1,3)){
      dataList$nparam = 1
    }
    if(model %in% c(2,4)){
      dataList$nparam = 2
    }
    if(incl_DN){
      dataList$nparam <- dataList$nparam + 1
    }
  }

  # Get model
  stan_directory <- system.file("executables",package="Denitrification")
  old_wd <- getwd()
  setwd(stan_directory)

  # Estimate
  StanFit <- rstan::stan("nn2_model.stan", data = dataList, iter = niter, chains = nChains, verbose = verbose, warmup = burnin)
  setwd(old_wd)

  # Rename parameters

  return(StanFit)
}


#' Plot and diagnose
#'
#' @param StanFit Standmodel output by \code{fitmod}
#' @param dataList dataList created by \code{create_dataList}
#' @param model Model number for fitting algorithm (see details bleow).
#' @param file filname to save the parameter estimates. Will not save if NULL.
#'
#' @details Model determines which model to estimate. 0 is the oxygen model (Eq. 2 from Nifong et al), 1 is the single station model without N consumption (DN base model; Eq. 3 from Nifong et al.), 2 is the single station model with N consumption (DN + Nconsume; Eq. 4 from Nifong et al.), the  3 being the two-station model without N consumption (DN base; Eq. 5 from Nifong et al.), and 4 being the two station model with N consumption (DN N consume; Eq. 6 from Nifong et al.).
#'
#' @examples
#' # Run model Eq. 5 from Nifong et al.
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, depth = 0.5588)
#' mod <- fitmod(dataList, model = 3)
#' plotmod(mod, dataList = dataList, model = 3, file = NULL)
#'
#' # Run model Eq. 6 from Nifong et al.
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.19097290833, depth = 0.5588)
#' mod2 <- fitmod(dataList, model = 4, verbose = FALSE)
#' plotmod(StanFit = mod2, dataList = dataList, model = 4)
#'
#' @export
plotmod <- function(StanFit, dataList, model = 3, file = NULL, incl_DN = TRUE){

  ##########################################
  # Results of parameters to return
  results <- rstan::summary(StanFit)$summary
  rows_sub <- c(grep( "estVal", rownames(results)), grep( "predVal", rownames(results)))
  results <- results[-rows_sub,]

  if(model == 0){
    params_names <- c("K", "GPP", "ER", "sigma2", "logPost")
  }

  if(model %in% c(1,3)){
    params_names <- c("K600", "sigma2", "logPost")
  }
  if(model %in% c(2,4)){
    params_names <- c("K600", "Nfix", "sigma2", "logPost")
  }

  if(incl_DN & model != 0){
    params_names <- c(params_names[1], "DN", params_names[2:length(params_names)])
  }

  rownames(results) <- params_names


  if(!is.null(file)){
    write.csv(results, file = paste0(file, ".csv"))
  }

  # Get list of returned objects
  returned_objects <- extract(StanFit)
  estVal <- returned_objects$estVal # Predicted mean
  predVal <- returned_objects$predVal # Posterior predictive

  ##########################################
  # Summarize posterior and posterior predictive
  row_names <- c("mean", "median",
                 "2.5%PI", "97.5%PI",
                 "5%PI", "95%PI",
                 "min", "max", "n")

  posterior_summary <- matrix(nrow = length(row_names), ncol = dim(estVal)[2])
  posterior_summary[1, ] <- colMeans(estVal)
  posterior_summary[2:6, ] <- apply(estVal, 2, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
  posterior_summary[7, ] <- apply(estVal, 2, min)
  posterior_summary[8, ] <- apply(estVal, 2, max)
  posterior_summary <-as.data.frame(posterior_summary)
  names(posterior_summary) <- names(estVal)
  row.names(posterior_summary) <- row_names

  predictive_summary <- matrix(nrow = length(row_names), ncol = dim(predVal)[2])
  predictive_summary[1, ] <- colMeans(predVal)
  predictive_summary[2:6, ] <- apply(predVal, 2, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
  predictive_summary[7, ] <- apply(predVal, 2, min)
  predictive_summary[8, ] <- apply(predVal, 2, max)
  predictive_summary <-as.data.frame(predictive_summary)
  names(predictive_summary) <- names(predVal)
  row.names(predictive_summary) <- row_names

  ##########################################
  # Plots

  # Plot mean with CI
  for(i in 1:(1+!is.null(file))){

    if(i == 2){
      png(filename = paste0(file, "CI.png"), height = 5, width = 7, units = "in", res = 300)
    }

    # Main plot
    plot(dataList$time ,rep(NA, length(dataList$time)), type="l",xlab="Time", ylab=ifelse(model > 0, "Nitrogen  (mg/L)", "Oxygen  (mg/L)"),
         ylim=c(min(c(as.numeric(posterior_summary[7, ]), dataList$observed)), max(c(as.numeric(posterior_summary[8, ]), dataList$observed))), lwd=3, col="black" , xaxt = "n", las = 1)

    # Axis
    x_pos <- dataList$time[c(FALSE, TRUE, FALSE, FALSE)]
    x_axt <- chron::hours(dataList$time[c(FALSE, TRUE, FALSE, FALSE)])
    x_axt <- ifelse(x_axt == 0, 24, x_axt)
    xlab <- ifelse(x_axt < 12, paste0(x_axt, "am"), paste0(x_axt - 12, "pm"))
    xlab <- ifelse(xlab == "12pm", "12am", xlab)
    xlab <- ifelse(xlab == "0pm", "12pm", xlab)
    axis(1,at=x_pos, labels=xlab,las=1)

    # Credible interval
    polygon(
      x = c(dataList$time, rev(dataList$time)),
      y = c(posterior_summary[3, ], rev(posterior_summary[4, ])),
      col = "grey80", border = NA) # 95% CI
    polygon( x = c(dataList$time, rev(dataList$time)),
             y = c(posterior_summary[5, ], rev(posterior_summary[6, ])),
             col = "grey60", border = NA) # 90% CI

    # Median
    lines( x = dataList$time, y = posterior_summary[2, ], lty = 1, lwd = 3, col = 1) # Median

    # Observed
    points(dataList$time,dataList$observed)


    if(i == 2){
      dev.off()
    }
  }


  # Plot posterior predictive
  for(i in 1:(1+!is.null(file))){
    if(i == 2){
      png(filename = paste0(file, "posterior_predictive.png"), height = 5, width = 6, units = "in", res = 300)
    }

    # Main plot
    plot(dataList$time,rep(NA, length(dataList$time)), type="l",xlab="Time", ylab= ifelse(model > 0, "Nitrogen  (mg/L)", "Oxygen  (mg/L)"),
         ylim=c(min(c(as.numeric(predictive_summary[3, ]), dataList$observed)), max(c(as.numeric(predictive_summary[4, ]), dataList$observed))), lwd=3, col="black", main = "Posterior Predictive Check", xaxt = "n", las = 1)

    # Axis
    x_pos <- dataList$time[c(FALSE, TRUE, FALSE, FALSE)]
    x_axt <- chron::hours(dataList$time[c(FALSE, TRUE, FALSE, FALSE)])
    x_axt <- ifelse(x_axt == 0, 24, x_axt)
    xlab <- ifelse(x_axt < 12, paste0(x_axt, "am"), paste0(x_axt - 12, "pm"))
    xlab <- ifelse(xlab == "12pm", "12am", xlab)
    xlab <- ifelse(xlab == "0pm", "12pm", xlab)
    axis(1,at=x_pos, labels=xlab,las=1)

    # Credible interval
    polygon(
      x = c(dataList$time, rev(dataList$time)),
      y = c(predictive_summary[3, ], rev(predictive_summary[4, ])),
      col = "grey80", border = NA) # 95% CI
    polygon( x = c(dataList$time, rev(dataList$time)),
             y = c(predictive_summary[5, ], rev(predictive_summary[6, ])),
             col = "grey60", border = NA) # 90% CI

    # Median
    lines( x = dataList$time, y = predictive_summary[2, ], lty = 1, lwd = 3, col = 1)

    # Observed
    points(dataList$time, dataList$observed)

    if(i == 2){
      dev.off()
    }
  }


  return(results)
}

#' Extract observed and estimated values
#'
#' @description Extracts the observed and estimated N2 values from a stan model fit.
#'
#' @param StanFit Standmodel output by \code{fitmod}
#' @param dataList dataList created by \code{create_dataList}
#' @param model Model number for fitting algorithm (see details bleow).
#' @param file filname to save the estimates. Will not save if NULL.
#'
#' @return a dataframe include the Temperature, equilibrium concentration of N2, observed N2, dT and estimated mean, median, 95% and 90% confidence interval (CI) and prediction interval (PI).
#'
#' @details Model determines which model to estimate. 0 is the oxygen model (Eq. 2 from Nifong et al), 1 is the single station model without N consumption (DN base model; Eq. 3 from Nifong et al.), 2 is the single station model with N consumption (DN + Nconsume; Eq. 4 from Nifong et al.), the  3 being the two-station model without N consumption (DN base; Eq. 5 from Nifong et al.), and 4 being the two station model with N consumption (DN N consume; Eq. 6 from Nifong et al.).
#'
#' @examples
#' # Run model Eq. 5 from Nifong et al.
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, depth = 0.5588)
#' mod <- fitmod(dataList, model = 3)
#' plotmod(mod, dataList = dataList, model = 3, file = NULL)
#' extract_values(mod, dataList = dataList, model = 3, file = NULL)
#'
#' # Run model Eq. 6 from Nifong et al.
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.19097290833, depth = 0.5588)
#' mod2 <- fitmod(dataList, model = 4, verbose = FALSE)
#' plotmod(StanFit = mod2, dataList = dataList, model = 4)
#' extract_values(mod, dataList = dataList, model = 3, file = NULL)
#'
#' @export
extract_values <- function(StanFit = mod3, dataList = dataList, model = 3, file = NULL){

  # Get list of returned objects
  returned_objects <- extract(StanFit)
  estVal <- returned_objects$estVal # Predicted mean
  predVal <- returned_objects$predVal # Posterior predictive

  ##########################################
  # Summarize posterior and posterior predictive
  row_names <- c("mean", "median",
                 "2.5%CI", "97.5%CI",
                 "5%CI", "95%CI","2.5%PI", "97.5%PI",
                 "5%PI", "95%PI")

  if(model == 0){
    valName <- "Est02Down_"
  }
  if(model %in% c(1,2)){
    valName <- "EstN2_"
  }
  if(model %in% c(3,4)){
    valName <- "EstN2Down_"
  }

  row_names <- as.character(sapply(row_names, function(x) paste0(valName, x)))

  posterior_summary <- matrix(nrow = length(row_names), ncol = dim(estVal)[2])
  posterior_summary[1, ] <- colMeans(estVal)
  posterior_summary[2:6, ] <- apply(estVal, 2, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
  posterior_summary[7:10, ] <- apply(predVal, 2, quantile, probs= c(0.025, 0.975, 0.25, 0.75))
  posterior_summary <-as.data.frame(posterior_summary)
  row.names(posterior_summary) <-  row_names
  posterior_summary <- t(posterior_summary)

  # Combine results with data
  if(model == 0){
    data_names <- c("TempUp", "TempDown", "O2EquilUp", "O2EquilDown", "O2Up", "O2Down")
  }
  if(model %in% c(1,2)){
    data_names <- c("Temp", "N2Equil", "N2", "dT")
  }
  if(model > 2){
    data_names <- c("TempUp", "TempDown", "N2EquilUp", "N2EquilDown", "N2Up", "N2Down")
  }
  data <- data.frame(dataList$data_obj)
  colnames(data) <- data_names
  results <- cbind(data, posterior_summary)

  if(!is.null(file)){
    write.csv(results, file = paste0(file, "results.csv"))
  }

  return(results)
}
