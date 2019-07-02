########################################
# Install necessary packages
install.packages("devtools") # Install devtool
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE) # Install stan
install.packages("hydroGOF")

# Install code
library(devtools)
devtools::install_github("grantdadams/Denitrification", auth_token = "4925b42ac46f1e0aefd671e9dc0c1cf1b3157017")

# Setup
library(rstan)
library(Denitrification)
library(hydroGOF)


########################################
# One station models
########################################
onestation <- read.csv("Ditch2onestat_DDmetrop.csv")
dataList1 <- create_dataList(onestation, model = 1, Kmean = 4.03, Ksd = 4.0, depth = 0.5588)


########################################
# Run model Eq. 3
mod1 <- fitmod(dataList1, model = 1)
rstan:: traceplot(mod1, pars = c("params", "sigma")) # Look at trace of K600, DN, sigma2
params1 <- plotmod(StanFit = mod1, dataList = dataList1, model = 1, file = "model1") # Plot it and get parameter values
params1

values1 <- extract_values(StanFit = mod1, dataList = dataList1, model = 1, file = "model1") # Get observed and predicted values
head(values1)

# Model validation
gof(sim=values1$EstN2_mean,obs=values1$N2)


########################################
# Run model Eq. 3 without DN
mod1a <- fitmod(dataList1, model = 1, incl_DN = FALSE)
rstan:: traceplot(mod1a, pars = c("params", "sigma")) # Look at trace of K600, sigma2
params1a <- plotmod(StanFit = mod1a, dataList = dataList1, model = 1, file = "model1a", incl_DN = FALSE) # Plot it and get parameter values
params1a

values1a <- extract_values(StanFit = mod1a, dataList = dataList1, model = 1, file = "model1a") # Get observed and predicted values
head(values1a)

# Model validation
gof(sim=values1a$EstN2_mean,obs=values1a$N2)


########################################
# Run model Eq. 4
mod2 <- fitmod(dataList1, model = 2)
rstan:: traceplot(mod2, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params2 <- plotmod(StanFit = mod2, dataList = dataList1, model = 2, file = "model2") # Plot it and get parameter values
params2

values2 <- extract_values(StanFit = mod2, dataList = dataList1, model = 2, file = "model2") # Get observed and predicted values
head(values2)

# Model validation
gof(sim=values2$EstN2_mean,obs=values2$N2)


########################################
# Run model Eq. 4 without DN
mod2a <- fitmod(dataList1, model = 2, incl_DN = FALSE)
rstan:: traceplot(mod2a, pars = c("params", "sigma")) # Look at trace of K600, Nfix, sigma2
params2a <- plotmod(StanFit = mod2a, dataList = dataList1, model = 2, file = "model2a", incl_DN = FALSE) # Plot it and get parameter values
params2a

values2a <- extract_values(StanFit = mod2a, dataList = dataList1, model = 2, file = "model2a") # Get observed and predicted values
head(values2a)

# Model validation
gof(sim=values2a$EstN2_mean,obs=values2a$N2)


########################################
# Plot it
tiff("onestation.tiff", width = 3.30, height = 6, units = 'in', res = 800)
par(mfrow=c(4,1), oma=c(5.5,4,1,1), mar=c(1,2,0.1,1))

# Create list of model objects
values_list <- list(values1a, values1, values2a, values2)
legend_text <- c("(a)", "(b)", "(c)", "(d)")

# Run through model objects
for(i in 1:length(values_list)){

  # -- Plot each run
  # Estimate
  plot(seq(1:length(values_list[[i]]$N2)),values_list[[i]]$EstN2_median,type="l",lty=1,lwd=2,col="black",ylim=c(10,14),xaxt="n",cex.axis=1.3,xlab="",ylab="")

  # Observed
  points(seq(1:length(values_list[[i]]$N2)),values_list[[i]]$N2, cex=1.25)

  # Saturation
  lines(seq(1:length(values_list[[i]]$N2)),values_list[[i]]$N2Equil, lty=2, lwd=2, col="gray40")

  # Legend
  legend("topleft", legend_text[i], bty = "n", cex = 1.25, adj = c(1,-.15))

  # Add side legend
  if(i == 2){
    mtext(expression(paste("N"[2] *"-N ", ( "mg"*" L"^-1), sep="")),side =2, line=1.5, outer=TRUE, cex=1)
  }
}

# Add bottom axis
axis(1,at=c(3,7,11,15,19,23), labels=c("4pm","8pm","12am","4am","8am","12pm"),las=2,cex.axis=1.3)
mtext(expression(paste("Time")),side =1, line=4, outer=TRUE, cex=1)


dev.off()



########################################
# Two station models
########################################
twostation <- read.csv("Ditch2Twostat_DDmetrop.csv")
dataList2 <- create_dataList(twostation, model = 3, Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, depth = 0.5588)

########################################
# Run model Eq. 5
mod3 <- fitmod(dataList2, model = 3)
rstan:: traceplot(mod3, pars = c("params", "sigma")) # Look at trace of K600, DN, sigma2
params3 <- plotmod(StanFit = mod3, dataList = dataList2, model = 3, file = "model3") # Plot it and get parameter values
params3

values3 <- extract_values(StanFit = mod3, dataList = dataList2, model = 3, file = "model3") # Get observed and predicted values

# Model validation
gof(sim=values3$EstN2Down_mean,obs=values3$N2Down)


########################################
# Run model Eq. 5 without DN
mod3a <- fitmod(dataList2, model = 3, incl_DN = FALSE)
rstan:: traceplot(mod3a, pars = c("params", "sigma")) # Look at trace of K600, sigma2
params3a <- plotmod(StanFit = mod3a, dataList = dataList2, model = 3, file = "model3a", incl_DN = FALSE) # Plot it and get parameter values
params3a

values3a <- extract_values(StanFit = mod3a, dataList = dataList2, model = 3, file = "model3a") # Get observed and predicted values
head(values3a)

# Model validation
gof(sim=values3a$EstN2Down_mean,obs=values3a$N2Down)


########################################
# Run model Eq. 6
mod4 <- fitmod(dataList2, model = 4)
rstan:: traceplot(mod4, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params4 <- plotmod(StanFit = mod4, dataList = dataList2, model = 4, file = "model4") # Plot it and get parameter values
params4

values4 <- extract_values(StanFit = mod4, dataList = dataList2, model = 4, file = "model4") # Get observed and predicted values
head(values4)

# Model validation
gof(sim=values4$EstN2Down_mean,obs=values4$N2Down)


########################################
# Run model Eq. 6 without DN
mod4a <- fitmod(dataList2, model = 4, incl_DN = FALSE)
rstan:: traceplot(mod4a, pars = c("params", "sigma")) # Look at trace of K600, sigma2
params4a <- plotmod(StanFit = mod4a, dataList = dataList2, model = 4, file = "model4a", incl_DN = FALSE) # Plot it and get parameter values
params4a

values4a <- extract_values(StanFit = mod4a, dataList = dataList2, model = 4, file = "model4a") # Get observed and predicted values
head(values4a)

# Model validation
gof(sim=values4a$EstN2Down_mean,obs=values4a$N2Down)


########################################
# Plot it
tiff("twostation.tiff", width = 3.30, height = 6, units = 'in', res = 800)
par(mfrow=c(4,1), oma=c(5.5,4,1,1), mar=c(1,2,0.1,1))

# Create list of model objects
values_list <- list(values3a, values3, values4a, values4)
legend_text <- c("(a)", "(b)", "(c)", "(d)")

# Run through model objects
for(i in 1:length(values_list)){

  # -- Plot each run
  # Estimate
  plot(seq(1:length(values_list[[i]]$N2Down)),values_list[[i]]$EstN2Down_median,type="l",lty=1,lwd=2,col="black",ylim=c(10,14),xaxt="n",cex.axis=1.3,xlab="",ylab="")

  # Observed
  points(seq(1:length(values_list[[i]]$N2Down)),values_list[[i]]$N2Down, cex=1.25)

  # Saturation
  lines(seq(1:length(values_list[[i]]$N2Down)),values_list[[i]]$N2EquilDown, lty=2, lwd=2, col="gray40")

  # Legend
  legend("topleft", legend_text[i], bty = "n", cex = 1.25, adj = c(1,-.15))

  # Add side legend
  if(i == 2){
    mtext(expression(paste("N"[2] *"-N ", ( "mg"*" L"^-1), sep="")),side =2, line=1.5, outer=TRUE, cex=1)
  }
}

# Add bottom axis
axis(1,at=c(3,7,11,15,19,23), labels=c("4pm","8pm","12am","4am","8am","12pm"),las=2,cex.axis=1.3)
mtext(expression(paste("Time")),side =1, line=4, outer=TRUE, cex=1)


dev.off()

