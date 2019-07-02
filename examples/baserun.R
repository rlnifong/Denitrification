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

# NOTE: Change the working directory for the data

########################################
# Oxygen model
oxy_data <- read.csv("D7ASUMetDataJune2018Day2.csv") #The file name needs to be updated
dataList0 <- create_dataList(oxy_data, model = 0, up="up1", down="down1", tt=0.18611125, lag_divisor = 0.010947721, PPFDstart = 14, depth = 0.53594)

# Run model Eq. 2
mod0 <- fitmod(dataList0, model = 0)
rstan::traceplot(mod0, pars = c("params", "sigma")) # Look at trace of DN, K600, sigma2
params0 <- plotmod(StanFit = mod0, dataList = dataList0, model = 0, file = "model0") # Plot it and get parameter values
params0

values0 <- extract_values(StanFit = mod0, dataList = dataList0, model = 0, file = "model0") # Get observed and predicted values
head(values0)

# Model validation
gof(sim=values0$Est02Down_mean,obs=values0$O2Down)


########################################
# One station N2 models
onestation <- read.csv("Ditch2onestat_DDmetrop.csv")
dataList1 <- create_dataList(onestation, model = 1, Kmean = 4.03, Ksd = 4.0, PPFDstart = 14, depth = 0.5588)

# Run model Eq. 3
mod1 <- fitmod(dataList1, model = 1)
rstan::traceplot(mod1, pars = c("params", "sigma")) # Look at trace of K600, DN, sigma2
params1 <- plotmod(StanFit = mod1, dataList = dataList1, model = 1, file = "model1") # Plot it and get parameter values
params1

values1 <- extract_values(StanFit = mod1, dataList = dataList1, model = 1, file = "model1") # Get observed and predicted values
head(values1)

# Model validation
gof(sim=values1$EstN2_mean,obs=values1$N2)



# Run model Eq. 4
mod2 <- fitmod(dataList1, model = 2)
rstan::traceplot(mod2, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params2 <- plotmod(StanFit = mod2, dataList = dataList1, model = 2, file = "model2") # Plot it and get parameter values
params2

values2 <- extract_values(StanFit = mod2, dataList = dataList1, model = 2, file = "model2") # Get observed and predicted values
head(values2)

# Model validation
gof(sim=values2$EstN2_mean,obs=values2$N2)



########################################
# Two station N2 models
twostation <- read.csv("Ditch2Twostat_DDmetrop.csv")
dataList2 <- create_dataList(twostation, model = 3, Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt=0.1909720833, lag_divisor = 0.063657361, PPFDstart = 14, depth = 0.5588)

# Run model Eq. 5
mod3 <- fitmod(dataList2, model = 3)
rstan::traceplot(mod3, pars = c("params", "sigma")) # Look at trace of K600, DN, sigma2
params3 <- plotmod(StanFit = mod3, dataList = dataList2, model = 3, file = "model3") # Plot it and get parameter values
params3

values3 <- extract_values(StanFit = mod3, dataList = dataList2, model = 3, file = "model3") # Get observed and predicted values
head(values3)

# Model validation
gof(sim=values3$EstN2Down_mean,obs=values3$N2Down)



# Run model Eq. 6
mod4 <- fitmod(dataList2, model = 4)
rstan::traceplot(mod4, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params4 <- plotmod(StanFit = mod4, dataList = dataList2, model = 4, file = "model4") # Plot it and get parameter values
params4

values4 <- extract_values(StanFit = mod4, dataList = dataList2, model = 4, file = "model4") # Get observed and predicted values
head(values4)

# Model validation
gof(sim=values4$EstN2Down_mean,obs=values4$N2Down)

