# Setup
library(rstan)
library(Denitrification)

########################################
# One station N2 models
onestation <- read.csv("data/Ditch2onestat_DDmetrop.csv")

# 12 to 1 model
onestation12to1 <- onestation[14:15,]
dataList12to1 <- create_dataList(onestation12to1, model = 2, Kmean = 4.03, Ksd = 4.0, PPFDstart = 0, nhrs = 2, depth = 0.5588)

# Run model Eq. 4 for 12 to 1
mod12to1 <- fitmod(dataList12to1, model = 2)
rstan::traceplot(mod12to1, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params12to1 <- plotmod(StanFit = mod12to1, dataList = dataList12to1, model = 2, file = "model2") # Plot it and get parameter values
params12to1

values12to1 <- extract_values(StanFit = mod12to1, dataList = dataList12to1, model = 2, file = "model2") # Get observed and predicted values
head(values12to1)


# 12 to 2 model
onestation12to2 <- onestation[14:16,]
dataList12to2 <- create_dataList(onestation12to2, model = 2, Kmean = 4.03, Ksd = 4.0, PPFDstart = 0, nhrs = 3, depth = 0.5588)


# Check likelihood
paramsVec = c(1,1,1,1) #  K600, DN, Nfix, sigma2
model_likelihood(dataList = dataList12to2, paramsVec = paramsVec, model = 2, incl_DN = TRUE)

# Run model Eq. 4 for 12 to 2
mod12to2 <- fitmod(dataList12to2, model = 2)
rstan::traceplot(mod12to2, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params12to2 <- plotmod(StanFit = mod12to2, dataList = dataList12to2, model = 2, file = "model2") # Plot it and get parameter values
params12to2

values12to2 <- extract_values(StanFit = mod12to2, dataList = dataList12to2, model = 2, file = "model2") # Get observed and predicted values
head(values12to2)
