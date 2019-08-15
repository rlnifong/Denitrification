D2onestationN <- read.csv("data/ASUS2DDnight2.csv")
D2onestationN
D2onestationD<-read.csv("data/ASUS2DDDay.csv")
D2onestationD

dataListD2N <- create_dataList(D2onestationN, model = 1, Kmean = 0.279583, Ksd = 1.0, PPFDstart = 23, nhrs =13, depth = 0.5588)
dataListD2N$PPFDtotal
dataListD2D <- create_dataList(D2onestationD, model = 1, Kmean = 0.279583, Ksd = 1.0, PPFDstart = 11, nhrs=13, depth = 0.5588)
dataListD2D$PPFDtotal


mod2D <- fitmod(dataListD2D, model = 2, control = list(adapt_delta = 0.9) )
rstan::traceplot(mod2D, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params2D <- plotmod(StanFit = mod2D, dataList = dataListD2D, model = 2, file = "model2") # Plot it and get parameter values
model_likelihood(dataList = dataListD2D, paramsVec = params2D[,1], model = 2, incl_DN = TRUE)


# FOR NIGHT HERE
mod2N <- fitmod(dataListD2N, model = 2)
rstan::traceplot(mod2N, pars = c("params", "sigma")) # Look at trace of K600, DN, Nfix, sigma2
params2N <- plotmod(StanFit = mod2N, dataList = dataListD2N, model = 2, file = "model2") # Plot it and get parameter values
params2N
