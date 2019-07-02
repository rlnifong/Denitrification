functions {

  // Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  real KcorN2(real temp, real K600) {
    return K600/(600/(1615-(temp*92.15)+(2.349*temp^2)-(0.0240*temp^3)))^-0.5;
  }


  // Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  real KcorO2(real temp, real K600) {
    return K600/(600/(1568-(temp*86.04)+(2.142*temp^2)-(0.0216*temp^3)))^-0.5;
  }

}
data {
  // Model specifications
  int<lower = 0> mod;
  int<lower = 1> nobs;
  int<lower = 1> ncol;
  int<lower = 1> nparam;
  int incl_DN;

  // Data
  real data_obj[nobs, ncol];

  // For models 0 data_obj is the following
  // -- Column 1 = tempup (celcius)
  // -- Column 2 = tempdown (celcius)
  // -- Column 3 = o2equilup (g m^-3)
  // -- Column 4 = o2equildown (g m^-3)
  // -- Column 5 = o2up (g m^-3)
  // -- Column 6 = O2down (g m^-3)

  // For models 1 and 2 data_obj is the following
  // -- Column 1 = temp (celcius)
  // -- Column 2 = n2equil (g m^-3)
  // -- Column 3 = n2 (g m^-3)
  // -- Column 4 = delta t (days^-1)

  // For models 3 and 4 data_obj is the following
  // -- Column 1 = tempup (celcius)
  // -- Column 2 = tempdown (celcius)
  // -- Column 3 = n2equilup (g m^-3)
  // -- Column 4 = n2equildown (g m^-3)
  // -- Column 5 = n2up (g m^-3)
  // -- Column 6 = n2down (g m^-3)

  // Fixed values
  real z; // Water depth (m)
  real tt; // Travel time between stations
  int<lower = 0> lag; // Time lag (if mod < 3 it is 0)

  real PPFD[nobs + lag]; // Photosynthetic photon flux density (PPFD/light) (mol m-2 s-1)
  real PPFDtotal; // Daily total of PPFD (mol m-2 s-1 d-1)

  // Prior specifications
  real Kmean;
  real Ksd;
}
parameters {
  // -- For 02 models
  // params[1] = K = Gas exchange rates (d-1)
  // params[2] = GPP = Gross primary production  (g m-2 d-1)
  // params[3] = ER = Ecosystem respiration  (g m-2 d-1)

  // -- For N2 models
  // params[1] = K600 = Gas exchange rates (d-1)
  // params[2] = DN = Denitrification rate (g m-2 d-1)
  // params[3] = Nfix = N2 fixation rate (g m-2 d-1)
  real params[nparam];

  real<lower=0> sigma; // Variance
}
model {
    // --------------------------------------
  // Objects
  vector[nobs] estVal;

  // --------------------------------------
  // Priors
  sigma ~ normal(0,10);

  // -- O2 Model
  if(mod == 0){
    params[1] ~ normal(5, 10); // K
    params[2] ~ normal(-5, 10); // GPP
    params[3] ~ normal(5, 5); // ER
  }

  // -- N2 Models
  if(mod > 0){
    params[1] ~ normal(Kmean, Ksd);

    // Remove DN
    if(incl_DN == 1){
      params[2] ~ normal(1, 5);
    }

    // Additional prior for model 2 and 4
    if(mod == 2){
      params[2+incl_DN] ~ normal(-1, 5);
    }
    if(mod == 4){
      params[2+incl_DN] ~ normal(-1, 5);
    }
  }

  // --------------------------------------
  // Models
  // -- Two-station oxygen models (Models 0)
  if(mod == 0){
  // Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5] + params[3] * tt / z ;
      estVal[i] += params[2] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      estVal[i] +=  KcorO2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;
      estVal[i] /=  (1 + KcorO2(data_obj[i, 1], params[1]) * tt / 2);
    }

    // Likelihood
    data_obj[,6] ~ normal( estVal, sigma);
  }


  // -- One-station models (Models 1 and 2)
  if(mod != 0){
  if(mod < 3){
    // Initialize
    estVal[1] = data_obj[1, 3];

    // Model loop
    for(i in 2:nobs){
      estVal[i] = estVal[i - 1];

      // Remove DN
      if(incl_DN == 1){
        estVal[i] += params[2] * data_obj[i, 4] / z ;
      }

      estVal[i] +=  KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] * (data_obj[i, 2] + data_obj[i - 1, 2] - estVal[i - 1])/2;

      // Add N consumption (DN + Nconsume) for model 2
      if(mod == 2){
         estVal[i] += params[2+incl_DN] / z * PPFD[i] / PPFDtotal;
      }

      estVal[i] /=  (1 + KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] / 2);
    }

    // Likelihood
    data_obj[2:,3] ~ normal( estVal[2:], sigma);
  }
  }


  // -- Two-station models (Models 3 and 4)
  if(mod > 2){
    // Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5] ;

      // Remove DN
      if(incl_DN == 1){
        estVal[i] += params[2] * tt / z;
      }

      estVal[i] +=  KcorN2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;

      // Add N consumption (DN + Nconsume) for model 4
      if(mod == 4){
         estVal[i] += params[2+incl_DN] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      }

      estVal[i] /=  (1 + KcorN2(data_obj[i, 1], params[1]) * tt / 2);
    }

    // Likelihood
    data_obj[,6] ~ normal( estVal, sigma);
  }

}
generated quantities{
    // --------------------------------------
  // Objects
  vector[nobs] estVal;
  vector[nobs] predVal;

  // --------------------------------------
  // Models

  // Two-station O2 model
  if(mod == 0){
  // Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5] + params[3] * tt / z ;
      estVal[i] += params[2] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      estVal[i] +=  KcorO2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;
      estVal[i] /=  (1 + KcorO2(data_obj[i, 1], params[1]) * tt / 2);
    }

    // Posterior predictive
    for(i in 1:nobs){
      predVal[i] = normal_rng( estVal[i], sigma);
    }
  }

  // -- One-station models (Models 1 and 2)
  if(mod != 0){
  if(mod < 3){
   // Initialize
    estVal[1] = data_obj[1, 3];

    // Model loop
    for(i in 2:nobs){
      estVal[i] = estVal[i - 1];

      // Remove DN
      if(incl_DN == 1){
        estVal[i] += params[2] * data_obj[i, 4] / z ;
      }

      estVal[i] +=  KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] * (data_obj[i, 2] + data_obj[i - 1, 2] - estVal[i - 1])/2;

      // Add N consumption (DN + Nconsume) for model 2
      if(mod == 2){
         estVal[i] += params[2+incl_DN] / z * PPFD[i] / PPFDtotal;
      }

      estVal[i] /=  (1 + KcorN2(data_obj[i, 1], params[1]) * data_obj[i, 4] / 2);
    }

    // Posterior predictive
    predVal[1] =  estVal[1];
    for(i in 2:nobs){
      predVal[i] = normal_rng( estVal[i], sigma);
    }
  }
  }


  // -- Two-station models (Models 3 and 4)
  if(mod > 2){
    // Model loop
    for(i in 1:nobs){
      estVal[i] = data_obj[i, 5] ;

      // Remove DN
      if(incl_DN == 1){
        estVal[i] += params[2] * tt / z;
      }

      estVal[i] +=  KcorN2(data_obj[i, 1], params[1]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;

      // Add N consumption (DN + Nconsume) for model 4
      if(mod == 4){
         estVal[i] += params[2+incl_DN] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      }

      estVal[i] /=  (1 + KcorN2(data_obj[i, 1], params[1]) * tt / 2);
    }

    // Posterior predictive
    for(i in 1:nobs){
      predVal[i] = normal_rng( estVal[i], sigma);
    }
  }
}
