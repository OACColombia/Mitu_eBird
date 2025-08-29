library(tidyverse)
library(dclone);
library(mcmcplots)


# Time series

ebd_eff |>
  ungroup() |>
 filter(category %in% c("Ultrafine",
                         "Fine",
                         "Coarse")) |> # it was a captive individual
  filter(!any(observation_count == "NA"),
              scientific_name %in% c("Tinamus major", #max count 1
                                "Tinamus guttatus",
                                "Crypturellus casiquiare",
                                "Crypturellus cinereus",
                                "Crypturellus duidae", # max count 2
                                "Crypturellus soui",
                                "Crypturellus undulatus",
                                "Crypturellus variegatus",
                                "Nothocrax urumutum",
                                "Mitu tomentosum",
                                "Crax alector",
                                "Penelope jacquacu",
                                "Psophia crepitans")) |>
  group_by(scientific_name, observation_date, locality_id) |>
  count() |>
  group_by(scientific_name) |>
  count()

# perhaps only C. soui, C. variegatus, P. jacquacu, T. guttatus
mitu_huntSS <-
  ebd_eff |>
  ungroup() |>
  filter(category %in% c("Ultrafine",
                         "Fine",
                         "Coarse")) |> # it was a captive individual
  group_by(checklist_id) |>
  filter(!any(observation_count == "NA"),
         scientific_name %in% c("Tinamus guttatus",
                                "Crypturellus soui",
                                "Crypturellus variegatus",
                                "Penelope jacquacu")) |>
  group_by(scientific_name, month) |>
  mutate(Count = max(observation_count)) |>
  mutate(Date = floor_date(observation_date, "quarter")) |>
  dplyr::select(Date, scientific_name, Count) |>
  group_by(scientific_name, Date) |>
  summarise(Count = round(max(Count),0)) |>
  group_by(scientific_name) |>
  arrange(Date) |>
  as.data.frame()

mitu_huntSS |>
  ggplot(aes(x = Date, y = Count, color = scientific_name)) +
  facet_wrap(~factor(scientific_name,
                    levels = c("Tinamus guttatus",
                               "Crypturellus soui",
                               "Crypturellus variegatus",
                               "Penelope jacquacu")),
             scales = "free_y", ncol = 1) +
  geom_segment(aes(y = 0,
                   yend = Count),
               alpha = 0.5)+
  geom_point(alpha = 0.6) +
  labs(x = "Observation date",
       y = "Counts in a quarter of year")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "italic"),
        strip.background = element_rect(colour = "black", fill = NA),
        panel.border = element_rect(colour = "black", fill = NA))


# prepare the data for datacloning
prepare_dc_data_list <- function(df) {

  # Create full date sequence
  full_dates <- seq(min(df$Date), max(df$Date), by = "quarter") #adjust time steps

  # Expand grid of all combinations
  full_grid <- expand.grid(Date = full_dates,
                           scientific_name = unique(df$scientific_name),
                           stringsAsFactors = FALSE)

  # Merge and sort
  df_full <- full_grid |>
    left_join(df, by = c("Date", "scientific_name")) |>
    arrange(scientific_name, Date)

  # Extract month of the year for grouping
  df_full <- df_full |>
    mutate(Date = floor_date(Date, "month"))

  # Nest by Species
  nested <- df_full |>
    group_by(scientific_name) |>
    summarise(ts = list(Count), .groups = "drop") |>
    mutate(name = scientific_name)

  # Create named list of time series
  Y1_list <- set_names(nested$ts, nested$name)

  # Filter out time series that are all NA
  Y1_list <- Y1_list[!map_lgl(Y1_list, ~ all(is.na(.x)))]

  # Vector of time series lengths
  Tvec <- map_int(Y1_list, length)

  return(list(
    Y1 = Y1_list,
    Tvec = Tvec
  ))
}

dc_data <- prepare_dc_data_list(mitu_huntSS)
str(dc_data)

#GSS
StochGSS.dc <- function(){

  # Priors on model parameters. Priors are DC1 in Lele et al (2007)
  a1 ~ dnorm(0,1);   # constant, the population growth rate.
  c1 ~ dunif(-1,1);      # constant, the density dependence parameter.
  sig1 ~ dlnorm(-0.5,10); #variance parameter of stochastic  environment (process noise) in the system
  stovar1 <- 1/pow(sig1,2)
  tau1 ~ dunif(0,1); # detection probability (scaling factor that adjust expected counts to imperfect detection or measurement error)

  for(k in 1:K){
    # Simulate trajectory that depends on the previous
    mean_X1[1,k] <- a1/(1-c1) # Expected value of the first realization of the process
    # this is drawn from the stationary distribution of the process
    # Equation 14 (main text) and  A.4 in Appendix of Dennis et al 2006
    Varno1[k] <- pow(sig1,2)/(1-pow(c1,2)) #. Equation A.5 in Appendix of Dennis et al 2006

    # Updating the state: Stochastic process for all time steps
    X1[1,k]~dnorm(mean_X1[1,k], 1/Varno1[k]); #first estimation of population

    #iteration of the GSS model in the data
    for (t in 2:qp1) {
      mean_X1[t,k] <- a1 + c1 * X1[(t - 1),k]
      X1[t,k] ~ dnorm(mean_X1[t,k], stovar1) # Et is included here since a+cX(t-1) + Et ~ Normal(a+cX(t-1),sigma^2)
    }

    # Updating the observations, from the counts under Poisson observation error
    # incorporating detection probability directly into the Poisson mean
    for (t in 1:qp1) {
      Y1[t,k] ~ dpois(exp(X1[t,k])*tau1);
    }
  }
}

guess.calc <- function(Yobs,Tvec){

  T.t <-Tvec-Tvec[1]; #  For calculations, time starts at zero.
  q <- length(Yobs)-1;      #  Number of time series transitions, q.
  qp1 <- q+1;              #  q+1 gets used a lot, too.
  S.t <- T.t[2:qp1]-T.t[1:q];  #  Time intervals.
  Ybar <- mean(Yobs);
  Yvar <- sum((Yobs-Ybar)*(Yobs-Ybar))/q;
  mu1 <- Ybar;

  # Kludge an initial value for theta based on mean of Y(t+s) given Y(t).
  th1<- -mean(log(abs((Yobs[2:qp1]-mu1)/(Yobs[1:q]-mu1)))/S.t);
  bsq1<- 2*th1*Yvar/(1+2*th1);         # Moment estimate using stationary
  tsq1<- bsq1;                         #   variance, with betasq=tausq.

  #three 0's
  three0s <- sum(c(th1,bsq1,tsq1))
  if(three0s==0|is.na(three0s)){th1 <- 0.5;bsq1 <- 0.09; tsq1 <- 0.23;}


  out1 <- c(th1,bsq1,tsq1);
  if(sum(out1<1e-7)>=1){out1 <- c(0.5,0.09,0.23)}
  out <- c(mu1,out1);
  return(abs(out))

}

guess.calc2.0<- function(TimeAndNs){

  newmat <- TimeAndNs
  isnas <- sum(is.na(TimeAndNs))

  if(isnas >= 1){

    isnaind <- which(is.na(TimeAndNs[,2]), arr.ind=TRUE)
    newmat <- TimeAndNs[-isnaind,]
    newmat[,1] <- newmat[,1] - newmat[1,1]

  }

  init.guess <- guess.calc(Yobs = log(newmat[,2]), Tvec=newmat[,1])

  mu1  <- init.guess[1]
  th1  <- init.guess[2]
  bsq1 <- init.guess[3]
  sigsq1<- ((1-exp(-2*th1))*bsq1)/(2*th1)

  out <- c(mu=mu1, theta=th1, sigmasq = sigsq1)
  return(out)
}

Kalman.pred.fn <- function() {
  # Priors on model parameters: they are on the real line.
  parms ~ dmnorm(MuPost,PrecPost)
  a1 <- parms[1]
  c1 <- parms[2]
  sig1 <- parms[3]
  tau1 <- parms[4]
  stovar1 <- 1/pow(sig1,2)

  # Likelihood
  mean_X1[1] <- a1/(1-c1) # Expected value of the first realization of the process
  # this is drawn from the stationary distribution of the process
  # Equation 14 (main text) and  A.4 in Appendix of Dennis et al 2006
  Varno1 <- pow(sig1,2)/(1-pow(c1,2)) #. Equation A.5 in Appendix of Dennis et al 2006

  # Updating the state: Stochastic process for all time steps
  X1[1]~dnorm(mean_X1[1], 1/Varno1); #first estimation of population
  N[1] <- exp(X1[1])
  #iteration of the GSS model in the data
  for (t in 2:qp1) {
    mean_X1[t] <- a1 + c1 * X1[(t - 1)]
    X1[t] ~ dnorm(mean_X1[t], stovar1) # Et is included here since a+cX(t-1) + Et ~ Normal(a+cX(t-1),sigma^2)
    Y1[(t-1)] ~ dpois(exp(X1[t])*tau1)
    N[t] <- exp(X1[t])
  }
}

# generate range of dates for vector of time
dates_range <- seq(min(mitu_huntSS$Date), max(mitu_huntSS$Date),
                   by = "quarter")

# Tinamous guttatus fit ####

ts.4guess  <- dc_data$Y1[[4]]
tvec4guess  <- 1:length(ts.4guess)
onets4guess <- cbind(tvec4guess, ts.4guess)
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)

datalistGSS.dc <- list(K = 1,
                       qp1 = length(ts.4guess),
                       Y1 = dcdim(array(ts.4guess,
                                        dim = c(length(ts.4guess),1))))

dcrun.GSS.Tg <- dc.fit(data = datalistGSS.dc,
                    params = c("a1", "c1", "sig1", "tau1"), # previous attempt with tau1
                    model = StochGSS.dc,
                    n.clones = c(1,5,10),
                    multiply = "K",
                    unchanged = c("qp1"),
                    n.chains = 3,
                    n.adapt = 50000,
                    n.update = 100,
                    thin = 10,
                    n.iter = 100000)

saveRDS(dcrun.GSS.Tg, "PopDynamicFit_MaxCountQuarter_Tingut.RDS")

coef(dcrun.GSS.Tg)

# Predict
data4kalman1 <- list(qp1 = as.numeric(dc_data$Tvec[4]),
                    Y1 = array(dc_data$Y1[[4]],
                               dim = c(as.numeric(dc_data$Tvec[4]))),
                    MuPost = coef(dcrun.GSS.Tg),
                    PrecPost = solve(vcov(dcrun.GSS.Tg)))

BH_DC_Pred1 = jags.fit(data=data4kalman1,
                      params=c("N"),
                      model=Kalman.pred.fn)

# extract predictions and IQR around them
pred1 <- as.data.frame(t(mcmcapply(BH_DC_Pred1, quantile, c(0.25, 0.5, 0.75))))

popdyn1 <- as.data.frame(cbind(dates_range,pred1,dc_data$Y1[[4]]))
# modify names
names(popdyn1) <- c("Date", "Lower", "Estimated", "Upper", "Observed")

TingutCount <- popdyn1 |>
  pivot_longer(cols = c(Lower, Estimated, Upper, Observed),
               names_to = "Abundance",
               values_to = "Count") |>
  ggplot(aes(x = Date,
             y = Count,
             color = factor(Abundance,
                            levels = c("Observed",
                                       "Upper",
                                       "Estimated",
                                       "Lower")))) +
  geom_line(aes(linetype = factor(Abundance,
                                  levels = c("Observed",
                                             "Upper",
                                             "Estimated",
                                             "Lower")))) +
  geom_point(aes(shape = factor(Abundance,
                                levels = c("Observed",
                                           "Upper",
                                           "Estimated",
                                           "Lower")))) +
  labs(title = expression(italic("Tinamus guttatus")~"- Pamiwã"),
       x = "Time (quarter of year)",
       y = "Abundance",
       color = "",
       linetype = "",
       shape = "") +
  scale_linetype_manual(values = c(NA,"dashed","solid","dashed"))+
  scale_shape_manual(values = c(21, NA,NA,NA)) +
  scale_color_manual(values = c("blue","darkgray","black","darkgray")) +
#  scale_y_continuous(limits = c(0,7))+
  theme_classic() +
  theme(legend.position = "bottom")
TingutCount

# Crypturellus soui fit ####
str(dc_data) # 1
ts.4guess  <- dc_data$Y1[[1]]
tvec4guess  <- 1:length(ts.4guess)
onets4guess <- cbind(tvec4guess, ts.4guess)
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)

datalistGSS.dc <- list(K = 1,
                       qp1 = length(ts.4guess),
                       Y1 = dcdim(array(ts.4guess,
                                        dim = c(length(ts.4guess),1))))

dcrun.GSS.Cs <- dc.fit(data = datalistGSS.dc,
                       params = c("a1", "c1", "sig1", "tau1"), # previous attempt with tau1
                       model = StochGSS.dc,
                       n.clones = c(1,5,10),
                       multiply = "K",
                       unchanged = c("qp1"),
                       n.chains = 3,
                       n.adapt = 50000,
                       n.update = 100,
                       thin = 10,
                       n.iter = 100000)

saveRDS(dcrun.GSS.Cs, "PopDynamicFit_MaxCountQuarter_Crysou.RDS")

coef(dcrun.GSS.Cs)

# Predict
data4kalman2 <- list(qp1 = as.numeric(dc_data$Tvec[1]),
                     Y1 = array(dc_data$Y1[[1]],
                                dim = c(as.numeric(dc_data$Tvec[1]))),
                     MuPost = coef(dcrun.GSS.Cs),
                     PrecPost = solve(vcov(dcrun.GSS.Cs)))

BH_DC_Pred2 = jags.fit(data=data4kalman2,
                       params=c("N"),
                       model=Kalman.pred.fn)

# extract predictions and IQR around them
pred2 <- as.data.frame(t(mcmcapply(BH_DC_Pred2, quantile, c(0.25, 0.5, 0.75))))

popdyn2 <- as.data.frame(cbind(dates_range,pred2,dc_data$Y1[[1]]))
# modify names
names(popdyn2) <- c("Date", "Lower", "Estimated", "Upper", "Observed")

CrysouCount <- popdyn2 |>
  pivot_longer(cols = c(Lower, Estimated, Upper, Observed),
               names_to = "Abundance",
               values_to = "Count") |>
  ggplot(aes(x = Date,
             y = Count,
             color = factor(Abundance,
                            levels = c("Observed",
                                       "Upper",
                                       "Estimated",
                                       "Lower")))) +
  geom_line(aes(linetype = factor(Abundance,
                                  levels = c("Observed",
                                             "Upper",
                                             "Estimated",
                                             "Lower")))) +
  geom_point(aes(shape = factor(Abundance,
                                levels = c("Observed",
                                           "Upper",
                                           "Estimated",
                                           "Lower")))) +
  labs(title = expression(italic("Crypturellus soui")~"- Pamiwã"),
       x = "Time (quarter of year)",
       y = "Abundance",
       color = "",
       linetype = "",
       shape = "") +
  scale_linetype_manual(values = c(NA,"dashed","solid","dashed"))+
  scale_shape_manual(values = c(21, NA,NA,NA)) +
  scale_color_manual(values = c("blue","darkgray","black","darkgray")) +
  scale_y_continuous(limits = c(0,7))+
  theme_classic() +
  theme(legend.position = "bottom")
CrysouCount

# Crypturellus variegatus fit ####
str(dc_data) # 2
ts.4guess  <- dc_data$Y1[[2]]
tvec4guess  <- 1:length(ts.4guess)
onets4guess <- cbind(tvec4guess, ts.4guess)
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)

datalistGSS.dc <- list(K = 1,
                       qp1 = length(ts.4guess),
                       Y1 = dcdim(array(ts.4guess,
                                        dim = c(length(ts.4guess),1))))

dcrun.GSS.Cv <- dc.fit(data = datalistGSS.dc,
                       params = c("a1", "c1", "sig1", "tau1"), # previous attempt with tau1
                       model = StochGSS.dc,
                       n.clones = c(1,5,10),
                       multiply = "K",
                       unchanged = c("qp1"),
                       n.chains = 3,
                       n.adapt = 50000,
                       n.update = 100,
                       thin = 10,
                       n.iter = 100000)

saveRDS(dcrun.GSS.Cv, "PopDynamicFit_MaxCountQuarter_Cryvar.RDS")

coef(dcrun.GSS.Cv)

# Predict
data4kalman3 <- list(qp1 = as.numeric(dc_data$Tvec[2]),
                     Y1 = array(dc_data$Y1[[2]],
                                dim = c(as.numeric(dc_data$Tvec[2]))),
                     MuPost = coef(dcrun.GSS.Cv),
                     PrecPost = solve(vcov(dcrun.GSS.Cv)))

BH_DC_Pred3 = jags.fit(data=data4kalman3,
                       params=c("N"),
                       model=Kalman.pred.fn)

# extract predictions and IQR around them
pred3 <- as.data.frame(t(mcmcapply(BH_DC_Pred3, quantile, c(0.25, 0.5, 0.75))))

popdyn3 <- as.data.frame(cbind(dates_range,pred3,dc_data$Y1[[2]]))
# modify names
names(popdyn3) <- c("Date", "Lower", "Estimated", "Upper", "Observed")

CryvarCount <- popdyn3 |>
  pivot_longer(cols = c(Lower, Estimated, Upper, Observed),
               names_to = "Abundance",
               values_to = "Count") |>
  ggplot(aes(x = Date,
             y = Count,
             color = factor(Abundance,
                            levels = c("Observed",
                                       "Upper",
                                       "Estimated",
                                       "Lower")))) +
  geom_line(aes(linetype = factor(Abundance,
                                  levels = c("Observed",
                                             "Upper",
                                             "Estimated",
                                             "Lower")))) +
  geom_point(aes(shape = factor(Abundance,
                                levels = c("Observed",
                                           "Upper",
                                           "Estimated",
                                           "Lower")))) +
  labs(title = expression(italic("Crypturellus variegatus")~"- Pamiwã"),
       x = "Time (quarter of year)",
       y = "Abundance",
       color = "",
       linetype = "",
       shape = "") +
  scale_linetype_manual(values = c(NA,"dashed","solid","dashed"))+
  scale_shape_manual(values = c(21, NA,NA,NA)) +
  scale_color_manual(values = c("blue","darkgray","black","darkgray")) +
  scale_y_continuous(limits = c(0,7))+
  theme_classic() +
  theme(legend.position = "bottom")
CryvarCount

# Penelope jacquacu fit ####
str(dc_data) # 3
ts.4guess  <- dc_data$Y1[[3]]
tvec4guess  <- 1:length(ts.4guess)
onets4guess <- cbind(tvec4guess, ts.4guess)
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)

datalistGSS.dc <- list(K = 1,
                       qp1 = length(ts.4guess),
                       Y1 = dcdim(array(ts.4guess,
                                        dim = c(length(ts.4guess),1))))

dcrun.GSS.Pj <- dc.fit(data = datalistGSS.dc,
                       params = c("a1", "c1", "sig1", "tau1"), # previous attempt with tau1
                       model = StochGSS.dc,
                       n.clones = c(1,5,10),
                       multiply = "K",
                       unchanged = c("qp1"),
                       n.chains = 3,
                       n.adapt = 50000,
                       n.update = 100,
                       thin = 10,
                       n.iter = 100000)

saveRDS(dcrun.GSS.Pj, "PopDynamicFit_MaxCountQuarter_Penjac.RDS")

coef(dcrun.GSS.Pj)

# Predict
data4kalman4 <- list(qp1 = as.numeric(dc_data$Tvec[3]),
                     Y1 = array(dc_data$Y1[[3]],
                                dim = c(as.numeric(dc_data$Tvec[3]))),
                     MuPost = coef(dcrun.GSS.Pj),
                     PrecPost = solve(vcov(dcrun.GSS.Pj)))

BH_DC_Pred4 = jags.fit(data=data4kalman4,
                       params=c("N"),
                       model=Kalman.pred.fn)

# extract predictions and IQR around them
pred4 <- as.data.frame(t(mcmcapply(BH_DC_Pred4, quantile, c(0.25, 0.5, 0.75))))

popdyn4 <- as.data.frame(cbind(dates_range,pred4,dc_data$Y1[[3]]))
# modify names
names(popdyn4) <- c("Date", "Lower", "Estimated", "Upper", "Observed")

PenjacCount <- popdyn4 |>
  pivot_longer(cols = c(Lower, Estimated, Upper, Observed),
               names_to = "Abundance",
               values_to = "Count") |>
  ggplot(aes(x = Date,
             y = Count,
             color = factor(Abundance,
                            levels = c("Observed",
                                       "Upper",
                                       "Estimated",
                                       "Lower")))) +
  geom_line(aes(linetype = factor(Abundance,
                                  levels = c("Observed",
                                             "Upper",
                                             "Estimated",
                                             "Lower")))) +
  geom_point(aes(shape = factor(Abundance,
                                levels = c("Observed",
                                           "Upper",
                                           "Estimated",
                                           "Lower")))) +
  labs(title = expression(italic("Penelope jacquacu")~"- Pamiwã"),
       x = "Time (quarter of year)",
       y = "Abundance",
       color = "",
       linetype = "",
       shape = "") +
  scale_linetype_manual(values = c(NA,"dashed","solid","dashed"))+
  scale_shape_manual(values = c(21, NA,NA,NA)) +
  scale_color_manual(values = c("blue","darkgray","black","darkgray")) +
  scale_y_continuous(limits = c(0,8))+
  theme_classic() +
  theme(legend.position = "bottom")
PenjacCount

## Probability of persistence ####

randmvn <- function(n, mu.vec, cov.mat){

  # Save the length of the mean vector of the multivariate normal distribution to sample
  p         <- length(mu.vec);
  # The Cholesky decomposition
  #(factorization of a real symmetric positive-definite sqr matriz)
  Tau       <- chol(cov.mat, pivot=TRUE);
  # generate normal deviates outside loop
  Zmat      <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);

  # empty matrix
  out       <- matrix(0,nrow=p,ncol=n);
  # iterate
  for(i in 1:n){
    Z       <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
  }

  return(out)
}


Game.sim <- function(prediction.t = 10,
                   theta.list=
                     list(
                       p = 1,
                       muvec = log(round(popdyn1$Estimated,0)[prediction.t+10]),
                       sigmasq = coef(dcrun.GSS.Tg)[3],
                       c=coef(dcrun.GSS.Tg)[2],
                       a=coef(dcrun.GSS.Tg)[1]),
                   xo = log(round(popdyn1$Estimated,0)[prediction.t]),
                   nsteps=14*4, #14 years is two generations \times 12 months
                   rnames = "Tingut"){

  p <- theta.list$p
  muvec <- theta.list$muvec
  sigmasq <- theta.list$sigmasq
  c <- theta.list$c

  a  <- theta.list$a
  sigma <- (1/sigmasq)

  X <- matrix(0,nrow=p,ncol=nsteps)
  X[1:p,1] <- randmvn(n=1,mu.vec=a+c%*%xo, cov.mat=sigma)

  for(i in 2:nsteps){
    ## Process error
    X[1:p,i] <- randmvn(n=1,mu.vec=a + c%*%X[1:p,(i-1)], cov.mat=sigma)
  }

  X <- round(exp(cbind(xo,X)),0)
  row.names(X) <- rnames
  colnames(X) <- 0:nsteps
  return(X)
}

Game.sim()

sims.list.Tg <- list()
steps = length(popdyn1$Estimated)
lsim.Tg = 14*4
nsims <- 1000

for(i in 1:steps){
  sims.list.Tg[[i]] <- matrix(0,nrow=nsims,ncol=(lsim.Tg)+1)
}

for(i in 1:nsims){
  for(j in 10:steps){ #from the position 10 in the time series
    ith.sim <- Game.sim(prediction.t = j)
    sims.list.Tg[[j]][i, ] <- as.numeric(ith.sim)
  }
}

saveRDS(sims.list.Tg, "Simulations_Tingut_10k.RDS")

coef(dcrun.GSS.Tg)

n.crit.Tg <- 4 #50% of the mean value

pers.prop.Tg <- matrix(NA, nrow=steps, ncol=1)

for(i in 10:steps){
pers.prop.Tg[i,]<- (1 - sum(sims.list.Tg[[i]][,2:lsim.Tg]<n.crit.Tg) / (lsim.Tg*nsims))
}

popdyn1.phi <- cbind(popdyn1,pers.prop.Tg)

ggplot(popdyn1.phi,
       aes(x = Date, y = pers.prop.Tm))+
  geom_line()+
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "Probability of persistance",
       title = expression(italic("Tinamus guttatus")~"- Pamiwã"))+
  theme_classic()

# Penelope jacquacu

Game.sim()

sims.list.Pj <- list()
steps = length(popdyn4$Estimated)
lsim.Pj = 16*4
nsims <- 1000

for(i in 1:steps){
  sims.list.Pj[[i]] <- matrix(0,nrow=nsims,ncol=(lsim.Pj)+1)
}

for(i in 1:nsims){
  for(j in 10:steps){ #from the position 10 in the time series
    ith.sim <- Game.sim(prediction.t = j,
                        theta.list=
                          list(
                            p = 1,
                            muvec = log(round(popdyn4$Estimated,0)[prediction.t+10]),
                            sigmasq = coef(dcrun.GSS.Pj)[3],
                            c=coef(dcrun.GSS.Pj)[2],
                            a=coef(dcrun.GSS.Pj)[1]),
                        xo = log(round(popdyn4$Estimated,0)[prediction.t]),
                        nsteps=16*4, #16 years is two generations \times 4 times per year
                        rnames = "Penjac")
    sims.list.Pj[[j]][i, ] <- as.numeric(ith.sim)
  }
}

saveRDS(sims.list.Pj, "Simulations_Penjac_10k.RDS")

coef(dcrun.GSS.Pj)

n.crit.Pj <- 4 #50% of the mean value

pers.prop.Pj <- matrix(NA, nrow=steps, ncol=1)

for(i in 10:steps){
  pers.prop.Pj[i,]<- (1 - sum(sims.list.Pj[[i]][,2:lsim.Pj]<n.crit.Pj) / (lsim.Pj*nsims))
}

popdyn4.phi <- cbind(popdyn4,pers.prop.Pj)

ggplot(popdyn4.phi,
       aes(x = Date, y = pers.prop.Tm))+
  geom_line()+
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "Probability of persistance",
       title = expression(italic("Penelope jacquacu")~"- Pamiwã"))+
  theme_classic()
