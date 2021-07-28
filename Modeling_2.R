library(mgcv)
library(tidyverse)
library(lubridate)
library(PresenceAbsence)
library(mgcViz)

vetula_freq_model_2 <- gam(data = na.omit(bah_data_2), 
                         (Abundance != 0) ~ Exp + Current + s(Btemp, bs = "tp") + 
                          s(nMonth, bs = "cc") + s(nDate, bs = "tp") + s(Btime, bs = "tp") +
                          s(Start, bs = "cc") + AverageDepth + Habitat +
                           te(lat, lon, bs = "tp", by = vetula_breeding) 
                         ,family = "binomial")
AIC(vetula_freq_model_2)

#AIC optimized: 8122.528
#AIC - Visibility: NA
#AIC - Exp: 8171.915
#AIC - Current: 8142.777
#AIC - Btemp: 8136.134
#AIC - nMonth: 8151.698
#AIC - nDate: 8248.008
#AIC - BTime: 8211.306
#AIC - Start: 8226.489
#AIC - moon_phase: NA
#AIC - AverageDepth: 8145.251
#AIC - Habitat: 8131.789
#AIC - lat/long: 8636.79
#AIC - by = vetula: 8137.504


#some summary information for the model
summary(vetula_freq_model_2)
plot(getViz(vetula_freq_model_2))
gam.check(vetula_freq_model_2)


#produces lists of fit metrics for each fold of cross validation
r_square_bin <- list()
rmse_bin <- list()
avg_error_bin <- list()
kappa_list <- list()
AUC_list <- list()
specificity_list <- list()
sensitivity_list <- list()
PCC_list <- list()
for(i in unique(bah_data_2$k_fold)){
  kth_gam <- gam(data = filter(bah_data_2, k_fold != i), 
                 (Abundance != 0) ~ Exp + Current + s(Btemp, bs = "tp") + 
                   s(nMonth, bs = "cc") + s(nDate, bs = "tp") + s(Btime, bs = "tp") +
                   s(Start, bs = "cc") + AverageDepth + Habitat +
                   te(lat, lon, bs = "tp", by = vetula_breeding) 
                 ,family = "binomial")
  kth_predictions <- predict.gam(kth_gam, newdata = filter(bah_data_2, k_fold == i), type = "response")
  kth_actual <- as.numeric(filter(bah_data_2, k_fold == i)$Abundance[!is.na(kth_predictions)]!=0)
  kth_predictions <- kth_predictions[!is.na(kth_predictions)]
  r_square_bin <- append(r_square_bin, r_square(as.numeric(kth_predictions),kth_actual))
  rmse_bin <- append(rmse_bin, RMSE(as.numeric(kth_predictions), kth_actual))
  avg_error_bin <- append(avg_error_bin, avg_error(as.numeric(kth_predictions), kth_actual))
  presence_absence_pred <- tibble(id = seq_along(kth_predictions),
                                  obs = kth_actual,
                                  pred = as.numeric(kth_predictions))
  pressence_absence_metrics <- presence.absence.accuracy(presence_absence_pred)
  kappa_list <- append(kappa_list, pressence_absence_metrics$Kappa)
  AUC_list <- append(AUC_list, pressence_absence_metrics$AUC)
  specificity_list <- append(specificity_list, pressence_absence_metrics$specificity)
  sensitivity_list <- append(sensitivity_list, pressence_absence_metrics$sensitivity)
  PCC_list <- append(PCC_list, pressence_absence_metrics$PCC)
}

#the model has high specificity but a low sensitivity
#ie the model is pretty good at saying where the fish isnt, but is not good at knowing where the fish will be 
#---------------------------------------- negative binomial model -------------#


vetula_nb_model_2 <- gam(data = na.omit(filter(bah_data_2, Abundance != 4)), 
                       SFMA_to_num(Abundance) ~Visibility + Exp + Current + s(Btemp, bs = "tp") + 
                         s(nDate, bs = "tp") + s(Btime, bs = "tp") +
                         s(Start, bs = "cc")  + Habitat +
                         te(lat, lon, bs = "tp")
                       ,family = "nb")
AIC(vetula_nb_model_2)
summary(vetula_nb_model_2)
gam.check(vetula_nb_model_2)
plot(getViz(vetula_nb_model_2))

#AIC optimized: 16031.22
#AIC - Visibility: 16046.54
#AIC - Exp: 16046.6
#AIC - Current: 16035.41
#AIC - Btemp: 16090.05
#AIC - nMonth: NA
#AIC - nDate:16219.91
#AIC - BTime: 16094.14
#AIC - Start: 16081.36
#AIC - moon_phase: NA
#AIC - AverageDepth: NA
#AIC - Habitat: 16057.14
#AIC - lat/long: 16442.87
#AIC - by = vetula: 16044.04

#k fold cross validation
r_square_nbin <- list()
rmse_nbin <- list()
avg_error_nbin <- list()
cats_percent_nbin <- list()

for(i in unique(bah_data_2$k_fold)){
  kth_gam <- gam(data = filter(bah_data_2, k_fold != i, Abundance != 4), 
                 SFMA_to_num(Abundance) ~Visibility + Exp + Current + s(Btemp, bs = "tp") + 
                   s(nDate, bs = "tp") + s(Btime, bs = "tp") +
                   s(Start, bs = "cc")  + Habitat +
                   te(lat, lon, bs = "tp")
                 ,family = "nb")
  kth_predictions <- predict.gam(kth_gam, newdata = filter(bah_data_2, k_fold == i,Abundance != 4), type = "response")
  kth_actual <-SFMA_to_num(filter(bah_data_2, k_fold == i, Abundance != 4)$Abundance[!is.na(kth_predictions)])
  kth_predictions <- as.numeric(kth_predictions[!is.na(kth_predictions)])
  r_square_nbin <- append(r_square_nbin, r_square(as.numeric(kth_predictions),kth_actual))
  rmse_nbin <- append(rmse_nbin, RMSE(as.numeric(kth_predictions), kth_actual))
  avg_error_nbin <- append(avg_error_nbin, avg_error(as.numeric(kth_predictions), kth_actual))
  
  kth_predictions <- num_to_SFMA(predict.gam(kth_gam, newdata = filter(bah_data_2, k_fold == i,Abundance != 4), type = "response"))
  kth_actual <- filter(bah_data_2, k_fold == i)$Abundance
  kth_actual <- kth_actual[!is.na(kth_predictions)]
  kth_predictions <- kth_predictions[!is.na(kth_predictions)]
  cats_percent_nbin <- append(cats_percent_nbin, mean(kth_actual == kth_predictions, na.rm = TRUE))
}

#------------------------------ multinomial model and comparison model -----------------------------------#

vetula_multinom_model <- gam( data = filter(bah_data_2, Abundance != 4 ),
                              list(Abundance~s(nDate, bs = "tp") + te(lon,lat, bs = "tp"),
                                   ~s(nDate, bs = "tp") + te(lon, lat, bs = "tp"),
                                   ~s(nDate, bs = "tp") + te(lon, lat, bs = "tp")),
                              family = multinom(K = 3))
summary(vetula_multinom_model)
plot(getViz(vetula_multinom_model))

#cross validation
percent_guess_multinom <- list()

for(i in unique(bah_data_2$k_fold)){
  kth_gam <- gam(data = filter(bah_data_2, k_fold != i, Abundance != 4), 
                 list(Abundance~s(nDate, bs = "tp") + te(lon,lat, bs = "tp"),
                      ~s(nDate, bs = "tp") + te(lon, lat, bs = "tp"),
                      ~s(nDate, bs = "tp") + te(lon, lat, bs = "tp")),
                 family = multinom(K = 3))
  kth_predictions <- predict.gam(kth_gam, newdata = filter(bah_data_2, k_fold == i,Abundance != 4), type = "response")
  kth_predictions <- apply(kth_predictions, 1, which.max)
  kth_predictions <- kth_predictions[kth_predictions != "integer(0)"]
  kth_actual <- filter(bah_data_2, !is.na(nDate), !is.na(lon), k_fold == i)$Abundance
  kth_actual <- kth_actual[kth_actual != 4]
  percent_guess_multinom <- append(percent_guess_multinom, mean(kth_actual == (unlist(kth_predictions) -1), na.rm = TRUE))
}


#nb categorical comparison model
percent_guess_nbin <- list()

for(i in unique(bah_data_2$k_fold)){
  kth_gam <- gam(data = filter(bah_data_2, k_fold != i, Abundance != 4), 
                 SFMA_to_num(Abundance) ~ s(nDate, bs = "tp") + te(lon, lat, bs = "tp"),
                 family = "nb")
  kth_predictions <- num_to_SFMA(predict.gam(kth_gam, newdata = filter(bah_data_2, k_fold == i,Abundance != 4), type = "response"))
  kth_predictions <- kth_predictions[!is.na(kth_predictions)]
  kth_actual <- filter(bah_data_2, !is.na(nDate), !is.na(lon), k_fold == i)$Abundance
  kth_actual <- kth_actual[kth_actual != 4]
  percent_guess_nbin <- append(percent_guess_nbin, mean(kth_actual == kth_predictions, na.rm = TRUE))
}

#------------------------------ population Trends with time -------------------------#

vetula_abundance_decline <- gam(data = filter(bah_data_2, Abundance != 4), 
                         SFMA_to_num(Abundance) ~Visibility + Exp + Current + s(Btemp, bs = "tp") + 
                           nDate + s(Btime, bs = "tp") +
                           s(Start, bs = "cc")  + Habitat +
                           te(lat, lon, bs = "tp") 
                         ,family = "nb")
summary(vetula_abundance_decline)
#highly significant decrease in abundance over time


vetula_year_factor <- gam(data = filter(bah_data_2, Abundance != 4), 
                   SFMA_to_num(Abundance) ~Visibility + Exp + Current + s(Btemp, bs = "tp") + 
                     Year + s(Btime, bs = "tp") +
                     s(Start, bs = "cc")  + Habitat +
                     te(lat, lon, bs = "tp") 
                   ,family = "nb")


year_abund_trend <- bah_data_2 %>%
  group_by(year(Date))%>%
  summarise(estimate = mean(SFMA_to_num(Abundance)), CI = CI_from_SFMA(Abundance, include_observational = TRUE))%>%
  na.omit()


year_standardized_sim <- tibble(Exp = 'E',
                                nMonth = .3,
                                moon_phase = 2,
                                Btemp = 80,
                                Btime = 50,
                                Start = 12,
                                Visibility = 4,
                                Current = 2,
                                Habitat = 1,
                                Year = c(1994:2021),
                                lat = 24.88667,
                                lon = -77.885,
                                vetula_breeding = FALSE,
                                AverageDepth = 5)

model_abundance_trend <- predict.gam(vetula_year_factor,newdata = year_standardized_sim,
                                          se.fit = TRUE)
abundance_trends <- tibble(year = 1994:2021,
                      abundance = exp(as.vector(model_abundance_trend$fit)),
                      min = exp(as.vector(model_abundance_trend$fit) - 1.96 *as.vector(model_abundance_trend$se.fit)),
                      max = exp(as.vector(model_abundance_trend$fit) + 1.96 *as.vector(model_abundance_trend$se.fit)))


#used for plotting both standardized and actual data, needs rework
#year_trends <- tibble(Abundance = c(year_abund_trend$estimate,as.vector(corrected_year_abund_trend$fit)),
#                      year = c(year_abund_trend$`year(Date)`, 1994:2021),
#                      CImin = c(year_abund_trend$estimate/(year_abund_trend$CI+1), as.vector(corrected_year_abund_trend$fit)-1.96*as.vector(corrected_year_abund_trend$se.fit)),
#                      CImax = c(year_abund_trend$estimate*(year_abund_trend$CI+1), as.vector(corrected_year_abund_trend$fit)+1.96*as.vector(corrected_year_abund_trend$se.fit)),
#                      source = c(rep("unstandardized", times = length(year_abund_trend$`year(Date)`)), rep("standardized", times = length(corrected_year_abund_trend$fit))))

ggplot(data = filter(abundance_trends, year != 2021), mapping = aes(y = abundance, x = year)) +
  geom_point() +
  geom_errorbar(aes(ymax = min, ymin = max)) +
  labs(title = "Standardized Abundance by Year", subtitle = "with 95% confidence intervals")

#frequency plot
vetula_year_factor_freq <- gam(data = filter(bah_data_2, Abundance != 4), 
                               (Abundance != 0) ~ Exp + Current + s(Btemp, bs = "tp") + 
                                 s(nMonth, bs = "cc") + Year + s(Btime, bs = "tp") +
                                 s(Start, bs = "cc") + s(moon_phase, bs = "cc") + AverageDepth + Habitat +
                                 te(lat, lon, bs = "tp", by = vetula_breeding) 
                               ,family = "binomial")

freq_standardized <- predict.gam(vetula_year_factor_freq,newdata = year_standardized_sim,
                                 se.fit = TRUE)

freq_trends <- tibble(year = 1994:2021,
                      frequency = exp(as.vector(freq_standardized$fit))/(1 +exp(freq_standardized$fit)),
                      max = exp(as.vector(freq_standardized$fit) + 1.96*as.vector(freq_standardized$se.fit))/(1 + exp(freq_standardized$fit+1.96*as.vector(freq_standardized$se.fit))),
                      min = exp(as.vector(freq_standardized$fit) - 1.96*as.vector(freq_standardized$se.fit))/(1 + exp(freq_standardized$fit-1.96*as.vector(freq_standardized$se.fit))))

ggplot(data = filter(freq_trends, year != 2021), mapping = aes(x = year, y = frequency)) + 
  geom_point() +
  geom_errorbar(aes(ymin = min, ymax = max)) + 
  labs( title = "Standardized Observation Frequency by Year", subtitle = "with 95% confidence intervals")


#------------------------- Simulated Data ----------------------------#


#Balistes vetula mean ~2
df1<-data.frame(x<- runif(10000, min = -3, max = 3))
#df1$x[df1$x <= 0] <- 0
df1$y<-rnbinom(10000,mu=exp(df1$x),size=10)
#size and mu for variance relations
df1$c<- num_to_SFMA(df1$y)
table(df1$c)

sim_gam1 <- gam(data = df1, y ~ x, family = "nb")
sim_gam2 <- gam(data = df1, SFMA_to_num(c) ~ x, family = "nb")
summary(sim_gam1)
summary(sim_gam2)

#RMSE

test_predictions_sim1 <- predict.gam(object = sim_gam1, type = "response")
rmse1 = (sum((test_predictions_sim1-df1$y)^2)/length(df1$y))^1/2
avg_error1 = mean(test_predictions_sim1-df1$y)

test_predictions_sim2 <- predict.gam(object = sim_gam2, type = "response")
rmse2 = (sum((test_predictions_sim2-df1$y)^2)/length(df1$y))^1/2
avg_error2 = mean(test_predictions_sim2-df1$y)


#Coefficients, significance test
#We can interpret the negative binomial regression coefficient as follows: for a one unit change in the predictor variable, the difference in the logs of expected counts of the response variable is expected to change by the respective regression coefficient, given the other predictor variables in the model are held constant.

mean(SFMA_to_num(df1$c))
mean(predict(sim_gam2, type = "response"))
mean(df1$y)

plot(sim_gam2)
ggplot(data = df1, mapping = aes(x = x, y = y)) + geom_smooth()

#------------------------------------- Residual plots(DHarma) ---------------#

library(DHARMa)

#Function to simulate DHARMa residuals from a negative binomial GAM
simulateNegBinGam <- function(modfit, nsims=250, offsetval=1){
  muval = predict(modfit, type = "response")*offsetval  #Get the mean with offset
  nObs = length(muval)  
  thetaval = modfit$family$getTheta(trans=TRUE)  #Get theta not log transformed
  sim = replicate(nsims,rnbinom(nObs,mu=muval, size=thetaval))  #Simulate negative binomial data
  sim
}

simvals<-simulateNegBinGam(vetula_nb_model_2 ,nsim=1000)  #Make simulated values from neg bin
simulation_output_nb = createDHARMa(simulatedResponse = simvals, 
                                   observedResponse = SFMA_to_num(na.omit(filter(bah_data_2, Abundance != 4))$Abundance) , 
                                   fittedPredictedResponse = predict(vetula_nb_model_2,type="response"))
plot(simulation_output_nb)
testDispersion(vetula_freq_model_2)

#function to simulate DHAMa residuals from binomial GAM
simulateBinGam <- function(modfit, nsims = 250, offsetval = 1){
  
}
