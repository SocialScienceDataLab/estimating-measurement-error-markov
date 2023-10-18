setwd ("/Users/Panko001/Downloads")

# install.packages("seqHMM")
library(seqHMM)
library(foreign)
library(tidyr)
library(dplyr)


data<-read.spss("data_hmm_final.sav", to.data.frame = T)
data<- data %>%
  mutate(fb_use_survey_cat5 = recode(fb_use_survey_cat5,"daily"="2.daily", "several times a day"="1.several times a day", "weekly"="3.weekly", "monthly"="4.monthly", 
                                     "less than monthly"= "5.less than monthly "),
         fb_use_tracking_cat5 = recode(fb_use_tracking_cat5,"daily"="2.daily", "several times a day"="1.several times a day", "weekly"="3.weekly", "monthly"="4.monthly", 
                                    "less than monthly"= "5.less than monthly "))

data_fb_survey<-data[,c(1,2,4)]

data_fb_tracking<-data[,c(1,2,5)]

data_fb_wide_survey<-spread(data_fb_survey, wave, fb_use_survey_cat5) %>%
  rename(fb_survey1="1", fb_survey2="2", fb_survey3="3")

head(data_fb_wide_survey)

data_fb_wide_tracking<-spread(data_fb_tracking, wave, fb_use_tracking_cat5) %>%
  rename(fb_tracking1="1", fb_tracking2="2", fb_tracking3="3")

head(data_fb_wide_tracking)


data_fb_wide<-inner_join(data_fb_wide_survey, data_fb_wide_tracking)
data_fb_seq_survey<-seqdef(data_fb_wide_survey[,2:4], start=1)
data_fb_seq_tracking<-seqdef(data_fb_wide_tracking[,2:4], start=1)

results<-NA
for (i in 2:7) {
  model_class <- build_hmm(observations = data_fb_seq_survey, n_states = i)
  model_fit<-fit_model(model_class,, control_em = list(restart = list(times = 5)))
  results<-append(results,BIC(model_fit$model))
}
plot(results[-1], type="b", xlab = "Number of classes", ylab = "BIC", main = "Survey")

model_4class<- build_hmm(observations = data_fb_seq_survey, n_states = 4)

model4_fit <- fit_model(model_4class, control_em = list(restart = list(times = 5)))
model4_fit

model4_fit$model$emission_probs
model4_fit$model$initial_probs

results_tracking<-NA
for (i in 2:7) {
  model_class_tracking <- build_hmm(observations = data_fb_seq_tracking, n_states = i)
  model_fit_tracking<-fit_model(model_class_tracking,control_em = list(restart = list(times = 5)))
  results_tracking<-append(results_tracking,BIC(model_fit_tracking$model))
}
plot(results_tracking[-1], type="b", xlab = "Number of classes", ylab = "BIC", main = "Tracking")

model_4class_tracking<- build_hmm(observations = data_fb_seq_tracking, n_states = 4)

model4_fit_tracking <- fit_model(model_4class_tracking,control_em = list(restart = list(times = 5)))
model4_fit_tracking
model4_fit_tracking$model$emission_probs
model4_fit_tracking$model$initial_probs

data_fb_seq_combo <- list(data_fb_seq_survey, data_fb_seq_tracking)

results_combo<-NA
for (i in 2:7) {
  model_class_combo <- build_hmm(observations = data_fb_seq_combo, n_states = i)
  model_fit_combo<-fit_model(model_class_combo, control_em = list(restart = list(times =5)))
  results_combo<-append(results_combo,BIC(model_fit_combo$model))
}

plot(results_combo[-1], type="b", xlab = "Number of classes", ylab = "BIC", main = "survey and tracking")

model_4class_combo<- build_hmm(observations = data_fb_seq_combo, n_states = 4, channel_names = c("survey", "tracking"))

model4_fit_combo <- fit_model(model_4class_combo, control_em = list(restart = list(times = 5)))
model4_fit_combo

model_3class_combo<- build_hmm(observations = data_fb_seq_combo, n_states = 3, channel_names = c("survey", "tracking"))

model3_fit_combo <- fit_model(model_3class_combo, control_em = list(restart = list(times = 5)))
model3_fit_combo

# include covariates
### Starting values for emission matrices
emiss_survey <- matrix(c(
  0.5023,	0.3438,	0.1305,	0.0206,	0.0028,
  0.44,	0.4649,	0.0901,	0.00425,	0.00075,
  0,	0.0158,	0.3739,	0.31077,	0.29953),
nrow = 3, ncol = 5, byrow = TRUE
)

emiss_survey2 <- matrix(c(
  0.5023,	0.3438,	0.1305,	0.0206,	0.0028,
  0.44,	0.4649,	0.0901,	0.00425,	0.00075,
  0,	0.0158,	0.3739,	0.31077,	0.29953),
  nrow = 3, ncol = 5, byrow = TRUE
)

emiss_tracking <- matrix(c(
  0.8200,  0.0000,    0.151,    0.0157,  0.0133,
  0.01012,  0.0903,    0.244,    0.1334,  0.52218,
  0.00624,  0.0529,    0.197,    0.1411,  0.60276),
  nrow = 3, ncol = 5, byrow = TRUE
)

emiss_tracking2 <- matrix(c(
  0.8200,  0.0000,    0.151,    0.0157,  0.0133,
  0.01012,  0.0903,    0.244,    0.1334,  0.52218,
  0.00624,  0.0529,    0.197,    0.1411,  0.60276),
  nrow = 3, ncol = 5, byrow = TRUE
)

emiss_list<-list(list(emiss_survey,emiss_survey2),list(emiss_tracking,emiss_tracking2))
### Starting values for transition matrix
trans <- matrix(
  c(
    1.00000, 0.00000, 0.00000,
      0.06393, 0.93498, 0.00109,
     0.00494, 0.00357, 0.99149
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

trans2 <- matrix(
  c(
    1.00000, 0.00000, 0.00000,
    0.06393, 0.93498, 0.00109,
    0.00494, 0.00357, 0.99149
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

trans_list<- list(trans,trans2)
### Starting values for initial state probabilities
inits <- c(0.297,0.405, 0.298)
inits2 <- c(0.297,0.405, 0.298)

inits_list<- list(inits, inits2)
data_covariates<- as.data.frame(data[complete.cases(data$device),])

data_fb_survey_cov<-data_covariates[,c(1,2,4)]

data_fb_tracking_cov<-data_covariates[,c(1,2,5)]

data_fb_wide_survey_cov<-spread(data_fb_survey_cov, wave, fb_use_survey_cat5) %>%
  rename(fb_survey1="1", fb_survey2="2", fb_survey3="3")

data_fb_wide_tracking_cov<-spread(data_fb_tracking_cov, wave, fb_use_tracking_cat5) %>%
  rename(fb_tracking1="1", fb_tracking2="2", fb_tracking3="3")

data_fb_wide_cov<-inner_join(data_fb_wide_survey_cov, data_fb_wide_tracking_cov)
data_fb_seq_survey_cov<-seqdef(data_fb_wide_survey_cov[,2:4], start=1)
data_fb_seq_tracking_cov<-seqdef(data_fb_wide_tracking_cov[,2:4], start=1)
data_fb_seq_combo_cov <-list(data_fb_seq_survey_cov,data_fb_seq_tracking_cov)

data_covariates_final<- data_covariates[!duplicated(data_covariates$new_id),]

model_3class_combo_cov <- build_mhmm(
  observations = data_fb_seq_combo_cov,
  transition_probs = trans_list,
  emission_probs = emiss_list,
  initial_probs = inits_list,
  formula = ~ device+1, data = data_covariates_final
)


model3_fit_combo_cov <- fit_model(model_3class_combo_cov, control_em = list(restart = list(times = 5)))
model_results<-model3_fit_combo_cov$model
summary(model_results)
