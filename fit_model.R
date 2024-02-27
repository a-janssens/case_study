
# Fit spatial model on LRTI data ------------------------------------------

library(INLA)
library(tidyverse)
library(sf)

# Working directory and others 
cd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cd)
getwd()

# dir_data_raw <- "../../github-repo-data2/data/downloaded_from_server/"
# dir_data_clean <- "../../github-repo-data2/data/clean/"
# dir_output <- "../output_model_fit_case/"

# Source useful scripts 
source("functions_inla.R")
source("functions_geo.R")

# Load the cleaned data 
data_inla <- read.csv("data/data_lrti_2019.csv") %>% dplyr::select(-X)
map_fl <- st_read("data/", layer  = "map_fl")


# Prepare for modelling ----

## Data aggregated according to communities 
data_inla_model_aggr <- data_inla %>% 
  group_by(nis_code, ID, sex, age_category, increased_compensation) %>% 
  summarize(ncases = sum(ncases),
            npatients = sum(npatients)) %>% 
  mutate_at(c("sex","age_category","increased_compensation"), factor) %>% 
  mutate(age_category = relevel(age_category,"(16,50]")) %>% 
  mutate(ID2 = ID)


## Data inclduing practice_id variable
data_inla_model <- data_inla %>% 
  mutate_at(c("sex","age_category","increased_compensation"), factor) %>% 
  mutate(age_category = relevel(age_category,"(16,50]")) %>% 
  mutate(ID2 = ID)


# Spatial models ----------------------------------------------------------

## Define prior of precision hyperparameters ----

# PC prior 
prec.prior = list(prec = list(prior = "pc.prec", param = c(1,0.01)))

## Model with BESAG + PRACTICE ----
# - Model we selected to simulate from 
formula_besag_practice <- ncases ~ sex + age_category + increased_compensation +
  f(ID, model = "besag", graph = paste0("data/","map.gra"), scale.model = T,
    hyper = prec.prior) +
  f(practice_id, model="iid",
    hyper = prec.prior)

s <- Sys.time()
model_besag_practice <- inla(formula_besag_practice,
                             data = data_inla_model,
                             family="binomial",
                             Ntrials = npatients,
                             control.compute = list(dic = T, waic = T, cpo = T, config = T),
                             control.predictor = list(link = 1))
Sys.time() - s # 25''

# Summary model output 
summary(model_besag_practice)
model_besag_practice$mode$mode.status


# Plot Random Effects -----------------------------------------------------

re_besag_practice <- get_random_effects_inla(model_besag_practice, OR= T)
# pdf(file = paste0(dir_output,"02_random_effect_besag_practice.pdf"))
plot_random_effects(re_besag_practice, plot_map = map_fl)
# dev.off()

# ### Save data of figures gy
# write.csv(re_besag_practice[[1]], "../results_data/S3_data_figure3a.csv")
# write.csv(re_besag_practice[[2]], "../results_data/S3_data_figure3b.csv")


# BESAG + PRACTICE save estimates after model selection script ---------------------------

# Plot posterior of intercept 
plot(model_besag_practice$marginals.fixed[[1]], type = "l", 
     xlab = expression(alpha), ylab = "density")

# Get effect, save and show
fe_inla <- get_fixed_effects_inla(model_besag_practice,OR=T)
# write.csv(fe_inla, file = paste0(dir_output,"02_fixed_OR.csv"))
# pdf(file = paste0(dir_output,"02_fixed_OR.pdf"))
plot_fixed_effects(fe_inla, OR = T)+
  theme_bw()
# dev.off()

re_inla <- get_random_effects_inla(model_besag_practice,OR=T)
# for(df in names(re_inla)){
  # write.csv(re_inla[[df]],file = paste0(dir_output,"02_",df,"_OR.csv"))
# }
# write.csv(fe_inla, file = paste0(dir_output,"/02_fixed_OR.csv"))
# pdf(file = paste0(dir_output,"02_random_OR.pdf"))
plot_random_effects(re_inla,plot_map = map_fl)
# dev.off()

# Logit level estimates
fe_inla <- get_fixed_effects_inla(model_besag_practice,OR=F)
# write.csv(fe_inla, file = paste0(dir_output,"02_fixed.csv"))
# pdf(file = paste0(dir_output,"02_fixed.pdf"))
plot_fixed_effects(fe_inla, OR = F)
# dev.off()

re_inla <- get_random_effects_inla(model_besag_practice, OR = F)
# for(df in names(re_inla)){
#   write.csv(re_inla[[df]],file = paste0(dir_output,"02_",df,".csv"))
# }
# pdf(file = paste0(dir_output,"02_random.pdf"))
# pdf(file = "~/Desktop/random_logit.pdf")
plot_random_effects(re_inla,plot_map = map_fl)
# dev.off()

# Precisions of random effects 
precisions <- summary(model_besag_practice)$hyperpar
# write.csv(precisions, file = paste0(dir_output,"02_precisions.csv"))

# Model fit statistics 
model_fit_statistics <- get_model_fit_statistics_inla(model_besag_practice)
# write.csv(model_fit_statistics, file = paste0(dir_output,"02_model_fit_statistics.csv"))

# Fitted values 
fitted_values <- model_besag_practice$summary.fitted.values
# write.csv(fitted_values, file = paste0(dir_output,"02_fitted_values.csv"))
