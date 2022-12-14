
# Load packages ----
library(tidyverse)
library(lcmm)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")

run_models = FALSE

# Create folder for LCMM model outputs -----
dir.create(here::here("01_data", "lcmm_models"), recursive = TRUE,
           showWarnings = FALSE)

# Load data
resource = read_rds(
  here::here("01_data", paste0("resource_", timestamp, ".rds")))

# Create dataset for LCMM model ----
## Filter for emergency beddays prior to index admission
resource_pre = resource %>% 
  filter(period < 0) %>% 
  filter(resource_type == "emergency", measurement == "beddays") %>% 
  pivot_wider(names_from = c("measurement", "resource_type")) %>% 
  filter(n_days == period_days)

## Calculate total beddays during period, assign numerical ID
resource_pre = resource_pre %>% 
  group_by(PatientID) %>% 
  mutate(total_beddays = sum(beddays_emergency),
         PatientID_num = cur_group_id()) %>%
  ungroup()

## Dataset containing patients with at least one admission 
resource_pre_nonzero = resource_pre %>% 
  filter(total_beddays > 0) %>% 
  select(PatientID, PatientID_num, period, beddays_emergency)

## Save datasets ----
write_rds(resource_pre, 
          here::here("01_data", paste0("resource_pre_", timestamp, ".rds")))

write_rds(resource_pre_nonzero, 
          here::here("01_data", paste0("resource_pre_nonzero_", timestamp, ".rds")))

# Model parameters ----
## Settings ----
ng_max = 5                 # Maximum number of clusters
max_iter = 5000            # Max iterations
subject = "PatientID_num"
nproc = 2                  # Number of cores for parallel computation

## Model equations ----
link_function = "7-equi-splines"


if(run_models){
  # Run model on single cluster ----
  lcmm_model_1 = lcmm(fixed = beddays_emergency ~ 1 + period + I(period^2) + I(period^3),
                      random = ~ 1 + period,
                      ng = 1,
                      data = resource_pre_nonzero,
                      subject = subject,
                      maxiter = max_iter,
                      link = "7-equi-splines",
                      nproc = nproc)
  
  ## Save model ----
  write_rds(lcmm_model_1, here::here("01_data", "lcmm_models",
                                     paste0("lcmm_", timestamp, "_1.rds")))
  
  
  # Run models on 2 or more clusters ---- 
  c(2:ng_max) %>% 
    walk(function(ng){
      
      lcmm_model = lcmm(fixed = beddays_emergency ~ 1 + period + I(period^2) + I(period^3),
                        mixture = ~ 1 + period + I(period^2) + I(period^3),
                        random = ~ 1 + period,
                        ng = ng,
                        B = lcmm_model_1,
                        data = resource_pre_nonzero,
                        subject = subject,
                        maxiter = max_iter,
                        link = "7-equi-splines",
                        nproc = nproc)
      
      # Save model ----
      write_rds(lcmm_model, here::here("01_data", "lcmm_models",
                                       paste0("lcmm_", timestamp, "_", ng, ".rds")))
      
    })
}
