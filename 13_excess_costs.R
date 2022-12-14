
# Load packages ----
library(tidyverse)
library(finalfit)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "09_excess_cost"), recursive = TRUE,
           showWarnings = FALSE)

# Bootstrap samples
n_bootstrap = 500

# Load datasets ----
survivor_index = read_rds(
  here::here("01_data", paste0("survivor_index_", timestamp, ".rds")))

cluster_labeled = read_rds(
  here::here("01_data", paste0("cluster_labeled_", timestamp, ".rds"))) 

resource_emerg_2yr_prior = read_rds(
  here::here("01_data", paste0("resource_emerg_2yr_prior_", timestamp, ".rds"))) 

resource = read_rds(
  here::here("01_data", paste0("resource_", timestamp, ".rds")))

# Join to cohort cluster and prior emergency data ----
survivor_index = survivor_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID")  %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))

# Determine annualised baseline (year -2 to -0.5) and followup (first 6 months) costs -----
resource_prepost = resource %>%
  filter(measurement == "cost") %>%
  group_by(PatientID, period) %>% 
  summarise(
    n_days = first(n_days),
    cost = sum(value)
  ) %>% 
  ungroup() %>% 
  mutate(
    interval = case_when(
      (period >= -24) & (period < -6) ~ "baseline",
      (period > 0) & (period <= 6) ~ "followup",
      TRUE ~ NA_character_
    )
  ) %>% 
  drop_na(interval) %>% 
  group_by(PatientID, interval) %>% 
  summarise(
    days = sum(n_days),
    cost = sum(cost)
  ) %>% 
  pivot_wider(names_from = interval, values_from = c(days, cost)) %>% 
  mutate(
    # Calucluate annualised costs
    annualised_baseline = cost_baseline/days_baseline*365.25,
    annualised_followup = cost_followup/days_followup*365.25,
    annualised_excess = annualised_followup - annualised_baseline
  )


# Filter for patients with at least 6 months potential follow-up -------
resource_prepost = resource_prepost %>%
  left_join(
    survivor_index %>%
      mutate(
        potential_followup_6months = if_else((discharge_date + days(180)) <= date_censor,
                                             "Yes", NA_character_)
      ) %>% 
      select(PatientID, potential_followup_6months),
    by = "PatientID"
  ) %>% 
  filter(potential_followup_6months == "Yes")
  


## Stratifying variables ----
var_statifaction = c(
  "overall", "age.factor", "sex",  
  "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
  "class_factor", "prior_emergency_beddays_factor"
)

# Extract labels
vlabel = extract_variable_label(survivor_index)


# Table of stratified baseline followup and excess costs 
tbl_excess_cost = var_statifaction %>% 
  map(function(group_var){
    
    resource_prepost %>% 
      left_join(survivor_index %>% 
                  select(PatientID, grouping = all_of(group_var)), by = "PatientID") %>% 
      group_by(grouping) %>% 
      summarise(
        statistic = c("mean", "ci_lower", "ci_upper"),
        baseline = Hmisc::smean.cl.boot(annualised_baseline, B = n_bootstrap),
        followup = Hmisc::smean.cl.boot(annualised_followup, B = n_bootstrap),
        excess   = Hmisc::smean.cl.boot(annualised_excess, B = n_bootstrap)
      ) %>% 
      pivot_longer(cols = c(baseline, followup, excess), names_to = "cost_type") %>% 
      pivot_wider(names_from = statistic) %>% 
      group_by(grouping, cost_type) %>% 
      summarise(
        output = paste0(scales::comma_format()(mean), " (",
                        scales::comma_format()(ci_lower), "; ",
                        scales::comma_format()(ci_upper), ")")
      ) %>% 
      pivot_wider(names_from = "cost_type", values_from = "output") %>% 
      mutate(group_var = group_var,
             group_label = vlabel[group_var]) %>% 
      relocate(group_var, group_label, group_level = grouping,
               `Baseline (£/yr)` = baseline,
               `Follow-up (£/yr)` = followup,
               `Excess (£/yr)` = excess)
    
  }) %>% 
  bind_rows()

## Save table ----
write_csv(tbl_excess_cost,
          here::here("02_outputs", "09_excess_cost", "tbl_excess_cost.csv"))

