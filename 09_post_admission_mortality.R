
# Load packages ----
library(tidyverse)
library(finalfit)
library(survival)
library(survminer)
library(gtsummary)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for outputs -----
dir.create(here::here("02_outputs", "05_postadm_mortality"), recursive = TRUE,
           showWarnings = FALSE)

# Create folder rds plots -----
dir.create(here::here("03_plots_rds", "postadm_mortality"), recursive = TRUE,
           showWarnings = FALSE)

# Load datasets ----
cohort_index = read_rds(
  here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))

cluster_labeled = read_rds(
  here::here("01_data", paste0("cluster_labeled_", timestamp, ".rds"))) 

resource_emerg_2yr_prior = read_rds(
  here::here("01_data", paste0("resource_emerg_2yr_prior_", timestamp, ".rds"))) 


## Join to cohort cluster and prior emergency data ----
cohort_index = cohort_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID")

# Calculate days to mortality and mortality flag -----
cohort_index = cohort_index %>% 
  mutate(
    mort_status = case_when(
      date_of_death <= censor_date ~ 1, # Death event
      TRUE ~ 0),                        # Alive and censored
    mort_days = (censor_date - admission_date) %>% as.numeric()
  )

# Extract lables -----
vlabel = extract_variable_label(cohort_index)

# Create overall group ----
cohort_index = cohort_index %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))

# Post discharge mortality Kaplan Meier ----------------

## Stratifying variables ----
var_statifaction = c(
  "overall", "age.factor", "sex",  
  "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
  "class_factor", "prior_emergency_beddays_factor"
)

## KM parameters ----
ftime   = cohort_index$mort_days
fstatus = cohort_index$mort_status


## Create Kaplan Meire plots ----
var_statifaction %>% 
  walk(function(var_strat){
    
    # Vector of group variable
    group = cohort_index %>% pull(all_of(var_strat))
    
    if(var_strat == "class_factor" | var_strat == "prior_emergency_beddays_factor"){
      legend_rows = 2
    } else {legend_rows = 1}
    
    # Create K-M plot
    plot_km = ff_km_plot(ftime, fstatus, group, legend_rows = legend_rows)
    
    # Save plot as jpeg
    ggsave(filename = here::here("02_outputs", "05_postadm_mortality",
                                 paste0("plot_postadm_km_", var_strat, ".jpeg")),
           plot = plot_km,
           width = 7, height = 5, dpi = 600)
    
    # Save plot as rds
    write_rds(plot_km,
              here::here("03_plots_rds", "postadm_mortality", 
                         paste0("plot_postadm_km_", var_strat, ".rds")))
    
  })


## Create Kaplan Meire table ----
tbl_postadm_km = var_statifaction %>% 
  map(function(var_strat){
    
    # Vector of group variable
    group = cohort_index %>% pull(all_of(var_strat))
    
    # Temp tibble ----
    data_temp = tibble(
      ftime = cohort_index$mort_days,
      fstatus = cohort_index$mort_status,
      group = group
    )
    
    # Survfit model ----
    survfit = surv_fit(Surv(ftime, fstatus) ~ group, data = data_temp)
    
    # Summary of survfit -----
    summary_survfit = summary(survfit, times = c(0, 30, 60, 90, 180, 365))
    
    if(var_strat == "overall"){summary_survfit$strata = "Overall"}
    
    # Output table ------
    tbl_output = tibble(
      level = summary_survfit$strata %>% str_replace("^group=", ""),
      time = summary_survfit$time,
      n_risk = summary_survfit$n.risk,
      cumulative_event = cumsum(summary_survfit$n.event),
      cumulative_censor = cumsum(summary_survfit$n.censor),
      surv = summary_survfit$surv,
      surv_lower = summary_survfit$lower,
      surv_upper = summary_survfit$upper
    ) %>% 
      mutate(var = var_strat,
             var_desc = vlabel[var_strat]
    ) %>% 
      relocate(var, var_desc)
  }) %>% 
  bind_rows()


## Save table ----
write_csv(tbl_postadm_km,
          here::here("02_outputs", "05_postadm_mortality", "tbl_postadm_km.csv"))


