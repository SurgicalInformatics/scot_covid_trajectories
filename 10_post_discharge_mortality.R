
# Load packages ----
library(tidyverse)
library(finalfit)
library(survival)
library(survminer)
library(scales)
library(broom)
library(gtsummary)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for outputs -----
dir.create(here::here("02_outputs", "06_postdis_mortality"), recursive = TRUE,
           showWarnings = FALSE)

# Create folder rds plots -----
dir.create(here::here("03_plots_rds", "postdis_mortality"), recursive = TRUE,
           showWarnings = FALSE)

# Load datasets ----
survivor_index = read_rds(
  here::here("01_data", paste0("survivor_index_", timestamp, ".rds")))

cohort_admissions = read_rds(
  here::here("01_data", paste0("cohort_admissions_", timestamp, ".rds")))

cluster_labeled = read_rds(
  here::here("01_data", paste0("cluster_labeled_", timestamp, ".rds"))) 

resource_emerg_2yr_prior = read_rds(
  here::here("01_data", paste0("resource_emerg_2yr_prior_", timestamp, ".rds"))) 

## Join to cohort cluster and prior emergency data ----
survivor_index = survivor_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID")

# Extract lables -----
vlabel = extract_variable_label(survivor_index)

# Calculate days to mortality and mortality flag -----
survivor_index = survivor_index %>% 
  mutate(
    mort_status = case_when(
      date_of_death <= censor_date ~ 1, # Death event
      TRUE ~ 0),                        # Alive and censored
    mort_days = (censor_date - discharge_date) %>% as.numeric()
  )

## Create overall group ----
survivor_index = survivor_index %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))

## Variable labels ----
var_labels = extract_variable_label(survivor_index)

# Main cuase of death ---------------------------
tbl_postdismort_cause_of_death = survivor_index %>% 
  filter(date_of_death <= date_censor) %>% 
  mutate(underlying_cause_of_death.factor = underlying_cause_of_death.factor %>% 
           fct_collapse(Other = "(Missing)")) %>% 
  count(underlying_cause_of_death.factor) %>% 
  mutate(percent = (n/sum(n)) %>%
           scales::percent(0.1))

## Save table ----
write_csv(tbl_postdismort_cause_of_death,
          here::here("02_outputs", "06_postdis_mortality", "tbl_postdismort_cause_of_death.csv"))


# Time to death/censor  ------------------
tbl_postdismort_time_to_death = survivor_index %>%
  mutate(mort_status = mort_status %>% 
           factor() %>% 
           fct_recode("Censored" = "0", "Died" = "1")) %>% 
  summary_factorlist(dependent = "mort_status", 
                     explanatory = "mort_days",
                     p = FALSE,
                     cont = "median",
                     add_col_totals = TRUE,
                     column = TRUE,
                     add_dependent_label  = FALSE)

## Save table ----
write_csv(tbl_postdismort_time_to_death,
          here::here("02_outputs", "06_postdis_mortality", "tbl_postdismort_time_to_death.csv"))



# Post-discharge death location --------------
survivor_index = survivor_index %>% 
  left_join(
    survivor_index %>%
      filter(mort_status == 1) %>%
      transmute(PatientID) %>% 
      left_join(
        cohort_admissions %>% 
          filter(cis_index > 0, discharge_transfer_to.factor == "Died") %>%
          group_by(PatientID) %>% 
          summarise(postdis_mort_location = "Hospital"),
        by = "PatientID"
      ) %>% 
      replace_na(list(postdis_mort_location = "Community")),
    by = "PatientID"
  )

## Table of location of death ----
tbl_postdismort_location = survivor_index %>% 
  filter(mort_status == 1) %>% 
  count(postdis_mort_location) %>% 
  mutate(percentage = (n/sum(n)) %>% scales::percent(0.1))

## Save table ----
write_csv(tbl_postdismort_location,
          here::here("02_outputs", "06_postdis_mortality", "tbl_postdismort_location.csv"))



# Cox regression (vaccination as time varying covariate) --------------------

## Plot settings ----
breaks = c(0.03, 0.1, 0.3, 1, 3, 10)
remove_ref = TRUE
column_space = c(-0.5, -0.1, 0.2)
table_text_size = 4
title_text_size = 16

## Filter out missing data ----
survivor_filtered = survivor_index %>% 
  filter(ethnicity.factor != "(Missing)", !is.na(simd2020_sc_quintile)) %>% 
  mutate(ethnicity.factor = ethnicity.factor %>% 
           fct_drop()
  )

## Baseline covariates ----
survivor_baseline = survivor_filtered %>% 
  select(PatientID, mort_days, mort_status, age.factor, 
         sex, ethnicity.factor, simd2020_sc_quintile, wave, n_comorb_charl.factor,
         any_icu, wave, class_factor, prior_emergency_beddays_factor
  )

## Time varying covariate - Vaccination -----
survivor_tvc = survivor_filtered %>% 
  transmute(PatientID,
            vacc_day = (vacc_date_1 - discharge_date) %>% as.numeric() + 21,
            vacc_status = if_else(is.na(vacc_day), 0, 1)
  ) %>% 
  drop_na()

## Model data input ----
survivor_model = tmerge(survivor_baseline, survivor_baseline,
                        id = PatientID,
                        mort_status = event(mort_days, mort_status)) %>%
  tmerge(survivor_tvc, id = PatientID, vacc_status = tdc(vacc_day, vacc_status)) %>% 
  replace_na(list(vacc_status = 0)) %>% 
  ff_relabel(var_labels) %>% 
  mutate(vacc_status = vacc_status %>% 
           factor() %>% 
           fct_recode("No" = "0", "Yes" = "1") %>% 
           ff_label("Vaccinated"))

## Predictor list ----
predictor_list = 
  list(
    cluster = c(
      "age.factor", "sex", "ethnicity.factor", "simd2020_sc_quintile", 
      "n_comorb_charl.factor", "any_icu", "wave", "vacc_status",
      "class_factor"),
    
    prioremerg = c(
      "age.factor", "sex", "ethnicity.factor", "simd2020_sc_quintile", 
      "n_comorb_charl.factor", "any_icu", "wave", "vacc_status",
      "prior_emergency_beddays_factor")
  )


## Run cox regression for each predictor list ----
walk2(predictor_list, names(predictor_list), 
      function(predictors, pred_label){
        
        # Dependent variable -----
        dependent = "Surv(tstart, tstop, mort_status)"
        dependent_label = "Post-discharge mortality"
        
        # Run logisitc regression ------
        model = survivor_model %>%
          coxphmulti(dependent, predictors)
        
        # Table of odds ratios -----------
        tbl_model_or = model %>%
          tidy(exponentiate = TRUE)
        
        write_csv(tbl_model_or,
                  here::here("02_outputs", "06_postdis_mortality", 
                             paste0("tbl_postdismort_coef_", pred_label, ".csv")))
        
        # Table of model metrics ----------
        tbl_model_metrics = model %>% 
          glance()
        
        write_csv(tbl_model_metrics,
                  here::here("02_outputs", "06_postdis_mortality", 
                             paste0("tbl_postdismort_metrics_", pred_label, ".csv")))
        
        # HR plot ----
        plot_model_hr = survivor_model %>% 
          hr_plot(dependent, predictors,
                  coxfit = model,
                  breaks = breaks,
                  remove_ref = remove_ref,
                  column_space = column_space,
                  table_text_size = table_text_size,
                  title_text_size = title_text_size,
                  dependent_label = dependent_label)
        
        ggsave(filename = here::here("02_outputs", "06_postdis_mortality", 
                                     paste0("plot_postdismort_hr_", pred_label, ".jpeg")),
               plot = plot_model_hr,
               width = 10, height = 6, dpi = 600)
        
        write_rds(plot_model_hr,
                  here::here("03_plots_rds", "postdis_mortality", 
                             paste0("plot_postdismort_hr_", pred_label, ".rds")))
        
      })


# Post discharge mortality Kaplan Meier ----------------

## Stratifying variables ----
var_statifaction = c(
  "overall", "age.factor", "sex",  
  "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
  "class_factor", "prior_emergency_beddays_factor"
)

## KM parameters ----
ftime   = survivor_index$mort_days
fstatus = survivor_index$mort_status


## Create Kaplan Meire plots ----
var_statifaction %>% 
  walk(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    if(var_strat == "class_factor" | var_strat == "prior_emergency_beddays_factor"){
      legend_rows = 2
    } else {legend_rows = 1}
    
    # Create K-M plot
    plot_km = ff_km_plot(ftime, fstatus, group, legend_rows = legend_rows)
    
    # Save plot as jpeg
    ggsave(filename = here::here("02_outputs", "06_postdis_mortality",
                                 paste0("plot_postdis_km_", var_strat, ".jpeg")),
           plot = plot_km,
           width = 7, height = 5, dpi = 600)
    
    # Save plot as rds
    write_rds(plot_km,
              here::here("03_plots_rds", "postdis_mortality", 
                         paste0("plot_postdis_km_", var_strat, ".rds")))
    
  })

# Create table of mortality over time ----
tbl_postdismort_mortality = as.list(var_statifaction) %>% 
  map(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    # Create cumulative incidence table
    tab_km_postdis = ff_km_table(ftime, fstatus, group) %>%
      as_tibble() %>% 
      rename(level = 1) %>%  # rename first column
      mutate(variable = var_strat,
             variable_desc = var_labels[var_strat]) %>% 
      relocate(variable, variable_desc)
    
    return(tab_km_postdis)
    
  }) %>% 
  bind_rows() %>% 
  drop_na()

# Create table of mortality over time ----
tbl_postdismort_mortality = as.list(var_statifaction) %>% 
  map(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    # Create cumulative incidence table
    tab_km_postdis = ff_km_table(ftime, fstatus, group) %>%
      as_tibble() %>% 
      rename(level = 1) %>%  # rename first column
      mutate(variable = var_strat,
             variable_desc = var_labels[var_strat]) %>% 
      relocate(variable, variable_desc)
    
    return(tab_km_postdis)
    
  }) %>% 
  bind_rows() %>% 
  drop_na()



## Save table ----
write_csv(tbl_postdismort_mortality,
          here::here("02_outputs", "06_postdis_mortality", "tbl_postdismort_mortality.csv"))



## Create Kaplan Meire table ----
tbl_postdis_km = var_statifaction %>% 
  map(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    # Temp tibble ----
    data_temp = tibble(
      ftime = survivor_index$mort_days,
      fstatus = survivor_index$mort_status,
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
write_csv(tbl_postdis_km,
          here::here("02_outputs", "06_postdis_mortality", "tbl_postdis_km.csv"))



