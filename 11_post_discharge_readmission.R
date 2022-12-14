
# Load packages ----
library(tidyverse)
library(finalfit)
library(survival)
library(cmprsk)
library(survminer)
library(scales)
library(broom)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "07_emerg_readmission"), recursive = TRUE,
           showWarnings = FALSE)

# Create folder rds plots -----
dir.create(here::here("03_plots_rds", "emerg_readmission"), recursive = TRUE,
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

# Calculate and join 1st emergency post-discharge readmission -------
survivor_index = survivor_index %>% 
  left_join(
    cohort_admissions %>%
      group_by(PatientID) %>% 
      filter(cis_index > 0, any_emergency == "Yes") %>%
      slice_head() %>% 
      select(PatientID, readmission_date = admission_date,
             readmission_condition = main_condition.factor),
    by = "PatientID"
  )

# Calculate days to readmission and readmission flag -----
survivor_index = survivor_index %>% 
  mutate(
    readm_status = case_when(
      readmission_date <= censor_date ~ 1, # Emergency readmission
      date_of_death == censor_date ~ 2,    # Death as competing risk
      TRUE ~ 0),                           # Alive and censored
    readm_days = case_when(
      readmission_date == discharge_date ~ 0.5,
      TRUE ~ (pmin(readmission_date, censor_date, na.rm = TRUE) -
                    discharge_date) %>% as.numeric()
    )
  )

## Create overall group ----
survivor_index = survivor_index %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))

## Variable labels ----
var_labels = extract_variable_label(survivor_index)

# Main condition for readmission ---------------------------
tbl_emergreadm_main_condition = survivor_index %>% 
  filter(!is.na(readmission_condition)) %>% 
  count(readmission_condition) %>%
  mutate(percent = (n/sum(n)) %>%
           scales::percent(0.1))

## Save table ----
write_csv(tbl_emergreadm_main_condition,
          here::here("02_outputs", "07_emerg_readmission", "tbl_emergreadm_main_condition.csv"))


# Time to readmission censor  ------------------
tbl_emergreadm_time_to_readmission = survivor_index %>%
  mutate(readm_status = readm_status %>% 
           factor() %>% 
           fct_recode("Censored" = "0", "Readmission" = "1", "Died" = "2")) %>% 
  summary_factorlist(dependent = "readm_status", 
                     explanatory = "readm_days",
                     p = FALSE,
                     cont = "median",
                     add_col_totals = TRUE,
                     column = TRUE,
                     add_dependent_label  = FALSE)

## Save table ----
write_csv(tbl_emergreadm_time_to_readmission,
          here::here("02_outputs", "07_emerg_readmission", "tbl_emergreadm_time_to_readmission.csv"))


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
  select(PatientID, readm_days, readm_status, age.factor, 
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
                    readm_status = event(readm_days, readm_status)) %>%
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
        dependent = "Surv(tstart, tstop, ifelse(readm_status == 1, 1, 0))"
        dependent_label = "Post-discharge emergency readmission"
        
        # Run logisitc regression ------
        model = survivor_model %>%
          coxphmulti(dependent, predictors)
        
        # Table of odds ratios -----------
        tbl_model_or = model %>%
          tidy(exponentiate = TRUE)
        
        write_csv(tbl_model_or,
                  here::here("02_outputs", "07_emerg_readmission", 
                             paste0("tbl_emergreadm_hr_", pred_label, ".csv")))
        
        # Table of model metrics ----------
        tbl_model_metrics = model %>% 
          glance()
        
        write_csv(tbl_model_metrics,
                  here::here("02_outputs", "07_emerg_readmission", 
                             paste0("tbl_emergreadm_metrics_", pred_label, ".csv")))
        
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
        
        ggsave(filename = here::here("02_outputs", "07_emerg_readmission", 
                                     paste0("plot_emergreadm_hr_", pred_label, ".jpeg")),
               plot = plot_model_hr,
               width = 10, height = 6, dpi = 600)
        
        write_rds(plot_model_hr,
                  here::here("03_plots_rds", "emerg_readmission", 
                             paste0("plot_emergreadm_hr_", pred_label, ".rds")))
        
      })


# Emergency readmission cumulative incidence ----------------

## Stratifying variables ----
cuminc_statifaction = c(
  "overall", "age.factor", "sex",  
  "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
  "class_factor", "prior_emergency_beddays_factor"
)

## Cumulative incidence table parameters ------
time_points = c(30,60,90,180,270,365)
accuracy = 0.1
alpha = 0.05
time_suffix = "days"

## cuminc parameters -----
ftime = survivor_index$readm_days
fstatus = survivor_index$readm_status

## Create cumulative incidence plots ----
cuminc_statifaction %>% 
  walk(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    if(var_strat == "class_factor" | var_strat == "prior_emergency_beddays_factor"){
      legend.row = 2
    } else {legend.row = 1}
    
    # Create cumulative incidence plot
    plot_ci = ff_ci_plot(
      ftime = ftime, fstatus = fstatus, group = group,
      cencode = 0, eventcode = 1,
      title = NULL, subtitle = NULL, 
      x_max = 360, x_breaks = 60, y_max = 80, y_breaks = 10, 
      fontsize = 3, font.main = 9, plot_confint = TRUE,
      legend.row = legend.row
    )
    
    # Save plot as jpeg
    ggsave(filename = here::here("02_outputs", "07_emerg_readmission",
                                 paste0("plot_emergreadm_cuminc_", var_strat, ".jpeg")),
           plot = plot_ci,
           width = 7, height = 5, dpi = 600)
    
    # Save plot as rds
    write_rds(plot_ci,
              here::here("03_plots_rds", "emerg_readmission", 
                         paste0("plot_emergreadm_cuminc_", var_strat, ".rds")))
    
  })



tbl_emergreadm_cumulative_incidence = as.list(cuminc_statifaction) %>% 
  map(function(var_strat){
    
    # Vector of group variable
    group = survivor_index %>% pull(all_of(var_strat))
    
    # Create cumulative incidence table
    tab_ci_readm = ff_ci_table(ftime = ftime, fstatus = fstatus, group = group,
                               cencode = 0, eventcode = 1, time_points = time_points,
                               accuracy = accuracy, alpha = alpha, time_suffix = time_suffix) %>% 
      rename(level = 1) %>%  # rename first column
      mutate(variable = var_strat,
             variable_desc = var_labels[var_strat])  %>% 
      relocate(variable, variable_desc)
    
    return(tab_ci_readm)
    
  }) %>% 
  bind_rows()

## Save table ----
write_csv(tbl_emergreadm_cumulative_incidence,
          here::here("02_outputs", "07_emerg_readmission", "tbl_emergreadm_cumulative_incidence.csv"))
