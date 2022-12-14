
# Load packages ----
library(tidyverse)
library(finalfit)
library(broom)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "04_index_mortality"), recursive = TRUE,
           showWarnings = FALSE)

# Create folder rds plots -----
dir.create(here::here("03_plots_rds", "index_mortality"), recursive = TRUE,
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

## Summary table ----------------
# Cohort summary table -----
vars_covariate = c(
  # Prior healthcare use
  "class_factor",
  "prior_emergency_beddays",
  "prior_emergency_beddays_factor",
  
  # Demographics
  "age", "age.factor", "sex", "simd2020_sc_quintile", 
  "ethnicity.factor", "covid_source", "wave", "vacc_status_index",
  
  # Comorbidity
  "n_comorb_charl.factor",
  "ami", "chf", "pvd", "cevd", "dementia", "copd", "rheumd", "pud", 
  "hp", "renal", "HIV",
  "liver.combined", "diab.combined", "cancer.combined",
  
  "n_comorb_PIS.factor",
  
  # Index COVID-19 Admission
  "length_of_stay", "length_of_stay.factor", "any_icu", "icu_days",
  "any_tracheo", "tracheo", "any_endo_tube", "endo_tube",
  "any_haemofilt", "haemofilt", "any_imv", "imv", 
  "any_iv_va", "iv_va"
  
)

## Unstratified ------
tbl_indexmort_cohort_summary = cohort_index %>%
  summary_factorlist(dependent = "mort_in_hosp", 
                     explanatory = vars_covariate,
                     p = FALSE,
                     cont = "median",
                     add_col_totals = TRUE,
                     column = TRUE,
                     add_dependent_label  = FALSE) %>%
  ff_remove_ref2() %>% 
  ff_remove_low_counts()

write_csv(tbl_indexmort_cohort_summary,
          here::here("02_outputs", "04_index_mortality", "tbl_indexmort_cohort_summary.csv"))



# Cause of death --------------------------
tbl_indexmor_cause_of_death = cohort_index %>% 
  filter(mort_in_hosp == "Non-survivors") %>% 
  mutate(
    underlying_cause_of_death.factor = underlying_cause_of_death.factor %>% 
      fct_explicit_na() %>% 
      fct_collapse(Other = c("Other", "(Missing)"))
  ) %>% 
  count(underlying_cause_of_death.factor, name = "non_survivors") %>% 
  mutate(underlying_cause_of_death.factor = if_else(
    non_survivors < 6, "Other", underlying_cause_of_death.factor %>% as.character())) %>% 
  group_by(underlying_cause_of_death.factor) %>% 
  summarise(
    non_survivors = sum(non_survivors)
  ) %>% arrange(desc(non_survivors)) %>%
  mutate(percent = (non_survivors/sum(non_survivors)) %>%
           scales::percent(0.1))

write_csv(tbl_indexmor_cause_of_death,
          here::here("02_outputs", "04_index_mortality", "tbl_indexmort_cause_of_death.csv"))


# Index admission mortality logistic regression ----------

dependent = "mort_in_hosp"
dependent_label = "Index-admission mortality"

## Plot settings ----
breaks = c(0.03, 0.1, 0.3, 1, 3, 10)
remove_ref = TRUE
column_space = c(-0.5, -0.1, 0.2)
table_text_size = 4
title_text_size = 16

predictor_list = 
  list(
    cluster = c(
      "age.factor", "sex", "ethnicity.factor", "simd2020_sc_quintile", 
      "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
      "class_factor"),
    
    prioremerg = c(
      "age.factor", "sex", "ethnicity.factor", "simd2020_sc_quintile", 
      "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
      "prior_emergency_beddays_factor")
  )

# Drop missing ethnicity
cohort_filtered = cohort_index %>% 
  filter(ethnicity.factor != "(Missing)") %>% 
  mutate(ethnicity.factor = ethnicity.factor %>% 
           fct_drop())


walk2(predictor_list, names(predictor_list), 
      function(predictors, pred_label){
      
        # Run logisitc regression ------
        model = cohort_filtered %>%
          glmmulti(dependent, predictors)
        
        # Table of odds ratios -----------
        tbl_model_or = model %>%
          tidy(exponentiate = TRUE)
        
        write_csv(tbl_model_or,
                  here::here("02_outputs", "04_index_mortality", 
                  paste0("tbl_indexmort_coef_", pred_label, ".csv")))
        
        # Table of model metrics ----------
        tbl_model_metrics = model %>% 
          glance()
        
        write_csv(tbl_model_metrics,
                  here::here("02_outputs", "04_index_mortality", 
                             paste0("tbl_indexmort_metrics_", pred_label, ".csv")))
        
        # OR plot ----
        plot_model_or = cohort_filtered %>% 
          or_plot(dependent, predictors,
                  glmfit = model,
                  breaks = breaks,
                  remove_ref = remove_ref,
                  column_space = column_space,
                  table_text_size = table_text_size,
                  title_text_size = title_text_size,
                  dependent_label = dependent_label)
        
        ggsave(filename = here::here("02_outputs", "04_index_mortality", 
                                     paste0("plot_indexmort_or_", pred_label, ".jpeg")),
               plot = plot_model_or,
               width = 10, height = 6, dpi = 600)
        
        write_rds(plot_model_or,
                  here::here("03_plots_rds", "index_mortality", 
                             paste0("plot_indexmort_or_", pred_label, ".rds")))
        
      })

