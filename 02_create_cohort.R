
# Load packages ----
library(lubridate)
library(tidyverse)
library(finalfit)

# Source timestamp -----
source("timestamp.R")

# Create folder for outputs -----
dir.create(here::here("02_outputs", "01_data_flowchart"), recursive = TRUE,
           showWarnings = FALSE)

# Load datasets ----
linked = read_rds(
  here::here("01_data", paste0("linked_", timestamp, ".rds")))

smr00 = read_rds(
  here::here("01_data", paste0("smr00_", timestamp, ".rds")))

# Create cohort and inclusion flowchart from linked SMR01 ----
## Inclusion criteria ----
# Flag for record inclusion
linked = linked %>% 
  mutate(
    recordID = row_number(),
    discharge_after_start_of_lookback = discharge_date > lookback_start_date,
    admission_prior_to_censor_date = admission_date < censor_date
  )

# Flag for patient inclusion
linked = linked %>%
  group_by(PatientID) %>% 
  mutate(
    covid_admission = any(cis_index == 0),
    index_discharge_before_cutoff = discharge_date[cis_index == 0] <= date_index_discharge_cutoff,
    age_18_and_over = age[cis_index == 0] >= 18,
    survived_index_admission = mort_in_hosp == "Survivors"
  ) %>%
  ungroup() %>%
  replace_na(list(covid_admission = FALSE,
                  index_discharge_before_cutoff = FALSE,
                  age_18_and_over = FALSE,
                  survived_index_admission = FALSE))

inclusion = linked %>% 
  transmute(
    recordID,
    PatientID,
    c0 = TRUE,
    c1 = c0 & covid_admission,
    c2 = c1 & age_18_and_over,
    c3 = c2 & index_discharge_before_cutoff,
    c4 = c3 & discharge_after_start_of_lookback & admission_prior_to_censor_date,
    c5 = c4 & survived_index_admission,
  ) 

# Create flowchart ----
flowchart = inclusion %>% 
  summarise(
    across(.cols = -c(recordID, PatientID),
           .fns = list(
             n_patient = ~n_distinct(PatientID[.x]),
             n_records = ~sum(.x)
             )
           )
  ) %>% 
  pivot_longer(cols = everything(), 
               names_to = c("crit", "count_type"),
               names_pattern = "^(c\\d+)_([[:alpha:]\\_]+)$",
               values_to = "value") %>% 
  pivot_wider(names_from = "count_type") %>% 
  mutate(
    criteria = case_when(
      crit == "c0" ~ "Scottish Safe Haven extract",
      crit == "c1" ~ "COVID admission: hospitalised with COVID-19/tested positive during hospitalisation",
      crit == "c2" ~ "Aged 18 or above during index COVID admission",
      crit == "c3" ~ paste0("Index discharge prior to cut-off of ", date_index_discharge_cutoff),
      crit == "c4" ~ "Admission stay during study period (2 years prior to index admission upto censor date)",
      crit == "c5" ~ "Discharged alive from COVID admission"
    )
  )

# Construct cohort datasets ----
## All admissions within study period for cohort
cohort_admissions = linked %>% 
  left_join(inclusion, by = c("recordID", "PatientID")) %>% 
  filter(c4)

## One-row per patient for cohort, drop unused levels 
cohort_index = cohort_admissions %>% 
  filter(cis_index == 0) %>%
  mutate(across(where(is.factor), fct_drop))
  

## One-row per patient for survivors
survivor_index = cohort_index %>% 
  filter(c5)


# Create cohort and inclusion flowchart from linked SMR01 ----
## Create inclusion flags ----
smr00 = smr00 %>%
  select(PatientID, clinic_date) %>% 
  left_join(
    cohort_index %>% 
      select(PatientID, admission_date, discharge_date,
             lookback_start_date, censor_date,
             covid_index_admission = c4,
             covid_index_survivor = c5),
    by = "PatientID"
  ) %>% 
  mutate(
    smr00_recordID = row_number(),
    clinic_after_start_of_lookback = clinic_date > lookback_start_date,
    clinic_prior_to_censor_date = clinic_date < censor_date
  ) %>% 
  replace_na(
    list(covid_index_admission = FALSE,
         covid_index_survivor = FALSE,
         clinic_after_start_of_lookback = FALSE,
         clinic_prior_to_censor_date = FALSE)
    )

## Calculate inlcusion steps -----
inclusion_smr00 = smr00 %>% 
  transmute(
    smr00_recordID,
    PatientID,
    c0 = TRUE,
    c1 = c0 & covid_index_admission,
    c2 = c1 & clinic_after_start_of_lookback & clinic_prior_to_censor_date,
    c3 = c2 & covid_index_survivor
  )

## Create flowchart
flowchart_smr00 = inclusion_smr00 %>% 
  summarise(
    across(.cols = -c(smr00_recordID, PatientID),
           .fns = list(
             n_patient = ~n_distinct(PatientID[.x]),
             n_records = ~sum(.x)
           )
    )
  ) %>% 
  pivot_longer(cols = everything(), 
               names_to = c("crit", "count_type"),
               names_pattern = "^(c\\d+)_([[:alpha:]\\_]+)$",
               values_to = "value") %>% 
  pivot_wider(names_from = "count_type") %>% 
  mutate(
    criteria = case_when(
      crit == "c0" ~ "Scottish Safe Haven extract",
      crit == "c1" ~ "Patient with index COVID-19 admission meeting study inclusion",
      crit == "c2" ~ "Outpatient clinic within study period (2 years prior to index admission upto censor date)",
      crit == "c3" ~ "Discharged alive from COVID admission"
    )
  )

# Filter SMR00 records for patients meeting inclusion and within study period (c2) ----
cohort_outpatient = smr00 %>% 
  left_join(inclusion_smr00, by = c("smr00_recordID", "PatientID")) %>% 
  filter(c2)


# Save datasets ----
write_rds(survivor_index,
          file = here::here("01_data", paste0("survivor_index_", timestamp, ".rds")))

write_rds(cohort_index,
          file = here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))

write_rds(cohort_admissions,
          file = here::here("01_data", paste0("cohort_admissions_", timestamp, ".rds")))

write_rds(cohort_outpatient,
          file = here::here("01_data", paste0("cohort_outpatient_", timestamp, ".rds")))


# Save flowcharts ----
write_csv(flowchart,
          here::here("02_outputs", "01_data_flowchart", "smr01_cohort_flowchart.csv"))

write_csv(flowchart_smr00,
          here::here("02_outputs", "01_data_flowchart", "smr00_cohort_flowchart.csv"))
  