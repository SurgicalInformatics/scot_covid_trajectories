
# Load library
library(lubridate)
library(tidyverse)
library(finalfit)

# Source parameters and functions
source("timestamp.R")
source("resource_parameters.R")


## Load datasets ----
cohort_index      = read_rds(here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))
cohort_admissions  = read_rds(here::here("01_data", paste0("cohort_admissions_", timestamp, ".rds")))
cohort_outpatient = read_rds(here::here("01_data", paste0("cohort_outpatient_", timestamp, ".rds")))


# Table with one row /patient/30-day period during look-back and follow-up ----
resource = cohort_index %>% 
  transmute(
    PatientID,
    mort_in_hosp,
    lookback_days = as.numeric((admission_date - days(1)) - lookback_start_date),
    followup_days = as.numeric(censor_date - (discharge_date)),
    lookback_period = lookback_days %/% period_days + 1,
    lookback_remainder = lookback_days %% period_days,
    followup_period = followup_days %/% period_days + 1,
    followup_remainder = followup_days %% period_days
  ) %>% 
  group_by(PatientID) %>% 
  summarise(
    period = c(-lookback_period:followup_period),
    n_days = case_when(
      period == min(period) ~ first(lookback_remainder),
      period == max(period) ~ first(followup_remainder),
      TRUE ~ period_days
    )
  ) %>%
  ungroup() %>% 
  filter(period != 0, n_days != 0, period >= -24)


# Calculate admission beddays and costs by type during period ----
admission_resource = cohort_admissions %>% 
  group_by(
    PatientID, cis_index, admission_date, discharge_date, any_emergency,
    inpatient_daycase_identifier,
    lookback_start_date, censor_date,
    admission_date_index, discharge_date_index,
  ) %>%
  summarise(
    date = seq(admission_date, discharge_date, by = "day")
  ) %>% 
  ungroup() %>%
  mutate(
    period = case_when(
      date < admission_date_index ~ 
        ceiling(as.numeric(date - admission_date_index)/period_days),
      date > discharge_date_index ~
        floor(as.numeric(date - discharge_date_index)/period_days),
      TRUE ~ 0
    ),
    beddays = case_when(
      date == admission_date ~ 0.5,
      date == discharge_date ~ 0.5,
      TRUE ~ 1
    ),
    type = case_when(
      (inpatient_daycase_identifier == "D") & (admission_date == date) ~ "daycase",
      any_emergency == "Yes" ~ "emergency",
      TRUE ~ "elective"
    ),
    cost = case_when(
      type == "daycase" ~ cost_daycase,
      TRUE ~ beddays * cost_inpatient
    )
  ) %>%
  pivot_longer(cols = c(beddays, cost), names_to = "metric") %>% 
  group_by(PatientID, period, type, metric) %>% 
  summarise(
    value = sum(value, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  pivot_wider(names_from = c(metric, type))

# Calculate outpatient appointments and costs ----
outpatient_resource = cohort_outpatient %>% 
  left_join(
    cohort_index %>% 
      select(PatientID, lookback_start_date, censor_date,
             admission_date_index, discharge_date_index),
    by = "PatientID"
  ) %>% 
  mutate(
    period = case_when(
      clinic_date < admission_date_index ~ 
        floor(as.numeric(clinic_date - admission_date_index)/period_days),
      clinic_date > discharge_date_index ~
        ceiling(as.numeric(clinic_date - discharge_date_index)/period_days),
      TRUE ~ 0
    )
  ) %>% 
  group_by(PatientID, period) %>% 
  summarise(
    appointments_outpatient = n(),
    cost_outpatient = cost_outpatient*appointments_outpatient
  ) %>% 
  ungroup()
  

# Join admission and outpatient resource use ----
resource = resource %>% 
  left_join(admission_resource, by = c("PatientID", "period")) %>% 
  left_join(outpatient_resource, by = c("PatientID", "period")) %>% 
  pivot_longer(cols = -c(PatientID, period, n_days),
               names_to = c("measurement", "resource_type"),
               names_pattern = "([[:alpha:]]+)_([[:alpha:]]+)") %>% 
  replace_na(list(value = 0))


# Save resource dataset -----
write_rds(resource,
          file = here::here("01_data", paste0("resource_", timestamp, ".rds")))


