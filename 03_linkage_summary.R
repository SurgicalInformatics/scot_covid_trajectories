
# Load packages ----
library(tidyverse)
library(lubridate)

# Source timestamp -----
source("timestamp.R")
source("filepath_linked_data.R")

# Load cohort index admission ----

linked = read_rds(
  here::here("01_data", paste0("linked_", timestamp, ".rds")))

cohort_index = read_rds(
  here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))

vaccinations = read_rds(
  here::here(filepath_linked_data, paste0("Vaccinations_", timestamp, ".rds")))

ecoss = read_rds(
  here::here(filepath_linked_data, paste0("COVID_TESTING_", timestamp, ".rds")))

sicsag = read_rds(
  here::here(filepath_linked_data, paste0("SICSAG_episode_", timestamp, ".rds")))

prescribing = read_rds(
  here::here(filepath_linked_data, paste0("PIS_", timestamp, ".rds")))

deaths = read_rds(
  here::here(filepath_linked_data, paste0("Deaths_", timestamp, ".rds")))


# ECOSS linkage ----

ecoss = ecoss %>% 
  mutate(recordID = row_number(),
         specimen_date = ymd(specimen_date))

ecoss_covid_admission = ecoss %>%
  left_join(
    linked %>% 
      select(PatientID, admission_date, discharge_date),
    by = "PatientID"
  ) %>% 
  mutate(flag_covid_admission =
           (specimen_date >= admission_date - days(14)) & 
           (specimen_date <= discharge_date)
  ) %>% 
  group_by(recordID) %>% 
  summarise(
    flag_covid_admission = any(flag_covid_admission)
  )

ecoss = ecoss %>%
  left_join(ecoss_covid_admission, by = "recordID") %>% 
  replace_na(list(flag_covid_admission = FALSE))

linkage_summary_ecoss = ecoss %>% 
  summarise(
    extract_records = n(),
    extract_patients = n_distinct(PatientID),
    linked_covid_admission_records = sum(flag_covid_admission),
    linked_covid_admission_patients = n_distinct(PatientID[flag_covid_admission])
  )

# Save linkage summary ----
write_csv(linkage_summary_ecoss,
          here::here("02_outputs", "01_data_flowchart", "linkage_summary_ecoss.csv"))

# Prescribing dataset ----

prescribing = prescribing %>% 
  mutate(prescribed_date = ymd(prescribed_date)) %>% 
  left_join(cohort_index %>% select(PatientID, admission_date), by = "PatientID") %>% 
  mutate(
    flag_1_year_prior = prescribed_date >= (admission_date - days(365)) &
      (prescribed_date < admission_date)
  ) 

linkage_summary_pis = prescribing %>% 
  summarise(
    extract_records = n(),
    extract_patients = n_distinct(PatientID),
    linked_covid_admission_1yr_records = sum(flag_1_year_prior, na.rm = TRUE),
    linked_covid_admission_1yr_patients = n_distinct(PatientID[flag_1_year_prior])
  )

## Save linkage summary ----
write_csv(linkage_summary_pis,
          here::here("02_outputs", "01_data_flowchart", "linkage_summary_pis.csv"))

# SICSAG dataset ----

sicsag = sicsag %>% 
  mutate(
    recordID = row_number(),
    admit_unit = ymd(admit_unit))

sicsag_covid_admission = sicsag %>% 
  left_join(
    cohort_index %>% 
      select(PatientID, admission_date, discharge_date), by = "PatientID") %>% 
  mutate(flag_covid_admission_icu = (admit_unit >= admission_date) & 
           (admit_unit <= discharge_date)) %>% 
  group_by(recordID) %>% 
  summarise(
    flag_covid_admission_icu = any(flag_covid_admission_icu)
  ) 

sicsag = sicsag %>%
  left_join(sicsag_covid_admission, by = "recordID") %>% 
  replace_na(list(flag_covid_admission_icu = FALSE))

linkage_summary_sicsag = sicsag %>% 
  summarise(
    extract_records = n(),
    extract_patients = n_distinct(PatientID),
    linked_covid_admission_records = sum(flag_covid_admission_icu),
    linked_covid_admission_patients = n_distinct(PatientID[flag_covid_admission_icu])
  )

## Save linkage summary ----
write_csv(linkage_summary_sicsag,
          here::here("02_outputs", "01_data_flowchart", "linkage_summary_sicsag.csv"))

# Vaccination dataset ----
vaccinations = vaccinations %>% 
  mutate(vacc_occurence_time = ymd(vacc_occurence_time),
         recordID = row_number()) %>% 
  group_by(PatientID) %>% 
  arrange(PatientID, vacc_occurence_time) %>% 
  mutate(vacc_index = row_number()) %>% 
  ungroup()

vaccinations_admission = vaccinations %>% 
  left_join(
    cohort_index %>% 
      select(PatientID, admission_date), by = "PatientID") %>% 
  mutate(
    flag_vacc_admission = case_when(
      vacc_product_name == "COVID-19 Vaccine (Not Administered)" ~ FALSE,
      vacc_occurence_time <= ymd("2020-05-01") ~ FALSE,
      (vacc_index == 1) & 
        (vacc_occurence_time + days(21) < admission_date) &
        (vacc_occurence_time + days(270) > admission_date) ~ TRUE,
      (vacc_index > 1) & 
        (vacc_occurence_time + days(14) < admission_date) &
        (vacc_occurence_time + days(270) > admission_date) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>% 
  group_by(recordID) %>% 
  summarise(flag_vacc_admission = any(flag_vacc_admission))
           
vaccinations = vaccinations %>%
  left_join(vaccinations_admission, by = "recordID") %>% 
  replace_na(list(flag_vacc_admission = FALSE))

linkage_summary_vaccinations = vaccinations %>% 
  summarise(
    extract_records = n(),
    extract_patients = n_distinct(PatientID),
    linked_covid_admission_records = sum(flag_vacc_admission),
    linked_covid_admission_patients = n_distinct(PatientID[flag_vacc_admission])
  )

## Save linkage summary ----
write_csv(linkage_summary_vaccinations,
          here::here("02_outputs", "01_data_flowchart", "linkage_summary_vaccinations.csv"))

# Deaths dataset ----

linkage_summary_deaths = deaths %>% 
  mutate(date_of_death = ymd(date_of_death)) %>% 
  left_join(cohort_index %>% 
              transmute(PatientID, flag_death = TRUE), by = "PatientID") %>% 
  summarise(
    extract_records = n(),
    linked_covid_admission_records = sum(flag_death, na.rm = TRUE),
  )

## Save linkage summary ----
write_csv(linkage_summary_deaths,
          here::here("02_outputs", "01_data_flowchart", "linkage_summary_deaths.csv"))
           
           
           
           