
# Load packages ----
library(lubridate)
library(tidyverse)
library(finalfit)

# Source timestamp -----
source("timestamp.R")
source("filepath_linked_data.R")

# Create folder for outputs -----
dir.create(here::here("01_data"), recursive = TRUE, showWarnings = FALSE)

# Load datasets ----
linked = read_rds(here::here(filepath_linked_data,
                             paste0("linked_", timestamp, ".rds")))

vaccination = read_rds(here::here(filepath_linked_data,
                                  paste0("Vaccinations_", timestamp, ".rds")))

smr00 = read_rds(here::here(filepath_linked_data,
                                  paste0("smr00_", timestamp, ".rds")))

## Collapse/recode factors ----
linked = linked %>% 
  mutate(
    n_comorb_PIS = if_else(is.na(n_comorb_PIS), 0, n_comorb_PIS),
    n_comorb_PIS.factor = case_when(
      n_comorb_PIS < 3 ~ "0-2",
      n_comorb_PIS < 6 ~ "3-5",
      n_comorb_PIS >= 6 ~ "6+",
    ) %>% 
      factor() %>% 
      ff_label("Prescriptions"),
    
    covid_source = covid_source %>% 
      fct_collapse(Community = c("Community", "Probable community"),
                   Nosocomial = c("Nosocomial", "Probable nonsocomial")),
    
    ethnicity.factor = ethnicity.factor %>% 
      fct_collapse(White = "White", "(Missing)" = "(Missing)", other_level = "Other ethnicity") %>% 
      fct_relevel("(Missing)", after = Inf) %>% 
      ff_label("Ethnicity"),
    
    # Length of stay ----
    length_of_stay.factor = case_when(
      length_of_stay < 2.5 ~ "0.5-2",
      length_of_stay < 8 ~ "3-7",
      length_of_stay < 22 ~ "8-21",
      TRUE ~ "22+") %>%
      factor() %>%
      fct_relevel("0.5-2", "3-7", "8-21", "22+") %>% 
      ff_label("Length of stay (days)")
  )

# Covid Wave ----
linked = linked  %>% 
  left_join(
    linked %>% 
      filter(cis_index == 0) %>%
      mutate(wave_date = coalesce(specimen_date, admission_date),
             wave = case_when(
               wave_date < as.Date("2020-09-01") ~ "Wave 1",
               wave_date < as.Date("2021-05-01") ~ "Wave 2",
               TRUE                              ~ "Wave 3") %>%
               factor() %>% 
               ff_label("COVID-19 wave")) %>% 
      select(PatientID, wave),
    by = "PatientID"
  )


# ICU treatment, make NA if 0 ----
linked = linked %>% 
  mutate(across(.cols = c(imv, tracheo, endo_tube, haemofilt, iv_va),
                .fns = ~na_if(.x, 0))
  )

# Index admission survival ----
linked = linked %>%
  left_join(
    linked %>% 
      filter(cis_index == 0) %>% 
      mutate(
        mort_in_hosp = case_when(
          discharge_transfer_to.factor == "Died" ~ "Non-survivors",
          discharge_date == date_of_death ~ "Non-survivors",
          TRUE ~ "Survivors"
        ) %>% 
          factor(levels = c("Survivors", "Non-survivors")) %>% 
          ff_label("Index admission survival")) %>% 
      select(PatientID, mort_in_hosp),
    by = "PatientID"
  )


# Lookback start date and censor date ----
linked = linked %>%
  group_by(PatientID) %>% 
  mutate(
    admission_date_index = admission_date[cis_index == 0],
    discharge_date_index = discharge_date[cis_index == 0]
  ) %>% 
  ungroup() %>% 
  mutate(
    lookback_start_date = admission_date_index - days(365*2),
    censor_date = pmin(date_of_death, date_censor, na.rm = TRUE),
    study_period = case_when(
      cis_index < 0 ~ "Lookback",
      cis_index == 0 ~ "Index",
      cis_index > 0 ~ "Follow-up",
      TRUE ~ NA_character_
    )
  )


# Link vaccination -----
## Clean vaccination data ----
vaccination = vaccination %>% 
  filter(vacc_product_name != "COVID-19 Vaccine (Not Administered)") %>%
  mutate(vacc_occurence_time = ymd(vacc_occurence_time)) %>%
  filter(vacc_occurence_time > ymd("2020-05-01")) %>% 
  arrange(vacc_occurence_time) %>% 
  group_by(PatientID) %>% 
  mutate(
    vacc_index = row_number()) %>% 
  ungroup()

## Vaccinated status for index admission ----
linked = linked %>% 
  left_join(
    vaccination %>% 
      left_join(
        linked %>% 
          filter(cis_index == 0) %>% 
          select(PatientID, admission_date),
        by = "PatientID"
      ) %>% 
      mutate(
        vacc_status_adm = case_when(
          (vacc_index == 1) & 
            (vacc_occurence_time + days(21) < admission_date) &
            (vacc_occurence_time + days(270) > admission_date) ~ TRUE,
          (vacc_index > 1) & 
            (vacc_occurence_time + days(14) < admission_date) &
            (vacc_occurence_time + days(270) > admission_date) ~ TRUE,
          TRUE ~ FALSE
        )
      ) %>% 
      group_by(PatientID) %>% 
      summarise(
        vacc_status_index = if_else(any(vacc_status_adm), "Yes", "No")
      ),
    by = "PatientID"
  ) %>% 
  replace_na(list(vacc_status_index = "No")) %>% 
  mutate(
    vacc_status_index = vacc_status_index %>% 
      factor() %>%
      fct_relevel("No") %>% 
      ff_label("Vaccinated")
  )
 
## Append vaccination dates ----
linked = linked %>% 
  left_join(
    vaccination %>%
      select(PatientID, vacc_occurence_time, vacc_index) %>% 
      pivot_wider(names_from = vacc_index, values_from = vacc_occurence_time,
                  names_prefix = "vacc_date_"),
    by = "PatientID"
  )

# Save cleaned dataset ----
write_rds(linked,
          file = here::here("01_data", paste0("linked_", timestamp, ".rds")))

# Clean SMR00 ----

smr00 = smr00 %>% 
  mutate(clinic_date = ymd(clinic_date))

# Save cleaned SMR00 dataset ----
write_rds(smr00,
          file = here::here("01_data", paste0("smr00_", timestamp, ".rds")))