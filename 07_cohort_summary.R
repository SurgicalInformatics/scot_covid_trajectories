
# Load packages ----
library(tidyverse)
library(finalfit)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "03_cohort_summary"), recursive = TRUE,
           showWarnings = FALSE)

# Load datasets ----
cohort_index = read_rds(
  here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))

cluster_assignment = read_rds(
  here::here("01_data", paste0("cluster_assignment_", timestamp, ".rds")))

resource = read_rds(
  here::here("01_data", paste0("resource_", timestamp, ".rds")))

# Join cluster assignments ----

## Lookup table for cluster labels ----
cluster_labels = tribble(
  ~class,      ~class_factor,
  0,           "C1 - No emergency admissions",
  1,           "C3 - Recently high",
  2,           "C4 - Persistently high",
  3,           "C2 - Minimal admissions"
)

## Filter for chosen number of clusters and join with labels
cluster_labeled = cluster_assignment %>% 
  filter(n_clusters == n_clusters_selected) %>%
  left_join(cluster_labels, by = "class") %>% 
  mutate(class_factor = class_factor %>% 
           factor() %>% 
           ff_label("Cluster")) %>% 
  select(PatientID, class_factor)

write_rds(cluster_labeled,
          here::here("01_data", paste0("cluster_labeled_", timestamp, ".rds")))


# Calculate total emergency bed-days 2 years prior to index admission
resource_emerg_2yr_prior = resource %>%
  filter(period < 0, resource_type == "emergency", measurement == "beddays") %>% 
  group_by(PatientID) %>% 
  summarise(
    prior_emergency_beddays = sum(value)
  ) %>% 
  ungroup() %>% 
  left_join(cluster_labeled) %>% 
  mutate(
    prior_emergency_beddays = prior_emergency_beddays %>% 
      ff_label("Emergency bed-days 2 years prior to index admission"),
    prior_emergency_beddays_factor = case_when(
      str_starts(class_factor, "C1") ~ "None",
      prior_emergency_beddays <= 7  ~ "0.5 - 7",
      prior_emergency_beddays <= 21 ~ "7.5 - 21",
      TRUE ~ "21+"
    ) %>% 
      factor() %>% 
      fct_relevel("None", "0.5 - 7", "7.5 - 21", "21+") %>% 
      ff_label("Prior emergency bed-days")
  ) %>% 
  select(-class_factor)

write_rds(resource_emerg_2yr_prior,
          here::here("01_data", paste0("resource_emerg_2yr_prior_", timestamp, ".rds")))

## Joing to cohort data 
cohort_index = cohort_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID")


# Cohort summary table -----
vars_covariate = c(
  # Prior healthcare use
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
  "any_iv_va", "iv_va", "mort_in_hosp"
  
)

## Unstratified ------
tbl_cohort_summary_unstratified = cohort_index %>%
  summary_factorlist(dependent = NULL, 
                     explanatory = vars_covariate,
                     p = FALSE,
                     cont = "median",
                     add_col_totals = TRUE,
                     column = TRUE,
                     add_dependent_label  = FALSE) %>%
  ff_remove_ref2() %>% 
  ff_remove_low_counts()

write_csv(tbl_cohort_summary_unstratified,
          here::here("02_outputs", "03_cohort_summary", "tbl_cohort_summary_unstratified.csv"))



## By Cluster -----
tbl_cohort_summary_by_cluster = cohort_index %>%
  summary_factorlist(dependent = "class_factor", 
                     explanatory = vars_covariate,
                     p = TRUE,
                     cont = "median",
                     add_col_totals = TRUE,
                     column = TRUE,
                     add_dependent_label  = FALSE) %>%
  ff_remove_ref2() %>% 
  ff_remove_low_counts()

write_csv(tbl_cohort_summary_by_cluster,
          here::here("02_outputs", "03_cohort_summary", "tbl_cohort_summary_by_cluster.csv"))


# Weekly index admission date ----
tbl_index_admission_date = cohort_index %>% 
  mutate(admission_month = admission_date %>% floor_date("month")) %>% 
  count(class_factor, admission_month) %>% 
  mutate(n = if_else(n < 6, NA_integer_, n))

plot_index_admission_date = tbl_index_admission_date %>%
  ggplot(aes(admission_month, n)) +
  geom_col() +
  facet_wrap(~class_factor, scales = "free_y") +
  scale_x_date(date_labels = "%b %Y") +
  labs(
    x = NULL, y = "Monthly count",
    caption = "Counts \u2264 5 have been redacted"
  )

ggsave(filename = here::here("02_outputs", "03_cohort_summary", "plot_index_admission_date.jpeg"),
       plot = plot_index_admission_date,
       width = 6, height = 6, dpi = 600)



# Discharge destination for survivors -----
tbl_survivor_discharge_destination = cohort_index %>%
  filter(mort_in_hosp == "Survivors") %>% 
  count(discharge_transfer_to.factor) %>% 
  mutate(perc = (n/sum(n)) %>% scales::percent(0.1))

write_csv(tbl_survivor_discharge_destination,
          here::here("02_outputs", "03_cohort_summary", "tbl_survivor_discharge_destination.csv"))
