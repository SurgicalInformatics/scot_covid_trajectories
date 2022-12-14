
# Load packages ----
library(cowplot)
library(gridExtra)

# Create folder for outputs -----
dir.create(here::here("02_outputs", "10_combined_plots"), recursive = TRUE,
           showWarnings = FALSE)


# Post-discharge Kaplan Meier ------------------
## Directory ----
dir_KM = here::here("03_plots_rds","postadm_mortality")

## Load KM plots --------
### A
plot_km_overall = read_rds(
  here::here(dir_KM, "plot_postadm_km_overall.rds"))

### B
plot_km_sex = read_rds(
  here::here(dir_KM, "plot_postadm_km_sex.rds"))

### C
plot_km_class = read_rds(
  here::here(dir_KM, "plot_postadm_km_class_factor.rds"))

### D
plot_km_prior = read_rds(
  here::here(dir_KM, "plot_postadm_km_prior_emergency_beddays_factor.rds")) 

### E
plot_km_age = read_rds(
  here::here(dir_KM, "plot_postadm_km_age.factor.rds"))

### F
plot_km_charlson = read_rds(
  here::here(dir_KM, "plot_postadm_km_n_comorb_charl.factor.rds"))

### G
plot_km_icu = read_rds(
  here::here(dir_KM, "plot_postadm_km_any_icu.rds"))

### H
plot_km_vacc = read_rds(
  here::here(dir_KM, "plot_postadm_km_vacc_status_index.rds"))

## Combine plots -------
plot_KM_combined = ggdraw() +
  draw_plot(plot_km_overall,   x = 0.04, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_km_sex,       x = 0.54, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_km_class,     x = 0.04, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_km_prior,     x = 0.54, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_km_age,       x = 0.04, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_km_charlson,  x = 0.54, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_km_icu,       x = 0.04, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot(plot_km_vacc,      x = 0.54, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot_label(label = LETTERS[1:8],
                  x = rep(c(0, 0.5), 4),
                  y = rep(c(1, 0.75, 0.5, 0.25), each = 2),
                  size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_postadm_KM_combined.jpeg"),
          plot = plot_KM_combined,
          base_height = 3,
          base_width = 6,
          ncol = 2,
          nrow = 4,
          dpi = 600)







# Post-discharge Kaplan Meier ------------------
## Directory ----
dir_KM = here::here("03_plots_rds","postdis_mortality")

## Load KM plots --------
### A
plot_km_overall = read_rds(
  here::here(dir_KM, "plot_postdis_km_overall.rds"))

### B
plot_km_sex = read_rds(
  here::here(dir_KM, "plot_postdis_km_sex.rds"))

### C
plot_km_class = read_rds(
  here::here(dir_KM, "plot_postdis_km_class_factor.rds"))

### D
plot_km_prior = read_rds(
  here::here(dir_KM, "plot_postdis_km_prior_emergency_beddays_factor.rds")) 

### E
plot_km_age = read_rds(
  here::here(dir_KM, "plot_postdis_km_age.factor.rds"))

### F
plot_km_charlson = read_rds(
  here::here(dir_KM, "plot_postdis_km_n_comorb_charl.factor.rds"))

### G
plot_km_icu = read_rds(
  here::here(dir_KM, "plot_postdis_km_any_icu.rds"))

### H
plot_km_vacc = read_rds(
  here::here(dir_KM, "plot_postdis_km_vacc_status_index.rds"))

## Combine plots -------
plot_KM_combined = ggdraw() +
  draw_plot(plot_km_overall,   x = 0.04, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_km_sex,       x = 0.54, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_km_class,     x = 0.04, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_km_prior,     x = 0.54, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_km_age,       x = 0.04, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_km_charlson,  x = 0.54, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_km_icu,       x = 0.04, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot(plot_km_vacc,      x = 0.54, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot_label(label = LETTERS[1:8],
                  x = rep(c(0, 0.5), 4),
                  y = rep(c(1, 0.75, 0.5, 0.25), each = 2),
                  size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_postdis_KM_combined.jpeg"),
          plot = plot_KM_combined,
          base_height = 3,
          base_width = 6,
          ncol = 2,
          nrow = 4,
          dpi = 600)



# Cumulative incidence ----------------
## Directory ----
dir_cuminc = here::here("03_plots_rds", "emerg_readmission")

## Load KM plots --------
### A
plot_cuminc_overall = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_overall.rds"))

### B
plot_cuminc_sex = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_sex.rds"))

### C
plot_cuminc_class = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_class_factor.rds"))

### D
plot_cuminc_prior = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_prior_emergency_beddays_factor.rds")) 

### E
plot_cuminc_age = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_age.factor.rds"))

### F
plot_cuminc_charlson = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_n_comorb_charl.factor.rds"))

### G
plot_cuminc_icu = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_any_icu.rds"))

### H
plot_cuminc_vacc = read_rds(
  here::here(dir_cuminc, "plot_emergreadm_cuminc_vacc_status_index.rds"))

## Combine plots -------
plot_cuminc_combined = ggdraw() +
  draw_plot(plot_cuminc_overall,   x = 0.04, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_sex,       x = 0.54, y = 0.75, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_class,     x = 0.04, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_prior,     x = 0.54, y = 0.50, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_age,       x = 0.04, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_charlson,  x = 0.54, y = 0.25, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_icu,       x = 0.04, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot(plot_cuminc_vacc,      x = 0.54, y = 0.00, width = 0.46, height = 0.25) +
  draw_plot_label(label = LETTERS[1:8],
                  x = rep(c(0, 0.5), 4),
                  y = rep(c(1, 0.75, 0.5, 0.25), each = 2),
                  size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_cuminc_combined.jpeg"),
          plot = plot_cuminc_combined,
          base_height = 3,
          base_width = 6,
          ncol = 2,
          nrow = 4,
          dpi = 600)


# Forest plot (cluster) -------------
## Load plots -----
### A
plot_forest_index_mort_cluster = read_rds(
  here::here("03_plots_rds", "index_mortality", "plot_indexmort_or_cluster.rds"))

### B
plot_forest_postdis_readm_cluster = read_rds(
  here::here("03_plots_rds", "emerg_readmission", "plot_emergreadm_hr_cluster.rds"))

### C
plot_forest_postdis_mort_cluster = read_rds(
  here::here("03_plots_rds", "postdis_mortality", "plot_postdismort_hr_cluster.rds"))

## Combine plots -------
plot_forest_cluster = ggdraw() +
  draw_plot(plot_forest_index_mort_cluster,    x = 0.03, y = 0.667, width = 0.97, height = 0.333) +
  draw_plot(plot_forest_postdis_readm_cluster, x = 0.03, y = 0.333, width = 0.97, height = 0.333) +
  draw_plot(plot_forest_postdis_mort_cluster,  x = 0.03, y = 0.000, width = 0.97, height = 0.333) +
  draw_plot_label(label = LETTERS[1:3],
                  x = rep(0, 3),
                  y = c(1, 0.667, 0.333), size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_forest_cluster.jpeg"),
          plot = plot_forest_cluster,
          base_height = 5,
          base_width = 10.5,
          ncol = 1,
          nrow = 3,
          dpi = 600)

# Forest plot (prior emergency bed-days) ----------
## Load plots -----
### A
plot_forest_index_mort_prioremerg = read_rds(
  here::here("03_plots_rds", "index_mortality", "plot_indexmort_or_prioremerg.rds"))

### B
plot_forest_postdis_readm_prioremerg = read_rds(
  here::here("03_plots_rds", "emerg_readmission", "plot_emergreadm_hr_prioremerg.rds"))

### C
plot_forest_postdis_mort_prioremerg = read_rds(
  here::here("03_plots_rds", "postdis_mortality", "plot_postdismort_hr_prioremerg.rds"))

## Combine plots -------
plot_forest_prioremerg = ggdraw() +
  draw_plot(plot_forest_index_mort_prioremerg,    x = 0.03, y = 0.667, width = 0.97, height = 0.333) +
  draw_plot(plot_forest_postdis_readm_prioremerg, x = 0.03, y = 0.333, width = 0.97, height = 0.333) +
  draw_plot(plot_forest_postdis_mort_prioremerg,  x = 0.03, y = 0.000, width = 0.97, height = 0.333) +
  draw_plot_label(label = LETTERS[1:3],
                  x = rep(0, 3),
                  y = c(1, 0.667, 0.333), size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_forest_prioremerg.jpeg"),
          plot = plot_forest_prioremerg,
          base_height = 5,
          base_width = 10.5,
          ncol = 1,
          nrow = 3,
          dpi = 600)



# Resource use -----------------------
## Load plots -----
### A
plot_resource_class = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_class_factor.rds"))

### B
plot_resource_piroemerg = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_prior_emergency_beddays_factor.rds"))

### C
plot_resource_age = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_age.factor.rds"))

### D
plot_resource_charlson = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_n_comorb_charl.factor.rds"))

### E
plot_resource_icu = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_any_icu.rds"))

### F
plot_resource_vacc = read_rds(
  here::here("03_plots_rds/resource_plots/plot_resource_cost_vacc_status_index.rds"))


## Combine plots -------
plot_resource_combined = ggdraw() +
  draw_plot(plot_resource_class,     x = 0.03, y = 0.667, width = 0.47, height = 0.33) +
  draw_plot(plot_resource_piroemerg, x = 0.53, y = 0.667, width = 0.47, height = 0.33) +
  draw_plot(plot_resource_age,       x = 0.03, y = 0.333, width = 0.47, height = 0.33) +
  draw_plot(plot_resource_charlson,  x = 0.53, y = 0.333, width = 0.47, height = 0.33) +
  draw_plot(plot_resource_icu,       x = 0.03, y = 0.000, width = 0.47, height = 0.33) +
  draw_plot(plot_resource_vacc,      x = 0.53, y = 0.000, width = 0.47, height = 0.33) +
  draw_plot_label(label = LETTERS[1:6],
                  x = rep(c(0, 0.5), 3),
                  y = rep(c(1, 0.667, 0.333), each = 2),
                  size = 20)

save_plot(here::here("02_outputs", "10_combined_plots", "plot_resource_combined.jpeg"),
          plot = plot_resource_combined,
          base_height = 5,
          base_width = 6,
          ncol = 2,
          nrow = 3,
          dpi = 600)


