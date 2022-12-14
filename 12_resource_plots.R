
# Load packages ----
library(tidyverse)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")
source("myfunctions.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "08_resource_plots"), recursive = TRUE,
           showWarnings = FALSE)

# Create folder rds plots -----
dir.create(here::here("03_plots_rds", "resource_plots"), recursive = TRUE,
           showWarnings = FALSE)

# Bootstrap samples
n_bootstrap = 500

# Load datasets ----
cohort_index = read_rds(
  here::here("01_data", paste0("cohort_index_", timestamp, ".rds")))

survivor_index = read_rds(
  here::here("01_data", paste0("survivor_index_", timestamp, ".rds")))

cluster_labeled = read_rds(
  here::here("01_data", paste0("cluster_labeled_", timestamp, ".rds"))) 

resource_emerg_2yr_prior = read_rds(
  here::here("01_data", paste0("resource_emerg_2yr_prior_", timestamp, ".rds"))) 

resource = read_rds(
  here::here("01_data", paste0("resource_", timestamp, ".rds")))

# Join to cohort cluster and prior emergency data ----
cohort_index = cohort_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID") %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))

survivor_index = survivor_index %>% 
  left_join(cluster_labeled, by = "PatientID") %>% 
  left_join(resource_emerg_2yr_prior, by = "PatientID")  %>% 
  mutate(overall = "Overall" %>% 
           ff_label(" "))
  
# Extract labels
vlabel = extract_variable_label(cohort_index)

# Costs by resource type ----
resource_cost_type = resource %>% 
  filter(measurement == "cost") %>% 
  pivot_wider(names_from = "measurement")

# Total cost ----
resource_cost = resource_cost_type %>% 
  group_by(PatientID, period) %>% 
  summarise(n_days = first(n_days),
            cost = sum(cost)) %>% 
  ungroup()

# Plot parameters ---------------
min_n_prop = 0.3
colour_scheme = c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FE6100", "#000000")


# Plot of costs by cluster, resource type and index survival -----------------
resource_summary_type_cluster_surv = resource_cost_type %>% 
  left_join(cohort_index %>% 
              select(PatientID, mort_in_hosp, class_factor),
            by = "PatientID") %>% 
  group_by(period, resource_type, mort_in_hosp, class_factor) %>% 
  summarise(n = n(),
            mean_cost = ci_weighted_mean(cost, n_days, R = n_bootstrap)) %>% 
  unnest(mean_cost)  %>% 
  group_by(resource_type, mort_in_hosp, class_factor) %>% 
  filter(n >= min_n_prop*n[period == 1]) %>% 
  ungroup() %>% 
  mutate(year = period*30/365.25,
         resource_type = if_else(resource_type == "daycase", "Day case",
                            str_to_sentence(resource_type)) %>% 
           factor() %>% fct_relevel("Emergency")
  )

plot_resource_type_cluster_surv = resource_summary_type_cluster_surv %>% 
  filter(period < 0) %>% 
  ggplot(aes(x = year, y = est, ymin = est.L, ymax = est.U,
             colour = mort_in_hosp, fill = mort_in_hosp)) +
  facet_grid(resource_type ~ class_factor, scales = "free_y") +
  geom_line() +
  geom_ribbon(alpha = 0.1, linetype = "dotted") +
  geom_line(data = resource_summary_type_cluster_surv %>% filter(period > 0)) +
  geom_ribbon(alpha = 0.1, linetype = "dotted",
              data = resource_summary_type_cluster_surv %>% filter(period > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  labs(x = "Years relative to index COVID-19 admission",
       y = "Hospital cost (2019 \u00a3 per patient per 30 days)",
       fill = vlabel["mort_in_hosp"],
       colour = vlabel["mort_in_hosp"]) +
  scale_x_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = colour_scheme) +
  scale_colour_manual(values = colour_scheme)

# Save plot as .jpeg
ggsave(here::here("02_outputs", "08_resource_plots",
                  paste0("plot_resource_type_cluster_surv.jpeg")),
       plot = plot_resource_type_cluster_surv,
       width = 11, height = 8, dpi = 600)

# Save plot as .rds
write_rds(plot_resource_type_cluster_surv,
          here::here("03_plots_rds", "resource_plots", 
                     paste0("plot_resource_type_cluster_surv.rds")))

# Plot costs for survivors -------------------------------------------------

## Stratifying variables ----
var_statifaction = c(
  "overall", "age.factor", "sex",  
  "n_comorb_charl.factor", "any_icu", "wave", "vacc_status_index",
  "class_factor", "prior_emergency_beddays_factor"
)

var_statifaction %>% 
  walk(function(group_var){
    
    
    if(group_var == "class_factor" | group_var == "prior_emergency_beddays_factor"){
      legend_rows = 2
    } else {legend_rows = 1}
    
    # Calculate mean cost by period and group
    resource_summary = survivor_index %>%
      select(PatientID, grouping = all_of(group_var)) %>% 
      left_join(resource_cost, by = "PatientID") %>% 
      group_by(period, grouping) %>% 
      summarise(
        n = n(),
        mean_cost = ci_weighted_mean(cost, n_days, R = n_bootstrap)
      ) %>% 
      unnest(mean_cost) 
    
    # Filter out row if sample is smaller than min_n_prop*n of period 1
    resource_summary = resource_summary  %>% 
      group_by(grouping) %>% 
      filter(n >= min_n_prop*n[period == 1]) %>% 
      ungroup() %>% 
      mutate(year = period*30/365.25)
    
    # Extract labels
    vlabel = extract_variable_label(resource_summary)
    
    # Create resource plot
    plot_resource = resource_summary %>%
      filter(period < 0) %>% 
      ggplot(aes(x = year, y = est, ymin = est.L, ymax = est.U,
                 colour = grouping, fill = grouping)) +
      geom_line() +
      geom_ribbon(alpha = 0.1, linetype = "dotted") +
      geom_line(data = resource_summary %>% filter(period > 0)) +
      geom_ribbon(alpha = 0.1, linetype = "dotted",
                  data = resource_summary %>% filter(period > 0)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_bw() +
      theme(legend.position = "bottom",
            text = element_text(size = 12),
            axis.text = element_text(size = 12)) +
      guides(fill = guide_legend(nrow = legend_rows),
            colour = guide_legend(nrow = legend_rows)) +
      labs(x = "\nYears relative to index COVID-19 admission",
           y = "Hospital cost (2019 \u00a3 per patient per 30 days)\n",
           fill = vlabel["grouping"],
           colour = vlabel["grouping"]) +
      scale_x_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 0.5)) +
      scale_y_continuous(limits = c(0, NA)) +
      scale_fill_manual(values = colour_scheme) +
      scale_colour_manual(values = colour_scheme)
    
    # Save plot as .jpeg
    ggsave(here::here("02_outputs", "08_resource_plots",
                      paste0("plot_resource_cost_", group_var, ".jpeg")),
           plot = plot_resource,
           width = 7, height = 6, dpi = 600)
    
    # Save plot as .rds
    write_rds(plot_resource,
              here::here("03_plots_rds", "resource_plots", 
                         paste0("plot_resource_cost_", group_var, ".rds")))
    
  })




