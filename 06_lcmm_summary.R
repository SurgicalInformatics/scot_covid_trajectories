
# Load packages ----
library(tidyverse)
library(lcmm)

# Source parameters ----
source("timestamp.R")
source("resource_parameters.R")

# Create folder for LCMM model outputs -----
dir.create(here::here("02_outputs", "02_lcmm_models"), recursive = TRUE,
           showWarnings = FALSE)

# Plot settings ----
colour_scheme = c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FE6100", "#000000")
theme_set(theme_bw())

# Load resource datasets ----
resource_pre = read_rds(here::here("01_data", paste0("resource_pre_", timestamp, ".rds")))
resource_pre_nonzero = read_rds(here::here("01_data", paste0("resource_pre_nonzero_", timestamp, ".rds")))

# Load LCMM models ----
lcmm_models = list.files(here::here("01_data", "lcmm_models"),
                         pattern = "lcmm_\\d{4}-\\d{2}-\\d{2}_\\d+.rds",
                         full.names = TRUE) %>% 
  map(function(lcmm_file){read_rds(lcmm_file)})


# Assigned clusters ----
## Assigned cluster dataset ----
cluster_assignment = full_join(
  
  resource_pre %>%
    group_by(PatientID, PatientID_num) %>%
    summarise(n_clusters = 1:length(lcmm_models)) %>% 
    ungroup(),
  
  lcmm_models %>% 
    map(function(lcmm_model){
      lcmm_model$pprob %>% 
        as_tibble() %>%
        mutate(n_clusters = lcmm_model$ng) %>% 
        select(n_clusters, PatientID_num, class)
    }) %>% 
    bind_rows()
) %>% 
  replace_na(list(class = 0))

## Save dataset ----
write_rds(cluster_assignment,
          here::here("01_data", paste0("cluster_assignment_", timestamp, ".rds")))

# Predicted trajectories ----
## Create new data to fit ----
newdata = data.frame(period = seq(-24, -1, length = 100))

## Table of predicted trajectories ----
tbl_predicted_trajectories = lcmm_models %>% 
  map(function(lcmm_model){
    
    resource_predicted = predictY(lcmm_model, newdata = newdata, 
                                  var.time = "period", draws = TRUE, ndraws = 20) %>% 
      do.call(cbind, .) %>% 
      pivot_longer(cols = -period) %>% 
      mutate(n_clusters = lcmm_model$ng) %>% 
      relocate(n_clusters) %>% 
      mutate(
        metric = case_when(
          str_starts(name, "pred.Ypred_50") ~ "y",
          str_starts(name, "pred.Ypred_2.5") ~ "y_lower",
          str_starts(name, "pred.Ypred_97.5") ~ "y_upper",
        ),
        class = if_else(is.na(str_extract(name, "class(\\d+)$")),
                        "1", str_extract(name, "\\d+$"))
      ) %>%
      select(-name) %>% 
      pivot_wider(names_from = "metric")
    
  }) %>% 
  bind_rows()

write_csv(tbl_predicted_trajectories,
          here::here("02_outputs", "02_lcmm_models", "tbl_predicted_trajectories.csv"))

## Plot of predicted trajectories ----
plot_predicted_trajectories = tbl_predicted_trajectories %>% 
  ggplot(aes(x = period, y = y, ymin = y_lower, ymax = y_upper,
             fill = class, colour = class)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, linetype = "dotted") +
  facet_wrap(~ n_clusters) +
  scale_color_manual(values = colour_scheme) +
  scale_fill_manual(values = colour_scheme) +
  labs(y = "Predicted emergency bed-days per 30 day period",
       x = "30-day periods prior to COVID-19 index admission",
       fill = "Cluster", colour = "Cluster"
  )

ggsave(filename = here::here("02_outputs", "02_lcmm_models", "plot_predicted_trajectories.jpeg"),
       plot = plot_predicted_trajectories,
       width = 7, height = 6, dpi = 600)

# Posterior cluster-membership probability ----
## Table of cluster-membership posterior probability by number of clusters ----
tbl_posterior_probability = lcmm_models %>% 
  map(function(lcmm_model) {
    pprob = lcmm_model$pprob %>% 
      as_tibble() %>% 
      group_by(class) %>% 
      summarise(across(starts_with("prob"), list(mean = mean, sd = sd), .names = "{.col}.{.fn}")) %>% 
      mutate(n_clusters = lcmm_model$ng) %>% 
      pivot_longer(cols = -c(class, n_clusters),
                   names_pattern = "prob(\\d).(.*)",
                   names_to = c("class_prob", "stat"))
  }) %>% 
  bind_rows() %>% 
  arrange(n_clusters, class) %>% 
  relocate(n_clusters, class) %>%
  mutate(class = factor(class),
         n_clusters = factor(n_clusters)) %>% 
  pivot_wider(names_from = "stat") %>% 
  filter(class == class_prob)

write_csv(tbl_posterior_probability,
          here::here("02_outputs", "02_lcmm_models", "tbl_posterior_probability.csv"))

## Plot poster probability ----
plot_posterior_probability = tbl_posterior_probability %>% 
  ggplot(aes(x = n_clusters, y = mean, ymin = mean - sd, ymax = mean + sd, 
             fill = class)) +
  geom_bar(stat = "identity", colour = "black",
           position = "dodge") +
  geom_errorbar(width = 0.2,
                position = position_dodge(0.9)) +
  labs(x = "Total number of clusters",
       y = "Posterior probability",
       fill = "Cluster") +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, by = 0.25)) +
  scale_fill_manual(values = colour_scheme)

ggsave(filename = here::here("02_outputs", "02_lcmm_models", "plot_posterior_probability.jpeg"),
       plot = plot_posterior_probability,
       width = 6, height = 6, dpi = 600)

# Cluster proportion membership ----
# Table of cluster membership counts and proportions ----
tbl_cluster_membership = lcmm_models %>% 
  map(function(lcmm_model) {
  class_count = lcmm_model$pprob %>% 
    as_tibble() %>% 
    count(class) %>% 
    mutate(prop = n/sum(n),
           n_clusters = lcmm_model$ng) %>% 
    relocate(n_clusters)
}) %>% 
  bind_rows() %>% 
  mutate(class = class %>% factor(),
         n_clusters = n_clusters %>% factor())

write_csv(tbl_cluster_membership,
          here::here("02_outputs", "02_lcmm_models", "tbl_cluster_membership.csv"))

# Plot of cluster membership counts and proportions ----
plot_cluster_membership = tbl_cluster_membership %>% 
  ggplot(aes(x = n_clusters, y = prop, fill = class)) +
  geom_bar(stat = "identity", colour = "black",
           position = "dodge") +
  labs(x = "Total number of clusters",
       y = "Cluster membership",
       fill = "Cluster") +
  scale_fill_manual(values = colour_scheme) +
  scale_y_continuous(label = scales::percent_format())

ggsave(filename = here::here("02_outputs", "02_lcmm_models", "plot_cluster_membership.jpeg"),
       plot = plot_cluster_membership,
       width = 7, height = 6, dpi = 600)

# Model fit ----
# Table of AIC, BIC and log-likelihood by number of clusters ----
tbl_model_fit = lcmm_models %>% 
  map(function(lcmm_model) {
    model_fit = tibble(
      n_clusters = lcmm_model$ng,
      BIC = lcmm_model$BIC,
      AIC = lcmm_model$AIC,
      loglik = lcmm_model$loglik)
  }) %>% 
  bind_rows() %>% 
  pivot_longer(cols = -n_clusters, names_to = "metric") %>% 
  mutate(metric = if_else(metric == "loglik", "Log-likelihood", metric))

write_csv(tbl_model_fit,
          here::here("02_outputs", "02_lcmm_models", "tbl_model_fit.csv"))


# Plot of AIC, BIC and log-likelihood by number of clusters ----
plot_model_fit = tbl_model_fit %>% 
  ggplot(aes(x = n_clusters, y = value)) +
  geom_point() + geom_line() +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(x = "Total number of clusters", y = NULL)

ggsave(filename = here::here("02_outputs", "02_lcmm_models", "plot_model_fit.jpeg"),
       plot = plot_model_fit,
       width = 6, height = 6, dpi = 600)
      