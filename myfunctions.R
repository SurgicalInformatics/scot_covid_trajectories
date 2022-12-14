
# ff_ci_plot: creates a cumulative incidence plot with risk and cumulative event tables ----
ff_ci_plot = function(
  # Cumulative incidence parameters
  ftime,
  fstatus,
  group = NULL,
  cencode = 0,
  eventcode = 1,
  factor_levels = NULL,
  time_points = c(35, 360),
  alpha = 0.05,
  plot_confint = FALSE,
  add_risk_table = FALSE,
  
  # Plot parameters
  x_min = 0,
  x_max = 360,
  x_breaks = 30,
  y_min = 0,
  y_max = 100,
  y_breaks = 10,
  colour_scheme = c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FE6100", "#000000"),
  x_lab = "Days",
  y_lab = "Cumulative incidence (%)",
  title = NULL,
  subtitle = NULL,
  group_lab = NULL,
  legend.margin = margin(-5, 0,-5,0),
  legend.row = 1,
  
  # Risk/cumulative event table parameters
  tables.xlim = c(x_min, x_max),
  tables.xbreak = x_breaks,
  tables.y.text = FALSE,
  fontsize = 3.5,
  font.main = 10,
  
  # Combined plot parameters
  rel_heights = c(4.5,1,1)
  
) {
  
  # If no group specified
  if(is.null(group)){
    group = rep("All", length(fstatus))
    group_lab = ""
  } else if (is.null(group_lab)){
    group_lab = attr(group, "label")
  }
  
  # Create tibble from data
  data_temp = tibble(ftime = ftime,
                     fstatus = fstatus,
                     group = group)
  
  # If factor levels not specified
  if(is.null(factor_levels)){
    # If group is not a factor, make a factor
    if (class(group) != "factor"){
      group = factor(group)
    }
    factor_levels = levels(group)
  } else { # Else, factor group with specified factor levels
    group = factor(group, levels = factor_levels)
  }
  
  
  
  
  # Cumulative incidence fit
  ci_fit =
    cuminc(
      ftime = ftime,
      fstatus = fstatus,
      group = group,
      cencode = cencode
    )
  
  # z-score for confidence intervals
  z = qnorm(1-alpha/2)
  
  # Format plot data
  ci_plotdat =
    ci_fit %>% 
    list_modify("Tests" = NULL) %>% 
    map_df(`[`, c("time", "est", "var"), .id = "id") %>%
    filter(id %in% paste0(unique(group), " ", eventcode)) %>% 
    mutate(group = sub(paste0(" ", eventcode, "$"), "", id) %>% 
             factor(levels = factor_levels),
           ci_lower = est^exp(-z*sqrt(var)/ (est*log(est))),
           ci_upper = est^exp( z*sqrt(var)/ (est*log(est))))
  
  
  # Plot cumulative incidence
  ci_plot = ci_plotdat %>% 
    ggplot(aes(x = time, y = est*100, colour = group, fill = group)) +
    geom_step(lwd = 1) +
    theme_bw() +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    scale_x_continuous(breaks = seq(x_min, x_max, x_breaks)) +
    scale_y_continuous(breaks = seq(y_min, y_max, y_breaks)) +
    scale_colour_manual(values = colour_scheme) +
    labs(x = x_lab,
         y = y_lab,
         title = title,
         subtitle = subtitle,
         colour = group_lab,
         fill = group_lab
    ) +
  theme(legend.position = "bottom",
        legend.margin = legend.margin,
        legend.box.margin = margin(0,0,5,0),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)
        ) +
    guides(colour = guide_legend(byrow = TRUE,
                                 nrow = legend.row)
             ) +
    annotate("text", x = 0, y = 0.95*y_max, hjust = 0,
             label = ifelse(is.null(ci_fit$Tests[1,2]),
                            "",
                            paste0(
                              "P-value ",
                              ifelse(ci_fit$Tests[1,2] < 0.001,
                                     "<0.001",
                                     paste0("= ", round(ci_fit$Tests[1,2],3)))
                            )))
  
  
  if(plot_confint == TRUE){
    # Add confidence intervals
    ci_plot = ci_plot + 
      geom_ribbon(aes(ymin = ci_lower*100, ymax= ci_upper*100), alpha = 0.1, linetype = "dotted") +
      scale_fill_manual(values = colour_scheme) +
      labs(fill = group_lab)
  }
  
  if(add_risk_table == TRUE){
  
  # Get number at risk from composite composite of all competing risks
  survfit_composite = surv_fit(
    Surv(ftime, ifelse(fstatus != cencode, 1, 0)) ~ group,
    data = data_temp
  )
  
  ci_tab_composite = ggsurvplot(
    fit = survfit_composite,
    risk.table = TRUE,
    cumevents = TRUE,
    tables.y.text = tables.y.text,
    fontsize = fontsize,
    tables.theme = theme_survminer(font.main = font.main),
    palette = colour_scheme,
    title = "Test",
    xlim = tables.xlim,
    break.time.by = tables.xbreak
  )
  
  # Get number of events of event of interest
  survfit_cause_specific = surv_fit(
    Surv(ftime, ifelse(fstatus == eventcode, 1, 0)) ~ group,
    data = data_temp
  )
  
  ci_tab_cause_specific = ggsurvplot(
    fit = survfit_cause_specific,
    risk.table = TRUE,
    cumevents = TRUE,
    tables.y.text = tables.y.text,
    fontsize = fontsize,
    tables.theme = theme_survminer(font.main = font.main),
    palette = colour_scheme,
    title = "Test",
    xlim = tables.xlim,
    break.time.by = tables.xbreak
  )
  
  # Combine cumulative incidence plot with risk and cumulative event tables
  ci_combined = plot_grid(
    ci_plot,
    ci_tab_composite$table + theme_cleantable(),
    ci_tab_cause_specific$cumevents + theme_cleantable(),
    nrow = 3,
    rel_heights = rel_heights,
    align = "v",
    axis = "b"
  )
  
  # Output plot
  ci_combined
  
  } else{
    
    ci_plot
  }
}

# ff_ci_table: creates a cumulative incidence table at given time points ----
ff_ci_table = function(ftime, fstatus, group = NULL, cencode = 0, eventcode = 1, time_points,
                       accuracy = 0.1, alpha = 0.05, time_suffix = "days", group_lab = NULL
) {
  
  if(is.null(group)){
    group = rep("All", length(fstatus))
    group_lab = " "
  } else if (is.null(group_lab)){
    group_lab = attr(group, "label")
  }
  
  
  ci_fit = 
    cuminc(
      ftime = ftime,
      fstatus = fstatus,
      group = group,
      cencode = cencode
    )
  
  tp_list = ci_fit %>% 
    timepoints(times = time_points)
  
  tbl_est = tp_list$est
  tbl_var = tp_list$var
  
  est = tbl_est %>% 
    as_tibble(rownames = "id") %>% 
    pivot_longer(-id, names_to = "time", values_to = "est") %>% 
    left_join(
      tbl_var %>% 
        as_tibble(rownames = "id") %>% 
        pivot_longer(-id, names_to = "time", values_to = "var"),
      by = c("id", "time")
    )
  
  z = qnorm(1-alpha/2)
  
  est = est %>%
    mutate(ci_lower = est^exp(-z*sqrt(var)/ (est*log(est))),
           ci_upper = est^exp( z*sqrt(var)/ (est*log(est))),
           value_output = paste0(percent(est, accuracy), " (", 
                                 percent(ci_lower, accuracy), ", ", 
                                 percent(ci_upper, accuracy), ")"),
           time_output = paste0(time, " ", time_suffix))
  
  est = est %>%
    filter(id %in% paste0(unique(group), " ", eventcode) | 
             id %in% paste0("All", " ", eventcode)) %>% 
    mutate(group = gsub(paste0(" ", eventcode), "", id)) %>% 
    select(group, time_output, value_output) %>% 
    pivot_wider(names_from = time_output, values_from = value_output) %>% 
    rename(!!quo_name(group_lab) := group)
  
  # Output
  est
}

# ff_km_plot: creates Kaplan-Meier plot ----
ff_km_plot = function(ftime, fstatus, group = NULL,
                      factor_levels = NULL, group_lab = NULL,
                      x_min = 0, x_max = 360, x_breaks = 60,
                      xlab = "Time (days)",
                      legend_rows = 1, plot_risk_table = FALSE){
  
  # If no group specified
  if(is.null(group)){
    group = rep("Overall", length(fstatus))
    group_lab = ""
  } else if (is.null(group_lab)){
    group_lab = attr(group, "label")
  }
  
  # Create tibble from data
  data_temp = tibble(ftime = ftime,
                     fstatus = fstatus,
                     group = group)
  
  # If factor levels not specified
  if(is.null(factor_levels)){
    # If group is not a factor, make a factor
    if (class(group) != "factor"){
      group = factor(group)
    }
    factor_levels = levels(group)
  } else { # Else, factor group with specified factor levels
    group = factor(group, levels = factor_levels)
  }
  
  data_temp = tibble(
    ftime = ftime,
    fstatus = fstatus,
    group = group
  )
  
  survfit = surv_fit(
    Surv(ftime, fstatus) ~ group,
    data = data_temp
  )
  
  plot_surv = ggsurvplot(
    fit = survfit,
    risk.table = TRUE,
    cumevents = TRUE,
    tables.y.text = FALSE,
    censor = FALSE,
    xlab = xlab,
    legend = "bottom",
    legend.labs = data_temp$group %>% levels(),
    legend.title = group_lab,
    pval = TRUE,
    fontsize = 4,
    tables.theme = theme_survminer(font.main = 12, font.legend = 12),
    palette = c("#648FFF", "#DC267F", "#FFB000", "#785EF0", "#FE6100", "#000000"),
    xlim = c(x_min, x_max),
    break.time.by = x_breaks,
    conf.int = TRUE,
    ggtheme = theme_survminer(font.legend = 12,
                              font.main = 12)
  ) +
    guides(colour = guide_legend(byrow = TRUE,
                                 nrow = legend_rows))
  if(plot_risk_table){
  
  km_output = plot_grid(
    plot_surv$plot,
    plot_surv$table + theme_cleantable(),
    plot_surv$cumevents + theme_cleantable(),
    nrow = 3,
    rel_heights = c(4,1.4,1.4),
    align = "v",
    axis = "b"
  )
  
  } else {
    
    km_output = plot_surv$plot
    
  }
  
  return(km_output)
  
}


# ff_km_table: Creates K-M table ----
ff_km_table = function(ftime, fstatus, group, factor_levels = NULL, group_lab = NULL,
                       times = c(30, 60, 90, 180, 365), time_suffix = "days"){
  
  data_temp = tibble(
    ftime = ftime,
    fstatus = fstatus,
    group = group
  )
  
  # If no group specified
  if(is.null(group)){
    group = rep("Overall", length(fstatus))
    group_lab = ""
  } else if (is.null(group_lab)){
    group_lab = attr(group, "label")
  }
  
  # If factor levels not specified
  if(is.null(factor_levels)){
    # If group is not a factor, make a factor
    if (class(group) != "factor"){
      group = factor(group)
    }
    factor_levels = levels(group)
  } else { # Else, factor group with specified factor levels
    group = factor(group, levels = factor_levels)
  }
  
  survfit = 
    surv_fit(
      Surv(ftime, fstatus) ~ group,
      data = data_temp
    )
  
  style_percent_1dp = function(x) {style_percent(x, digits = 1)}
  
  tbl_survfit(
    survfit,
    times=times,
    reverse=TRUE,
    label_header = paste0("{time} ", time_suffix),
    estimate_fun = style_percent_1dp
  ) 
  
}

# ff_remove_low_counts: Obfuscates low counts in finalfit summary table ----
# .data = output of summary_factorlist()
ff_remove_low_counts = function (.data, threshold = 5, 
                                 replace_value = paste0("<", threshold),
                                 ignore = c("label", "levels", "p")){
  if (!any(names(.data) == "label"))
    stop("summary_factorlist() must include: add_dependent_label = FALSE")
  .data %>% 
    dplyr::mutate(across(-dplyr::any_of(ignore),
                         function(.){
                           value_count = as.numeric(stringr::str_extract(., "[:digit:]+"))
                           value_perc = as.numeric(stringr::str_extract(., "(?<=\\().+?(?=\\))"))
                           replace_perc = "-"
                           
                           dplyr::case_when(!levels %in% c("Mean (SD)", "Median (IQR)") &
                                              value_count < threshold ~
                                              paste0(replace_value,
                                                     " (",
                                                     replace_perc,
                                                     ")"),
                                            TRUE ~ .)
                         }))
}

# ff_valid_tests: create summary tables with p values for valid tests only ----

ff_valid_tests = function(.data, dependent, vars, add_dependent_label = TRUE){
  
  explanatory = vars[which(vars != dependent)]
  
  # Create summary table without tests
  summary_pre_tests = .data %>% 
    summary_factorlist(dependent = dependent, 
                       explanatory = explanatory,
                       cont = "median",
                       include_row_missing_col = TRUE,
                       add_row_totals = TRUE,
                       add_col_totals = TRUE,
                       column = TRUE,
                       add_dependent_label  = add_dependent_label,
                       p = FALSE)
  
  # Column 1 (variable labels)
  col_1 = summary_pre_tests %>% select(1)
  
  # Identify explanatory variables producing bad tests
  var_index_bad = summary_pre_tests %>% 
    mutate(var_index = if_else((col_1 == "" | col_1 == "Total N (%)"), 0, 1),
           var_index = cumsum(var_index)) %>%
    filter(if_any(5:(ncol(.)-1), ~ . == "NA (NA to NA)")) %>% 
    pull(var_index)
  
  p_test = replicate(length(explanatory), 1)
  p_test[var_index_bad] = 0
  
  # Create summary table
  map2_df(explanatory, p_test, ~ if(.y ==1){
    summary_factorlist(.data, dependent, explanatory = .x,   
                       cont = "median",
                       include_row_missing_col = TRUE,
                       add_row_totals = TRUE,
                       add_col_totals = TRUE,
                       column = TRUE,
                       add_dependent_label  = add_dependent_label,
                       p = TRUE)
  } else {
    summary_factorlist(.data, dependent, explanatory = .x,
                       cont = "median",
                       include_row_missing_col = TRUE,
                       add_row_totals = TRUE,
                       add_col_totals = TRUE,
                       column = TRUE,
                       add_dependent_label  = add_dependent_label,
                       p = FALSE)
  }) %>% 
    replace_na(list(p = "-")) %>% 
    mutate(across(everything(), ~replace(. , . == "NA (NA to NA)", "-"))) %>% 
    filter(!((.[[1]] == "Total N (%)") & (row_number() != 1)))
}


ci_weighted_mean = function(x, w, probs = c(0.025, 0.975), R = 100){
  
  estimate = weighted.mean(x,w)
  
  est_ci = 1:R %>% 
    map(function(iter){
      i = sample(1:length(x), replace = TRUE)
      est_ci_iter = sum(x[i]*w[i])/sum(w[i])
      est_ci_iter
    }) %>% 
    unlist() %>% 
    quantile(prob = probs)
  
  result = tibble(est = estimate,
                  est.L = est_ci[1],
                  est.U = est_ci[2])
  return(result)
}


ff_remove_ref2 = function (.data, only_binary = TRUE) 
{
  if (!any(names(.data) == "label")) 
    stop("finalfit function must include: add_dependent_label = FALSE")
  df.out = .data %>% dplyr::mutate(label = ifelse(label == 
                                                    "", NA, label)) %>% tidyr::fill(label) %>% dplyr::group_by(label)
  if (only_binary) {
    df.out = df.out %>% dplyr::filter(levels %in% c("Mean (SD)", "Median (IQR)") |
                                        dplyr::row_number() != 1 |
                                        dplyr::n() > 2 |
                                        label %in% c("Total N (%)"))
  }
  else {
    df.out = df.out %>% dplyr::filter(levels %in% c("Mean (SD)", 
                                                    "Median (IQR)") | dplyr::row_number() != 1)
  }
  df.out %>% as.data.frame() %>% rm_duplicate_labels()
}

