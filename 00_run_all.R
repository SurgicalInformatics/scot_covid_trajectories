

# Run all scripts ----

source("01_cleaning.R")
rm(list = ls(all.names = TRUE))

source("02_create_cohort.R")
rm(list = ls(all.names = TRUE))

source("03_linkage_summary.R")
rm(list = ls(all.names = TRUE))

source("04_calculate_resource_use.R")
rm(list = ls(all.names = TRUE))

source("05_run_lcmm_models.R")
rm(list = ls(all.names = TRUE))

source("06_lcmm_summary.R")
rm(list = ls(all.names = TRUE))

source("07_cohort_summary.R")
rm(list = ls(all.names = TRUE))

source("08_index_admission_mortality.R")
rm(list = ls(all.names = TRUE))

source("09_post_admission_mortality.R")
rm(list = ls(all.names = TRUE))

source("10_post_discharge_mortality.R")
rm(list = ls(all.names = TRUE))

source("11_post_discharge_readmission.R")
rm(list = ls(all.names = TRUE))

source("12_resource_plots.R")
rm(list = ls(all.names = TRUE))

source("13_excess_costs.R")
rm(list = ls(all.names = TRUE))

source("14_publication_plots.R")
rm(list = ls(all.names = TRUE))