library(lubridate)

timestamp = "2022-01-17"
date_censor = as.Date(timestamp) %m-% weeks(8)
date_index_discharge_cutoff = date_censor %m-% weeks(4)