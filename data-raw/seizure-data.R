if (!require("here")) install.packages("here")

# clean data ----
market <- read.csv(here::here("data-raw","acciones.csv"))

# write data in correct format to data folder ----
usethis::use_data(market, overwrite = TRUE)
