#!/usr/bin/env Rscript
library(dplyr)

# Get the directory where this script is located
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(sub("^--file=", "", script_arg[1]))
} else {
  script_dir <- getwd()
}

# If data files don't exist here, try parent directory
if (!file.exists(file.path(script_dir, "archive/lap_times.csv"))) {
  script_dir <- file.path(script_dir, "..")
}
if (!file.exists(file.path(script_dir, "archive/lap_times.csv"))) {
  script_dir <- file.path(script_dir, "Time-Series-Analysis")
}
setwd(script_dir)

cat("Working directory:", getwd(), "\n")

args <- commandArgs(trailingOnly = TRUE)
driver_id <- if(length(args) >= 1) as.integer(args[1]) else NULL
race_id <- if(length(args) >= 2) as.integer(args[2]) else NULL

lap_times <- read.csv("archive/lap_times.csv")
pit_stops <- read.csv("archive/pit_stops.csv")

if(is.null(race_id)) {
  races <- lap_times %>% 
    group_by(raceId) %>% 
    summarise(n_drivers = n_distinct(driverId), max_laps = max(lap)) %>%
    filter(n_drivers >= 18, max_laps >= 50)
  
  if(!is.null(driver_id)) {
    races <- races %>% 
      filter(raceId %in% unique(lap_times$raceId[lap_times$driverId == driver_id]))
  }
  
  race_id <- races$raceId[nrow(races)]
}

df <- lap_times %>%
  filter(raceId == race_id) %>%
  select(driverId, lap, milliseconds) %>%
  rename(lap_number = lap, lap_time = milliseconds)

pit <- pit_stops %>%
  filter(raceId == race_id) %>%
  select(driverId, lap) %>%
  mutate(pitstop = TRUE)

AllDat <- df %>%
  left_join(pit, by = c("driverId", "lap_number" = "lap")) %>%
  mutate(pitstop = ifelse(is.na(pitstop), FALSE, pitstop)) %>%
  arrange(driverId, lap_number) %>%
  group_by(driverId) %>%
  mutate(pitstop_lagged = as.integer(lag(pitstop, default = FALSE)),
         first_lap = as.integer(lap_number == 1)) %>%
  ungroup() %>%
  select(driverId, first_lap, lap_number, pitstop, pitstop_lagged, lap_time)

write.csv(AllDat, "AllDat.csv", row.names = FALSE)
save(AllDat, file = "DataA1.RData")

cat("Race:", race_id, "| Drivers:", n_distinct(AllDat$driverId), 
    "| Laps:", max(AllDat$lap_number), "\n")
cat("Driver IDs:", paste(sort(unique(AllDat$driverId)), collapse = ", "), "\n")