# The scripts have to be called in the order specified.

rm(list=ls())

# Set defaults for ggplot -> must berun, it builds the colour palette referenced
# in plots!
source("latex/plotDefaults.R")

# Read raw data, construct data frame
source("analysis/prepare_individual_series.R")

# Achieve stationarity and build exploratory plots / data frames
source("analysis/make_series_stationary.R")

# Fit ts models and do diagnostics, build plots / data frames for report
source("analysis/model_identification.R")

# Predictions: Impute using different methods + build plots / data frames for report
source("analysis/missing_values.R")
