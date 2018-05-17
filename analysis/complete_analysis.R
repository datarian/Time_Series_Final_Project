# The scripts have to be called in the following order.

rm(list=ls())

# Read raw data, build dataframes for later.
source("analysis/prepare_individual_series.R")

# achieve stationarity and build plots / dataframes for the chapter Description
source("analysis/make_series_stationary.R")

# Fit ts models and do diagnostics
source("analysis/model_identification.R")

# Predictions: Impute using different methods
source("analysis/missing_values.R")
