source("analysis/prepare_individual_series.R")
library(dplyr)

rwl <- read.Dendro.toRWL("data/62167.txt", saveRWL=T)

library(utils)

rwl <- read.rwl("data/900+.rwl", format ="tucson")

summary(rwl)

library(ggplot2)

rwl_to_join <- tibble::rownames_to_column(rwl,var = "Jahr")
