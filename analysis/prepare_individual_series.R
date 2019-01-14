# This script parses the raw data file 900+.txt and builds a dataframe in 'rwl'
# format, a standard for dendrochronological data.

library(readtext)
library(stringr)
library(dplR)
library(itsmr)
library(MASS)


read.Dendro.toRWL <- function(filepath, saveRWL=T, format="tucson"){

file_rwl_name <- paste0(str_match(filepath,"(.*).txt")[2],".rwl")
file_contents <- readtext(filepath)

# Extract the individual data series from the unstructured text file
pattern <- "\\sNr\\. ([0-9]*)\\.\\d{1} Dat: ([0-9]{3,4}).*\\n.*\\n.*\\n{1,2}((?:[ ]{1,3}[0-9 -]*\\n){3,})"
parsed <- str_match_all(file_contents,pattern)
rm(file_contents) # Free up memory

# Extract the sample identification numbers and their corresponding dated year
chronology <- as.data.frame(parsed[[1]][,-1],stringsAsFactors = F)
rm(parsed) # Free up memory


ring_pattern <- "\\s{1,3}[0-9 -]{8,11}\\s+((([0-9]{1,})\\s*)+)\\\n"
# Check whether we have a single chronology or more, then assign the variables
# and extract the ring width sequences for individual samples.
if(dim(chronology)[2] == 1){
    chron_numbers <- as.numeric(chronology[1,1])
    chron_years <- as.numeric(chronology[2,1])
    rings <- str_match_all(chronology[3,1],ring_pattern)
} else {
    chron_numbers <- as.numeric(chronology[,1])
    chron_years <- as.numeric(chronology[,2])
    rings <- str_match_all(chronology[,3],ring_pattern)
}


# We end up with a list of matrices containing the group matches from regex.
# We need the 2nd column of each matrix which contains the data. Extract and
# coerce to a vector:
chron_rings <- list()
for (i in 1:length(chron_numbers) ) {
    chron_rings[[i]] <- scan(text = rings[[i]][,2])
}

# We now have three objects to construct the aggregated chronology dataframe:
# chron_numbers: the sample numbers (identifier)
# chron_years: the dated year of the samples (growth year of youngest ring measured)
# chron_rings: The sequence of yearly ring widths of the samples

# construct the sequence of years spanning the data. For the oldest year, we need
# the minimum of dated year - age of tree
oldest <- +Inf # initialize with a ridiculously high number
for(i in 1:length(chron_years)){
    currentsample <- chron_years[i] - (length(chron_rings[[i]])-1)
    oldest <- ifelse(currentsample < oldest, currentsample, oldest)
}
newest <- max(chron_years)
year <- seq(oldest,newest,by=1)


# Result dataframe. Structure: Columns represent one sample each,
# colname is the sample number. Years without tree ring widhts are coded as NA.
# Rownames are the years for which we have data.
rwl_df <- data.frame(row.names = year)

# iterate over the individual sample chronologies, fill space before and after
# with NA, append to rwl_df
for(i in 1:length(chron_numbers)){
    datedyear <- chron_years[i]
    fillbefore <- numeric()
    fillafter <- numeric()
    if(datedyear > oldest){
        fillbefore <- rep(NA,times=(datedyear-(length(chron_rings[[i]])-1)-oldest))
    }
    if(datedyear < newest){
        fillafter <- rep(NA,times=newest-datedyear)
    }
    rings <- chron_rings[[i]]
    series <- as.data.frame(c(fillbefore,rings,fillafter))
    colnames(series) <- chron_numbers[i]
    rwl_df <- cbind(rwl_df,series)
}

# cast to rwl:
rwl_df

if(saveRWL == T){
    dplR::write.rwl(rwl_df, fname = file_rwl_name, format = "tucson")
}

rwl_df
}
