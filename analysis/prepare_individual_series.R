####
#### Individual series
library(readtext)
library(stringr)

file_contents <- readtext("data/900+.txt")

# Extract the individual data series from the unstructured text file
pattern <- "\\sNr\\. ([0-9]*)\\.0 Dat: ([0-9]{3,4}).*\\n.*\\n.*\\n{1,2}((?:[ ]{1,3}[0-9 -]*\\n){3,})"
parsed <- str_match_all(file_contents,pattern)
rm(file_contents) # Free up memory

# Extract the sample numbers and their corresponding dated year
chronology <- as.data.frame(parsed[[1]][,-1],stringsAsFactors = F)
rm(parsed) # Free up memory
chron_numbers <- as.numeric(chronology[,1])
chron_years <- as.numeric(chronology[,2])

# Extract the ring width sequences for individual samples.
pat <- "\\s{1,3}[0-9 -]{8,11}\\s+((([0-9]{1,})\\s*)+)\\\n"
rings <- str_match_all(chronology[,3],pat)

# We end up with a list of matrices containing the group matches from regex.
# We need the 2nd column of each matrix which contains the data. Extract and
# coerce to a vector:
chron_rings <- list()
for (i in 1:nrow(chronology)) {
    chron_rings[[i]] <- scan(text = rings[[i]][,2])
}

# We now have three objects to construct the aggregated chronology dataframe:
# chron_numbers: the sample numbers
# chron_years: the dated years of the sample
# chron_rings: The sequence of yearly ring widths of the sample

# construct the sequence of years spanning the data. For the oldest year, we need
# the minimum of dated year - age of tree
oldest <- 10000
for(i in 1:length(chron_years)){
    currentsample <- chron_years[i]-(length(chron_rings[[i]])-1)
    oldest <- ifelse(currentsample < oldestring,currentsample,oldestring)
}

newest <- max(chron_years)
year <- seq(oldest,newest,by=1)

# Result dataframe. Structure: First column is the year, subsequent columns
# represent one sample each, colname is the sample number. Years without
# tree ring widhts are coded as NA
rwl_df <- as.data.frame(year); rownames(rwl_df) <- year

for(i in 1:length(chron_numbers)){
    datedyear <- chron_years[i]
    fillbefore <- numeric()
    fillafter <- numeric()
    if(datedyear > oldest){
        fillbefore <- rep(NA,times=(datedyear-(length(chron_rings[[i]])-1)-oldest))
    }
    if(datedyear<newest){
        fillafter <- rep(NA,times=newest-datedyear)
    }
    rings <- chron_rings[[i]]
    series <- as.data.frame(c(fillbefore,rings,fillafter))
    colnames(series) <- chron_numbers[i]
    rwl_df <- cbind(rwl_df,series)
}

# remove year column:
rwl_df <- rwl_df[,-1]

# Some prior analysis:
dplR::rwl.report(rwl_df)
dplR::summary.rwl(rwl_df)

rwi <- detrend(rwl_df,make.plot = T,method = "ModNegExp")

plot(rwl_df,plot.type="spag")
