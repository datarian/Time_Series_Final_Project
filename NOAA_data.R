# Time series project: Dendro data 

library(treeclim)
library(dplR) # package containing the Arizona dendro dataset
library(tseries)

# Literature
# [1] https://cran.r-project.org/web/packages/dplR/dplR.pdf
# [2] http://waldwachstum.wzw.tum.de/fileadmin/publications/2012_Pretzsch_dvffa.pdf
# [3] https://www.rdocumentation.org/packages/treeclim/versions/2.0.0
# https://www.ncdc.noaa.gov/cag/statewide/time-series/2/pcp/1/2/1895-2018?base_prd=true&firstbaseyear=1901&lastbaseyear=2000 (US weather data: Temp.,
# precipitation, drought index and many more)
# [4] https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1583&context=jmasm (Variance stabilizing power transformation for time series)

setwd("/Volumes/Samsung_T3/Master_Statistics_UNINE/Semester_4/Time_Series/Time_series_project2")

# wd <- getwd()
# 
# # NOAA-Daten
# path <- paste0(wd,"/Data")
# 
# file <- paste0(path,"/Average_Temperature")
# 
# filename <- paste0(file,"/2-tavg-1-1-1895-2018.csv")
# 
# data <- read.csv(file=filename, header=T, sep=",")


# -----------------------

wd <- getwd()

# NOAA-Daten
path <- paste0(wd,"/data/NOAA")


# List to store all the data
data_list <- list()

folderNames <- c("Average_Temperature", "Cooling_Degree_Days", "Heating_Degree_Days", "Maximum_Temperature", 
               "Minimum_Temperature", "Palmer_Drought_Severity_Index", "Palmer_Hydrological_Drought_Index",
               "Palmer_Modified_Drought_Index", "Palmer_Z-Index", "Percipation")

for (name in folderNames){
  pathToFolder <- paste0(paste0(path, "/"), name)
  print(pathToFolder)
  #print(list.files(pathToFile))
  #print("-------")
  
  files <- list.files(pathToFolder)
  
  # data frame to store the data
  numberOfROws <- 2017-1895+1
  
  df <- data.frame(rep(0,numberOfROws)) # data.frame with 1 columns and 125 rows
  
  for (file in files){
    path_tmp <- paste0(paste0(pathToFolder, "/"), file)
    print(path_tmp)
    
    data <- read.csv(path_tmp, header=F, sep=",")
    
    # Exlude text
    data <- data[-c(1,2,3,4,5),]
    
    # Data frame has numberOfRows rows
    colnames(data) <- c("Date", "Value", "Anomaly")
    data <- data[1:numberOfROws,]
    print(data)
    df <- cbind(df, data)
  }
  print(df)
  # change the order from 1, 10, 11, 12, 4, ..., 9 to 1, 2, ..., 12
  # throw away first column
  data_list[[file]] <- cbind(df[,2:4], df[,14:(13+(8*3))], df[,5:13])
}
  

data_list[[1]]

# Date Value Anomaly   Date Value Anomaly Date.1 Value.1 Anomaly.1 Date.2 Value.2 Anomaly.2 Date.3 Value.3 Anomaly.3 Date.4 Value.4 Anomaly.4 Date.5 Value.5 Anomaly.5 Date.6 Value.6 Anomaly.6 Date.7 Value.7 Anomaly.7
# 6   189501  40.4    -0.6 189502  43.5    -1.4 189503    50.0       0.0 189504    58.3       1.3 189505    65.8       0.5 189506    71.6      -3.3 189507    78.3      -1.6 189508    77.8      -0.2 189509    71.5      -0.7
# 7   189601  44.6     3.6 189602  45.3     0.4 189603    51.3       1.3 189604    55.2      -1.8 189605    65.8       0.5 189606    77.9       3.0 189607    78.5      -1.4 189608    77.3      -0.7 189609    71.5      -0.7



# create one single row for each variable
n <- (2017-1895+1)
n12 = 12*n
data_df <- data.frame(rep(0,n12))

for (i in 1:length(data_list)){
  # data frame to store the numeric values
  #data_df_tmp <- data.frame(rep(0,n))
  
  # data
  data <- data_list[[i]]
  
  # Extract the Columns with name 'Value'  
  data_values <- data[,seq(from=2, to=12*3, by=3)]
  print(data_values)
  
  # matrix
  data_mat <- matrix(0, ncol=ncol(data_values), nrow=nrow(data_values))
  
  for (col in 1:ncol(data_values)){
    #data_df_tmp <- cbind(data_df_tmp, as.numeric(as.character(data_values[,col])))
    data_mat[,col] <- as.numeric(as.character(data_values[,col]))
  }
  
  # print(data_mat)
  
  # Vectorize the data
  data_vec <- c(t(data_mat))
  # print(data_vec)
  data_df <- cbind(data_df, data_vec)
}

data_df <- data_df[,-1]
colnames(data_df) <- folderNames
data_df

# Average_Temperature Cooling_Degree_Days Heating_Degree_Days Maximum_Temperature Minimum_Temperature Palmer_Drought_Severity_Index Palmer_Hydrological_Drought_Index Palmer_Modified_Drought_Index Palmer_Z-Index
# 1                   40.4                   0                 508                49.0                31.8                         -0.79                              0.87                         -0.36          -1.41

data_df_ts <- ts(data_df, start=c(1895,1), freq=12)
data_df_ts

# Average_Temperature Cooling_Degree_Days Heating_Degree_Days Maximum_Temperature Minimum_Temperature Palmer_Drought_Severity_Index Palmer_Hydrological_Drought_Index Palmer_Modified_Drought_Index Palmer_Z-Index
# Jan 1895                40.4                   0                 508                49.0                31.8                         -0.79                              0.87                         -0.36     

########################################################################################
########################################################################################

# Dendro data - commented out

# # data [1], p.37: This data set includes ring-width measurements for ponderosa pine (Pinus ponderosa) increment cores collected at the Gus Pearson Natural Area (GPNA) 
# # in northern Arizona, USA. There are 58 series from 29 trees (2 cores per tree). Data are further described by Biondi and Qeadan (2008) and references therein.
# data(gp.rwl)
# 
# str(gp.rwl)
# head(gp.rwl)
# tail(gp.rwl) # 1990
# data <- as.matrix(gp.rwl)
# 
# start =  nrow(data) - (1990-1895) #421
# 
# data_sub_mat <- data[start:nrow(data), ]
# 
# # take the means
# means <- rowMeans(data_sub_mat, na.rm=T)
# 
# # plot the means
# plot(y=means, x=1:length(means), type="l")
# 
# # differenced time sereis
# data_diff <- diff(means)
# 
# plot(y=data_diff, x=1:length(data_diff), type="l")
# 
# # H0: stationary cannot be rejected. But data 
# kpss.test(data_diff)
# 
# # is this MA(2)?
# acf(data_diff)
# pacf(data_diff)
# 
# # Artificial MA(2) process, quite similar to dendro process
# set.seed(999)
# noise <- rnorm(400, 1, 1)
# n <- length(noise)
# x <- noise[3:n] - 0.5*noise[2:(n-1)] - 0.2*noise[1:(n-2)]
# plot(x, type="l")
# acf(x)
# pacf(x)
# 
# # -------------------------------------------------------------------
# 
# # Nice :))))))
# data(anos1)
# anos1
# str(anos1)
# 
# yrs <- as.numeric(rownames(anos1))
# # Plot rings with data of two radii from same individual tree
# res <- plotRings(yrs, anos1[,4], trwW = anos1[,5],
#                  species.name = "Cedrela odorata")
# 

