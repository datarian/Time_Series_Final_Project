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

#setwd("/Volumes/Samsung_T3/Master_Statistics_UNINE/Semester_4/Time_Series/Time_series_project2")

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

# --------------------


create_data_set <- function(data, toChoose, aggregationMethod){

  # INPUT:
  # - data: data.frame with monthly data
  # - toChoose: matrix, indicating, which values to choose: 1: Jan, 2: Feb, .... The rest of the column is filled with zeros.
  # - aggregationMethod: "sum", "mean"

  # OUTPUT:
  # aggregated yearly data

  year <- floor(time(data)) # extract the years
  data_cleaned <- cbind(data,year) # year is stored in the last column
  colnames(data_cleaned)  <- c(colnames(data), "Year")

  n <- floor(nrow(data)/12)*12
  data_cleaned <- as.data.frame(data_cleaned[1:n,]) # data set starts in Jan. With that it ends in Dez

  data_new <- data.frame(Year=unique(data_cleaned$Year)) # this data frame will be filled with values in the next for loop

  for (col in 1:ncol(data)){
  data_sub <- data_cleaned[toChoose[,col]!=0,c(col, ncol(data_cleaned))] # year is stored in the last column
  data_sub <- as.data.frame(data_sub)
  colnames(data_sub) <- colnames(data_cleaned)[c(col,ncol(data_cleaned))]

    if (aggregationMethod[col]=="mean"){
      agg_values <- aggregate(data_sub, by=list(data_sub$Year), FUN=mean)
    } else if (aggregationMethod[col]=="sum"){
      agg_values <- aggregate(data_sub, by=list(data_sub$Year), FUN=sum)
    }
    data_new <- cbind(data_new, agg_values[,2])
  }

  colnames(data_new) <- c("Year", colnames(data))
  return(data_new)

}

# --------------------

data_df_ts_sub <- data_df_ts[,1:6]

toChoose <- matrix(c(1:12, # Average Temperature
                      1:12, # Cooling Degree Days
                      1:12, # Heating Degree Days
                      1:12, # Maximum Temperature
                      1:12, # Minimum Temperature
                      1:12), # Palmer Drought Severity Index
                    nrow=12, ncol=ncol(data_df_ts_sub), byrow=F)

aggregationMethod <- c("mean", "mean", "mean", "mean", "mean", "mean")


# -------------

toChoose2 <- matrix(c(c(0,0,0,0,5:12), # Average Temperature
                     1:12, # Cooling Degree Days
                     1:12, # Heating Degree Days
                     1:12, # Maximum Temperature
                     1:12, # Minimum Temperature
                     1:12), # Palmer Drought Severity Index
                   nrow=12, ncol=ncol(data_df_ts_sub), byrow=F)

aggregationMethod2 <- c("mean", "mean", "mean", "mean", "mean", "sum")


data_aggr <- create_data_set(data_df_ts_sub, toChoose, aggregationMethod)
# Year Average_Temperature Cooling_Degree_Days Heating_Degree_Days Maximum_Temperature Minimum_Temperature Palmer_Drought_Severity_Index
# 1   1895            58.49167            206.2500            192.3333            71.41667            45.55000                   -0.18416667
# 2   1896            59.98333            215.5000            161.9167            72.62500            47.38333                    0.21833333


# ARIMA with exogeneous regressors

# Dendro data end in 1990
n <- 1990-1895 +1
data_aggr_sub <-  data_aggr[1:n,]
class(data_aggr_sub)
data_aggr_sub # fine

# Dendro
data(gp.rwl)
str(gp.rwl)
head(gp.rwl)
tail(gp.rwl) # 1990
data <- as.matrix(gp.rwl)

# starts im 1895
start <- 1895 - 1570 + 1
end <- 1990 -1570 +1
data_dendro_sub <- apply(data[start:end,],1,function(x) {mean(x, na.rm=T)})
# 1895      1896      1897      1898      1899      1900      1901      1902      1903      1904      1905      1906      1907
# 1.4122807 1.1721053 1.2643860 1.4264912 0.9003448 0.7215517 1.2863793 0.8929310 1.3768966 0.9870690 1.3189655 1.4232759 1.7353448
# 1908      1909      1910      1911      1912      1913      1914      1915      1916      1917      1918      1919      1920
# 2.0539655 1.7146552 1.4084483 1.5632759 1.2960345 1.0543103 1.3579310 1.1782759 1.5205172 1.3132759 0.9708621 1.4182759 1.1510345
# 1921      1922      1923      1924      1925      1926      1927      1928      1929      1930      1931      1932      1933
# 1.3660345 1.0737931 1.3918966 1.1610345 1.2122414 1.1486207 1.2131034 1.1927586 1.1812069 1.1465517 1.1362069 0.9908621 0.9063793
# 1934      1935      1936      1937      1938      1939      1940      1941      1942      1943      1944      1945      1946
# 0.8632759 0.9615517 0.9091379 1.0048276 0.8444828 0.7198276 0.8194828 0.9722414 1.0070690 0.9248276 0.7986207 0.7362069 0.7339655
# 1947      1948      1949      1950      1951      1952      1953      1954      1955      1956      1957      1958      1959
# 0.5487931 0.4372414 0.8994828 0.6348276 0.3446552 0.8746552 0.7518966 0.7272414 0.7398276 0.4224138 0.6770690 0.6286207 0.7568966
# 1960      1961      1962      1963      1964      1965      1966      1967      1968      1969      1970      1971      1972
# 0.5710345 0.5670690 0.6517241 0.2501724 0.5570690 0.6605172 0.6120690 0.6267241 0.5086207 0.5301724 0.4825862 0.2775862 0.5306897
# 1973      1974      1975      1976      1977      1978      1979      1980      1981      1982      1983      1984      1985
# 0.5129310 0.2670690 0.6039655 0.4917241 0.1751724 0.3793103 0.4131034 0.4343103 0.3396552 0.4151724 0.3553448 0.3060345 0.3979310
# 1986      1987      1988      1989      1990
# 0.4175862 0.4189655 0.4664912 0.2777193 0.5433333

plot(y=data_dendro_sub, x=1:length(data_dendro_sub), type="l")


model <- arima(data_dendro_sub, order=c(1,0,0), xreg=data_aggr_sub,
      seasonal=list(order=c(0,0,0), period=NA))

# Call:
#   arima(x = data_dendro_sub, order = c(1, 0, 0), seasonal = list(order = c(0,
#                                                                            0, 0), period = NA), xreg = data_aggr_sub)
#
# Coefficients:
#          ar1  intercept     Year  Average_Temperature  Cooling_Degree_Days  Heating_Degree_Days  Maximum_Temperature  Minimum_Temperature
#       0.3631    48.4113  -0.0134              -0.0224               0.0098              -0.0116              -0.1899              -0.1341
# s.e.  0.1104     8.6395   0.0013               0.0361               0.0040               0.0043               0.0608               0.0679
# Palmer_Drought_Severity_Index
#                             -0.0132
# s.e.                         0.0103
#
# sigma^2 estimated as 0.03021:  log likelihood = 31.69,  aic = -43.38

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

# aggregate data frame mtcars by cyl and vs, returning means
# for numeric variables

