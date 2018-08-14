# after downloading .csv files from Earth Engine:
# 	remove the (empty) .geo column
#	remove the brackets and other data surrounding the coordinates
#		in the centroid column
#	replace the comma separating the coordinates in the centroid
#		column with an "m"
# 	format: -44.123456m-20.123456
# 	put Time column in Number format to show all sig figs
# 		otherwise saving as .csv will save rounded values
#		resulting in multiple observations per polygon in some years
# In the future, don't have any processing steps that involve Excel

library(dplyr)
library(lubridate)

addcoords=function(dat){
	cents=as.character(dat$centroid)
	dat$long=as.vector(lapply(strsplit(cents,'m'),function(x){
			first=x[1]
			first
		}),mode='numeric')
	dat$lat=as.vector(lapply(strsplit(cents,'m'),function(x){
			second=x[2]
			second
		}),mode='numeric')
	dat = dat %>%
            mutate(DATETIME=as.Date.POSIXct(Time/1000, 
			origin="1970-01-01", tz='UTC')) %>%
            mutate(YEAR=year(DATETIME)) %>%
            mutate(MONTH=month(DATETIME)) %>%
            mutate(DOY=yday(DATETIME)) 
	#dat$ones=rep(1)
	dat
}

# just datetime
datetime=function(dat){
	dat = dat %>%
            mutate(DATETIME=as.Date.POSIXct(Time/1000, 
			origin="1970-01-01", tz='UTC')) %>%
            mutate(YEAR=year(DATETIME)) %>%
            mutate(MONTH=month(DATETIME)) %>%
            mutate(DOY=yday(DATETIME)) 
	#dat$ones=rep(1)
	dat
}

#just coords
getcoords=function(dat){
	cents=as.character(dat$centroid)
	dat$long=as.vector(lapply(strsplit(cents,'m'),function(x){
			first=x[1]
			first
		}),mode='numeric')
	dat$lat=as.vector(lapply(strsplit(cents,'m'),function(x){
			second=x[2]
			second
		}),mode='numeric')
	dat
}
