# This script prepares NIRv time series files from Earth Engine
#     for analysis and rotation detection
# Files are named by WRS path/row

setwd("C:\\Users\\Devin\\Documents\\RSvctrl")

library(sp)
library(dplyr)
library(lubridate)
library(zoo)
#library(changepoint)
library(lattice)
#library(reshape2)

source('addcoords.R')
source('rotviz.R')

#props=read.csv('properties_219072_foreuc.csv')
#props=read.csv('properties_218073.csv')
#props=read.csv('properties_218072-3.csv')
#props=read.csv('properties_217074.csv')
#props=read.csv('properties_219074.csv')
props=read.csv('properties_219073_2.csv')
#head(props)
props=props[,-which(names(props) %in% c('system.index','.geo'))]

#eucts=read.csv('LS219072_eucnotGer_11-27-17.csv')
#eucts=read.csv('LS218073_eucTS_1-08-18.csv')
#eucts=read.csv('LS218072_eucTS_2-12-18.csv')
#names(eucts)[which(names(eucts)=='lon')]='long'
#eucts1 = read.csv('LS217074_foundeucTS_region1_3-18-18.csv')
#eucts2 = read.csv('LS217074_foundeucTS_region2_3-18-18.csv')
#eucts3 = read.csv('LS217074_foundeucTS_region3_3-18-18.csv')
#identical(sort(names(eucts1)),sort(names(eucts2)))
#eucts1$.geo=rep(NA) # make sure you're consistent in Excel pre-processing
#nat1=nat1[,-which(names(nat1)=='system.index')]
#eucts=rbind(eucts1,eucts2,eucts3)
#eucts=read.csv('LS219074_eucTS_3-23-18.csv')
eucts=read.csv('LS219073_eucTS_3-4-18.csv')
if(is.element('system.index',names(eucts))==T){
  eucts=eucts[,-which(names(eucts)=='system.index')]
}
if(is.element('.geo',names(eucts))==T){
  eucts=eucts[,-which(names(eucts)=='.geo')]
}

head(eucts)
# 219073--lat was exported as a point, not a number
# paste back on from props--assume longs unique

props2=props
if(is.element("class",names(eucts))){
  props2=props2[,-which(names(props2)=='class')]}
if(is.element("sharedlat",names(eucts))){
  props2=props2[,-which(names(props2)=='sharedlat')]}
if(is.element("sharedlong",names(eucts))){
  props2=props2[,-which(names(props2)=='sharedlong')]}
#props2=getcoords(props2)
# for 218072: 
props2$lat=round(props2$lat,8) # rounding diffs between getcoords and exported long
props2$long=round(props2$long,8)
head(props2)

source('genprep_climprops_TS.R')
#eucts2=merge(eucts,props2,by='centroid')
#eucts$lat=round(eucts$lat,8)
eucts$long=round(eucts$long,8)
eucts=eucts[,-which(names(eucts)=='lat')]
#eucts2=merge(eucts,props2,by=c('lat','long'))
eucts2=merge(eucts,props2,by='long')
head(eucts2)
length(unique(eucts$long))
length(unique(eucts2$lat)) #same, good

#eucts2=genericprep(eucts2)
eucts2=genericprepll(eucts2)


#nat1=read.csv('LS219072_nat6_alldates_11-28-17.csv')
#head(nat1)
#nat2=read.csv('LS219072_nat3_alldates_11-29-17.csv')
#nat3=read.csv('LS219072_nat4_alldates_11-29-17.csv')
#nat4=read.csv('LS219072_nat5_alldates_11-29-17_part1.csv')
#nat5=read.csv('LS219072_nat5_alldates_11-29-17_part2.csv')
as.POSIXct(tail(nat4$Time)/1000,origin='1970-01-01')
# make sure you got the right number of significant digits

#nat1=read.csv('LS218073_nat1_alldates_1-08-18.csv')
#nat2=read.csv('LS218073_nat2_1-08-18.csv')
#nat3=read.csv('LS218073_nat3_alldates_1-08-18.csv')
#nat4=read.csv('LS218073_nat4_alldates_1-08-18.csv')

#nat1=read.csv('LS218072_nat1_alldates_2-12-18.csv')
#nat2=read.csv('LS218072_nat2_alldates_2-13-18.csv')
#nat3=read.csv('LS218072_nat3_alldates_2-26-18.csv')

#nat1=read.csv('LS217074_nat1_alldates_3-21-18.csv')
#nat2=read.csv('LS217074_nat2_alldates_3-21-18.csv')
#nat3=read.csv('LS217074_nat3_alldates_3-21-18.csv')

#nat1=read.csv('LS219074_nat1_alldates_2-28-18.csv')
#nat2=read.csv('LS219074_nat2_alldates_3-1-18.csv')
#nat3=read.csv('LS219074_nat3_alldates_3-2-18.csv')
#nat4=read.csv('LS219074_nat4_alldates_3-2-18.csv')


nat1=read.csv('LS219073_nat1_alldates_3-5-18.csv')
nat2=read.csv('LS219073_nat2_alldates_3-7-18.csv')
nat3=read.csv('LS219073_nat3_alldates_3-10-18.csv')
nat4=read.csv('LS219073_nat4_alldates_3-12-18.csv')
identical(sort(names(nat1)),sort(names(nat2)))
#nat1=nat1[,-which(names(nat1)=='system.index')]

nats=rbind(nat1,nat2,nat3,nat4)#,nat5)
if(is.element('system.index',names(nats))==T){
  nats=nats[,-which(names(nats)=='system.index')]
}
if(is.element('.geo',names(nats))==T){
  nats=nats[,-which(names(nats)=='.geo')]
}
# saveRDS(nats,'allnatobs_LS219072.Rds')
# saveRDS(nats,'allnatobs_LS218073.Rds')
# saveRDS(nats,'allnatobs_LS218072.Rds')
# saveRDS(nats,'allnatobs_LS217074.Rds')
# saveRDS(nats,'allnatobs_LS219074.Rds')
#saveRDS(nats,'allnatobs_LS219073.Rds')
#nats=readRDS('allnatobs_LS218072.Rds')
#nats=readRDS('allnatobs_LS217074.Rds')
#nats2=merge(nats,props2,by='centroid')
nats$lat=round(nats$lat,8)
nats$long=round(nats$long,8)
#nats2=merge(nats,props2,by='centroid')
# nats lat got rounded--is this a opening-in-Excel problem?
# I think so
nats2=merge(nats,props2,by=c('lat','long'))
#length(nats$NDVI)
#nats2=natprep(nats2)
nats2=natprepll(nats2)

nats2$sharedlat=round(nats2$sharedlat,8)
nats2$sharedlong=round(nats2$sharedlong,8)
eucts2$sharedlat=round(eucts2$sharedlat,8)
eucts2$sharedlong=round(eucts2$sharedlong,8)

source('TSprep_function_2-9-18.R')
eux=TSprep(eucts2,nats2)

source('tsviz.R')
tsvizn(eux)
length(unique(eux$ID)) 
length(unique(eux$ID[eux$ndates>70])) 

plot(lat~long,data=distinct(eux,ID,.keep_all=T))
points(lat~long,data=distinct(eux[eux$ndates<70,],ID,.keep_all=T),
       col=2)
points(lat~long,data=distinct(nats,lat,.keep_all=T),
       col=2) 
# the ones on the edges of the region (for 218073)
# Skip those, I guess
# all ok for 218072

eux2=eux[eux$ndates<100 & eux$ndates>70 & eux$nadj>=eux$ndates-5,]
eux2=ungroup(eux2)
length(unique(eux2$ID)) # 614 of the 674 stands have natveg
# I think the remaining ones overlap cen stands?
# 1006/1075 for 219074
# 1381/1551 for 219073


source('simple_rotationfinder_feb18.R')
foundrots=findrotations(eux2)

newrotvizadjt(foundrots)
#saveRDS(foundrots,'TS_rotations_LS218073.Rds')
#saveRDS(foundrots,'TS_rotations_LS219072.Rds')
#saveRDS(foundrots,'TS_rotations_LS218072.Rds')
#saveRDS(foundrots,'TS_rotations_LS217074.Rds')
length(unique(foundrots$ID[foundrots$badrots>0]))

foundrots=group_by(foundrots,ID) %>% mutate(maxrots=max(ezrot))
#saveRDS(foundrots,'TS_rotations_LS219074.Rds')
#saveRDS(foundrots,'TS_rotations_LS219073.Rds')


