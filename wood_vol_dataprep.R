
#library(sp)
library(dplyr)
library(lubridate)
library(zoo)
#library(changepoint)
library(lattice)

# Prep wood volume data for analysis 
togroup1=read.csv('wood_vols_rotations.csv')
length(unique(togroup1$ID)) # 707

# Split the stands with known biomass into test and training sets
# stratify sampling within projects (management units of 1-107 stands; usu < 20)
#set.seed(10)
#set.seed(150) 
set.seed(500)
#set.seed(3)
testIDs=group_by(togroup1,project) %>% sample_n(2,replace=T)
length(unique(testIDs$ID)) # 115
names(togroup1)
 
togroup1=togroup1[,order(names(togroup1))] # make column name order replicable
# in case matrix multiplication with coefficients is necessary later (it isn't)

### VI-only version:
# variables: wood volume, stand ID, rotation length (months) and transformations,
# peak age = age at which NIRv within 200 units of max NIRv occurs
# mean and median of pre- and post-peak NIRv
# mean and median of NIRv up to age 2
# mnat = mean NIRv of native vegetation over all observed dates
togroup=togroup1[,which(names(togroup1) %in% c(
  'Vol', 'ID','rotlen','rotlenrt','rotlenln','peakage','maxVI','meanVI','medVI',
  'prediff','prepeakmean','postpeakmean','prepeakmed','postpeakmed',
  'meanage2','medage2', #'rotlensq',
  'slope','north','MAP','MAT','precip_wetmo','precip_drymo','elev','mnat'))]

### VI-only version with no environmental data:
ntogroup=togroup1[,which(names(togroup1) %in% c(
  'Vol', 'ID','rotlen','rotlenrt','rotlenln','peakage','maxVI','meanVI','medVI',
  'prediff','prepeakmean','postpeakmean',#'rotlensq',
  'prepeakmed','postpeakmed','meanage2','medage2'))]
# don't allow rotation length squared term because it doesn't make biological sense 
#   for longer rotations

### adjVI-only version:
# adjVI = difference between eucalyptus and native NIRv on a given date
atogroup=togroup1[,which(names(togroup1) %in% c(
  'Vol','ID','rotlen','rotlenrt','rotlenln','adjpeakage','maxadj','meanadj','medadj',
  'preadiff','preamean','postamean','preamed','postamed','prexmed','prexmean',
  'meanaage2','medaage2','meannat','mednat',#'rotlensq',
  'slope','north','MAP','MAT','precip_wetmo','precip_drymo','elev','mnat'))]

# retain only the variables to be used in the regressions
traindat=togroup[-which(togroup$ID %in% testIDs$ID),]
testdat=togroup[which(togroup$ID %in% testIDs$ID),]

traindat3=traindat[,-which(names(traindat)=='ID')]
testdat3=testdat[,-which(names(testdat)=='ID')]

atraindat=atogroup[-which(atogroup$ID %in% testIDs$ID),]
atestdat=atogroup[which(atogroup$ID %in% testIDs$ID),]

atraindat3=atraindat[,-which(names(atraindat)=='ID')]
atestdat3=atestdat[,-which(names(atestdat)=='ID')]

ntraindat=ntogroup[-which(ntogroup$ID %in% testIDs$ID),]
ntestdat=ntogroup[which(ntogroup$ID %in% testIDs$ID),]

ntraindat3=ntraindat[,-which(names(ntraindat)=='ID')]
ntestdat3=ntestdat[,-which(names(ntestdat)=='ID')]

alldat3=rbind(traindat3,testdat3)
aalldat3=rbind(atraindat3,atestdat3)
nalldat3=rbind(ntraindat3,ntestdat3)

#### Scaling doesn't affect the results; skip it
