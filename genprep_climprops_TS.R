# This function rearranges data exported from Earth Engine
#   so that the 50th, 75th, and 90th percentile NIRv values from each
#   year form a single "NIRv" column
# Note that MAP, MAT, precip_wetmo, and precip_drymo are from WorldClim
#   and are replaced by averaging the values from CHIRPS and NCEP
#   over the period of analysis for each area
#   in the actual analysis; precip is a less-useful 15-year CHIRPs average
# sat indicates which Landsat took the observation (averaged for euc mosaics)
source('addcoords.R')

genericprep=function(pct){
columnnames50=c('areaHa','aspect','centroid','elev','precip',
	'slope','soil','class','NIRv_p50','NIRv_50','Time_50',
	'NDVI_50','Path_50','Row_50','sat_50',
	'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')
columnnames75=c('areaHa','aspect','centroid','elev','precip',
	'slope','soil','class','NIRv_p75','NIRv_75','Time_75',
	'NDVI_75','Path_75','Row_75','sat_75',
	'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')
columnnames90=c('areaHa','aspect','centroid','elev','precip',
	'slope','soil','class','NIRv_p90','NIRv_90','Time_90',
	'NDVI_90','Path_90','Row_90','sat_90',
	'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')

pct50=pct[,which(is.element(names(pct),columnnames50))]
pct50$percentile=rep(50)
pct50=pct50[,sort(names(pct50))]
#head(pct50)
names(pct50)=c('areaHa','aspect','centroid','class','elev',#'ID',
  'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
	'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
	'slope','soil','Time')
pct75=pct[,which(is.element(names(pct),columnnames75))]
pct75$percentile=rep(75)
pct75=pct75[,sort(names(pct75))]
names(pct75)=c('areaHa','aspect','centroid','class','elev',#'ID',
  'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
  'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
  'slope','soil','Time')
pct90=pct[,which(is.element(names(pct),columnnames90))]
pct90$percentile=rep(90)
pct90=pct90[,sort(names(pct90))]
names(pct90)=c('areaHa','aspect','centroid','class','elev',#'ID',
  'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
  'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
  'slope','soil','Time')
pct2=rbind(pct50,pct75,pct90)
pct2=addcoords(pct2)
pct2$ID=paste(pct2$class,round(abs(pct2$lat),3)*1000,
	round(abs(pct2$long),3)*1000,sep='')
pct2$ID=as.numeric(pct2$ID)
pct2$soilclass=as.vector(lapply(strsplit(
		as.character(pct2$soil),'[[:lower:]]'),
	function(x){
		first=x[1]
		first
	}),mode='character')
pct2$soilclass=factor(pct2$soilclass,levels=sort(unique(pct2$soilclass)))
pct2$north=cos(pct2$aspect*pi/180)
pct2$elevbin=round(pct2$elev,-2) # -1 previously: too precise
pct2$sat2=round(pct2$sat,0)
pct2$sat2[pct2$sat2==6 & pct2$sat<6]=5
pct2$sat2[pct2$sat2==6 & pct2$sat>6]=7
pct2
}

genericprepll=function(pct){ # ll for lat/long
  columnnames50=c('areaHa','aspect','lat','long','elev','precip',
                  'slope','soil','class','NIRv_p50','NIRv_50','Time_50',
                  'NDVI_50','Path_50','Row_50','sat_50',
                  'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')
  columnnames75=c('areaHa','aspect','lat','long','elev','precip',
                  'slope','soil','class','NIRv_p75','NIRv_75','Time_75',
                  'NDVI_75','Path_75','Row_75','sat_75',
                  'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')
  columnnames90=c('areaHa','aspect','lat','long','elev','precip',
                  'slope','soil','class','NIRv_p90','NIRv_90','Time_90',
                  'NDVI_90','Path_90','Row_90','sat_90',
                  'MAP','MAT','precip_wetmo','precip_drymo','sharedlat','sharedlong')
  
  pct50=pct[,which(is.element(names(pct),columnnames50))]
  pct50$percentile=rep(50)
  pct50=pct50[,sort(names(pct50))]
  #head(pct50)
  names(pct50)=c('areaHa','aspect','class','elev',#'ID',
                 'lat','long',
                 'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
                 'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
                 'slope','soil','Time')
  pct75=pct[,which(is.element(names(pct),columnnames75))]
  pct75$percentile=rep(75)
  pct75=pct75[,sort(names(pct75))]
  names(pct75)=c('areaHa','aspect','class','elev',#'ID',
                 'lat','long',
                 'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
                 'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
                 'slope','soil','Time')
  pct90=pct[,which(is.element(names(pct),columnnames90))]
  pct90$percentile=rep(90)
  pct90=pct90[,sort(names(pct90))]
  names(pct90)=c('areaHa','aspect','class','elev',#'ID',
                 'lat','long',
                 'MAP','MAT','NDVI','NIRv','NIRv_p','Path','percentile','precip',
                 'precip_drymo','precip_wetmo','Row','sat','sharedlat','sharedlong',
                 'slope','soil','Time')
  pct2=rbind(pct50,pct75,pct90)
  pct2=datetime(pct2)
  pct2$ID=paste(pct2$class,round(abs(pct2$lat),3)*1000,
                round(abs(pct2$long),3)*1000,sep='')
  pct2$ID=as.numeric(pct2$ID)
  pct2$soilclass=as.vector(lapply(strsplit(
    as.character(pct2$soil),'[[:lower:]]'),
    function(x){
      first=x[1]
      first
    }),mode='character')
  pct2$soilclass=factor(pct2$soilclass,levels=sort(unique(pct2$soilclass)))
  pct2$north=cos(pct2$aspect*pi/180)
  pct2$elevbin=round(pct2$elev,-2) # -1 previously: too precise
  pct2$sat2=round(pct2$sat,0)
  pct2$sat2[pct2$sat2==6 & pct2$sat<6]=5
  pct2$sat2[pct2$sat2==6 & pct2$sat>6]=7
  pct2
}

natprep=function(pct2){
  pct2=addcoords(pct2)
  pct2$ID=paste(pct2$class,round(abs(pct2$lat),3)*1000,
                round(abs(pct2$long),3)*1000,sep='')
  pct2$ID=as.numeric(pct2$ID)
  pct2$soilclass=as.vector(lapply(strsplit(
    as.character(pct2$soil),'[[:lower:]]'),
    function(x){
      first=x[1]
      first
    }),mode='character')
  pct2$soilclass=factor(pct2$soilclass,levels=sort(unique(pct2$soilclass)))
  pct2$north=cos(pct2$aspect*pi/180)
  pct2$elevbin=round(pct2$elev,-2) 
  pct2$sat2=round(pct2$sat,0)
  pct2$sat2[pct2$sat2==6 & pct2$sat<6]=5
  pct2$sat2[pct2$sat2==6 & pct2$sat>6]=7
  pct2 
}
natprepll=function(pct2){
  pct2=datetime(pct2)
  pct2$ID=paste(pct2$class,round(abs(pct2$lat),3)*1000,
                round(abs(pct2$long),3)*1000,sep='')
  pct2$ID=as.numeric(pct2$ID)
  pct2$soilclass=as.vector(lapply(strsplit(
    as.character(pct2$soil),'[[:lower:]]'),
    function(x){
      first=x[1]
      first
    }),mode='character')
  pct2$soilclass=factor(pct2$soilclass,levels=sort(unique(pct2$soilclass)))
  pct2$north=cos(pct2$aspect*pi/180)
  pct2$elevbin=round(pct2$elev,-2) 
  pct2$sat2=round(pct2$sat,0)
  pct2$sat2[pct2$sat2==6 & pct2$sat<6]=5
  pct2$sat2[pct2$sat2==6 & pct2$sat>6]=7
  pct2 
}