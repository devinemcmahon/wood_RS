# This script matches NIRv observations from Eucalyptus stands
#   with the corresponding native vegetation NIRv observations

TSprep=function(eucTS,natTS){
  natTS$pctok=natTS$clearpix*(30*30*10^-2)/natTS$areaHa
  # Proportion of native vegetation area with pixels flagged as not cloudy
  natTS=natTS[natTS$areaHa>0,]
  natTS$NIRv[natTS$NIRv<0]=NA # do this before calculating stats on NIRv
  
  natsum=group_by(natTS,ID,YEAR) %>% summarise(yrVI=mean(NIRv,na.rm=T),
                                               yrobs=n())
  natsum=group_by(natsum,ID) %>% 
    mutate(sdyrVI=sd(yrVI),meanyrVI=mean(yrVI[yrobs>5]),normVI=yrVI/meanyrVI, 
           maxnorm=max(normVI[yrobs>5]),minnorm=min(normVI[yrobs>5]))
  # avoid high values due to years (especially at end of time series)
  #   with just a few observations
  oknats=unique(natsum$ID[natsum$maxnorm<1.5 & natsum$minnorm>0.5 ])
  length(unique(natTS$ID))
  length(oknats) # 97% of original stands for Cenibra (91% for Cen in 217074)
  
  nv=natTS[which(natTS$ID %in% oknats),]
  
  nv$NIRv[nv$NIRv<0]=NA
  # throw out extremely high NIRv values (top 0.01%)
  nv4sig=quantile(nv$NIRv,0.9999,na.rm=T)
  nv$NIRv[nv$NIRv>nv4sig]=NA
  
  eucTS$clearpix=rep(NA)
  eucTS$pctok=rep(100) # all pixels represented in composites
  alreadyinboth=names(nv)[which(names(nv) %in% names(eucTS))]
  justnv=names(nv)[-which(names(nv) %in% names(eucTS))] 
  eucTS2=eucTS
  if(length(justnv)>0){eucTS2[,((ncol(eucTS)+1):(ncol(eucTS)+length(justnv)))]=rep(NA)
  names(eucTS2)=c(names(eucTS),justnv)}
  justeuc=names(eucTS2)[-which(names(eucTS2) %in% names(nv))]
  nv2=nv
  nv2[,((ncol(nv)+1):(ncol(nv)+length(justeuc)))]=rep(NA)
  names(nv2)=c(names(nv),justeuc)
  
  eucTS2=ungroup(eucTS2)
  eucTS2$ID=as.numeric(as.character(eucTS2$ID))
  all=rbind(eucTS2,nv2)
  
  pctgr=all[all$ID %in% oknats | all$ID %in% eucTS$ID,] %>% 
    group_by(sharedlat,sharedlong,soil,elevbin,DATETIME) %>% 
    mutate(datepix=sum(clearpix[class==2],na.rm=T),
           natVI=sum(NIRv[class==2]*clearpix[class==2],na.rm=T)/datepix,
           adjVI=NIRv-natVI)
  # Each eucalyptus NIRv observation is compared to the mean of all
  #   native vegetation patches within the same cluster of stands
  #   and with the same elevation, soil type, and date of observation
  # Weight native vegetation NIRv by number of clear pixels 
  #   in that patch on that date
  # class == 2 for native vegetation; class == 1, 8, 9, or 10 for eucalyptus
  
  pctgr=group_by(pctgr,sharedlat,sharedlong,soil,elevbin) %>% 
    mutate(natarea=sum(areaHa[class==2 & duplicated(ID)==F]))
  
  pctgr=ungroup(pctgr)
  eux = distinct(pctgr[pctgr$class!=2,],ID,DATETIME,.keep_all=T)
  eux=group_by(eux,ID) %>% mutate(nadj=sum(!is.na(adjVI)),ndates=n(),
                                  pctadj=nadj*100/ndates)
  return(eux)
}