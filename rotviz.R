# This script includes several functions for visualizing 
#   NIRv time series and rotation detection


# Deprecated; applied to previous versions of rotation finder
##########
rotviz=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$NIRv, order.by=sub$DATETIME)
    plot(ts,ylim=c(500,4000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    points(NIRv~DATETIME,data=sub,col=rotation)
    abline(v=sub$lastplant,lty=3)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotviz2=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$NIRv, order.by=sub$DATETIME)
    plot(ts,ylim=c(500,4000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    points(NIRv~DATETIME,data=sub,col=rotation)
    lines(natVI~DATETIME,data=sub,col='darkorange')
    abline(v=sub$lastplant,lty=3)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotvizadj=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    points(adjVI~DATETIME,data=sub,col=rotation)
    abline(h=0)
    abline(h=-200,lty=2)
    abline(h=200,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotvizadj2=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    abline(v=sub$DATETIME[sub$start==1],col=3)
    abline(v=sub$DATETIME[sub$end==1],col=2)
    abline(v=sub$DATETIME[sub$putstart==1],col=3,lty=2)
    abline(v=sub$DATETIME[sub$putend==1],col=2,lty=2)
    #abline(v=sub$DATETIME[sub$type2==1],col=4,lty=2)
    points(adjVI~DATETIME,data=sub,col=rotation,pch=5)
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0)
    abline(h=-300,lty=2)
    abline(h=300,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}
##########

# Applicable to new rotation finder used in article
# ezrot is the name of the column designating into which rotation
#   a particular observation falls
# This replaces the "rotation" column from the old (less "easy") 
#   rotation finding algorithm
newrotvizadj=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    abline(v=sub$DATETIME[sub$start==1],col=3)
    abline(v=sub$DATETIME[sub$end==1],col=2)
    abline(v=sub$DATETIME[sub$putstart==1],col=3,lty=2)
    abline(v=sub$DATETIME[sub$putend==1],col=2,lty=2)
    abline(v=sub$DATETIME[sub$latestart==1],col=5)
    #points(adjVI~DATETIME,data=sub,col=rotation,pch=5)
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0)
    abline(h=-300,lty=2)
    abline(h=300,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}
newrotvizadjt=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    abline(v=sub$DATETIME[sub$start==1],col=3)
    abline(v=sub$DATETIME[sub$end==1],col=2)
    abline(v=sub$DATETIME[sub$putstart==1],col=3,lty=2)
    abline(v=sub$DATETIME[sub$putend==1],col=2,lty=2)
    abline(v=sub$DATETIME[sub$latestart==1],col=5)
    #points(adjVI~DATETIME,data=sub,col=rotation,pch=5)
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0)
    #abline(h=-300,lty=2)
    #abline(h=300,lty=2)
    abline(h=-1*sub$thresh,lty=2)
    abline(h=sub$thresh,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

newrotviz=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$NIRv, order.by=sub$DATETIME)
    plot(ts,ylim=c(500,4000),main=unique(sub$ID))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    points(NIRv~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    lines(natVI~DATETIME,data=sub,col='orange')
    abline(h=0)
    abline(h=-300,lty=2)
    abline(h=300,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotvizadjbio=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),#main=unique(sub$ID))
         main=paste(unique(sub$project),
                    unique(sub$date_pl),unique(sub$age)))
    #	abline(v=sub$peakdate[!is.na(sub$rotation)],lty=2)
    abline(v=sub$DATETIME[sub$start==1],col=3)
    abline(v=sub$DATETIME[sub$end==1],col=2)
    #abline(v=sub$DATETIME[sub$putstart==1],col=3,lty=2)
    #abline(v=sub$DATETIME[sub$putend==1],col=2,lty=2)
    abline(v=sub$DATETIME[sub$latestart==1],col=5)
    abline(v=sub$date_pl,col=6)
    abline(v=sub$date_med,col=6)
    #points(adjVI~DATETIME,data=sub,col=rotation,pch=5)
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0)
    abline(h=-300,lty=2)
    #abline(h=300,lty=2)
    #abline(h=-1*sub$thresh,lty=2)
    abline(h=sub$thresh,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotvizadjlite=function(dfr,lob=1){
  par(mar=c(2,2,2,1),mfrow=c(5,1))
  for(i in seq(lob,lob+4)){
    sub=dfr[dfr$ID==unique(dfr$ID)[i],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1000,3000),
         main=paste(sub$co[1],sub$soil[1],sub$elevbin[1],
                    round(sub$sharedlat[1],2),
                    round(sub$sharedlong[1],2),sep=' '))
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0,lty=2)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,1))
}

rotvizlite1=function(dfr,lob=1){
  par(mfrow=c(1,1),mar=c(4,4,2,1))
  sub=dfr[dfr$ID==unique(dfr$ID)[lob],]
    ts=zoo(x=sub$adjVI, order.by=sub$DATETIME)
    plot(ts,ylim=c(-1200,3000),
         main=paste(sub$co[1],sub$soil[1],sub$elevbin[1],
                    round(sub$sharedlat[1],2),
                    round(sub$sharedlong[1],2),sep=' '),
         xlab='Year (3 observations per year)',
         ylab='NIRv_euc')
    points(adjVI~DATETIME,data=sub[sub$badrots==0,],
           col=ezrot,pch=1)
    abline(h=0,lty=2)
}

rotvizlite1_VI=function(dfr,lob=1){
  par(mfrow=c(1,1),mar=c(4,4,2,1))
  sub=dfr[dfr$ID==unique(dfr$ID)[lob],]
  ts=zoo(x=I(sub$NIRv/1000), order.by=sub$DATETIME)
  plot(ts, xlab='Year (3 observations per year)',
       ylab='NIRv',las=1,ylim=c(0,4))
  points(I(NIRv/1000)~DATETIME,data=sub[sub$ezrot>0,],
         col=ezrot,pch=1)
  #arrows(x0=unique(sub$peakdate[sub$ezrot==1]),y0=0,
  #       x1=unique(sub$peakdate[sub$ezrot==1]),
  #       y1=max(sub$NIRv[sub$ezrot==1]),col=4,lwd=2,
  #       angle=90,code=2)
  abline(h=0,lty=2)
} 
