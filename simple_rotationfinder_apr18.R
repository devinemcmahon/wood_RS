# This script segments a time series of NIRv values into rotations
## !diagnostics off

myrollsum=function(x){sum(x,na.rm=T)}

findrotations = function(dfr,leavecols=FALSE){
  # Assign unique date to each observation to order by date
  
  # preserve old DATETIME values
  dfr$oldDATETIME=dfr$DATETIME
  dfr$IDDATE=paste(dfr$ID,dfr$DATETIME,sep='-')
  dfr=dfr[order(dfr$DATETIME),] 
  
  # when two observations have same date, add a day to the second
  dfr$DATETIME[which(duplicated(dfr$IDDATE))]=
    dfr$DATETIME[which(duplicated(dfr$IDDATE))]+days(1)
  # do it again in case all 3 observations for some year come from same date
  dfr$IDDATE=paste(dfr$ID,dfr$DATETIME,sep='-')
  dfr$DATETIME[which(duplicated(dfr$IDDATE))]=
    dfr$DATETIME[which(duplicated(dfr$IDDATE))]+days(1)
  dfr$IDDATE=paste(dfr$ID,dfr$DATETIME,sep='-') #leave the column to check
  dfr=dfr[order(dfr$DATETIME),]
  
  dfr=group_by(dfr,ID) %>% mutate(meanadj=mean(adjVI,na.rm=T),
                                  q80adj=quantile(adjVI,0.8,na.rm=T))
  dfr$thresh=rep(300)
  dfr$thresh[dfr$meanadj<300]=200
  
  # Create the necessary columns
  dfr$pos=rep(0)
  dfr$pos[dfr$adjVI>0]=1
  dfr$gt=rep(0)
  dfr$gt[dfr$adjVI>dfr$thresh]=1 
  dfr$ggt=as.numeric(dfr$q80adj<dfr$thresh*2)
  # if adjVI is always small, void the rule
  #   that ggt must be exceeded to start a rotation
  #   by setting default to 1 rather than 0
  dfr$ggt[dfr$adjVI>dfr$thresh*2]=1
  dfr$lt=rep(0)
  #dfr$lt[dfr$adjVI< -1*dfr$thresh]=1
  dfr$lt[dfr$adjVI< -300]=1
  dfr$ones=rep(1)
  
  dfr=group_by(dfr,ID) %>%
    mutate(lag1gt=lag(gt,1),lead1gt=lead(gt,1),
           lag1pos=lag(pos,1),
           lag1lt=lag(lt,1),lead1lt=lead(lt,1),
           lag1vi=lag(adjVI,1),lead1vi=lead(adjVI,1),
           lag2pos=lag(pos,2),lag2gt=lag(gt,2),lag2lt=lag(lt,2)
    )
  
  # Replace missing values by repeating prior values
  dfr$pos[is.na(dfr$adjVI)]=dfr$lag1pos[is.na(dfr$adjVI)]
  dfr$gt[is.na(dfr$adjVI)]=dfr$lag1gt[is.na(dfr$adjVI)]
  dfr$lt[is.na(dfr$adjVI)]=dfr$lag1lt[is.na(dfr$adjVI)]
  dfr$lead1gt[is.na(dfr$lead1vi)]=dfr$gt[is.na(dfr$lead1vi)]
  dfr$lead1lt[is.na(dfr$lead1vi)]=dfr$lt[is.na(dfr$lead1vi)]
  dfr$lag1pos[is.na(dfr$lag1vi)]=dfr$lag2pos[is.na(dfr$lag1vi)]
  dfr$lag1gt[is.na(dfr$lag1vi)]=dfr$lag2gt[is.na(dfr$lag1vi)]
  dfr$lag1lt[is.na(dfr$lag1vi)]=dfr$lag2lt[is.na(dfr$lag1vi)]
  
  # Create the necessary derived columns
  dfr= group_by(dfr,ID) %>%
    mutate(next15gt=rollsum(gt,k=15,fill=0,align='left'),
           # using rollsum prevents new rotation from starting
           #    within 5 years of end of time series
           last12gt=rollsum(gt,k=12,fill=0,align='right'),
           last12pos=rollsum(pos,k=12,fill=0,align='right'),
           # or ending within 4 years of the start
           next9gt=rollsum(gt,k=9,fill=0,align='left'),
           next15ggt=rollsum(ggt,k=15,fill=0,align='left'),
           last8gt=rollsum(gt,k=8,fill=0,align='right'),
           #next15lt=rollsum(lt,k=15,fill=0,align='left'),
           #next10lt=rollsum(lt,k=10,fill=0,align='left'),
           next15lt=rollapply(lt,width=15, align='left',
                              FUN=myrollsum,partial=T),
           next10lt=rollapply(lt,width=10, align='left',
                              FUN=myrollsum,partial=T),
           next6lt=rollapply(lt,width=6, align='left',
                              FUN=myrollsum,partial=T),
           next3lt=rollapply(lt,width=3, align='left',
                              FUN=myrollsum,partial=T),
           # myrollsum allows this to apply at end of time series, too
           last4lt=rollsum(lt,k=4,fill=0,align='right'),
           next15sum=rollapply(adjVI,width=15,align='left',
                               FUN=myrollsum,partial=T),
           next15ones=rollapply(ones,width=15,align='left',
                                FUN=myrollsum,partial=T),
           next18ones=rollapply(ones,width=18,align='left',
                                FUN=myrollsum,partial=T),
           next6sum=rollapply(adjVI,width=6,align='left',
                              FUN=myrollsum,partial=T),
           next7to15mean=(next15sum-next6sum)/9, #change this to address NAs? 
           next3sum=rollapply(adjVI,width=3,align='left',
                              FUN=myrollsum,partial=T),
           next4to6mean=(next6sum-next3sum)/3,
           next3mean=next3sum/3,
           next2to3mean=(next3sum-adjVI)/2,
           ctr5min=rollapply(adjVI,width=5,align='center',fill=NA,partial=T,
                             FUN=function(x){min(x,na.rm=T)}), 
           next12min=rollapply(adjVI,width=12,align='left',fill=NA,partial=T,
                               FUN=function(x){min(x,na.rm=T)}),
           next15min=rollapply(adjVI,width=15,align='left',fill=NA,partial=T,
                               FUN=function(x){min(x,na.rm=T)})
    )
  
  
  # All the putative starts 
  # Leave these identified while algorithm under development
  dfr$putstart=rep(0)
  dfr$putstart[dfr$gt==0 & # at start of rotation, euc VI <= native VI
                 dfr$next9gt >= 2 &  # 2 years to get to gt (euc > nat)
                 dfr$next15gt>=8 & # then stay gt for at least 3 years
                 dfr$adjVI==dfr$next12min & # and don't drop back down
                 dfr$next4to6mean>dfr$next3mean & 
                 # don't start rotation before a bump and decrease in low VI
                 #dfr$next7to15mean>400 & # peak VI can't be < nat+1 sd (was 600)
                 #dfr$next7to15mean>dfr$adjVI+600 &
                 # VI must increase after start
                 (dfr$next15ggt>=3 | dfr$lt==1) & # new 4-7-18
                 # increase threshold for calling a rotation
                 # unless VI has dipped really low
                 # lt part added 4-10-18
                 # effectively, require VI to increase by >600
                 #    after start
                 dfr$adjVI==dfr$ctr5min]=1 
  # start rotation at lowest part of clearcut period
  
  # or close to it, if a long stretch of lt precedes start
  dfr$putstart[dfr$gt==0 & 
                 dfr$next9gt >= 1 &  
                 dfr$next15gt>=7 &
                 dfr$adjVI==dfr$next12min & 
                 dfr$last4lt>=1 &
                 #dfr$next7to15mean>300 &
                 #dfr$next7to15mean>400 &
                 #dfr$next7to15mean>dfr$adjVI+600 &
                 dfr$next4to6mean>dfr$next3mean & 
                 #dfr$adjVI-dfr$ctr5min < 100]=1
                 dfr$next15ggt>=3 & # new 4-7-18
                 dfr$adjVI-dfr$ctr5min < dfr$thresh/3]=1
  
  # next15gt prevents declaring a start within 5 years of end of time series
  
  # dfr$putstart[dfr$meanadj<300 & dfr$lt==1 & dfr$leaddiff < -1*]

  # Flag the starts of rotations not completed by end of time series
  dfr$latestart=rep(0)
  dfr$latestart[dfr$gt==0 & dfr$next15ones<15 &
                  dfr$next9gt >= 2 &  # still requires 3 years of observations
                  dfr$adjVI==dfr$next15min & 
                  dfr$next4to6mean>dfr$next3mean & 
                  dfr$adjVI==dfr$ctr5min]=1
  
  # Pick the actual starts 
  dfr$start=dfr$putstart

  # Not a rotation start if adjVI doesn't increase after putative planting
  dfr=mutate(dfr,next10start=rollapply(start,width=10,align='left',
                                       FUN=myrollsum,partial=T))
  
  # Don't start a rotation too early if there's a long pre-start lag
  # Use next2to3mean to avoid identifying a rotation start prior to
  #    short bump in VI (caused by noise or competing vegetation)
  dfr$start[dfr$next10start>1 & 
            #  abs(dfr$next2to3mean-dfr$adjVI)<200]=0
              abs(dfr$next2to3mean-dfr$adjVI)<dfr$thresh*2/3]=0
  dfr=mutate(dfr,last10start=rollapply(start,width=10,align='right',
                                       FUN=myrollsum,partial=T))
  
  # Don't start a new rotation less than 4 years after the previous
  dfr$start[dfr$last10start>1]=0 
  
  # Update the last10start column
  dfr=mutate(dfr,last10start=rollapply(start,width=10,align='right',
                                       FUN=myrollsum,partial=T),
             next6start=rollapply(start,width=6,align='left',
                                  FUN=myrollsum,partial=T),
             next18start=rollapply(start,width=18,align='left',
                                  FUN=myrollsum,partial=T),
             sumstart1=cumsum(start))
  
  # Once starts are defined, define ends
  dfr$putend=rep(0)

  # A rotation ends when VI drops to or below native veg levels
  dfr$putend[dfr$gt==1 & dfr$lead1gt==0 &
             dfr$last12gt>=8 & dfr$last8gt>=5 & # after a stable period
             dfr$sumstart1>0 & # after a rotation has started
             #(dfr$next4to6mean < 300 | # followed by a period of low VI
               (dfr$next4to6mean < dfr$thresh |
                  dfr$next6start>0 | # or the start of a new rotation
                dfr$next18ones<18)]=1 # or the end of the time series
                # end of time series previously defined by next 15, but leave room
                # for gap between end and next start
  
  # (added 1-23-18)
  # Or when VI hovering around native veg levels drops way off
  dfr$putend[dfr$pos==1 & dfr$next3lt>=1 & dfr$last4lt==0 &
               dfr$last12pos>=8 & dfr$last8gt <=3 & # changed from 5 to 3
               dfr$sumstart1>0 & #(dfr$next4to6mean < 300 | 
                  ((dfr$next4to6mean < dfr$thresh) |
                  dfr$next6start>0 | dfr$next18ones<18)]=1
  
  
  # Pick the actual ends
  dfr$end=dfr$putend
  
  # Don't end a rotation less than 4 years after it starts
  dfr$end[dfr$last10start>0]=0 
  #dfr$toosoon=rep(0)
  #dfr$toosoon[dfr$putend==1 & dfr$last10start>0]=1
  # no need to flag these; extra starts should be removed below

  # If there are several possible ends in a row, choose the one 
  #   closest to a sharp decline
  
  #dfr=mutate(dfr, next10end=rollapply(end,width=10,align='left',
  #                                    FUN=myrollsum,partial=T))
  #dfr$end[dfr$next10end>1 & dfr$next10lt==0 & dfr$next15lt>=1]=0 #& dfr$next3mean>0

  dfr=mutate(dfr,end3lt = end*next3lt, start3lt= start*next3lt,
             #next15end=rollapply(end,width=15,align='left',
             #                      FUN=myrollsum,partial=T)
             next12end3lt=rollapply(end3lt,width=12,align='left',
                                    FUN=myrollsum,partial=T), # was 15
             next12start3lt=rollapply(start3lt,width=12,align='left',
                                     FUN=myrollsum,partial=T)) 
             # this wasn't working before because of a typo 
             # (both next12s were end3lt, not start3lt)
             # fixed 2-6-18
  #dfr$end[dfr$next15end>1 & dfr$next3lt==0 & dfr$next6lt>=1]=0
  dfr$end[dfr$next12end3lt>dfr$next3lt*2 & # pick end closest to the drop-off
             (dfr$next12end3lt-dfr$next3lt) >dfr$next12start3lt]=0 
  # unless there's a start before the end closest to a drop-off
  # the *2 (added 2-7-18) is so next12end3lt - next3lt > next3lt
  # i.e. there are more observations below -300 in the next 12 than the next 3
  # it's really the depth of the decline, not the number lt, that's important
  # this seems to be working ok, though
  

  # Turn off warnings before taking mins of sets of numbers that do not
  #   exist for all IDs
  oldw <- getOption("warn")
  options(warn = -1)
  # Remove starts with no ends for 20 years (and 20 years left in time series)
  dfr = mutate(dfr,next60end=rollapply(end,width=60,align='left',
                                        FUN=myrollsum,partial=T))
  dfr$start[dfr$next60end==0 & dfr$DATETIME<as.Date('1996-01-01')]=0
  
  # Can't have two starts in a row with no ends
  dfr = mutate(dfr, #lagend=lag(end,1,default=0),
               sumstart=cumsum(start),sumend1=cumsum(end),
               #badflag=sum(max(sumend1-sumstart)>0) +
               #   sum(max(sumstart-sumend1)>1), # 2-6-18
                 # specify whether too many starts or too many ends
               badstart=sum(max(sumstart-sumend1)>1),
               #next10toosoon=rollapply(toosoon,width=10,align='left',
               #                         FUN=myrollsum,partial=T),
               next18start=rollapply(start,width=18,align='left',
                                     FUN=myrollsum,partial=T),
               next18end=rollapply(end,width=18,align='left',
                                    FUN=myrollsum,partial=T)) 
  # remove starts followed shortly by another start
  dfr$start[dfr$badstart>0 & dfr$next18start>dfr$next18end+1]=0

  # leave in starts that correspond to big decreases
  #     (lt constraints added 2-6-18)
  dfr = mutate(dfr, 
               sumstart=cumsum(start),
               badstart=sum(max(sumstart-sumend1)>1),
               firstbadstart=min(DATETIME[sumstart>sumend1+1]))
  dfr$start[dfr$badstart>0 & dfr$DATETIME==dfr$firstbadstart &
              dfr$last4lt==0 & dfr$next3lt==0]=0
  
  # can't have two ends in a row
  # usually want to pick the later of two possible ends
  # unless there's a clearcut in between
  
  dfr = mutate(dfr, #last15end=rollapply(end,width=15,align='right', #was 15 then 12
                    #               FUN=myrollsum,partial=T),
               #last15start=rollapply(start,width=15,align='right',
               #                     FUN=myrollsum,partial=T),
               #last30end=rollapply(end,width=30,align='right',
               #                     FUN=myrollsum,partial=T),
               #last30start=rollapply(start,width=30,align='right', #was 40
               #                     FUN=myrollsum,partial=T),
               last10start=rollapply(start,width=10,align='right',
                                    FUN=myrollsum,partial=T),
               next15end=rollapply(end,width=15,align='left', 
                                    FUN=myrollsum,partial=T),
               next15start=rollapply(start,width=15,align='left',
                                    FUN=myrollsum,partial=T),
               next5end=rollapply(end,width=5,align='left', 
                                   FUN=myrollsum,partial=T),
               next5start=rollapply(start,width=5,align='left',
                                     FUN=myrollsum,partial=T),
               next3end=rollapply(end,width=3,align='left', 
                                  FUN=myrollsum,partial=T),
               next3start=rollapply(start,width=3,align='left',
                                    FUN=myrollsum,partial=T))
  
  dfr$end[dfr$next15end>1 & dfr$next15start==0 & dfr$next3lt==0]=0
  dfr$end[dfr$next5end>1 & dfr$next5start==0 & dfr$next3lt==0]=0
  dfr$end[dfr$next3end>1 & dfr$next3start==0 & dfr$lead1lt==0]=0
  
  # Force an end before a clearcut and second start for which 
  #   end of previous rotation was not detected
  dfr=mutate(dfr,next3start=rollapply(start,width=3,align='left',
                                      FUN=myrollsum,partial=T),
             sumstart=cumsum(start),
             badstart=sum(max(sumstart-sumend1)>1),
             firstbadstart=min(DATETIME[sumstart>sumend1+1]))
  dfr$end[dfr$next3start==1 & dfr$firstbadstart-dfr$DATETIME < 400 & 
            dfr$firstbadstart-dfr$DATETIME > 0 & dfr$next3lt>0 &
            dfr$pos==1]=1
  
  # can't have ends with no starts
  # do this after all edits of starts are done (2-6-18)
  #dfr$end[dfr$sumstart==0]=0
   
  # Update columns
  # As of 2-6-18, remake sumend1 here
  dfr = mutate(dfr, sumstart=cumsum(start),sumend1=cumsum(end),
               badend=sum(max(sumend1-sumstart)>0),
               firstbadend=min(DATETIME[sumend1>sumstart]))
  
  # Clean up remaining two-ends-in-a-row
  dfr$end[dfr$badend>0 & dfr$DATETIME==dfr$firstbadend]=0
  # Do it again in case of multiple violations
  dfr = mutate(dfr, sumend1=cumsum(end),
               badend=sum(max(sumend1-sumstart)>0),
               firstbadend=min(DATETIME[sumend1>sumstart]))
    dfr$end[dfr$badend>0 & dfr$DATETIME==dfr$firstbadend]=0
  
  # As of 2-6-18: remove extra starts again
  dfr = mutate(dfr, badstart=sum(max(sumstart-sumend1)>1),
               firstbadstart=min(DATETIME[sumstart>sumend1+1]))
    
  # Preferentially remove starts that probably aren't real 
  #   Those that don't correspond to a large VI decrease
  dfr$start[dfr$badstart>0 & dfr$DATETIME==dfr$firstbadstart &
              dfr$last4lt==0 & dfr$next3lt==0]=0
  #   And those that are followed by long low-VI periods
  #     (should be just one of these per stand)
  #dfr$start[dfr$badflag>0 & dfr$DATETIME<dfr$firstbadstart &
  #            dfr$next15sum<(dfr$thresh*15)]=0 
  # i.e. next 15 observations not gt on average
  # Should the too-early starts be removed earlier? Nah
  # And again
  dfr = mutate(dfr, badstart=sum(max(sumstart-sumend1)>1),
               firstbadstart=min(DATETIME[sumstart>sumend1+1]))
  # remove clearcut constraint: those starts should be preceded by ends
  dfr$start[dfr$badstart>0 & dfr$DATETIME==dfr$firstbadstart]=0
  # If there really is a new start--evaluate how often that happens
  #   and maybe make a rule to force an end before that (done above)
  
  
  # Can't have ends with no starts
  dfr = mutate(dfr, sumstart=cumsum(start))
  dfr$end[dfr$sumstart==0]=0
  
  # Remove extra latestarts
  dfr = mutate(dfr, last10start=rollapply(start,width=10,align='right', 
                                        FUN=myrollsum,partial=T),
               next15end=rollapply(end,width=15,align='left', 
                                        FUN=myrollsum,partial=T))
               
  dfr$latestart[dfr$last10start>=1]=0
  dfr$latestart[dfr$next15end>=1]=0
  dfr=mutate(dfr,last15late=rollapply(latestart,width=15,align='right',
                                      FUN=myrollsum,partial=T))
  dfr$latestart[dfr$last15late>1]=0
  
 
  # Label the rotations (including the end date as part of the rotation)
  dfr = mutate(dfr, lagend=lag(end,1,default=0),
               sumstart=cumsum(start),sumend=cumsum(lagend),
               badrots=sum(max(sumend-sumstart)>0) +
                 sum(max(sumstart-sumend)>1),
               ezrot=(sumstart-sumend)*sumstart)
  
  # New 2-16-18: remove rotations with tiny VI spread?
  dfr=group_by(dfr,ID,ezrot) %>% 
    mutate(rotspread=max(adjVI,na.rm=T)-min(adjVI,na.rm=T))
  dfr$start[dfr$rotspread<400]=0
  dfr$end[dfr$rotspread<400]=0
  
  dfr = group_by(dfr,ID) %>%
    mutate(lagend=lag(end,1,default=0),
               sumstart=cumsum(start),sumend=cumsum(lagend),
               badrots=sum(max(sumend-sumstart)>0) +
                 sum(max(sumstart-sumend)>1),
               ezrot=(sumstart-sumend)*sumstart)
  
  # Turn warnings back on 
  options(warn = oldw)
  
  if(leavecols==FALSE){
    dfr=dfr[,-which(names(dfr) %in% 
            c('pos','gt','lt','ones','lag1gt','lead1gt','lag1pos','lag1lt',
              'lead1lt','lag1vi','lead1vi','lag2pos','lag2gt','lag2lt',
              'next15gt','last12gt','last12pos','next9gt','last8gt',
              'next15lt','next10lt','next6lt','next3lt','last4lt',
              'next15sum','next15ones','next18ones','next6sum','next7to15mean',
              'next3sum','next4to6mean','next3mean','next2to3mean','ctr5min',
              'next12min','next15min','last10start','next6start','next18start',
              'sumstart1','end3lt','next12end3lt','next15start3lt','last15end',
              'last15start','last30end','last30start','last10start','next15end',
              'badstart','badend','firstbadstart','firstbadend','next60end',
              'start3lt','next10start','next12start3lt','next5end',
              'next5start','next3end','next3start','sumend1','last15late'))]
  }
  dfr=ungroup(dfr)
  return(dfr)
}
