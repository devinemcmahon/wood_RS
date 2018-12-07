# Non-lasso (best subsets) approach to identify a simple statistical model of 
#   wood volume, based on the aspects of the NIRv time series for the rotation
# use with data from wood_vol_dataprep
source('wood_vol_dataprep.R')

library(leaps)
# The full best subset approach for two types of models: 
# based solely on plantation NIRv (strongest biophysical basis and most useful
#   for predicting within-site trends, used in final analysis)
# based on plantation NIRv and site properties
# Also tested methods based on difference between plantation and native NIRv
#   on a given date; these were sensitive to variation in timing of native vs
#   plantation observations, especially when multiple satellites contributed
#   to each composite image, and to variation in native NIRv
#   which were not clearly related to wood volume; these tests are not shown
leap1=regsubsets(Vol~.,data=traindat3,method='exhaustive',nvmax=20, really.big=T)
nleap1=regsubsets(Vol~.,data=ntraindat3,method='exhaustive',nvmax=20, really.big=T)
summary(leap1)$bic
which.min(summary(leap1)$bic) 
which.min(summary(nleap1)$bic) 

names(which(summary(leap1)$which[1,]==T))

# A function to find the optimal number of predictors, based on
#   predefined training and test datasets
bestlmfinder=function(bestsuboutput,traindata,testdata){
  num_vbls=length(summary(bestsuboutput)$which[,1])
  lmlist=vector('list',num_vbls)
  predlist=rep(NA,num_vbls)
  for(i in seq(num_vbls)){
    vbls=names(which(summary(bestsuboutput)$which[i,]==T))
    thislm=lm(Vol~.,data=traindata[,which(names(traindata) %in% c(vbls,'Vol'))])
    lmlist[[i]]=thislm
    preds=predict(thislm,newdata=testdata[,which(names(testdata) %in% c(vbls,'Vol'))])
    rmse=sqrt(sum((preds-testdata$Vol)^2,na.rm=T)/sum(!is.na(preds)))
    predlist[i]=rmse
  }
    return(list(lmlist,predlist,which(predlist==min(predlist))))
} 

## Apply this function to the regsubsets output generated above
bestleap1=bestlmfinder(leap1,traindat3,testdat3)
bestleap1[[3]] # 5
plot(bestleap1[[2]])

nbestleap1=bestlmfinder(nleap1,ntraindat3,ntestdat3)
nbestleap1[[3]] # 3
plot(nbestleap1[[2]]) # higher rmse than the other two

summary(bestleap1[[1]][[bestleap1[[3]]]]) 
# MAP, medage2, precip_drymo, precip_wetmo, rotlenln
# r-sq 0.43
summary(nbestleap1[[1]][[nbestleap1[[3]]]]) 
# maxVI, medage2, rotlen
# r-sq 0.35

# can we do nearly as well with fewer predictors?
summary(bestleap1[[1]][[3]]) # R-sq 0.35, precip_drymo, prediff, and rotlen
summary(nbestleap1[[1]][[3]]) # R-sq 0.35, maxVI, medage2, rotlen


1-sum((predict(bestleap1[[1]][[5]])-traindat3$Vol)^2)/
          sum((traindat3$Vol-mean(traindat3$Vol))^2) 
# .43
 1-sum((predict(bestleap1[[1]][[5]],newdata=testdat3)-testdat3$Vol)^2)/
          sum((testdat3$Vol-mean(testdat3$Vol))^2) 
# 0.44: overfit to test data?
 
plot(predict(bestleap1[[1]][[3]])~traindat3$Vol)
abline(0,1) # not horrible
plot(predict(bestleap1[[1]][[3]],newdata=testdat3)~testdat3$Vol)
abline(0,1)
plot(resid(bestleap1[[1]][[3]])~predict(bestleap1[[1]][[3]])) # blob, ok
plot(resid(bestleap1[[1]][[3]])~traindat3$rotlen)
# looks good for this, mnat, 
# some higher resids at medium-high precip_wetmo, nothing crazy
# slight decrease with MAT?

# try a cross-validation method for selecting optimal number of predictors
bestlmfinderxval=function(data,seed,nvbls){
  set.seed(seed)
  data$random=runif(length(data$Vol))
  #for(i in seq(round(length(data$Vol/20),0))){
  predmat=matrix(nrow=10,ncol=nvbls)
  for(i in seq(10)){
      indx=i/10
      traindata=data[-which(data$random<indx & data$random>=(indx-0.1)),]
      testdata=data[which(data$random<indx & data$random>=(indx-0.1)),]
      traindata=traindata[,-which(names(traindata)=='random')]
      testdata=testdata[,-which(names(testdata)=='random')]
      bestsuboutput=regsubsets(Vol~.,data=traindata,method='exhaustive',
                               nvmax=nvbls, really.big=T) 
      for(j in seq(nvbls)){
          vbls=names(which(summary(bestsuboutput)$which[j,]==T))
          thislm=lm(Vol~.,data=traindata[,which(names(traindata) %in% c(vbls,'Vol'))])
          preds=predict(thislm,newdata=testdata[,which(names(testdata) %in% 
                                                         c(vbls,'Vol'))])
          rmse=sqrt(sum((preds-testdata$Vol)^2,na.rm=T)/sum(!is.na(preds)))
          predmat[i,j]=rmse
      }
  }
  predsums=apply(predmat,2,sum)
  return(predsums)
}

# Initially performed on alldat3, not traindat

numpreds=bestlmfinderxval(traindat3,10,20)
which.min(numpreds) # 18--> 12 w/o rotlensq (8 for alldat)
plot(numpreds) #but anything over 7 (7 of more) is similar
nnumpreds=bestlmfinderxval(ntraindat3,10,14) # there are only 15 variables
which.min(nnumpreds) # 10 (10)
plot(nnumpreds) # 6 ok, too; 10 better

# 8 predictors
lm8=lm(Vol~.,data=traindat3[,which(names(traindat3) %in% 
                                   c(names(coef(leap1,8)),'Vol'))])
nlm8=lm(Vol~.,data=ntraindat3[,which(names(ntraindat3) %in% 
                                     c(names(coef(nleap1,8)),'Vol'))])

# How different is a model with only three predictors?
lm3=lm(Vol~.,data=traindat3[,which(names(traindat3) %in% 
                                     c(names(coef(leap1,3)),'Vol'))])
nlm3=lm(Vol~.,data=ntraindat3[,which(names(ntraindat3) %in% 
                                       c(names(coef(nleap1,3)),'Vol'))])

# Coefficient of determination with 3 or 6 predictors
1-sum((predict(lm3,newdata=traindat3)-traindat3$Vol)^2)/
  sum((traindat3$Vol-mean(traindat3$Vol))^2)  
1-sum((predict(lm3,newdata=testdat3)-testdat3$Vol)^2)/
  sum((testdat3$Vol-mean(testdat3$Vol))^2)  # 32% vs 35 for train
1-sum((predict(nlm3,newdata=ntestdat3)-ntestdat3$Vol)^2)/
  sum((ntestdat3$Vol-mean(ntestdat3$Vol))^2) # 35% for train or test

1-sum((predict(lm8,newdata=traindat3)-traindat3$Vol)^2)/
  sum((traindat3$Vol-mean(traindat3$Vol))^2)  # 40% for test, 48 for train
1-sum((predict(nlm8,newdata=ntestdat3[!is.na(ntestdat3$postpeakmed),])-
         ntestdat3$Vol[!is.na(ntestdat3$postpeakmed)])^2)/
  sum((ntestdat3$Vol[!is.na(ntestdat3$postpeakmed)]-
         mean(ntestdat3$Vol[!is.na(ntestdat3$postpeakmed)]))^2) 
# 33% test, 40% train


# Plots to visually assess the accuracy of the different models
#   and the differences between them
#   including when extrapolated over the whole dataset
# Can skip this
###################
plot(predict(lm8)~Vol,data=traindat3)
abline(0,1)
points(predict(nlm8)~Vol,data=ntraindat3,col=3)

plot(predict(lm8,newdata=mostrotsums)~vol1,data=mostrotsums)
abline(0,1)
plot(predict(nlm8,newdata=mostrotsums)~vol1,data=mostrotsums)
abline(0,1) # no good

plot(predict(lm3)~Vol,data=traindat3)
abline(0,1)
points(predict(nlm3)~Vol,data=ntraindat3,col=3)

plot(predict(lm8,newdata=testdat3)~Vol,data=testdat3)
abline(0,1)
points(predict(nlm8,newdata=ntestdat3)~Vol,data=ntestdat3,col=3)

mostrotsums=readRDS('cos+219072+218073_okrotsums_morevbls.Rds')

plot(predict(lm8,newdata=mostrotsums)~I(rotlen/12),data=mostrotsums)
points(predict(nlm8,newdata=mostrotsums)~I(rotlen/12),data=mostrotsums,col=3)

points(Vol~I(rotlen/12),data=bios3,col=6)
hist(I(mostrotsums$rotlen/12),xlab='Rotation length (years)')
quantile(mostrotsums$rotlen,0.9,na.rm=T)/12 # 90% < 11.5 years
# 95% < 13.5 years

plot(predict(lm3,newdata=mostrotsums)~I(rotlen/12),data=mostrotsums)
points(predict(nlm3,newdata=mostrotsums)~I(rotlen/12),data=mostrotsums,col=3)


boxplot(predict(lm8,newdata=mostrotsums)~ezrot,data=mostrotsums,
        varwidth=T,ylab='Predicted wood volume, lm8')
boxplot(predict(nlm8,newdata=mostrotsums)~ezrot,data=mostrotsums,
        varwidth=T,boxwex=0.4,border=3,add=T)

# how do the distributions of the training and extrapolation data compare?

boxplot(medage2~ezrot,data=mostrotsums,ylab='Median NIRv to age 2',
        xlim=c(0.5,6.5),varwidth=T)
boxplot(alldat3$medage2,at=6,add=T,border=4)



most3rots=mostrotsums[mostrotsums$maxrots==3,]

boxplot(predict(lm8,newdata=most3rots)~ezrot,data=most3rots,
        varwidth=T,ylab='Predicted wood volume, 8 variables',
        xlim=c(0.5,4.5),  ylim=c(0,700))
boxplot(predict(nlm8,newdata=most3rots)~ezrot,data=most3rots,
        varwidth=T,boxwex=0.4,border=3,add=T)
boxplot(alldat3$Vol,at=4,add=T,border=4)


boxplot(predict(lm3,newdata=most3rots)~ezrot,data=most3rots,
        varwidth=T,ylab='Predicted wood volume, 3 variables',
        xlim=c(0.5,4.5),ylim=c(0,700))
bboxplot(predict(nlm3,newdata=most3rots)~ezrot,data=most3rots,
        varwidth=T,boxwex=0.4,border=3,add=T)
boxplot(alldat3$Vol,at=4,add=T,border=4)
#legend('topright',legend=c('NIRv terms','delta NIRv terms',))


boxplot(I(rotlen/12)~ezrot,data=most3rots,xlim=c(0.5,4.5),
        varwidth=T,ylab='Rotation length (years)',las=1)
boxplot(I(alldat3$rotlen/12),at=4,add=T,border=4,las=1)


boxplot(medage2~ezrot,data=most3rots,ylab='Median NIRv to age 2',
        xlim=c(0.5,4.5),varwidth=T)
boxplot(alldat3$medage2,at=4,add=T,border=4)
# lower values for test data than extrapolation data
# does this just reflect the different geographical distribution?

points(predict(nlm8,newdata=ntestdat3)~Vol,data=ntestdat3,col=3)

inbios=mostrotsums[mostrotsums$ID %in% bios3$ID,]
boxplot(medage2~ezrot,data=inbios,ylab='Median deltaNIRv to age 2',
        xlim=c(0.5,5.5),varwidth=T)
boxplot(alldat3$medage2,at=5,add=T,border=4)

plot(medage2~rotstartdate,data=inbios,col=ezrot)
points(medage2~rotstartdate,data=bios3,pch=5,col=foundrot)

plot(medage2~rotstartdate,data=mostrotsums,col=ezrot)
points(medage2~rotstartdate,data=bios3,pch=5,col=7)
###################

# Best subsets analysis, but without withholding a test set
bestall=regsubsets(Vol~.,data=alldat3,method='exhaustive',
                     nvmax=20, really.big=T)
nbestall=regsubsets(Vol~.,data=nalldat3,method='exhaustive',
                    nvmax=20, really.big=T)


coef(leap1,8)
#   Too many terms; canceling occurring
coef(nleap1,8) 

lm8=lm(Vol~.,data=traindat3[,which(names(traindat3) %in% 
                                   c(names(coef(leap1,8)),'Vol'))])

summary(bestall)$bic # 7 is the best with all data
coef(bestall,7) 
# but has cancelling meanVI and medVI terms, and -MAP/+drymo/+wetmo
coef(bestall,4)
# elev, MAT, medage2, rotlensq 
# why MAT? cooler sites do have higher annual increment
# elev and MAT cancel, I think: cooler = higher, correlated
# go with 3 terms to avoid this effect?
# or use it to introduce nonlinearities?
# without rotlensq, it's MAP, wetmo, drymo, rotlen
# no good for trends
coef(bestall,3) 
# now MAT, medage2, rotlen
coef(nbestall,3) # maxVI, medage2, rotlen
# go with that
summary(abestall)$bic # 8 or 9 best fit; 
# substantial improvement with each term added 3 to 7
# but will those relationships hold over a larger data set?

alllm3=lm(Vol~maxVI+medage2+rotlen,data=alldat3)
summary(alllm3)

#saveRDS(alllm3,'vol_lm_3termsnonat.Rds') # final version
