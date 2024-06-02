
## Load libraries and auxillary functions
library(tidyverse)
library(MASS)
library(caret)
library(sensitivitymult)
library(optmatch)
library(glmnet)
library(kableExtra)
library(ggplot2)
library(xtable)
options("optmatch_max_problem_size" = Inf)

## Function for computing
## rank based Mahalanobis distance. Prevents an outlier from
## inflating the variance for a variable, thereby decreasing its importance.
## Also, the variances are not permitted to decrease as ties
## become more common, so that, for example, it is not more important
## to match on a rare binary variable than on a common binary variable
## z is a vector, length(z)=n, with z=1 for treated, z=0 for control
## X is a matrix with n rows containing variables in the distance
smahal=function(z,X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(z)
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-as.matrix(X[z==0,])
  Xt<-as.matrix(X[z==1,])
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  icov<-ginv(cv)
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
  out
}

## Function for adding a propensity score caliper to a distance matrix dmat
## calipersd is the caliper in terms of standard deviation of the logit propensity score
addcaliper=function(dmat,z,logitp,calipersd=.25,penalty=1000){
  # Pooled within group standard devation
  sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}


## Clean data and divide dat into matching columns and outcome columns
#dat <- dat %>% filter(faffect != "Moderate")
dat$treated <- ifelse(dat$faffect=='Not exposed', 0, 1)
dat=subset(dat,select=-c(hhnum,faffect))
matching_cols=c("treated", "perc_male_0to4", "perc_male_5to14", "perc_male_15to19", "perc_male_20to34", "perc_male_35to54", "perc_male_55up", "perc_female_0to4", "perc_female_5to14", "perc_female_15to19", "perc_female_20to34", "perc_female_35to54", "perc_female_55up", "nmale_none", "nmale_anyprimary", "nmale_anysecondary", "nfemale_none", "nfemale_anyprimary", "nfemale_anysecondary", "household_size", "quant.Assets.G9", "quant.Assets.G6", "quant.Assets.G5", "quant.Assets.G4", "quant.Assets.G3", "quant.Assets.G2", "quant.Assets.G1", "eslosp.Assets.G9", "eslosp.Assets.G6", "eslosp.Assets.G5", "eslosp.Assets.G4", "eslosp.Assets.G3", "eslosp.Assets.G2", "eslosp.Assets.G1", "femalehead", "agehead", "whousefa")
data_matching = subset(dat,select=matching_cols)
data_outcomes = subset(dat,select=!colnames(dat) %in% matching_cols)


## Run propensity score matching
datatemp=as.data.frame(data_matching)
datatemp=subset(datatemp,select=-c(perc_female_55up,femalehead,agehead))
rows_datamatching_na <- which(apply(datatemp,1,function(x) any(is.na(x))))
datatemp=datatemp[-rows_datamatching_na,]

## Propensity score model
set.seed(0)
best_lambda = cv.glmnet(x=model.matrix( ~.,subset(datatemp, select=-c(treated)) )[,-1], y=datatemp$treated, family='binomial', nfolds=15, alpha=0)$lambda.min
propscore.model = glmnet(x=subset(datatemp, select=-c(treated)), y=datatemp$treated, family='binomial', alpha = 0, lambda = best_lambda)
#propscore.model = glm(treated~., family='binomial', data=datatemp)
dmy=dummyVars('treated~.',data=datatemp)
dmy$lvls=dmy$lvls[-1]
Xmat=data.frame(predict(dmy,newdata=datatemp))
Xmatmahal=Xmat
datatemp$logit.ps=as.numeric(predict(propscore.model,newx=as.matrix(Xmat)))
## Use Hansen (2009)’s rule for removing subjects who lack overlap
logit.propscore=datatemp$logit.ps
pooled.sd.logit.propscore=sqrt(var(logit.propscore[datatemp$treated==1])/2+var
                               (logit.propscore[datatemp$treated==0])/2)
min.treated.logit.propscore=min(logit.propscore[datatemp$treated==1])
max.control.logit.propscore=max(logit.propscore[datatemp$treated==0])
## How many treated and control subjects lack overlap by Hansen's criterion
no.treated.lack.overlap=sum(logit.propscore[datatemp$treated==1]>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))
no.control.lack.overlap=sum(logit.propscore[datatemp$treated==0]<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore))
## If there are subjects who lack overlap, remove them from the datatemp dataset
datatemp.original=datatemp
datatemp.full=datatemp
Xmat.original=Xmat
Xmat.full=Xmat
which.remove=NA
if(no.treated.lack.overlap+no.control.lack.overlap>0){
  which.remove=which((logit.propscore>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))|(logit.propscore<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore)))
  datatemp=datatemp[-which.remove,]
  datatemp.full=rbind(datatemp,datatemp.original[which.remove,])
  Xmat=Xmat[-which.remove,]
  Xmat.full=rbind(Xmat,Xmat.original[which.remove,])
  Xmatmahal=Xmatmahal[-which.remove,]
  logit.propscore=c(logit.propscore[-which.remove],logit.propscore[which.remove])
}
## For the purposes of balance checking later, in datatemp.full, append
## the removed rows of datatemp to the end of datatemp
rownames(datatemp)=seq(1,nrow(datatemp),1)
Xmatmahal$logit.ps=datatemp$logit.ps
## Rank based Mahalanobis distance with caliper
distmat=smahal(datatemp$treated,as.matrix(Xmatmahal$logit.ps))
distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=.25)
rownames(distmat)=rownames(datatemp)[datatemp$treated==1]
colnames(distmat)=rownames(datatemp)[datatemp$treated==0]
## Match nocontrols.per.match to each treated unit
nocontrols.per.match=1
matchvec=pairmatch(distmat,controls=nocontrols.per.match,data=datatemp)
datatemp$matchvec=matchvec
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
## Have matchedset.index.numeric.full append 0 to matchedset.index.numeric for
## the removed subjects
if(no.control.lack.overlap+no.treated.lack.overlap==0){
  matchedset.index.numeric.full=matchedset.index.numeric
}
if(no.control.lack.overlap+no.treated.lack.overlap>0){
  matchedset.index.numeric.full=c(matchedset.index.numeric,rep(0,no.control.lack.overlap+no.treated.lack.overlap))
}
## Create a matrix saying which control units each treated unit is matched to
## Create vectors of the subject indices of the treatment units ordered by
## their matched set and corresponding control unit
treated.subject.index=rep(0,max(matchedset.index.numeric.full,na.rm=T))
matched.control.subject.index.mat=matrix(rep(0,nocontrols.per.match*length(treated.subject.index)),ncol=nocontrols.per.match)
for(i in 1:length(treated.subject.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  treated.temp.index=which(datatemp$treated[matched.set.temp]==1)
  treated.subject.index[i]=matched.set.temp[treated.temp.index]
  matched.control.subject.index.mat[i,]=matched.set.temp[-treated.temp.index]
}
matched.control.subject.index=matched.control.subject.index.mat
num_matched_sets=nrow(matched.control.subject.index.mat)
## Check balance
missing.mat=matrix(rep(0,ncol(Xmat.full)*nrow(Xmat.full)),ncol=ncol(Xmat.full)
)
Xmat.without.missing=Xmat.full
for(i in 1:ncol(Xmat.full)){
  Xmat.without.missing[missing.mat[,i]==1,i]=NA
}
## Also compute balance on logit propensity score
Xmat.without.missing$logit.ps=logit.propscore
treatedmat=Xmat.without.missing[datatemp.full$treated==1,];
## Standardized differences before matching
controlmat.before=Xmat.without.missing[datatemp.full$treated==0,];
controlmean.before=apply(controlmat.before,2,mean,na.rm=TRUE);
treatmean=apply(treatedmat,2,mean,na.rm=TRUE);
treatvar=apply(treatedmat,2,var,na.rm=TRUE);
controlvar=apply(controlmat.before,2,var,na.rm=TRUE);
stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);
## Standardized differences after matching
treatmat.after=Xmat.without.missing[treated.subject.index,]
controlmat.after=Xmat.without.missing[matched.control.subject.index,];
controlmean.after=apply(controlmat.after,2,mean,na.rm=TRUE);
treatmean.after=apply(treatmat.after,2,mean,na.rm=TRUE)
## Standardized differences after matching
stand.diff.after=(treatmean.after-controlmean.after)/sqrt((treatvar+controlvar)/2);
sd.bf=stand.diff.before[1:(length(stand.diff.before)-1)]
sd.af=stand.diff.after[1:(length(stand.diff.after)-1)]
standardized_difs=cbind(sd.bf,sd.af)
standardized_difs
as.data.frame(standardized_difs) %>%
  kable(format = "latex", 
        booktabs = TRUE,
        digits = 2,
        caption = "Standardized mean differences",
        label = "SMD") %>%
  kable_styling(latex_options = "hold_position") %>%
  save_kable("SMD.tex")

names(sd.bf) <- c("male.0.to.4","male.5.to.14","male.15.to.19",
                  "male.20.to.34","male.35.to.54","male.over.54","female.0.to.4",
                  "female.5.to.14","female.15.to.19","female.20.to.34","female.35.to.54",
                  "male.no.educ","male.primary.educ","male.secondary.educ",
                  "female.no.educ","female.primary.educ","female.secondary.educ",
                  "household.size","quant.livestock",
                  "quant.valuables","quant.domestic","quant.transport","quant.water",
                  "quant.agric.equip","quant.house","val.pre.livestock",
                  "val.pre.valuables","val.pre.domestic","val.pre.transport","val.pre.water",
                  "val.pre.agric.equip","val.pre.house","normal.flood.depth")

cov_order <- rev(c("male.0.to.4","male.5.to.14","male.15.to.19",
                   "male.20.to.34","male.35.to.54","male.over.54","female.0.to.4",
                   "female.5.to.14","female.15.to.19","female.20.to.34","female.35.to.54",
                   "male.no.educ","male.primary.educ","male.secondary.educ",
                   "female.no.educ","female.primary.educ","female.secondary.educ",
                   "household.size","quant.livestock",
                   "quant.valuables","quant.domestic","quant.transport","quant.water",
                   "quant.agric.equip","quant.house","val.pre.livestock",
                   "val.pre.valuables","val.pre.domestic","val.pre.transport","val.pre.water",
                   "val.pre.agric.equip","val.pre.house","normal.flood.depth"))

plot.dataframe=data.frame(SMD=c(abs(sd.af),abs(sd.bf)),Covariates=rep(names(sd.bf),2),type=c(rep("After Matching",length(names(sd.bf))),rep("Before Matching",length(names(sd.bf)))))
ggplot(plot.dataframe,aes(x=SMD,y=Covariates))+
  geom_point(size=3,aes(shape=factor(type,levels = c('Before Matching','After Matching'))))+
  scale_shape_manual(values =c(21,16))+
  scale_y_discrete(limits = cov_order)+
  geom_vline(xintercept=c(0,0.1),lty=2) +
  labs(x = "Absolute Standardized Mean Differences", y="Covariates") + theme_bw() +
  theme(axis.text.y=element_text(size=8),legend.title = element_blank(),legend.position="bottom")

## Update data_outcomes; extract treated and control indices
data_outcomes <- data_outcomes[-rows_datamatching_na,]
data_outcomes <- rbind(data_outcomes[-which.remove,],data_outcomes[which.remove,])
ix_subjects_treated=treated.subject.index
ix_subjects_control=matched.control.subject.index

## Get rid of outcomes where NA proportion is over 75%
outcomestoomanyna <- which( apply(data_outcomes,2,function(x) sum(is.na(x))/nrow(data_outcomes)) >= 0.75 )
data_outcomes <- subset(data_outcomes, select=-outcomestoomanyna)

names(data_outcomes) <- c('food.price.idx','val.post.livestock','val.post.valuables','val.post.domestic',
                          'val.post.transport','val.post.water','val.post.agric.equip','val.post.house',
                          'calories.pc','protein.pc',
                          'rice.consumed.pc','wheat.consumed.pc','pulses.consumed.pc','oil.consumed.pc',
                          'veg.consumed.pc','egg.consumed.pc','milk.consumed.pc','fruit.consumed.pc',
                          'fish.consumed.pc','meat.consumed.pc','cereals.consumed.pc',
                          'rice.budget','wheat.budget','pulses.budget','oil.budget',
                          'veg.budget','egg.budget','milk.budget','fruit.budget',
                          'fish.budget','meat.budget','cereals.budget',
                          'tot.expend', 'marriage.funerals.expend', 'tot.expend.incl.repairs',
                          'tot.expend.pc', 'food.expend.pc', 'share.food.expend',
                          'share.food.expend.incl.repairs', 'share.nonfood.expend',
                          'share.nonfood.expend.incl.repairs', 'nonfood.expend',
                          'nonfood.expend.incl.repairs', 'food.consumed.val', 'food.produced.val',
                          'food.received.val', 'food.purchased.val', 'credit.val', 'nonfood.expend.pc',
                          'food.produced.val.pc', 'food.received.val.pc', 'food.purchased.val.pc',
                          'food.credit.val.pc', 'adult.equivalent.tot.expend', 'calories.per.adult',
                          'protein.per.adult', 'icrs.fish.male', 'icrs.fish.female', 'amt.fish.male',
                          'amt.fish.female', 'amt.milk.meat.egg.male', 'amt.milk.meat.egg.female',
                          'ft.water', 'days.water', 'days.away.home', 'repair.cost',
                          'sanitary.latrine', 'tubewell.drinking', 'tubewell.cooking', 
                          'tubewell.washing', 'distance.drinking', 'distance.cooking',
                          'distance.washing', 'time.water.collection', 'liters.water.collection',
                          'any.electricity',
                          'any.ill.0-5yrs', 'any.ill.over18yrs', 'any.ill.6-18yrs',
                          'any.ill.gastrointest', 'any.ill.respiratory', 'any.ill.fever',
                          'any.ill.female', 'any.ill.male', 'days.ill',
                          'rd3.height.on.age', 'rd3.weight.on.age', 'rd3.weight.on.height',
                          'rd3.bmi',
                          'rd3.height.on.age.preschool', 'rd3.weight.on.age.preschool', 'rd3.weight.on.height.preschool',
                          'rd3.bmi.adolesc.women')

