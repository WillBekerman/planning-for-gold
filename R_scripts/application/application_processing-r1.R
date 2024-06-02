
## Assumes we're in ~/bangladesh_application/data/Data household level/RND1"
# setwd('./data/Data household level/RND1')

## Load libraries and auxillary functions
library(plyr)
library(tidyverse)

## Function to make list of datasets that have household-based information
get_hhdatalist <- function(keepnums=FALSE){
  
  all_files <- list.files(".")
  dta_files <- all_files %>%
    tibble() %>%
    rename(files = 1) %>%
    filter(stringr::str_detect(files, ".dta")) %>%
    pull()
  
  datalist <- vector("list", length(dta_files))
  
  for (ix in 1:length(dta_files)){
    filename <- dta_files[ix]
    if (keepnums) { datalist[[ix]] <- readstata13::read.dta13(filename,convert.factors=F) }
    if (!keepnums) { datalist[[ix]] <- readstata13::read.dta13(filename) }
    
  }
  
  names(datalist) <- dta_files
  
  datasets_no_hhnum <- numeric()
  for (datasetix in 1:length(datalist)){
    varnames=colnames(datalist[[datasetix]])
    if (!('hhnum' %in% varnames)) {
      datasets_no_hhnum <- c(datasets_no_hhnum, names(datalist)[datasetix])
    }
  }
  datasets_no_hhnum <- c(datasets_no_hhnum,'SECA1 (1).dta')
  
  smalldatalist <- datalist
  smalldatalist <- within(smalldatalist, rm(list=datasets_no_hhnum))
  
  smalldatalist
}


## Collect covariate data
## Gathering relevant demographic data
hhdatalist <- get_hhdatalist()
hhdatalist$SECA1.dta <- hhdatalist$SECA1.dta %>% select(hhnum,pcode,sex,reltin,agey,educa)
hhdatalist$hhsize.dta <- hhdatalist$hhsize.dta %>% select(hhnum,hhsizea)
hhdatalist$tothhexp.dta <- hhdatalist$tothhexp.dta %>% select(hhnum,fprinx)
hhdatalist$strata.dta <- hhdatalist$strata.dta %>% select(hhnum,faffect)
hhdatalist <- hhdatalist[c('SECA1.dta', 'hhsize.dta','tothhexp.dta','strata.dta')]
dat <- join_all(hhdatalist, by='hhnum', type='left')

## Note: some hhnums have no household head information  -- we put in NA
## unique(dat$hhnum)[which(!(unique(dat$hhnum) %in% (dat %>% filter(reltin=='H head'))$hhnum))]

## Add household head information to dat
head_dat <- dat %>% filter(reltin=='H head') %>% select(hhnum, sex, agey)
head_dat$femalehead <- head_dat$sex == 'Female'
head_dat$agehead <- head_dat$agey
dat <- merge(dat, head_dat[, c("hhnum", "femalehead", "agehead")], by="hhnum", all.x=T)

## Partition data by sex and age; get household breakdowns (percentages) by relevant groupings
sex_dat <- dat %>% select(hhnum, sex, agey, educa)
nmale_dat <- sex_dat %>% filter(sex=='Male') %>% group_by(hhnum) %>% summarise(nmale=length(sex))
nfemale_dat <- sex_dat %>% filter(sex=='Female') %>% group_by(hhnum) %>% summarise(nfemale=length(sex))
sex_dat <- merge(sex_dat, nmale_dat[, c("hhnum", "nmale")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_dat[, c("hhnum", "nfemale")], by="hhnum", all.x=T)

nmale_0to4_dat <- sex_dat %>% filter(sex=='Male',agey>=0,agey<=4) %>% group_by(hhnum) %>% summarise(nmale_0to4=length(sex))
nmale_5to14_dat <- sex_dat %>% filter(sex=='Male',agey>=5,agey<=14) %>% group_by(hhnum) %>% summarise(nmale_5to14=length(sex))
nmale_15to19_dat <- sex_dat %>% filter(sex=='Male',agey>=15,agey<=19) %>% group_by(hhnum) %>% summarise(nmale_15to19=length(sex))
nmale_20to34_dat <- sex_dat %>% filter(sex=='Male',agey>=20,agey<=34) %>% group_by(hhnum) %>% summarise(nmale_20to34=length(sex))
nmale_35to54_dat <- sex_dat %>% filter(sex=='Male',agey>=35,agey<=54) %>% group_by(hhnum) %>% summarise(nmale_35to54=length(sex))
nmale_55up_dat <- sex_dat %>% filter(sex=='Male',agey>=55) %>% group_by(hhnum) %>% summarise(nmale_55up=length(sex))
nfemale_0to4_dat <- sex_dat %>% filter(sex=='Female',agey>=0,agey<=4) %>% group_by(hhnum) %>% summarise(nfemale_0to4=length(sex))
nfemale_5to14_dat <- sex_dat %>% filter(sex=='Female',agey>=5,agey<=14) %>% group_by(hhnum) %>% summarise(nfemale_5to14=length(sex))
nfemale_15to19_dat <- sex_dat %>% filter(sex=='Female',agey>=15,agey<=19) %>% group_by(hhnum) %>% summarise(nfemale_15to19=length(sex))
nfemale_20to34_dat <- sex_dat %>% filter(sex=='Female',agey>=20,agey<=34) %>% group_by(hhnum) %>% summarise(nfemale_20to34=length(sex))
nfemale_35to54_dat <- sex_dat %>% filter(sex=='Female',agey>=35,agey<=54) %>% group_by(hhnum) %>% summarise(nfemale_35to54=length(sex))
nfemale_55up_dat <- sex_dat %>% filter(sex=='Female',agey>=55) %>% group_by(hhnum) %>% summarise(nfemale_55up=length(sex))

sex_dat <- merge(sex_dat, nmale_0to4_dat[, c("hhnum", "nmale_0to4")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nmale_5to14_dat[, c("hhnum", "nmale_5to14")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nmale_15to19_dat[, c("hhnum", "nmale_15to19")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nmale_20to34_dat[, c("hhnum", "nmale_20to34")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nmale_35to54_dat[, c("hhnum", "nmale_35to54")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nmale_55up_dat[, c("hhnum", "nmale_55up")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_0to4_dat[, c("hhnum", "nfemale_0to4")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_5to14_dat[, c("hhnum", "nfemale_5to14")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_15to19_dat[, c("hhnum", "nfemale_15to19")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_20to34_dat[, c("hhnum", "nfemale_20to34")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_35to54_dat[, c("hhnum", "nfemale_35to54")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, nfemale_55up_dat[, c("hhnum", "nfemale_55up")], by="hhnum", all.x=T)

sex_dat$nmale_0to4 <- with(sex_dat, ifelse(is.na(nmale_0to4), 0, nmale_0to4))
sex_dat$nmale_5to14 <- with(sex_dat, ifelse(is.na(nmale_5to14), 0, nmale_5to14))
sex_dat$nmale_15to19 <- with(sex_dat, ifelse(is.na(nmale_15to19), 0, nmale_15to19))
sex_dat$nmale_20to34 <- with(sex_dat, ifelse(is.na(nmale_20to34), 0, nmale_20to34))
sex_dat$nmale_35to54 <- with(sex_dat, ifelse(is.na(nmale_35to54), 0, nmale_35to54))
sex_dat$nmale_55up <- with(sex_dat, ifelse(is.na(nmale_55up), 0, nmale_55up))
sex_dat$nfemale_0to4 <- with(sex_dat, ifelse(is.na(nfemale_0to4), 0, nfemale_0to4))
sex_dat$nfemale_5to14 <- with(sex_dat, ifelse(is.na(nfemale_5to14), 0, nfemale_5to14))
sex_dat$nfemale_15to19 <- with(sex_dat, ifelse(is.na(nfemale_15to19), 0, nfemale_15to19))
sex_dat$nfemale_20to34 <- with(sex_dat, ifelse(is.na(nfemale_20to34), 0, nfemale_20to34))
sex_dat$nfemale_35to54 <- with(sex_dat, ifelse(is.na(nfemale_35to54), 0, nfemale_35to54))
sex_dat$nfemale_55up <- with(sex_dat, ifelse(is.na(nfemale_55up), 0, nfemale_55up))

sex_dat <- sex_dat %>% mutate(perc_male_0to4 = 100*nmale_0to4/nmale)
sex_dat <- sex_dat %>% mutate(perc_male_5to14 = 100*nmale_5to14/nmale)
sex_dat <- sex_dat %>% mutate(perc_male_15to19 = 100*nmale_15to19/nmale)
sex_dat <- sex_dat %>% mutate(perc_male_20to34 = 100*nmale_20to34/nmale)
sex_dat <- sex_dat %>% mutate(perc_male_35to54 = 100*nmale_35to54/nmale)
sex_dat <- sex_dat %>% mutate(perc_male_55up = 100*nmale_55up/nmale)
sex_dat <- sex_dat %>% mutate(perc_female_0to4 = 100*nfemale_0to4/nfemale)
sex_dat <- sex_dat %>% mutate(perc_female_5to14 = 100*nfemale_5to14/nfemale)
sex_dat <- sex_dat %>% mutate(perc_female_15to19 = 100*nfemale_15to19/nfemale)
sex_dat <- sex_dat %>% mutate(perc_female_20to34 = 100*nfemale_20to34/nfemale)
sex_dat <- sex_dat %>% mutate(perc_female_35to54 = 100*nfemale_35to54/nfemale)
sex_dat <- sex_dat %>% mutate(perc_female_55up = 100*nfemale_55up/nfemale)

## Note: some individuals have no education information  -- we put in NA
## sum(is.na(sex_dat$educa)) = 146 people who we have NAs for education
## These people are ignored when getting num people with none/anyprimary/anysecondary education

male_educ_sex_dat <- sex_dat %>% filter(sex=='Male') %>% select(hhnum, educa) %>% group_by(hhnum) %>% summarise(nmale_none=as.vector(table(educa))[1],nmale_anyprimary=sum(as.vector(table(educa))[2:18]),nmale_anysecondary=sum(as.vector(table(educa))[7:18]))
female_educ_sex_dat <- sex_dat %>% filter(sex=='Female') %>% select(hhnum, educa) %>% group_by(hhnum) %>% summarise(nfemale_none=as.vector(table(educa))[1],nfemale_anyprimary=sum(as.vector(table(educa))[2:18]),nfemale_anysecondary=sum(as.vector(table(educa))[7:18]))
sex_dat <- merge(sex_dat, male_educ_sex_dat[, c("hhnum", "nmale_none", "nmale_anyprimary", "nmale_anysecondary")], by="hhnum", all.x=T)
sex_dat <- merge(sex_dat, female_educ_sex_dat[, c("hhnum", "nfemale_none", "nfemale_anyprimary", "nfemale_anysecondary")], by="hhnum", all.x=T)

dat <- merge(dat, sex_dat[, c("hhnum", "perc_male_0to4", "perc_male_5to14", "perc_male_15to19", "perc_male_20to34", "perc_male_35to54", "perc_male_55up",
                              "perc_female_0to4", "perc_female_5to14", "perc_female_15to19", "perc_female_20to34", "perc_female_35to54", "perc_female_55up",
                              "nmale_none", "nmale_anyprimary", "nmale_anysecondary", "nfemale_none", "nfemale_anyprimary", "nfemale_anysecondary")], by="hhnum", all.x=T)
dat$household_size <- dat$hhsizea
dat <- subset(dat, select = -c(pcode, sex, reltin, agey, educa, hhsizea))
dat <- dat[!duplicated(dat$hhnum), ]

## Gathering relevant asset-related data
hhdatalist <- get_hhdatalist(keepnums=TRUE)
assetdat=hhdatalist$SECH1.dta %>% select(hhnum,Assets,quant,esval,eslosp)

## Groups of assets in accordance with data description
group1=1:3
group2=11:20
group3=31:34
group4=41:44
group5=51:55
group6=61:65
#group7=71
#group8=101:103
group9=301:313
assetdat <- assetdat %>% mutate(G1=Assets %in% group1,G2=Assets %in% group2,G3=Assets %in% group3,G4=Assets %in% group4,G5=Assets %in% group5,G6=Assets %in% group6,G9=Assets %in% group9)
assetdat$quant.Assets.G1=assetdat$quant.Assets.G2=assetdat$quant.Assets.G3=assetdat$quant.Assets.G4=assetdat$quant.Assets.G5=assetdat$quant.Assets.G6=assetdat$quant.Assets.G9=0
assetdat$esval.Assets.G1=assetdat$esval.Assets.G2=assetdat$esval.Assets.G3=assetdat$esval.Assets.G4=assetdat$esval.Assets.G5=assetdat$esval.Assets.G6=assetdat$esval.Assets.G9=0
assetdat$eslosp.Assets.G1=assetdat$eslosp.Assets.G2=assetdat$eslosp.Assets.G3=assetdat$eslosp.Assets.G4=assetdat$eslosp.Assets.G5=assetdat$eslosp.Assets.G6=assetdat$eslosp.Assets.G9=0

## Go through each class of assets for each household
for (hhnumval in unique(assetdat$hhnum)){
  hhnumdat=assetdat %>% filter(hhnum==hhnumval)
  hhnumdat=hhnumdat %>% mutate(denom=1+eslosp/100)
  assetdat_ix=which(assetdat$hhnum == hhnumval)
  
  assetdat$quant.Assets.G1[assetdat_ix]= hhnumdat%>%filter(G1==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G2[assetdat_ix]= hhnumdat%>%filter(G2==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G3[assetdat_ix]= hhnumdat%>%filter(G3==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G4[assetdat_ix]= hhnumdat%>%filter(G4==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G5[assetdat_ix]= hhnumdat%>%filter(G5==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G6[assetdat_ix]= hhnumdat%>%filter(G6==T)%>%select(quant)%>%sum
  #assetdat$quant.Assets.G7[assetdat_ix]= hhnumdat%>%filter(G7==T)%>%select(quant)%>%sum
  #assetdat$quant.Assets.G8[assetdat_ix]= hhnumdat%>%filter(G8==T)%>%select(quant)%>%sum
  assetdat$quant.Assets.G9[assetdat_ix]= hhnumdat%>%filter(G9==T)%>%select(quant)%>%sum
  
  assetdat$esval.Assets.G1[assetdat_ix]= hhnumdat%>%filter(G1==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G2[assetdat_ix]= hhnumdat%>%filter(G2==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G3[assetdat_ix]= hhnumdat%>%filter(G3==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G4[assetdat_ix]= hhnumdat%>%filter(G4==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G5[assetdat_ix]= hhnumdat%>%filter(G5==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G6[assetdat_ix]= hhnumdat%>%filter(G6==T)%>%select(esval)%>%sum
  #assetdat$esval.Assets.G7[assetdat_ix]= hhnumdat%>%filter(G7==T)%>%select(esval)%>%sum
  #assetdat$esval.Assets.G8[assetdat_ix]= hhnumdat%>%filter(G8==T)%>%select(esval)%>%sum
  assetdat$esval.Assets.G9[assetdat_ix]= hhnumdat%>%filter(G9==T)%>%select(esval)%>%sum
  
  assetdat$eslosp.Assets.G1[assetdat_ix]= sum( hhnumdat%>%filter(G1==T)%>%select(esval) / hhnumdat%>%filter(G1==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G2[assetdat_ix]= sum( hhnumdat%>%filter(G2==T)%>%select(esval) / hhnumdat%>%filter(G2==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G3[assetdat_ix]= sum( hhnumdat%>%filter(G3==T)%>%select(esval) / hhnumdat%>%filter(G3==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G4[assetdat_ix]= sum( hhnumdat%>%filter(G4==T)%>%select(esval) / hhnumdat%>%filter(G4==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G5[assetdat_ix]= sum( hhnumdat%>%filter(G5==T)%>%select(esval) / hhnumdat%>%filter(G5==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G6[assetdat_ix]= sum( hhnumdat%>%filter(G6==T)%>%select(esval) / hhnumdat%>%filter(G6==T)%>%select(denom),na.rm=T )
  #assetdat$eslosp.Assets.G7[assetdat_ix]= sum( hhnumdat%>%filter(G7==T)%>%select(esval) / hhnumdat%>%filter(G7==T)%>%select(denom),na.rm=T )
  #assetdat$eslosp.Assets.G8[assetdat_ix]= sum( hhnumdat%>%filter(G8==T)%>%select(esval) / hhnumdat%>%filter(G8==T)%>%select(denom),na.rm=T )
  assetdat$eslosp.Assets.G9[assetdat_ix]= sum( hhnumdat%>%filter(G9==T)%>%select(esval) / hhnumdat%>%filter(G9==T)%>%select(denom),na.rm=T )
}

assetdat <- subset(assetdat, select = -c(Assets, quant, esval, eslosp, G1,G2,G3,G4,G5,G6,G9))
assetdat <- assetdat[!duplicated(assetdat$hhnum), ]
dat<-merge(dat, assetdat, by='hhnum')

## Replace NAs with 0s; reorder some vars
femalehead=dat$femalehead
agehead=dat$agehead
dat=subset(dat,select=-c(femalehead,agehead))
dat[is.na(dat)] <- 0
dat$femalehead=femalehead
dat$agehead=agehead

## Also include 1997 flood info
hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist['secj3.dta']
hhdatalist$secj3.dta <- hhdatalist$secj3.dta %>% select(hhnum,whousefa)
dat <- join_all(list(dat,hhdatalist$secj3.dta), by='hhnum', type='left')


## Collect outcome data
## Start with food consumption and food security outcomes
hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('expf4.dta','tothhexp.dta')]
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% select(hhnum,fdgroup,qcons,vpurch,calcons,procons)
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% filter(fdgroup %in% 1:11)
hhdatalist$tothhexp.dta <- hhdatalist$tothhexp.dta %>% select(-c(hrepair,pccal,pcprot,hhrent)) # remove redundancies
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% group_by(hhnum,fdgroup) %>% mutate(tot_consumed=sum(qcons),
                                                                                    tot_spent=sum(vpurch),
                                                                                    tot_cals=sum(calcons),
                                                                                    tot_prot=sum(procons))
hhdatalist$expf4.dta <- subset(hhdatalist$expf4.dta, select = -c(qcons, vpurch, calcons, procons))
hhdatalist$expf4.dta <- join_all(list(hhdatalist$expf4.dta, subset(dat, select=c(hhnum, household_size))), by='hhnum', type='left')
hhdatalist$expf4.dta <- join_all(list(hhdatalist$expf4.dta, subset(hhdatalist$tothhexp.dta, select=c(hhnum, tothhexp))), by='hhnum', type='left')
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% mutate(consumed_capita=tot_consumed/household_size,
                                                        cals_capita=tot_cals/household_size,
                                                        prot_capita=tot_prot/household_size,
                                                        budget_share=100*tot_spent/tothhexp)
hhdatalist$expf4.dta <- subset(hhdatalist$expf4.dta, select = -c(tot_consumed, tot_spent, tot_cals, tot_prot, household_size, tothhexp))
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% distinct(hhnum, fdgroup, .keep_all = TRUE)
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>% group_by(hhnum) %>% mutate(cals_capita=sum(cals_capita,na.rm=T),
                                                                            prot_capita=sum(prot_capita,na.rm=T))
hhdatalist$expf4.dta <- hhdatalist$expf4.dta %>%
  group_by(hhnum, fdgroup) %>%
  mutate(
    consumed_col = paste0("consumed_", fdgroup),
    budget_col = paste0("budget_", fdgroup)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = consumed_col,
    values_from = consumed_capita,
    values_fill = 0
  ) %>%
  pivot_wider(
    names_from = budget_col,
    values_from = budget_share,
    values_fill = 0
  ) %>%
  select(-fdgroup) %>%
  group_by(hhnum) %>%
  summarize(
    cals_capita = first(cals_capita),
    prot_capita = first(prot_capita),
    across(starts_with("consumed_"), sum),
    across(starts_with("budget_"), sum)
  ) %>%
  ungroup()
hhdatalist$expf4.dta[is.na(hhdatalist$expf4.dta )] <- 0
smalldata <- join_all(hhdatalist, by='hhnum', type='left')
dat<-merge(dat, smalldata, by='hhnum')
dat <- subset(dat, select=-c(hhsizea,hhsizer,fprinx.y,quint,aequint))

## Gender discrimination outcomes
hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('SECA1.dta','secl1.dta','SECL2.dta')]
hhdatalist$SECA1.dta <- hhdatalist$SECA1.dta %>% select(hhnum,pcode,sex,agey)
hhdatalist$secl1.dta <- hhdatalist$secl1.dta %>% select(hhnum,Foodtype,ingrr,ingrw,ingrm,ingrf,ingrmm,ingreg)
hhdatalist$secl1.dta[is.na(hhdatalist$secl1.dta )] <- 0
hhdatalist$secl1.dta$ingrtotal <- rowSums(subset(hhdatalist$secl1.dta, select=-c(hhnum,Foodtype)))
hhdatalist$secl1.dta$propr <- hhdatalist$secl1.dta$ingrr/hhdatalist$secl1.dta$ingrtotal
hhdatalist$secl1.dta$propw <- hhdatalist$secl1.dta$ingrw/hhdatalist$secl1.dta$ingrtotal
hhdatalist$secl1.dta$propm <- hhdatalist$secl1.dta$ingrm/hhdatalist$secl1.dta$ingrtotal
hhdatalist$secl1.dta$propf <- hhdatalist$secl1.dta$ingrf/hhdatalist$secl1.dta$ingrtotal
hhdatalist$secl1.dta$propmm <- hhdatalist$secl1.dta$ingrmm/hhdatalist$secl1.dta$ingrtotal
hhdatalist$secl1.dta$propeg <- hhdatalist$secl1.dta$ingreg/hhdatalist$secl1.dta$ingrtotal

## Omit guests and focus on household members
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta %>% select(hhnum,Foodtype,pcode,ammeal)
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta %>%
  group_by(hhnum, Foodtype, pcode) %>%
  summarise(sum_ammeal = sum(ammeal))
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta %>%
  group_by(hhnum, Foodtype) %>%
  mutate(tot_ammeal_everyone = sum(sum_ammeal,na.rm=T))
hhdatalist$SECL2.dta$sum_ammeal[is.na(hhdatalist$SECL2.dta$sum_ammeal)] <- 0
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta %>% filter(pcode < 51)
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta[-which(is.na(hhdatalist$SECL2.dta$Foodtype)),]
hhdatalist$SECL2.dta <- hhdatalist$SECL2.dta[-which(hhdatalist$SECL2.dta$tot_ammeal_everyone==0),]

smalldata <- join_all(list(hhdatalist$SECL2.dta, hhdatalist$SECA1.dta), by=c('hhnum','pcode'), type='left')
smalldata <- smalldata[-c(which(apply(smalldata,1,function(x) any(is.na(x))))),]
smalldata <- join_all(list(smalldata, hhdatalist$secl1.dta), by=c('hhnum','Foodtype'), type='left')
smalldata <- smalldata[-c(which(apply(smalldata,1,function(x) any(is.na(x))))),]
smalldata <- smalldata %>%
  mutate(
    amt_ingrr = sum_ammeal * propr,
    amt_ingrw = sum_ammeal * propw,
    amt_ingrm = sum_ammeal * propm,
    amt_ingrf = sum_ammeal * propf,
    amt_ingrmm = sum_ammeal * propmm,
    amt_ingreg = sum_ammeal * propeg
  )
smalldata <- subset(smalldata, select=-c(sum_ammeal,tot_ammeal_everyone,ingrr,ingrw,ingrm,ingrf,ingrmm,ingreg,propr,propw,propm,propf,propmm,propeg))
smalldata <- smalldata %>%
  group_by(hhnum, pcode) %>%
  mutate(
    sum_amt_ingrr = sum(amt_ingrr),
    sum_amt_ingrw = sum(amt_ingrw),
    sum_amt_ingrm = sum(amt_ingrm),
    sum_amt_ingrf = sum(amt_ingrf),
    sum_amt_ingrmm = sum(amt_ingrmm),
    sum_amt_ingreg = sum(amt_ingreg)
  ) 
smalldata <- subset(smalldata, select=-c(Foodtype,amt_ingrr,amt_ingrw,amt_ingrm,amt_ingrf,amt_ingrmm,amt_ingreg))
smalldata <- smalldata %>% filter(agey > 4)
smalldata <- smalldata %>% distinct(hhnum, pcode, .keep_all = TRUE)
smalldata <- smalldata %>%
  group_by(hhnum) %>%
  mutate(
    # icrs_f = (sum_amt_ingrf/(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)) / (sum(sum_amt_ingrf)/sum(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)),
    # icrs_m = (sum_amt_ingrm/(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)) / (sum(sum_amt_ingrm)/sum(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)),
    # icrs_mm = (sum_amt_ingrmm/(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)) / (sum(sum_amt_ingrmm)/sum(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)),
    # icrs_eg = (sum_amt_ingreg/(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm)) / (sum(sum_amt_ingreg)/sum(sum_amt_ingrr+sum_amt_ingrw+sum_amt_ingrm+sum_amt_ingrf+sum_amt_ingrmm))
    
    icrs_f = (sum_amt_ingrf/sum(sum_amt_ingrf,na.rm=T)) / mean(sum_amt_ingrf/sum(sum_amt_ingrf,na.rm=T),na.rm=T),
    # icrs_m = (sum_amt_ingrm/sum(sum_amt_ingrm,na.rm=T)) / mean(sum_amt_ingrm/sum(sum_amt_ingrm,na.rm=T),na.rm=T),
    # icrs_mm = (sum_amt_ingrmm/sum(sum_amt_ingrmm,na.rm=T)) / mean(sum_amt_ingrmm/sum(sum_amt_ingrmm,na.rm=T),na.rm=T),
    # icrs_eg = (sum_amt_ingreg/sum(sum_amt_ingreg,na.rm=T)) / mean(sum_amt_ingreg/sum(sum_amt_ingreg,na.rm=T),na.rm=T),
    icrs_mmmeg = (sum_amt_ingrm+sum_amt_ingrmm+sum_amt_ingreg/sum(sum_amt_ingrm+sum_amt_ingrmm+sum_amt_ingreg,na.rm=T)) / mean(sum_amt_ingrm+sum_amt_ingrmm+sum_amt_ingreg/sum(sum_amt_ingrm+sum_amt_ingrmm+sum_amt_ingreg,na.rm=T),na.rm=T)
  ) %>%
  ungroup()

smalldata <- smalldata %>%
  group_by(hhnum) %>%
  summarise(
    icrs_f_male = mean(icrs_f[sex == "Male"],na.rm=T),
    icrs_f_female = mean(icrs_f[sex == "Female"],na.rm=T),
    # icrs_m_male = mean(icrs_m[sex == "Male"],na.rm=T),
    # icrs_m_female = mean(icrs_m[sex == "Female"],na.rm=T),
    # icrs_mm_male = mean(icrs_mm[sex == "Male"],na.rm=T),
    # icrs_mm_female = mean(icrs_mm[sex == "Female"],na.rm=T),
    # icrs_eg_male = mean(icrs_eg[sex == "Male"],na.rm=T),
    # icrs_eg_female = mean(icrs_eg[sex == "Female"],na.rm=T),
    icrs_mmmeg_male = mean(icrs_mmmeg[sex == "Male"],na.rm=T),
    icrs_mmmeg_female = mean(icrs_mmmeg[sex == "Female"],na.rm=T),
    ###
    # ingrr_male = mean(sum_amt_ingrr[sex == "Male"]),
    # ingrr_female = mean(sum_amt_ingrr[sex == "Female"]),
    # ingrw_male = mean(sum_amt_ingrw[sex == "Male"]),
    # ingrw_female = mean(sum_amt_ingrw[sex == "Female"]),
    # ingrr_plus_ingrw_male = ingrr_male+ingrw_male,
    # ingrr_plus_ingrw_female = ingrr_female+ingrw_female,
    
    # ingrm_male = mean(sum_amt_ingrm[sex == "Male"]),
    # ingrm_female = mean(sum_amt_ingrm[sex == "Female"]),
    ingrf_male = mean(sum_amt_ingrf[sex == "Male"]),
    ingrf_female = mean(sum_amt_ingrf[sex == "Female"]),
    # ingrmm_male = mean(sum_amt_ingrmm[sex == "Male"]),
    # ingrmm_female = mean(sum_amt_ingrmm[sex == "Female"]),
    # ingreg_male = mean(sum_amt_ingreg[sex == "Male"]),
    # ingreg_female = mean(sum_amt_ingreg[sex == "Female"]),
    # ingrm_plus_ingrmm_plus_ingreg_male = ingrm_male+ingrmm_male+ingreg_male,
    # ingrm_plus_ingrmm_plus_ingreg_female = ingrm_female+ingrmm_female+ingreg_female
    ingrm_plus_ingrmm_plus_ingreg_male = mean(sum_amt_ingrm[sex == "Male"])+mean(sum_amt_ingrmm[sex == "Male"])+mean(sum_amt_ingreg[sex == "Male"]),
    ingrm_plus_ingrmm_plus_ingreg_female = mean(sum_amt_ingrm[sex == "Female"])+mean(sum_amt_ingrmm[sex == "Female"])+mean(sum_amt_ingreg[sex == "Female"])
  ) %>%
  ungroup()

## NOTE: SOME HHNUMS MAY ONLY HAVE ONE GENDER REPRESENTED OR OTHER MISSINGNESS

dat=join_all(list(dat,smalldata), by=c('hhnum'), type='left')

## Illness outcomes
hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('secj3.dta','SECJ45.dta','SECJ78.dta')]

hhdatalist$secj3.dta <- hhdatalist$secj3.dta %>% select(hhnum,whousef,whoused,ahomed,exprep)
hhdatalist$secj3.dta$ftwater <- hhdatalist$secj3.dta$whousef
hhdatalist$secj3.dta$dayswater <- hhdatalist$secj3.dta$whoused
hhdatalist$secj3.dta$dayswater[which(is.na(hhdatalist$secj3.dta$dayswater))] <- 0
hhdatalist$secj3.dta$daysawayhome <- hhdatalist$secj3.dta$ahomed
hhdatalist$secj3.dta$daysawayhome[which(is.na(hhdatalist$secj3.dta$daysawayhome))] <- 0
hhdatalist$secj3.dta$repaircost <- hhdatalist$secj3.dta$exprep
hhdatalist$secj3.dta <- subset(hhdatalist$secj3.dta, select=-c(whousef,whoused,ahomed,exprep))

hhdatalist$SECJ45.dta <- hhdatalist$SECJ45.dta %>% select(hhnum,latrine,swater,swatera,swaterb,awater,awatera,awaterb)
hhdatalist$SECJ45.dta$sanitarylatrine <- hhdatalist$SECJ45.dta$latrine == 'Pacca' | hhdatalist$SECJ45.dta$latrine == 'Pacca(W'
hhdatalist$SECJ45.dta$tubewell_drinking <- hhdatalist$SECJ45.dta$swater == 1
hhdatalist$SECJ45.dta$tubewell_cooking <- hhdatalist$SECJ45.dta$swatera == 1
hhdatalist$SECJ45.dta$tubewell_washing <- hhdatalist$SECJ45.dta$swaterb == 'Tube wel'
hhdatalist$SECJ45.dta$dist_drinking <- hhdatalist$SECJ45.dta$awater
hhdatalist$SECJ45.dta$dist_cooking <- hhdatalist$SECJ45.dta$awatera
hhdatalist$SECJ45.dta$dist_washing <- hhdatalist$SECJ45.dta$awaterb
hhdatalist$SECJ45.dta <- subset(hhdatalist$SECJ45.dta, select=-c(latrine,swater,swatera,swaterb,awater,awatera,awaterb))

hhdatalist$SECJ78.dta <- hhdatalist$SECJ78.dta %>% group_by(hhnum) %>% mutate(prod1=twater*timtr/60, prod2=twaterb*timtrb/60, prod3=twaterc*timtrc/60)
hhdatalist$SECJ78.dta <- hhdatalist$SECJ78.dta %>% group_by(hhnum) %>% mutate(timewater=mean(c(prod1,prod2,prod3),na.rm=T))
hhdatalist$SECJ78.dta <- hhdatalist$SECJ78.dta %>% group_by(hhnum) %>% mutate(prod1=twater*amcol, prod2=twaterb*amcola, prod3=twaterc*amcolb)
hhdatalist$SECJ78.dta <- hhdatalist$SECJ78.dta %>% group_by(hhnum) %>% mutate(literswater=mean(c(prod1,prod2,prod3),na.rm=T))
hhdatalist$SECJ78.dta <- subset(hhdatalist$SECJ78.dta, select=-c(prod1,prod2,prod3))
hhdatalist$SECJ78.dta$electric <- ifelse(hhdatalist$SECJ78.dta$elesup=='No',0,1)
hhdatalist$SECJ78.dta <- subset(hhdatalist$SECJ78.dta, select=-c(midj78,twater,timtr,amcol,midj78a,twaterb,timtrb,amcola,midj78b,twaterc,timtrc,amcolb,elesup))

dat=join_all(list(dat,hhdatalist$secj3.dta), by=c('hhnum'), type='left')
dat=join_all(list(dat,hhdatalist$SECJ45.dta), by=c('hhnum'), type='left')
dat=join_all(list(dat,hhdatalist$SECJ78.dta), by=c('hhnum'), type='left')
repeatedhh=with(rle(dat$hhnum), values[lengths > 1])
dat=dat[-c(562,567,594),] # repeated HHs with less information from SECJ78.dta


hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('SECA1.dta','secp1.dta')]

hhdatalist$SECA1.dta <- hhdatalist$SECA1.dta %>% select(hhnum,pcode,agey,sex)
hhdatalist$secp1.dta <- hhdatalist$secp1.dta %>% select(hhnum,pcode,sick)
smalldata <- join_all(list(hhdatalist$SECA1.dta,hhdatalist$secp1.dta), by=c('hhnum','pcode'), type='left')
smalldata <- join_all(list(subset(dat,select=c(hhnum,household_size)),smalldata), by=c('hhnum'), type='left')
smalldata <- smalldata[-which(is.na(smalldata$sick)),]

## Calculate the proportion of individuals with ANY illness in different age groups
any_illness_data <- smalldata %>%
  mutate(age_group = case_when(
    agey <= 5 ~ "0-5yrs",
    # agey <= 10 ~ "6-10yrs",
    # agey <= 18 ~ "11-18yrs",
    agey <= 18 ~ "6-18yrs",
    TRUE ~ "over18yrs"
  )) %>%
  group_by(hhnum, age_group) %>%
  summarize(any_illness_prop = any(!is.na(sick))) %>%
  pivot_wider(names_from = age_group, values_from = any_illness_prop, names_prefix = "any_ill_") %>%
  ungroup()

## Calculate the proportion of individuals with ANY illness in different sex groups
any_illness_data_sex <- smalldata %>%
  group_by(hhnum, sex) %>%
  summarize(any_illness_prop_sex = any(!is.na(sick))) %>%
  pivot_wider(names_from = sex, values_from = any_illness_prop_sex, names_prefix = "any_ill_") %>%
  ungroup()

## Calculate the proportion of individuals with specific diseases
disease_data <- smalldata %>%
  group_by(hhnum, sick) %>%
  summarize(disease_incidence = ifelse(sum(!is.na(sick))>0,1,0) ) %>%
  pivot_wider(names_from = sick, values_from = disease_incidence, names_prefix = "ill_prop_") %>%
  ungroup()
disease_data <- disease_data %>% select(c("hhnum","ill_prop_1","ill_prop_5","ill_prop_6","ill_prop_7","ill_prop_8"))
disease_data[is.na(disease_data)]=0
disease_data <- disease_data %>%
  group_by(hhnum) %>%
  summarise(any_gastro=max(ill_prop_7,ill_prop_8),
            any_resp=max(ill_prop_5,ill_prop_6),
            any_fever=ill_prop_1) %>%
  ungroup()
smalldata <- left_join(any_illness_data, disease_data, by = "hhnum")

smalldata[is.na(smalldata)] <- 0
illness_columns <- grep("ill_", names(smalldata), value = TRUE)
dat <- left_join(dat, smalldata, by = "hhnum") %>%
  mutate_at(vars(all_of(illness_columns)), ~ ifelse(is.na(.), 0, .))
any_illness_data_sex[is.na(any_illness_data_sex)] <- 0
illness_columns <- grep("ill_", names(any_illness_data_sex), value = TRUE)
dat <- left_join(dat, any_illness_data_sex, by = "hhnum") %>%
  mutate_at(vars(all_of(illness_columns)), ~ ifelse(is.na(.), 0, .))



hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('SECA1.dta','secp1.dta')]

hhdatalist$SECA1.dta <- hhdatalist$SECA1.dta %>% select(hhnum,pcode,agey)
hhdatalist$secp1.dta <- hhdatalist$secp1.dta %>% select(hhnum,pcode,days)
smalldata <- join_all(list(hhdatalist$SECA1.dta,hhdatalist$secp1.dta), by=c('hhnum','pcode'), type='left')
smalldata <- join_all(list(subset(dat,select=c(hhnum,household_size)),smalldata), by=c('hhnum'), type='left')
smalldata$days[which(is.na(smalldata$days))] = 0
smalldata <- smalldata %>% group_by(hhnum,pcode) %>% mutate(daystot=sum(days))
smalldata$daystot[which(smalldata$daystot > 17)] = 17
smalldata <- smalldata %>%
  mutate(
    transdays = daystot
  )
smalldata <- smalldata %>%
  group_by(hhnum) %>%
  summarise(
    transdays_pc = mean(transdays)
  )
dat <- left_join(dat, smalldata, by = "hhnum")
