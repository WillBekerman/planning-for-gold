
## Assumes we're in ~/bangladesh_application/data/Data household level/Rnd3"
# setwd('../Rnd3')

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
  
  smalldatalist <- datalist
  smalldatalist <- within(smalldatalist, rm(list=datasets_no_hhnum))
  
  smalldatalist
}


hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('anthvar.dta')]
hhdatalist$anthvar.dta <- hhdatalist$anthvar.dta %>% select(hhnum,pcode,agey,sex,haz,waz,whz,bmi,waste,sevwaste,stunt,sevstunt,undwt,sevundwt,preg)

dat1 <- hhdatalist$anthvar.dta
dat1 <- dat1 %>% group_by(hhnum) %>% summarise(avg_haz_all=mean(haz,na.rm=T),avg_waz_all=mean(waz,na.rm=T),avg_whz_all=mean(whz,na.rm=T))
dat2 <- dat2 %>% group_by(hhnum) %>% summarise(bmi_all=mean(bmi,na.rm=T))
smalldata <- join_all(list(dat1,dat2), by=c('hhnum'), type='right')
colnames(smalldata) <- paste('rd3', colnames(smalldata), sep = '_')
colnames(smalldata)[1] <- 'hhnum'
dat <- join_all(list(dat,smalldata), by=c('hhnum'), type='left')


hhdatalist <- get_hhdatalist()
hhdatalist <- hhdatalist[c('anthvar.dta')]
hhdatalist$anthvar.dta <- hhdatalist$anthvar.dta %>% select(hhnum,pcode,agey,sex,haz,waz,whz,bmi,waste,sevwaste,stunt,sevstunt,undwt,sevundwt,preg)

preschooldat <- hhdatalist$anthvar.dta %>% filter(agey >= 0 & agey <= 5)
preschooldat <- preschooldat %>% group_by(hhnum) %>% summarise(avg_haz_prek=mean(haz,na.rm=T),avg_waz_prek=mean(waz,na.rm=T),avg_whz_prek=mean(whz,na.rm=T))
womendat <- hhdatalist$anthvar.dta %>% filter( (sex==2 & agey >= 10 & agey <= 18) | (sex==2 & agey > 18 & preg!=1) )
womendat <- womendat %>% group_by(hhnum) %>% summarise(bmi_all_w=mean(bmi,na.rm=T))
smalldata <- join_all(list(preschooldat,womendat), by=c('hhnum'), type='right')
colnames(smalldata) <- paste('rd3', colnames(smalldata), sep = '_')
colnames(smalldata)[1] <- 'hhnum'
dat <- join_all(list(dat,smalldata), by=c('hhnum'), type='left')
