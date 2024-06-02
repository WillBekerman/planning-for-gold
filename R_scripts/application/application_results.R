
## Load libraries and auxillary functions
library(tidyverse)
library(sensitivitymult)
library(DOS2)
library(glmnet)

## Computes sensitivity value for data.
## Inputs: mpdifs, a, alt.
get_sensitivity_val <- function(mpdifs,a=a,alt=alt) {
  senfx <- function(mpdifs,g,a,alt){
    senWilcox(mpdifs, gamma=g, alpha=a, alternative=alt)$pval-a
  }
  g<-1
  try(g<-uniroot(senfx,c(1,1e10),mpdifs=mpdifs,a=a,alt=alt)$root,silent=TRUE)
  kappa<-g/(1+g)
  return(kappa)
}

## Sensivitity Analysis for McNemar's Test Statistic
## Let D be the number of discordant pairs
## Tobs be the number of discordant pairs in which treated unit has a 1
## and Gamma be the sensitivity parameter exp(gamma)
sens.analysis.mcnemar=function(D,Tobs,Gamma){
  p.positive=Gamma/(1+Gamma);
  p.negative=1/(1+Gamma);
  lowerbound=1-pbinom(Tobs-1,D,p.negative);
  upperbound=1-pbinom(Tobs-1,D,p.positive);
  list(lowerbound=lowerbound,upperbound=upperbound);
}

get_sensitivity_val_binary <- function(dat,a=a,alt=alt) {
  ## Assumes dat is given as control,treated,control,treated...
  if (alt=='greater'){
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    senfx <- function(g,a,alt){
      sens.analysis.mcnemar(D,Tobs,g)$upperbound-a
    }
    g<-1
    try(g<-uniroot(senfx,c(1,1e10),a=a,alt=alt)$root,silent=TRUE)
    kappa<-g/(1+g)
    return(kappa)
  }
  if (alt=='less'){
    dat<-!dat
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    senfx <- function(g,a,alt){
      sens.analysis.mcnemar(D,Tobs,g)$upperbound-a
    }
    g<-1
    try(g<-uniroot(senfx,c(1,1e10),a=a,alt=alt)$root,silent=TRUE)
    kappa<-g/(1+g)
    return(kappa)
  }
  if (alt=='twosided'){
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    
    dat<-!dat
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    Dnew=sum(abs(mpdifs)==1)
    Tobsnew=sum(mpdifs==1)
    
    senfx <- function(g,a,alt){
      2*min(sens.analysis.mcnemar(D,Tobs,g)$upperbound,
            sens.analysis.mcnemar(Dnew,Tobsnew,g)$upperbound)-a
    }
    g<-1
    try(g<-uniroot(senfx,c(1,1e10),a=a,alt=alt)$root,silent=TRUE)
    kappa<-g/(1+g)
    return(kappa)
  }
}

sensitivity_analysis_binary <- function(dat,gamma=gamma,a=a,alt=alt) {
  ## Assumes dat is given as control,treated,control,treated...
  if (alt=='greater'){
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    p=sens.analysis.mcnemar(D,Tobs,gamma)$upperbound
  }
  if (alt=='less'){
    dat<-!dat
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    p=sens.analysis.mcnemar(D,Tobs,gamma)$upperbound
  }
  if(alt=='twosided'){
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    p=sens.analysis.mcnemar(D,Tobs,gamma)$upperbound
    
    dat<-!dat
    dat_filtered <- numeric()
    controlix <- seq(1,length(dat),2)
    treatedix <- controlix + 1
    for (set in 1:(length(dat)/2)){
      if (any(is.na(c(dat[controlix[set]],dat[treatedix[set]])))) { next }
      dat_filtered <- c(dat_filtered,dat[controlix[set]],dat[treatedix[set]])
    }
    mpdifs <- diff(dat_filtered)[seq(1,length(dat_filtered),2)]
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
    pnew=sens.analysis.mcnemar(D,Tobs,gamma)$upperbound
    
    return(2*min(p,pnew)) #twosided
  }
  return(p)
}

is_binary_outcome = function(x) { return(all(is.na(x) | x==0 | abs(x)==1 )) }


## Function for split-sampling framework to determine and test outcomes.
## Assume matched pairs framework.
## Inputs: data_outcomes, control_idx, treated_idx, planning_sample_prop, alpha,
## gamma, nboot.
## Outputs: outcomes_tested, outcomes_tested_oneshot, outcomes_rejected,
## outcomes_rejected_oneshot, outcomes_rejected_bonferroni,
## planning_sample_proportion, gamma.
run_method_matchedprs <- function(data_outcomes,control_idx,treated_idx,
                                  planning_sample_prop=0.10,alpha=0.05,
                                  gamma=1,nboot=100,seed=0){
  
  ## Set seed for reproducibility
  set.seed(seed)
  
  num_outcomes <- ncol(data_outcomes)
  num_matched_sets <- length(treated_idx)
  
  ## Create objects to return
  outcomes_tested <- outcomes_tested_oneshot <- outcomes_rejected <- outcomes_rejected_oneshot <- outcomes_rejected_bonf <- logical(length=num_outcomes)
  
  ## Divide data into planning and analysis samples
  planning_sample_size_sets <- floor(planning_sample_prop*num_matched_sets)
  planning_sample_size_subjects <- planning_sample_size_sets*2
  ix_subjects_treated_planning <- sample(x=treated_idx,size=planning_sample_size_sets)
  ix_subjects_control_planning <- control_idx[match(ix_subjects_treated_planning,treated_idx),]
  ix_subjects_treated_analysis <- treated_idx[-match(ix_subjects_treated_planning,treated_idx)]
  ix_subjects_control_analysis <- control_idx[match(ix_subjects_treated_analysis,treated_idx),]
  
  ## Some relevant parameters
  n1 <- planning_sample_size_sets
  n2 <- num_matched_sets-planning_sample_size_sets
  
  ## Get vector for alternatives to test using planning data + run outcome regressions
  alt_vec <- character(num_outcomes)
  binary_ix <- which( apply(data_outcomes,2,is_binary_outcome) )
  
  for (outcome in 1:num_outcomes) {
    y_plan <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])))) {
        next
      }
      y_plan <- c(y_plan,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])
    }
    mpdifs <- diff(y_plan)[seq(1,length(y_plan),2)]
    alt_vec[outcome] <- ifelse(mean(mpdifs)>=0,"greater","less")
  }
  
  
  ## BOOTSTRAP METHOD ##
  # Planning stage
  ixsetlist <- list()
  for (boot in 1:nboot){
    ixset <- sample(x=1:planning_sample_size_sets, size=planning_sample_size_sets, replace=T)
    ixsetlist[[boot]] <- c(as.matrix(ix_subjects_control_planning)[ixset,],ix_subjects_treated_planning[ixset])
  }
  sigma_hat_vec <- kappahatvec_whole <- numeric(length = num_outcomes)
  for (outcome in 1:num_outcomes){
    ## Compute sigmahat
    y <- as.matrix(data_outcomes)[,outcome]
    kappahatvec <- numeric()
    for (boot in 1:nboot){
      y_boot <- y[ixsetlist[[boot]]]
      mpdifs <- diff(y_boot)[seq(1,length(y_boot),2)]
      if (!(all(is.na(mpdifs))) & !is_binary_outcome(y_boot)){ kappahatvec <- c(kappahatvec, get_sensitivity_val(mpdifs[!is.na(mpdifs)],a=alpha,alt=alt_vec[outcome])) }
      if (!(all(is.na(mpdifs))) & is_binary_outcome(y_boot)){ kappahatvec <- c(kappahatvec, get_sensitivity_val_binary(y_boot,a=alpha,alt=alt_vec[outcome])) }
    }
    sigma_hat_vec[outcome] <- sqrt(var(kappahatvec)*n1)
    ## Compute kappahat on the whole planning sample
    y_boot <- numeric()
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])))) { next }
      y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])
    }
    if (!is_binary_outcome(y_boot)){ kappahatvec_whole[outcome] <- get_sensitivity_val(diff(y_boot)[seq(1,length(y_boot),2)], a=alpha, alt=alt_vec[outcome]) }
    if (is_binary_outcome(y_boot)){ kappahatvec_whole[outcome] <- get_sensitivity_val_binary(y_boot, a=alpha, alt=alt_vec[outcome]) }
  }
  aprimegrid <- alpha/(1:num_outcomes)
  tested_list <- vector("list", num_outcomes)
  difvec <- numeric(length = num_outcomes)
  for (aprime_idx in 1:num_outcomes){
    aprime <- aprimegrid[aprime_idx]
    our_tested_ix = logical(length = num_outcomes)
    for (outcome in 1:num_outcomes){
      sigma_hat <- sigma_hat_vec[outcome]
      kappahat <- kappahatvec_whole[outcome]
      sigma_r <- sigma_hat / sqrt(planning_sample_prop * (1 - planning_sample_prop))
      c1_hat <- -sqrt(4/3) * qnorm(1 - alpha) * sqrt(kappahat * (1 - kappahat))
      c2_hat <- -sqrt(4/3) * qnorm(1 - aprime) * sqrt(kappahat * (1 - kappahat))
      mu_r <- c1_hat / sqrt(planning_sample_prop) - c2_hat / sqrt(1 - planning_sample_prop)
      our_tested_ix[outcome] <- kappahat > gamma / (1 + gamma) + (mu_r - qnorm(1 - alpha) * sigma_r) / (sqrt(n1 + n2))
    }
    tested_list[[aprime_idx]] <- our_tested_ix
    difvec[aprime_idx] <- abs(alpha-aprime*sum(our_tested_ix, na.rm=TRUE))
  }
  aprime <- aprimegrid[which.min(difvec)]
  outcomes_tested <- tested_list[[which.min(difvec)]]
  outcomes_tested[is.na(outcomes_tested)] <- F
  if (sum(outcomes_tested) > 0){
    weights <- alpha/sum(outcomes_tested)
  }
  # Analysis stage
  for (outcome in 1:num_outcomes){
    if (outcomes_tested[outcome]){
      y_boot <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n2){
        if (any(is.na(c(y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])))) { next }
        y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])
      }
      if (!is_binary_outcome(y_boot)){reject <- senWilcox(d=diff(y_boot)[seq(1,length(y_boot),2)],gamma=gamma,alpha=alpha,alternative=alt_vec[outcome])$pval < alpha/sum(outcomes_tested)}#weights[outcome]
      if (is_binary_outcome(y_boot)){reject <- sensitivity_analysis_binary(dat=y_boot,gamma=gamma,a=alpha,alt=alt_vec[outcome]) < alpha/sum(outcomes_tested)}#weights[outcome]
      if (!is.na(reject)&reject) { outcomes_rejected[outcome] <- T }
    }
  }
  
  
  ## ONE-SHOT METHOD ##
  # Planning stage
  for (outcome in 1:num_outcomes){
    y_oneshot <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])))) { next }
      y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])
    }
    if (!is_binary_outcome(y_oneshot)){ outcomes_tested_oneshot[outcome] <- senWilcox(d=diff(y_oneshot)[seq(1,length(y_oneshot),2)],gamma=gamma,alpha=alpha,alternative=alt_vec[outcome])$pval < alpha}
    if (is_binary_outcome(y_oneshot)){ outcomes_tested_oneshot[outcome] <- sensitivity_analysis_binary(dat=y_oneshot,gamma=gamma,a=alpha,alt=alt_vec[outcome]) < alpha}
  }
  outcomes_tested_oneshot[is.na(outcomes_tested_oneshot)] <- F
  # Analysis stage
  for (outcome in 1:num_outcomes){
    if (outcomes_tested_oneshot[outcome]){
      y_oneshot <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n2){
        if (any(is.na(c(y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])))) { next }
        y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])
      }
      if (!is_binary_outcome(y_oneshot)){ reject <- senWilcox(d=diff(y_oneshot)[seq(1,length(y_oneshot),2)],gamma=gamma,alpha=alpha,alternative=alt_vec[outcome])$pval < alpha/sum(outcomes_tested_oneshot)}
      if (is_binary_outcome(y_oneshot)){ reject <- sensitivity_analysis_binary(dat=y_oneshot,gamma=gamma,a=alpha,alt=alt_vec[outcome]) < alpha/sum(outcomes_tested_oneshot)}
      if (!is.na(reject)&reject) { outcomes_rejected_oneshot[outcome] <- T }
    }
  }
  
  
  ## BONFERRONI CORRECTION ##
  for (outcome in 1:num_outcomes){
    y_bonf <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:num_matched_sets){
      if (any(is.na(c(y[as.matrix(control_idx)[set,]],y[treated_idx[set]])))) { next }
      y_bonf <- c(y_bonf,y[as.matrix(control_idx)[set,]],y[treated_idx[set]])
    }
    if (!is_binary_outcome(y_bonf)){ reject <- senWilcox(d=diff(y_bonf)[seq(1,length(y_bonf),2)],gamma=gamma,alpha=alpha,alternative='twosided')$pval < alpha/num_outcomes}
    if (is_binary_outcome(y_bonf)){ reject <- sensitivity_analysis_binary(dat=y_bonf,gamma=gamma,a=alpha,alt='twosided') < alpha/num_outcomes}
    if (!is.na(reject)&reject) { outcomes_rejected_bonf[outcome] <- T }
  }
  
  return(list(outcomes_tested=outcomes_tested,
              outcomes_tested_oneshot=outcomes_tested_oneshot,
              outcomes_rejected=outcomes_rejected,
              outcomes_rejected_oneshot=outcomes_rejected_oneshot,
              outcomes_rejected_bonferroni=outcomes_rejected_bonf,
              planning_sample_proportion=planning_sample_prop,
              gamma=gamma,
              num_outcomes_tested=sum(outcomes_tested),
              num_outcomes_tested_oneshot=sum(outcomes_tested_oneshot),
              rejected_boot=names(data_outcomes)[which(outcomes_rejected)],
              rejected_oneshot=names(data_outcomes)[which(outcomes_rejected_oneshot)],
              rejected_bonf=names(data_outcomes)[which(outcomes_rejected_bonf)]))
  
}


res1<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                            treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                            alpha=0.1,gamma=1+1e-5,nboot=250)

res1.25<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                               treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                               alpha=0.1,gamma=1.25,nboot=250)

res1.50<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                               treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                               alpha=0.1,gamma=1.5,nboot=250)

res1.75<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                               treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                               alpha=0.1,gamma=1.75,nboot=250)

res2<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                            treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                            alpha=0.1,gamma=2,nboot=250)

res2.5<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                              treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                              alpha=0.1,gamma=2.5,nboot=250)

res3<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                            treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                            alpha=0.1,gamma=3,nboot=250)

res4<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                            treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                            alpha=0.1,gamma=4,nboot=250)

res6<-run_method_matchedprs(data_outcomes,control_idx=ix_subjects_control,
                            treated_idx=ix_subjects_treated,planning_sample_prop=0.2,
                            alpha=0.1,gamma=6,nboot=250)



# Create a list of your results objects
results <- list(res1, res1.25, res1.50, res1.75, res2, res2.5, res3, res4, res6)
gammas <- c(1,1.25,1.5,1.75,2,2.5,3,4,6)

allvars <- res1$rejected_boot

ix=1
for (var in allvars){
  cat('var: ', allvars[ix],'\n\n')
  ix=ix+1
  
  ixgam=1
  for (gamma in gammas){
    if( var %in% results[[ixgam]]$rejected_bonf ) {cat('Gamma= ', gamma, ' Bonf');cat('\n')}
    if( var %in% results[[ixgam]]$rejected_oneshot ) {cat('Gamma= ', gamma, ' Naive');cat('\n')}
    if( var %in% results[[ixgam]]$rejected_boot ) {cat('Gamma= ', gamma, ' S-V');cat('\n')}
    
    ixgam=ixgam+1
  }
  cat('\n\n')
}

