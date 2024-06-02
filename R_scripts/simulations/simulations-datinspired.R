################################# SIMULATIONS #################################  
###############################################################################

# Load libraries and auxillary functions
library(tidyverse)
library(MASS)
library(caret)
library(scales)
library(cowplot)
library(sensitivitymult)
library(DOS2)
library(optmatch)
options("optmatch_max_problem_size" = Inf)

#' Computes sensitivity value for data.
#'
#' @param mpdifs Differences between matched pairs in data.
#' @param a Alpha level.
#' @return Sensitivity value.
get_sensitivity_val <- function(mpdifs,a=a) {
  senfx <- function(mpdifs,g,a){
    senWilcox(mpdifs, gamma=g, alpha=a)$pval - a
  }
  g <- 1
  try(g <- uniroot(senfx,c(1,1e10),mpdifs=mpdifs,a=a)$root,silent=TRUE)
  kappa <- g/(1+g)
  return(kappa)
}

#' Runs propensity score matching on data.
#'
#' @param data_matching .
#' @param data_outcomes .
#' @param treated Vector of treated indices.
#' @return .
run_PSM <- function(data_matching,data_outcomes,treated){
  ## Function for computing rank based Mahalanobis distance.
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
  addcaliper=function(dmat,z,logitp,calipersd=.01,penalty=1000){
    # Pooled within group standard devation
    sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
    adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
    adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
    dmat=dmat+adif*penalty
    dmat
  }
  datatemp=as.data.frame(cbind(data_matching,treated))
  names(datatemp)[ncol(datatemp)] <- 'treated'
  ## Propensity score model
  propscore.model=glm(treated ~ ., data = datatemp, family = "binomial")
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
  distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=.01)
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
  sd.bf=stand.diff.before
  sd.af=stand.diff.after
  standardized_difs=cbind(sd.bf,sd.af)
  if(no.control.lack.overlap+no.treated.lack.overlap>0){
    data_outcomes <- rbind(data_outcomes[-which.remove,],data_outcomes[which.remove,])
  }
  ix_subjects_treated=treated.subject.index
  ix_subjects_control=matched.control.subject.index
  return(list(standardized_difs=standardized_difs,
              num_matched_sets=num_matched_sets,
              data_outcomes=data_outcomes,
              ix_subjects_treated=ix_subjects_treated,
              ix_subjects_control=ix_subjects_control))
}




############################################################################

#' Runs simulation for split-sampling framework and conducts tests 
#'
#' @param num_subjects Total sample size
#' @param num_outcomes Total number of outcomes
#' @param true_theta Constant value of nonzero theta, greater than zero
#' @param planning_sample_prop Proportion of sample to use in planning stage
#' @param alpha Cutoff for statistical significance
#' @param uc_dgp Should we generate artificial unmeasured confounding in our data
#' @param gamma Level of unmeasured confounding, equal to one implies none
#' @param nsim Number of simulations
#' @param nboot Number of bootstraps
#' @return List of simulation metrics returned.

run_sim <- function(num_subjects=757, num_outcomes=93, true_theta=1,
                    planning_sample_prop=0.2, alpha=0.05,
                    uc_dgp=TRUE, gamma=1, nsim=1000, nboot=250,
                    num_effects=5,num_covariates=33){
  
  # create theta_vec: vector of our signals for each outcome
  ix_important <- sample(x=1:num_outcomes, size=num_effects) # setting true signals
  theta_vec <- rep(0,num_outcomes)
  theta_vec[ix_important] <- true_theta
  thetavecold=theta_vec
  
  # create objects to return after simulation is complete
  outcomes_tested_list <- outcomes_rejected_list <- vector("list",nsim)
  sensitivity_vec <- sensitivity_vec_old <- bonf_sensitivity_vec <- 
    type1error_vec <- specificity_vec <- type2error_vec <- 
    num_outcomes_tested_vec <- num_outcomes_tested_vec_old <- rep(NA,nsim)
  
  for (sim in 1:nsim){
    
    #browser()
    
    outcomes <- 1:num_outcomes
    covariates <- matrix(runif(num_subjects*num_covariates),nrow=num_subjects,ncol=num_covariates)
    #covariates <- matrix(truncnorm::rtruncnorm(num_subjects*num_covariates, a=0, b=1, mean=0.5, sd=.25),nrow=num_subjects,ncol=num_covariates)
    
    ## Execute appropriate DGP ##
    if (uc_dgp){
      alpha0grid <- seq(from = -10, to = 10, by = 0.1)
      alpha0res <- list()
      beta_values <- list()
      mu_values <- list()
      U_values <- list()
      trt_assignment_values <- list()
      desired_trt_ratio_rv <- 0.71
      for (alpha0val in alpha0grid) {
        alpha0 <- alpha0val#rnorm(n = num_subjects)#
        beta <- rnorm(n = num_outcomes*num_covariates, mean = 1, sd = 1)
        mu <- rnorm(n = num_covariates, mean = 0, sd = 1)
        U <- rnorm(n = num_subjects, sd = 1 + sin(3 * covariates[, 1]) / 2)
        Zvec <- numeric(length = num_subjects)
        for (i in 1:num_subjects) {
          Zvec[i] <- exp(alpha0 + as.numeric(t(covariates[i, ]) %*% mu) - log(gamma) * (U[i] > 0)) /
            (1 + exp(alpha0 + as.numeric(t(covariates[i, ]) %*% mu) - log(gamma) * (U[i] > 0)))
        }
        trt_assignment_vec <- rbinom(n = num_subjects, size = 1, prob = Zvec)
        alpha0res[[length(alpha0res) + 1]] <- abs(sum(trt_assignment_vec)/length(trt_assignment_vec) - desired_trt_ratio_rv)
        beta_values[[length(beta_values) + 1]] <- beta
        mu_values[[length(mu_values) + 1]] <- mu
        U_values[[length(U_values) + 1]] <- U
        trt_assignment_values[[length(trt_assignment_values) + 1]] <- trt_assignment_vec
      }
      min_index <- which.min(unlist(alpha0res))
      alpha0 <- alpha0grid[min_index]
      beta <- beta_values[[min_index]]
      beta <- matrix(beta,nrow=num_outcomes,ncol=num_covariates)
      mu <- mu_values[[min_index]]
      U <- U_values[[min_index]]
      trt_assignment_vec <- trt_assignment_values[[min_index]]
      y_control_big <- matrix(NA,nrow=num_subjects,ncol=num_outcomes)
      #browser()
      for (outcome in outcomes){
        if (outcome %in% ix_important){ y_control_outcome <- apply(covariates,1,function(x)as.numeric(t(beta[outcome,]) %*% x))+rnorm(n=num_subjects,sd=0.5) }
        else { y_control_outcome <- apply(covariates,1,function(x)as.numeric(t(beta[outcome,]) %*% x))+U+rnorm(n=num_subjects,sd=0.5) }
        y_control_big[,outcome] <- y_control_outcome
      }
      # # adjust null outcomes
      # for (i in 1:(2*r)){
      #   y_control_big[i,which(theta_vec==0)] <- qnorm(1-Uvec[i])+i/200+y_control_big[i,which(theta_vec==0)]
      # }
      matchres <- run_PSM(data_matching=covariates,data_outcomes=y_control_big,treated=trt_assignment_vec)
      standardized_difs=matchres$standardized_difs
      num_matched_sets=matchres$num_matched_sets
      data_outcomes=matchres$data_outcomes
      treated_idx=matchres$ix_subjects_treated
      control_idx=matchres$ix_subjects_control
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
    }  
    else if (!uc_dgp) {
      
      alpha0grid <- seq(from = -10, to = 10, by = 0.1)
      alpha0res <- list()
      beta_values <- list()
      mu_values <- list()
      U_values <- list()
      trt_assignment_values <- list()
      desired_trt_ratio_rv <- 0.71
      for (alpha0val in alpha0grid) {
        alpha0 <- alpha0val#rnorm(n = num_subjects)#
        beta <- rnorm(n = num_outcomes*num_covariates, mean = 1, sd = 1)
        mu <- rnorm(n = num_covariates, mean = 0, sd = 1)
        U <- rnorm(n = num_subjects, sd = 1 + sin(3 * covariates[, 1]) / 2)
        Zvec <- numeric(length = num_subjects)
        for (i in 1:num_subjects) {
          Zvec[i] <- exp(alpha0 + as.numeric(t(covariates[i, ]) %*% mu)) /
            (1 + exp(alpha0 + as.numeric(t(covariates[i, ]) %*% mu)))
        }
        trt_assignment_vec <- rbinom(n = num_subjects, size = 1, prob = Zvec)
        alpha0res[[length(alpha0res) + 1]] <- abs(sum(trt_assignment_vec)/length(trt_assignment_vec) - desired_trt_ratio_rv)
        beta_values[[length(beta_values) + 1]] <- beta
        mu_values[[length(mu_values) + 1]] <- mu
        U_values[[length(U_values) + 1]] <- U
        trt_assignment_values[[length(trt_assignment_values) + 1]] <- trt_assignment_vec
      }
      min_index <- which.min(unlist(alpha0res))
      alpha0 <- alpha0grid[min_index]
      beta <- beta_values[[min_index]]
      beta <- matrix(beta,nrow=num_outcomes,ncol=num_covariates)
      mu <- mu_values[[min_index]]
      U <- U_values[[min_index]]
      trt_assignment_vec <- trt_assignment_values[[min_index]]
      y_control_big <- matrix(NA,nrow=num_subjects,ncol=num_outcomes)
      #browser()
      for (outcome in outcomes){
        if (outcome %in% ix_important){ y_control_outcome <- apply(covariates,1,function(x)as.numeric(t(beta[outcome,]) %*% x))+rnorm(n=num_subjects,sd=0.5) }
        else { y_control_outcome <- apply(covariates,1,function(x)as.numeric(t(beta[outcome,]) %*% x))+U+rnorm(n=num_subjects,sd=0.5) }
        y_control_big[,outcome] <- y_control_outcome
      }
      # # adjust null outcomes
      # for (i in 1:(2*r)){
      #   y_control_big[i,which(theta_vec==0)] <- qnorm(1-Uvec[i])+i/200+y_control_big[i,which(theta_vec==0)]
      # }
      matchres <- run_PSM(data_matching=covariates,data_outcomes=y_control_big,treated=trt_assignment_vec)
      standardized_difs=matchres$standardized_difs
      num_matched_sets=matchres$num_matched_sets
      data_outcomes=matchres$data_outcomes
      treated_idx=matchres$ix_subjects_treated
      control_idx=matchres$ix_subjects_control
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
    }
    
    
    ## ONE-SHOT METHOD ##
    our_tested_ix_old=logical(length = length(outcomes))
    for (outcome in outcomes){
      y_oneshot <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n1){
        y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]]+theta_vec[outcome])
      }
      our_tested_ix_old[outcome] <- senWilcox(d=diff(y_oneshot)[seq(1,length(y_oneshot),2)],gamma=gamma,alpha=alpha,alternative='greater')$pval < alpha
    }
    num_outcomes_tested_vec_old[sim] <- sum(our_tested_ix_old)
    # create counters, vector for indices to reject
    true_positive_count <- false_positive_count <- true_negative_count <- 
      false_negative_count <- 0 
    our_rejected_ix <- rep(FALSE,length(outcomes))
    for (outcome in outcomes){
      true_effect <- theta_vec[outcome] != 0
      if (our_tested_ix_old[outcome]){ # if we choose to test for outcome
        y <- as.matrix(data_outcomes)[,outcome]
        y_oneshot <- numeric()
        for (set in 1:n2){
          y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]]+theta_vec[outcome])
        }
        reject <- senWilcox(d=diff(y_oneshot)[seq(1,length(y_oneshot),2)],gamma=gamma,alpha=alpha,alternative='greater')$pval < alpha/sum(our_tested_ix_old)
        # update metrics
        if(reject) our_rejected_ix[outcome] <- TRUE
        if(reject && true_effect) true_positive_count <- true_positive_count+1
        else if(reject && !true_effect) false_positive_count <- false_positive_count+1
        else if(!reject && true_effect) false_negative_count <- false_negative_count+1
        else true_negative_count <- true_negative_count+1
      }
      else { # if we choose to not test for outcome
        if(true_effect) false_negative_count <- false_negative_count+1
        else true_negative_count <- true_negative_count+1
      }
    }
    # compute metrics
    true_positive_rate_old <- true_positive_count/num_effects # sensitivity
    sensitivity_vec_old[sim] <- true_positive_rate_old
    
    
    ## Bootstrap Method ##
    our_tested_ix=logical(length = length(outcomes))
    ixsetlist=list()
    for (boot in 1:nboot){
      ixset <- sample(x=1:planning_sample_size_sets, size=planning_sample_size_sets, replace=T)
      ixsetlist[[boot]] <- c(as.matrix(ix_subjects_control_planning)[ixset,],ix_subjects_treated_planning[ixset])
    }
    
    sigma_hat_vec <- kappahatvec_whole <- numeric(length = length(outcomes))
    for (outcome in outcomes){
      ## Compute sigmahat
      y <- as.matrix(data_outcomes)[,outcome]
      kappahatvec <- numeric()
      for (boot in 1:nboot){
        y_boot <- y[ixsetlist[[boot]]]
        y_boot[seq(2,length(y_boot),2)] <- y_boot[seq(2,length(y_boot),2)]+theta_vec[outcome]
        mpdifs <- diff(y_boot)[seq(1,length(y_boot),2)]
        kappahatvec <- c(kappahatvec, get_sensitivity_val(mpdifs,a=alpha))
      }
      sigma_hat_vec[outcome] <- sqrt(var(kappahatvec)*n1)
      ## Compute kappahat on the whole planning sample
      y_boot <- numeric()
      for (set in 1:n1){
        y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]]+theta_vec[outcome])
      }
      kappahatvec_whole[outcome] <- get_sensitivity_val(diff(y_boot)[seq(1,length(y_boot),2)], a=alpha)
    }
    
    aprimegrid <- alpha/(1:num_outcomes)
    tested_list <- vector("list", length(outcomes))
    difvec <- numeric(length = length(outcomes))
    for (aprime_idx in 1:length(aprimegrid)){
      aprime <- aprimegrid[aprime_idx]
      our_tested_ix = logical(length = length(outcomes))
      for (outcome in outcomes){
        sigma_hat <- sigma_hat_vec[outcome]
        kappahat <- kappahatvec_whole[outcome]
        sigma_r <- sigma_hat / sqrt(planning_sample_prop * (1 - planning_sample_prop))
        c1_hat <- -sqrt(4/3) * qnorm(1 - alpha) * sqrt(kappahat * (1 - kappahat))
        c2_hat <- -sqrt(4/3) * qnorm(1 - aprime) * sqrt(kappahat * (1 - kappahat))
        mu_r <- c1_hat / sqrt(planning_sample_prop) - c2_hat / sqrt(1 - planning_sample_prop)
        our_tested_ix[outcome] <- kappahat > gamma / (1 + gamma) + (mu_r - qnorm(1 - alpha) * sigma_r) / (sqrt(n1 + n2))
      }
      tested_list[[aprime_idx]] <- our_tested_ix
      difvec[aprime_idx] <- abs(alpha-aprime*sum(our_tested_ix))
    }
    aprime <- aprimegrid[which.min(difvec)]
    our_tested_ix <- tested_list[[which.min(difvec)]]
    num_outcomes_tested_vec[sim] <- sum(our_tested_ix)
    cat('Num Outcomes Tested:', sum(our_tested_ix),'\n\n')
    # create counters, vector for indices to reject
    true_positive_count <- false_positive_count <- true_negative_count <- false_negative_count <- 0 
    our_rejected_ix <- rep(FALSE,length(outcomes))
    for (outcome in outcomes){
      true_effect <- theta_vec[outcome] != 0
      if (our_tested_ix[outcome]){ # if we choose to test for outcome
        y_boot <- numeric()
        y <- as.matrix(data_outcomes)[,outcome]
        for (set in 1:n2){
          y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]]+theta_vec[outcome])
        }
        reject <- senWilcox(d=diff(y_boot)[seq(1,length(y_boot),2)],gamma=gamma,alpha=alpha,alternative='greater')$pval < alpha/sum(our_tested_ix)
        # update metrics
        if(reject) our_rejected_ix[outcome] <- TRUE
        if(reject && true_effect) true_positive_count <- true_positive_count+1
        else if(reject && !true_effect) false_positive_count <- false_positive_count+1
        else if(!reject && true_effect) false_negative_count <- false_negative_count+1
        else true_negative_count <- true_negative_count+1
      }
      else { # if we choose to not test for outcome
        if(true_effect) false_negative_count <- false_negative_count+1
        else true_negative_count <- true_negative_count+1
      }
    }
    # compute metrics
    true_positive_rate <- true_positive_count/num_effects # sensitivity
    false_positive_rate <- false_positive_count/(length(theta_vec)-num_effects) # type i error
    true_negative_rate <- true_negative_count/(length(theta_vec)-num_effects) # specificity
    false_negative_rate <- false_negative_count/num_effects # type ii error
    # append metrics to lists, vectors
    outcomes_tested_list[[sim]] <- our_tested_ix
    outcomes_rejected_list[[sim]] <- our_rejected_ix
    sensitivity_vec[sim] <- true_positive_rate
    type1error_vec[sim] <- false_positive_rate
    specificity_vec[sim] <- true_negative_rate
    type2error_vec[sim] <- false_negative_rate
    
    
    ## BONFERRONI BASELINE METHOD ##
    bonf_optimal_ix <- rep(FALSE,length(outcomes))
    bonf_sensitivity_count <- 0 
    for (outcome in outcomes){
      true_effect <- theta_vec[outcome] != 0
      y_bonf <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:num_matched_sets){
        y_bonf <- c(y_bonf,y[as.matrix(control_idx)[set,]],y[treated_idx[set]]+theta_vec[outcome])
      }
      #browser()
      reject <- senWilcox(d=diff(y_bonf)[seq(1,length(y_bonf),2)],gamma=gamma,alpha=alpha,alternative='greater')$pval < alpha/num_outcomes
      if (reject) bonf_optimal_ix[outcome] <- TRUE
      if (reject && true_effect) bonf_sensitivity_count <- bonf_sensitivity_count+1
    }
    bonf_sensitivity_vec[sim] <- bonf_sensitivity_count/num_effects
    theta_vec=thetavecold
  }
  
  return(list(outcomes_tested=outcomes_tested_list, outcomes_rejected=outcomes_rejected_list,
              avg_sensitivity=mean(sensitivity_vec), avg_sensitivity_old=mean(sensitivity_vec_old),
              avg_type1error=mean(type1error_vec), avg_specificity=mean(specificity_vec),
              avg_type2error=mean(type2error_vec), bonf_sensitivity=mean(bonf_sensitivity_vec),
              avg_num_outcomes_tested_old=mean(num_outcomes_tested_vec_old),
              avg_num_outcomes_tested=mean(num_outcomes_tested_vec)))
}


############################################################################

# Simulation setting
gammas <- c(1,1.25,1.5,2,2.5,3)

true_theta_list=vector("list",length=length(gammas))
true_theta_list[[1]] <- seq(from=0.0,to=1.4,by=0.14)
true_theta_list[[2]] <- seq(from=0.2,to=1.6,by=0.14)
true_theta_list[[3]] <- seq(from=0.4,to=1.8,by=0.14)
true_theta_list[[4]] <- seq(from=0.8,to=2.2,by=0.14)
true_theta_list[[5]] <- seq(from=1.2,to=2.6,by=0.14)
true_theta_list[[6]] <- seq(from=1.6,to=3.0,by=0.14)

set.seed(0)
for (gamma_idx in 1:length(gammas)){
  
  gamma <- gammas[gamma_idx]
  true_thetas <- true_theta_list[[gamma_idx]]
  sensvec<-oneshotsensvec<-bonfsensvec<-numeric()
  
  for (true_theta in true_thetas){
    cat('\n')
    cat(paste(100*(which(true_thetas == true_theta)-1)/length(true_thetas), 'percent done'), '\n')
    sim_res=run_sim(true_theta=true_theta, gamma=gamma, uc_dgp = F)
    sensvec<-c(sensvec,sim_res$avg_sensitivity)
    oneshotsensvec<-c(oneshotsensvec,sim_res$avg_sensitivity_old)
    bonfsensvec<-c(bonfsensvec,sim_res$bonf_sensitivity)
  }
  
  data <- data.frame(
    true_thetas = true_thetas,
    Sensitivity = c(bonfsensvec, oneshotsensvec, sensvec),
    Method = factor(rep(c("Bonferroni", "Naive", "Sens-Val"), each = length(true_thetas)))
  )
  
  saveRDS(data, file = paste0("dat-NUC-df-","gamma",gamma,".rds"))
  
}


true_theta_list=vector("list",length=length(gammas))
true_theta_list[[1]] <- seq(from=0.0,to=1.4,by=0.14)+0.1
true_theta_list[[2]] <- seq(from=0.2,to=1.6,by=0.14)+0.1
true_theta_list[[3]] <- seq(from=0.4,to=1.8,by=0.14)+0.1
true_theta_list[[4]] <- seq(from=0.8,to=2.2,by=0.14)+0.1
true_theta_list[[5]] <- seq(from=1.2,to=2.6,by=0.14)+0.1
true_theta_list[[6]] <- seq(from=1.6,to=3.0,by=0.14)+0.1

set.seed(0)
for (gamma_idx in 1:length(gammas)){
  
  gamma <- gammas[gamma_idx]
  true_thetas <- true_theta_list[[gamma_idx]]
  sensvec<-oneshotsensvec<-bonfsensvec<-numeric()
  
  for (true_theta in true_thetas){
    cat('\n')
    cat(paste(100*(which(true_thetas == true_theta)-1)/length(true_thetas), 'percent done'), '\n')
    sim_res=run_sim(true_theta=true_theta, gamma=gamma, uc_dgp = T)
    sensvec<-c(sensvec,sim_res$avg_sensitivity)
    oneshotsensvec<-c(oneshotsensvec,sim_res$avg_sensitivity_old)
    bonfsensvec<-c(bonfsensvec,sim_res$bonf_sensitivity)
  }
  
  data <- data.frame(
    true_thetas = true_thetas,
    Sensitivity = c(bonfsensvec, oneshotsensvec, sensvec),
    Method = factor(rep(c("Bonferroni", "Naive", "Sens-Val"), each = length(true_thetas)))
  )
  
  saveRDS(data, file = paste0("dat-UC-df-","gamma",gamma,".rds"))
  
}

