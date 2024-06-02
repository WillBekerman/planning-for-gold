
## Load required libraries
library(DOS2)

## Computes sensitivity value for pair-matched data using Wilcoxon SRT.
## Inputs: mpdifs, a, alt.
## Outputs: kappa
get_sensitivity_val <- function(mpdifs,a=a,alt=alt) {
  senfx <- function(mpdifs,g,a,alt){
    senWilcox(mpdifs, gamma=g, alpha=a, alternative=alt)$pval-a
  }
  g<-1
  try(g<-uniroot(senfx,c(1,1e10),mpdifs=mpdifs,a=a,alt=alt)$root,silent=TRUE)
  kappa<-g/(1+g)
  return(kappa)
}

## Runs Sens-Val procedure on pair-matched data using Wilcoxon SRT.
## Inputs: data_outcomes, control_idx, treated_idx, alt_vec, 
## planning_sample_prop, alpha, gamma, nboot, seed.
## Outputs: outcomes_tested, outcomes_tested_oneshot, outcomes_rejected,
## outcomes_rejected_oneshot, outcomes_rejected_bonferroni,
## planning_sample_proportion, gamma.
sensval_pairs <- function(data_outcomes,control_idx,treated_idx,
                          alt_vec=NULL,planning_sample_prop=0.20,
                          alpha=0.05,gamma=1,nboot=250,seed=0){
  
  ## Set seed for reproducibility
  set.seed(seed)
  num_outcomes <- ncol(data_outcomes)
  num_matched_sets <- length(treated_idx)
  
  ## Create objects to return
  outcomes_tested <- outcomes_rejected <- logical(length=num_outcomes)
  
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
  
  ## Planning stage
  if (is.null(alt_vec)){
    ## Get vector for alternatives to test using planning data
    alt_vec <- character(num_outcomes)
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
  }
  ixsetlist <- list()
  for (boot in 1:nboot){
    ixset <- sample(x=1:planning_sample_size_sets, size=planning_sample_size_sets, replace=T)
    ixsetlist[[boot]] <- c(as.matrix(ix_subjects_control_planning)[ixset,],ix_subjects_treated_planning[ixset])
  }
  sigma_hat_vec <- kappahatvec_whole <- numeric(length = num_outcomes)
  for (outcome in 1:num_outcomes){
    ## Compute sigmahat via bootstrapping
    y <- as.matrix(data_outcomes)[,outcome]
    kappahatvec <- numeric()
    for (boot in 1:nboot){
      y_boot <- y[ixsetlist[[boot]]]
      mpdifs <- diff(y_boot)[seq(1,length(y_boot),2)]
      if (!(all(is.na(mpdifs)))){ kappahatvec <- c(kappahatvec, get_sensitivity_val(mpdifs[!is.na(mpdifs)],a=alpha,alt=alt_vec[outcome])) }
    }
    sigma_hat_vec[outcome] <- sqrt(var(kappahatvec)*n1)
    ## Compute kappahat on the whole planning sample
    y_boot <- numeric()
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])))) { next }
      y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])
    }
    kappahatvec_whole[outcome] <- get_sensitivity_val(diff(y_boot)[seq(1,length(y_boot),2)], a=alpha, alt=alt_vec[outcome])
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
  
  ## Analysis stage
  for (outcome in 1:num_outcomes){
    if (outcomes_tested[outcome]){
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n2){
        if (any(is.na(c(y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])))) { next }
        y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])
      }
      reject <- senWilcox(d=diff(y_boot)[seq(1,length(y_boot),2)],gamma=gamma,alpha=alpha,alternative=alt_vec[outcome])$pval < alpha/sum(outcomes_tested)
      if (!is.na(reject)&reject) { outcomes_rejected[outcome] <- T }
    }
  }
  
  return(list(outcomes_tested=names(data_outcomes)[which(outcomes_tested)],
              outcomes_rejected=names(data_outcomes)[which(outcomes_rejected)]))
}
