
## Run after application_matching.R

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

is_binary_outcome = function(x) { return(all(is.na(x) | x==0 | abs(x)==1 )) }


get_statistics <- function(data_outcomes,control_idx,treated_idx,
                           planning_sample_prop=0.10,alpha=0.05,
                           gamma=1,nboot=100,seed=0){
  
  ## Set seed for reproducibility
  set.seed(seed)
  
  num_outcomes <- ncol(data_outcomes)
  num_matched_sets <- length(treated_idx)
  
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
  binary_ix <- which( apply(data_outcomes,2,is_binary_outcome) )
  
  planning_list <- analysis_list <- whole_list <- vector("list",num_outcomes)
  names(planning_list) <- names(analysis_list) <- names(whole_list) <- names(data_outcomes)
  
  # Planning sample
  for (outcome in 1:num_outcomes){
    y_oneshot <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])))) { next }
      y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_planning)[set,]],y[ix_subjects_treated_planning[set]])
    }
    cix <- seq(1,length(y_oneshot),2)
    tix <- cix+1
    if (!is_binary_outcome(y_oneshot)){ planning_list[[outcome]]$control=summary(y_oneshot[cix]); planning_list[[outcome]]$treated=summary(y_oneshot[tix]) }
    if (is_binary_outcome(y_oneshot)){ planning_list[[outcome]]$control=proportions(table(y_oneshot[cix])); planning_list[[outcome]]$treated=proportions(table(y_oneshot[tix])) }
  }
  
  # Analysis sample
  for (outcome in 1:num_outcomes){
    y_oneshot <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:n2){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])))) { next }
      y_oneshot <- c(y_oneshot,y[as.matrix(ix_subjects_control_analysis)[set,]],y[ix_subjects_treated_analysis[set]])
    }
    cix <- seq(1,length(y_oneshot),2)
    tix <- cix+1
    if (!is_binary_outcome(y_oneshot)){ analysis_list[[outcome]]$control=summary(y_oneshot[cix]); analysis_list[[outcome]]$treated=summary(y_oneshot[tix]) }
    if (is_binary_outcome(y_oneshot)){ analysis_list[[outcome]]$control=proportions(table(y_oneshot[cix])); analysis_list[[outcome]]$treated=proportions(table(y_oneshot[tix])) }
  }
  
  # Whole sample
  for (outcome in 1:num_outcomes){
    y_oneshot <- numeric()
    y <- as.matrix(data_outcomes)[,outcome]
    for (set in 1:num_matched_sets){
      if (any(is.na(c(y[as.matrix(control_idx)[set,]],y[treated_idx[set]])))) { next }
      y_oneshot <- c(y_oneshot,y[as.matrix(control_idx)[set,]],y[treated_idx[set]])
    }
    if (!is_binary_outcome(y_oneshot)){ whole_list[[outcome]]=summary(y_oneshot) }
    if (is_binary_outcome(y_oneshot)){ whole_list[[outcome]]=proportions(table(y_oneshot)) }
  }
  
  return(list(planning_list=planning_list, analysis_list=analysis_list, whole_list=whole_list))
  
}

res<-get_statistics(data_outcomes,control_idx=ix_subjects_control,
                    treated_idx=ix_subjects_treated,planning_sample_prop=0.2)
num_outcomes_bin <- sum( apply(data_outcomes,2,is_binary_outcome) ); ix_outcomes_bin <- which( apply(data_outcomes,2,is_binary_outcome) )
num_outcomes_num <- ncol(data_outcomes) - num_outcomes_bin; ix_outcomes_num <- which( apply(data_outcomes,2,function(x)!is_binary_outcome(x)) )







### Planning Sample Summary Statistics

# Function to generate a data frame for each numeric outcome variable
generate_data_num <- function(varnum) {
  ix <- ix_outcomes_num[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 6),
    Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
    Control = as.numeric(res$planning_list[[ix]]$control),
    Treated = as.numeric(res$planning_list[[ix]]$treated)
  )
  return(df)
}

# Create a list of data frames for all numeric outcomes
planning_num_list <- setNames(lapply(1:length(ix_outcomes_num), generate_data_num), names(data_outcomes)[ix_outcomes_num])
planning_num <- do.call(rbind, planning_num_list)

# Set the number of digits for numerical columns and correct alignment
digits <- c(0, 2, 2)  # Assuming the 'Statistic' column does not need numerical precision
align <- "llrr"  # Column alignment for 'Variable', 'Statistic', 'Control', and 'Treated'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(planning_num)), planning_num$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & & \\\\\\hline", planning_num$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f & %.2f \\\\", planning_num$Statistic[rows], planning_num$Control[rows], planning_num$Treated[rows]))
}, headers, split(seq_len(nrow(planning_num)), planning_num$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Control", "Treated", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Key Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "planning_num_summarystats.tex"
writeLines(latex_table, con = output_file)



# Function to generate a data frame for each binary outcome variable
generate_data_bin <- function(varnum) {
  ix <- ix_outcomes_bin[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 2),
    Statistic = c("0", "1"),
    Control = as.numeric(res$planning_list[[ix]]$control),
    Treated = as.numeric(res$planning_list[[ix]]$treated)
  )
  return(df)
}

# Create a list of data frames for all binary outcomes
planning_bin_list <- setNames(lapply(1:length(ix_outcomes_bin), generate_data_bin), names(data_outcomes)[ix_outcomes_bin])
planning_bin <- do.call(rbind, planning_bin_list)

# Set the number of digits for numerical columns
# Note: No precision is needed for the 'Statistic' column as it's categorical, but place the '0' to keep the structure
digits <- c(0, 0, 2)  # For 'Statistic', 'Control', and 'Treated'
align <- "llrr"  # Column alignment for 'Variable', 'Statistic', 'Control', and 'Treated'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(planning_bin)), planning_bin$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & & \\\\\\hline", planning_bin$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f & %.2f \\\\", planning_bin$Statistic[rows], planning_bin$Control[rows], planning_bin$Treated[rows]))
}, headers, split(seq_len(nrow(planning_bin)), planning_bin$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Control", "Treated", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Binary Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "planning_bin_summarystats.tex"
writeLines(latex_table, con = output_file)




































### Analysis Sample Summary Statistics

# Function to generate a data frame for each numeric outcome variable
generate_data_num <- function(varnum) {
  ix <- ix_outcomes_num[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 6),
    Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
    Control = as.numeric(res$analysis_list[[ix]]$control),
    Treated = as.numeric(res$analysis_list[[ix]]$treated)
  )
  return(df)
}

# Create a list of data frames for all numeric outcomes
analysis_num_list <- setNames(lapply(1:length(ix_outcomes_num), generate_data_num), names(data_outcomes)[ix_outcomes_num])
analysis_num <- do.call(rbind, analysis_num_list)

# Set the number of digits for numerical columns and correct alignment
digits <- c(0, 2, 2)  # Assuming the 'Statistic' column does not need numerical precision
align <- "llrr"  # Column alignment for 'Variable', 'Statistic', 'Control', and 'Treated'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(analysis_num)), analysis_num$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & & \\\\\\hline", analysis_num$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f & %.2f \\\\", analysis_num$Statistic[rows], analysis_num$Control[rows], analysis_num$Treated[rows]))
}, headers, split(seq_len(nrow(analysis_num)), analysis_num$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Control", "Treated", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Key Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "analysis_num_summarystats.tex"
writeLines(latex_table, con = output_file)



# Function to generate a data frame for each binary outcome variable
generate_data_bin <- function(varnum) {
  ix <- ix_outcomes_bin[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 2),
    Statistic = c("0", "1"),
    Control = as.numeric(res$analysis_list[[ix]]$control),
    Treated = as.numeric(res$analysis_list[[ix]]$treated)
  )
  return(df)
}

# Create a list of data frames for all binary outcomes
analysis_bin_list <- setNames(lapply(1:length(ix_outcomes_bin), generate_data_bin), names(data_outcomes)[ix_outcomes_bin])
analysis_bin <- do.call(rbind, analysis_bin_list)

# Set the number of digits for numerical columns
# Note: No precision is needed for the 'Statistic' column as it's categorical, but place the '0' to keep the structure
digits <- c(0, 0, 2)  # For 'Statistic', 'Control', and 'Treated'
align <- "llrr"  # Column alignment for 'Variable', 'Statistic', 'Control', and 'Treated'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(analysis_bin)), analysis_bin$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & & \\\\\\hline", analysis_bin$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f & %.2f \\\\", analysis_bin$Statistic[rows], analysis_bin$Control[rows], analysis_bin$Treated[rows]))
}, headers, split(seq_len(nrow(analysis_bin)), analysis_bin$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Control", "Treated", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Binary Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "analysis_bin_summarystats.tex"
writeLines(latex_table, con = output_file)
































### Whole Sample Summary Statistics

# Function to generate a data frame for each numeric outcome variable
generate_data_num <- function(varnum) {
  ix <- ix_outcomes_num[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 6),
    Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
    Sample = as.numeric(res$whole_list[[ix]])
  )
  return(df)
}

# Create a list of data frames for all numeric outcomes
whole_num_list <- setNames(lapply(1:length(ix_outcomes_num), generate_data_num), names(data_outcomes)[ix_outcomes_num])
whole_num <- do.call(rbind, whole_num_list)

# Set the number of digits for numerical columns and correct alignment
digits <- c(0, 2)  # Assuming the 'Statistic' column does not need numerical precision
align <- "llr"  # Column alignment for 'Variable', 'Statistic', 'Sample'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(whole_num)), whole_num$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & \\\\\\hline", whole_num$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f \\\\", whole_num$Statistic[rows], whole_num$Sample[rows]))
}, headers, split(seq_len(nrow(whole_num)), whole_num$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Sample", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Key Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "whole_num_summarystats.tex"
writeLines(latex_table, con = output_file)



# Function to generate a data frame for each binary outcome variable
generate_data_bin <- function(varnum) {
  ix <- ix_outcomes_bin[varnum]
  df <- data.frame(
    Variable = rep(names(data_outcomes)[ix], 2),
    Statistic = c("0", "1"),
    Sample = as.numeric(res$whole_list[[ix]])
  )
  return(df)
}

# Create a list of data frames for all binary outcomes
whole_bin_list <- setNames(lapply(1:length(ix_outcomes_bin), generate_data_bin), names(data_outcomes)[ix_outcomes_bin])
whole_bin <- do.call(rbind, whole_bin_list)

# Set the number of digits for numerical columns
# Note: No precision is needed for the 'Statistic' column as it's categorical, but place the '0' to keep the structure
digits <- c(0, 2)  # For 'Statistic', 'Sample'
align <- "llr"  # Column alignment for 'Variable', 'Statistic', 'Sample'

# Group by 'Variable' and create headers in LaTeX format
headers <- tapply(seq_len(nrow(whole_bin)), whole_bin$Variable, function(rows) {
  c(sprintf("\\textbf{%s} & & \\\\\\hline", whole_bin$Variable[rows[1]]), rep("", length(rows) - 1))
})

# Combine headers with their respective rows
latex_content <- unlist(mapply(function(header, rows) {
  c(header, sprintf(" & %s & %.2f \\\\", whole_bin$Statistic[rows], whole_bin$Sample[rows]))
}, headers, split(seq_len(nrow(whole_bin)), whole_bin$Variable), SIMPLIFY = FALSE), use.names = FALSE)

# Construct the header row for the table
column_headers <- c("&", "Statistic", "Sample", "\\\\\\hline")

# Combine the column headers with the content
latex_table_content <- c(column_headers, latex_content)

# Print the table to get the LaTeX code
latex_table <- paste(
  "\\begin{longtable}{", align, "}\n",
  "\\caption{Summary Statistics for Binary Variables} \\\\\n",
  "\\hline\n",
  paste(latex_table_content, collapse = "\n"),
  "\\hline\n",
  "\\end{longtable}",
  sep = ""
)

# Save the LaTeX table to a file
output_file <- "whole_bin_summarystats.tex"
writeLines(latex_table, con = output_file)


