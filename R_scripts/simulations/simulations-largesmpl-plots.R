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

# Simulation setting
gammas <- c(1,1.25,1.5,2,2.5,3)

true_theta_list=vector("list",length=length(gammas))
true_theta_list[[1]] <- seq(from=0.0,to=1.4,by=0.14)
true_theta_list[[2]] <- seq(from=0.2,to=1.6,by=0.14)
true_theta_list[[3]] <- seq(from=0.2,to=1.6,by=0.14)
true_theta_list[[4]] <- seq(from=0.6,to=2.0,by=0.14)
true_theta_list[[5]] <- seq(from=0.8,to=2.2,by=0.14)
true_theta_list[[6]] <- seq(from=1.0,to=2.4,by=0.14)

plots_nuc <- list()  # Initialize an empty list to store plots
ncol <- 3  # Define the number of cols for your grid


files <- list.files(path=".")

set.seed(0)
for (gamma_idx in 1:length(gammas)){
  
  gamma <- gammas[gamma_idx]
  true_thetas <- true_theta_list[[gamma_idx]]
  show_ylabel <- (gamma_idx - 1) %% ncol == 0
  
  fileixs <- c( ((gamma_idx-1)*11+1):((gamma_idx-1)*11+11))
  files_gamma <- files[fileixs]
  thetavals <- numeric(length=11)
  for (fileix in 1:length(fileixs)){
    filename <- files_gamma[fileix]
    thetavals[fileix] <- as.numeric( strsplit(strsplit(filename, 'thetaxxx')[[1]][2], 'xxx.rds')[[1]] )
  }
  fileorder <- sort(thetavals, index.return=T)$ix
  data <- setNames(data.frame(matrix(ncol=3,nrow =0)), c("true_thetas", "Sensitivity", "Method"))
  for (fileix in 1:length(fileixs)){
    val <- thetavals[fileorder[fileix]]
    df <- readRDS(files_gamma[fileorder[fileix]])
    dfnew <- data.frame(true_thetas=rep(val,3),
                        Sensitivity=df$Sensitivity[1:3],
                        Method=c('Bonferroni','Naive','Sens-Val'))
    data <- rbind(data,dfnew)
  }
  
  p <- ggplot(data, aes(x = true_thetas, y = Sensitivity)) +
    geom_line(aes(linetype=Method), size=1) + scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
    labs(
      x = expression(tau),
      y = if(show_ylabel) "Average Power" else "",
      title = bquote(Gamma[con] ~ "=" ~ .(format(gamma, nsmall = 2)))
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 7)) +
    coord_cartesian(ylim = c(0, 1)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plots_nuc[[length(plots_nuc) + 1]] <- p  # Store the plot in the list
  
}

true_theta_list=vector("list",length=length(gammas))
true_theta_list[[1]] <- seq(from=0.0,to=1.4,by=0.14)+0.1
true_theta_list[[2]] <- seq(from=0.2,to=1.6,by=0.14)+0.1
true_theta_list[[3]] <- seq(from=0.2,to=1.6,by=0.14)+0.1
true_theta_list[[4]] <- seq(from=0.6,to=2.0,by=0.14)+0.1
true_theta_list[[5]] <- seq(from=0.8,to=2.2,by=0.14)+0.1
true_theta_list[[6]] <- seq(from=1.0,to=2.4,by=0.14)+0.1

plots_uc <- list()  # Initialize an empty list to store plots
ncol <- 3  # Define the number of cols for your grid

set.seed(0)
for (gamma_idx in 1:length(gammas)){
  
  gamma <- gammas[gamma_idx]
  true_thetas <- true_theta_list[[gamma_idx]]
  show_ylabel <- (gamma_idx - 1) %% ncol == 0
  
  fileixs <- c( ((gamma_idx-1)*11+1+66):((gamma_idx-1)*11+11+66) )
  files_gamma <- files[fileixs]
  thetavals <- numeric(length=11)
  for (fileix in 1:length(fileixs)){
    filename <- files_gamma[fileix]
    thetavals[fileix] <- as.numeric( strsplit(strsplit(filename, 'thetaxxx')[[1]][2], 'xxx.rds')[[1]] )
  }
  fileorder <- sort(thetavals, index.return=T)$ix
  data <- setNames(data.frame(matrix(ncol=3,nrow =0)), c("true_thetas", "Sensitivity", "Method"))
  for (fileix in 1:length(fileixs)){
    val <- thetavals[fileorder[fileix]]
    df <- readRDS(files_gamma[fileorder[fileix]])
    dfnew <- data.frame(true_thetas=rep(val,3),
                        Sensitivity=df$Sensitivity[1:3],
                        Method=c('Bonferroni','Naive','Sens-Val'))
    data <- rbind(data,dfnew)
  }
  
  p <- ggplot(data, aes(x = true_thetas, y = Sensitivity)) +
    geom_line(aes(linetype=Method), size=1) + scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
    labs(
      x = expression(tau),
      y = if(show_ylabel) "Average Power" else "",
      title = bquote(Gamma[con] ~ "=" ~ .(format(gamma, nsmall = 2)))
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 7)) +
    coord_cartesian(ylim = c(0, 1)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plots_uc[[length(plots_uc) + 1]] <- p  # Store the plot in the list
  
}

get_legend_plot <- ggplot(data, aes(x = true_thetas, y = Sensitivity)) +
  geom_line(aes(linetype=Method), size=1) + scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
  labs(
    x = expression(tau),
    y = "Average Power",
    title = bquote(Gamma ~ "=" ~ .(format(gamma, nsmall = 2)))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, 1))
all_plots_nuc <- plot_grid(plotlist = plots_nuc, nrow=2, ncol=3)
all_plots_uc <- plot_grid(plotlist = plots_uc, nrow=2, ncol=3)
legend <- get_legend(get_legend_plot +
                       guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "bottom"))
combined_plot <- plot_grid(NULL, all_plots_nuc, NULL, all_plots_uc, legend, ncol=1, rel_heights = c(.1,1,.1,1,.1), labels = c('', 'A', '', 'B', ''), vjust=-0.5)

combined_plot
# Save the combined plot
ggsave("combined_plots.pdf", combined_plot, device = "pdf", width = 12, height = 10)

