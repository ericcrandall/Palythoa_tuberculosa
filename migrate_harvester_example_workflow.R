# This script models a workflow for analyzing migrate-n output using the migrate_harvester_functions.R

## Dependencies
install.packages("tidyverse")
install.packages("Rmpfr")
install.packages("spatstat.explore")
install.packages("ggridges")

source("./migrate_harvester_functions.R")
library(tidyverse)
## starting with a set of 5 replicated runs, each containing the same set of model outputs, 
## use harvest_model_likelihoods to create a list of 5 data frames
repsDir <-"~/github/Palythoa_tuberculosa/migrate/run6/final/"

marg_like_list <- list()

for(r in 1:5){
  rep = paste0("rep",r)
  print(rep)
  marg_like_list[[rep]] <- harvest_model_likelihoods(modelsDir=file.path(repDir,rep))
}

## Use bf_calcs to perform model selection for each replicate...
marg_like_list %>% map(bf_calcs)
## OR better yet, collapse the likelihoods into a dataframe, 
## take the mean marginal likelihood across replicates 
## and perform model selection with bf_calcs on the summarized data frame.
like.df <-  marg_like_list %>% bind_rows() %>% group_by(model)
means <- like.df %>% summarize(bezier.corrected = mean(bezier.corrected))
final_model_selection <- bf_calcs(means)

# read in the posterior for the winning model
modelsDir <- "/Users/eric/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Documents/Research/Ptuberculosa/migrate/run6/rep1"
winningModel <- final_model_selection$model[which(final_model_selection$choice==1)]
posterior <- read_bayesallfile(Dir = file.path(modelsDir, winningModel), Nm_calc = T)

# read in a population key to translate numbers into names in figures
parameter_key <- make_parameter_key(Dir = "./", popkeyname = "popkey.csv", posterior = posterior)


##summarize the posterior for theta. 
##This will take up to 15-20 minutes depending on number of parameters, length of MCMC chains,
## and how many density points you specified with n
## the mean needs to be the same mean that you gave to migrate-n, 
## but upper and lower bounds can bracket whatever region of parameter space you'd like to focus on.
thetas <- summarize_posterior(posterior, parameter_type = "Theta", 
                              exponential_mean=0.001,lower.bound=0, upper.bound=0.1, n=2^14)

# this will calculate a table of statistics for all theta parameters
# (stats for each Theta parameter calculated by parameter_stats())
theta_stats <- posterior_stats(summarized_posterior = thetas, 
                               parameter_key = parameter_key,
                               parameter_type = "Theta")

# this will plot the theta posteriors as a ridge plot. 
# because colors aren't specified, it will pick them from the "turbo" palette for the viridis package
theta_plot <- plot_parameter_posteriors(summarized_posterior = thetas,
                                        parameter_key = parameter_key, 
                                        parameter_type = "Theta", x_limits = c(0,0.004), 
                                        param_colors = NULL)
