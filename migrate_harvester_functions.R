# Migrate Harvester Functions
# R Functions for reading and making sense of migrate-n output
# By Eric Crandall (ecrandall@psu.edu)
# Version 0.1 March 25, 2024

# Functions for reading in the bayesallfile and a user-created popkey.csv
#
#
read_bayesallfile <- function(Dir, Nm_calc = T){
  #this function reads in a Migrate-n bayesallfile from the specificed directory
  #and optionally calculates Nm from Theta and M parameters (needs migrants_per_gen() to do so)
 require(readr)
 require(dplyr)
 require(magrittr)
  
  if(Nm_calc){
    posterior <- read_table(file.path(Dir,"bayesallfile"), comment = "#") %>% 
      migrants.per.gen() %>% dplyr::select(starts_with(c("Locus","Steps","Replicate","Theta_", "M_", 
                                                         "Nm_","D_")))
  } else{
    posterior <- read_table(file.path(Dir,"bayesallfile"), comment = "#") %>% 
       dplyr::select(starts_with(c("Locus","Steps","Replicate","Theta_", "M_", 
                                                        "Nm_", "D_")))
  }
  return(posterior)
}


make_parameter_key <- function(Dir,popkeyname = "popkey.csv", posterior){
# this function reads in a comma-delimited population key (must be named "popkey.csv") with two columns
# column 1: Index, which gives the number of the population in the order specified in the data infile
# column 2: Pop, which is the name of the population you'd like to use in your figures
# It also takes a bayesallfile posterior data frame that has been loaded with read_bayesallfile()
# It then translates the popkey into a parameter key that can be used to label your parameter plots
  require(readr)
  require(dplyr)
  require(magrittr)
  require(stringr)
  
  popkey <- read_csv(file.path(Dir,popkeyname))
  names(popkey$Pop) <- popkey$Index
  #get all the parameter names from the posterior
  parameter_key <- names(posterior) %>% str_subset("Theta|M_|Nm_|D_")
  
  names(parameter_key) <- parameter_key %>% 
    str_replace("(\\d)_(\\d)", "\\1 to \\2") %>%
    str_replace("M_", "m/mu " ) %>% 
    str_replace("Theta_", "theta ") %>% 
    str_replace_all(popkey$Pop) %>% 
    str_replace("Nm_", "4Nem ") 
  
  return(parameter_key)
}


# Functions for harvesting model marginal likelihoods from multiple models and doing model selection
#
#
harvest_model_likelihoods <- function(modelsDir,
                                      outfileName = "outfile.txt",
                                      multilocus = T){
  # this function harvests model marginal likelihoods for models calculated by
  # the program migrate-n (Beerli & Felsenstein 2001).
  # It takes as input a directory comprising multiple migrate models, each in its own directory (named after the model), 
  # containing and outfile.txt for that model. 
  # This function will find the marginal likelihoods (Thermodynamic, bezier-corrected and harmonic) 
  # and harvest them into a data-frame.
  # NOTE! still need to transfer code to harvest single-locus outputs
  
  #initialize a data frame to take the values
  modelMarglikes <- data.frame(model=character(),
                               thermodynamic=numeric(),
                               bezier.corrected=numeric(), 
                               harmonic=numeric()) 
  # loop through directories in the working directory, each of which is name
  # after a different model
  for(i in list.dirs(modelsDir, full.names = F)[-1]){ #i<-"stepping.stone"
    modelDir<-file.path(modelsDir,i)
    print(modelDir)
    #scan in the outfile, separating at each newline
    outfile<-scan(file=file.path(modelDir,outfileName),what="character",sep="\n") 
    #find the line with the likelihoods on it and split on runs of spaces
    marglikeline <- grep("Scaling factor",outfile,value=F)-1
    marglikeline <- strsplit(outfile[marglikeline],
                             "\\s+", perl = T)[[1]][3:5]
    #  if(length(marglikeline)==0){next}
    marglikes <- c(i,marglikeline)
    
    modelMarglikes <- rbind(modelMarglikes,marglikes, deparse.level = 2)
  }
  names(modelMarglikes) <- c("model","thermodynamic","bezier.corrected","harmonic")
  modelMarglikes[2:4] <- sapply(modelMarglikes[2:4], as.numeric)
  return(modelMarglikes)
}

bf_calcs<-function(df,ml="bezier.corrected"){
  # This calculates log bayes factors on data frames output by
  # harvest.model.likelihoods(), following Johnson and Omland (2004)
  # You may choose the likelihood flavor with
  # ml = "bezier.corrected", "thermodynamic" or "harmonic"
  #df$thermodynamic <- as.numeric(df$thermodynamic)
  #df$bezier.corrected <- as.numeric(df$bezier.corrected)
  #df$harmonic <- as.numeric(df$harmonic)
  mlcol <- df[[ml]] 
  bmvalue <- mlcol[which.max(mlcol)]
  lbf <- 2*(mlcol-bmvalue)
  choice <- rank(-mlcol)
  modelprob <- exp(lbf/2)/sum(exp(lbf/2))
  dfall <- cbind(df,lbf,choice,modelprob)
  return(dfall)
}	

migrants_per_gen<-function(x){
  #a function for creating 4Nm vectors out of m and Theta vectors.
  m<-names(x)[which(grepl("M_",names(x)))] #names of m columns
  #theta<-names(x)[which(grepl("Theta_",names(x)))] #names of theta columns
  for(n in m){
    t<-paste("Theta",strsplit(n,split="_")[[1]][3],sep="_")
    x[,paste("Nm",strsplit(n,split="_")[[1]][2],strsplit(n,split="_")[[1]][3],sep="_")]<- x[,which(names(x)==n)]*x[,which(names(x)==t)] 
    #this hairy little statement makes a new column named "Nm_X_Y" and then fills it by multiplying the M_X_Y column by the Theta_Y column	
  }
  return(x)
}


# Functions for summarizing parameter estimates over multiple loci without overusing the prior
# The user need only call summarize_posterior(), which will call remove_prioro and sum_over_loci()
# Currently they only work for exponential priors

remove_prior <- function(densityd,prior,floor = 1e-10){
  # this function subtracts the prior from the posterior density for the y values of a single density column
  # and changes values less than zero to a small floor value
  # and then standardizes the result to sum to 1
  minusprior <- densityd - prior
  minusprior[minusprior <= 0] <- floor
  standardized <- minusprior / sum(minusprior)
  return(standardized)
}

sum_over_loci <- function(df,parameter){
  #this function takes a data frame of probability densities for many loci
  # that have had the prior removed, and also has a logged prior column called logPrior
  # as well as the name of a parameter (e.g. "Theta_1")
  # and sums the densities over loci.
  # Rmpfr package allows quadruple precision for calcs on very small numbers.
  require(Rmpfr)
  require(dplyr)
  require(magrittr)
  
  #log all values
  df2 <- df %>%  mutate(across(starts_with(c(parameter)),
                               .fns=log)) %>% 
    # convert the df to rowwise so that rows can be summed
    # and then sum across the row, including the prior
    rowwise() %>% 
    mutate(sum_param_prior = 
             sum(c_across(starts_with(c(parameter,"logPrior"))))) %>% 
    #convert back to a regular df
    ungroup()
  
  #need to convert to quadruple precision because 
  #these will exponentiate to very small numbers.
  sum_param_prior_exp <- exp(mpfr(df2$sum_param_prior, precBits = 128))
  # standardize by dividing by the sum of the column
  sum_param_prior_standardized <-
    sum_param_prior_exp/sum(sum_param_prior_exp)
  #drop the intermediate columns (no longer needed), change the standardized
  # output back to double precision so that it can be incorporated into the df
  # rename the summed column after the parameter
  df3 <- df2 %>% dplyr::select(-c(sum_param_prior)) %>%
    mutate(summed_over_loci =
             as.numeric(sum_param_prior_standardized)) %>% 
    rename_with(.fn = ~ paste(parameter), 
                .cols = summed_over_loci)
  return(df3)
}



summarize_posterior <- function(posterior, 
                                parameter_type = c("Theta","M","Nm","D"),
                                exponential_mean = 1,
                                lower.bound, upper.bound,
                                n = 2^10, floor = 1e-10){
  # this function takes a Migrate-n posterior "bayesallfile" as a dataframe
  # as well as one of the parameter types, 
  # and the mean of an exponential prior (currently only exponential priors supported)
  # The function will create densities for each parameter of the given type,
  # remove the prior from each, sum across loci, and re-add the prior (once)
  # upper.bound and lower.bound do not have to be the bounds of the prior
  # they can be used to focus on where most of the posterior density lies
  # n is passed to density() and is the number of points used for density calculation
  # I recommend using 2^10 to 2^15 points
  # floor is passed to remove_prior and should be a very small value used to replace negative values
  require(dplyr)
  require(magrittr)
  require(stringr)
  require(tidyr)
  require(purrr)
  
  parameters <- names(posterior) %>%
    str_subset(parameter_type)
  # create a tibble with the x values for all density plots 
  #  for a particular parameter type
  getx <- posterior %>% filter(Locus == 1) %>% 
    dplyr::select(parameters[1])
  
  # create a tibble with x values for a density plot
  #  of the chosen number of points
  dens <- tibble(.rows = n, x = density(getx[[1]],  n=n, 
                                        from =lower.bound, 
                                        to = upper.bound,
                                        bw = "nrd0")$x)
  print("calculating densities")
  # calculate densities for each parameter of a given type at each locus
  dens <- posterior %>% 
    dplyr::select(starts_with(c("Steps","Locus","rep",
                                paste0(parameter_type,"_")))) %>% 
    pivot_wider(names_from = "Locus", values_from = 
                  starts_with(paste0(parameter_type,"_")),
                names_sep = "-") %>% 
    dplyr::select(starts_with(paste0(parameter_type,"_"))) %>% 
    map_dfc(function(x) density(x, n = n, from = lower.bound,
                                to = upper.bound, 
                                bw = "nrd0")$y) %>%
    bind_cols(dens)
  
  # create, log and remove prior
  dens$prior <- dexp(dens$x, rate = 1/exponential_mean, 
                     log = F)
  
  dens$logPrior <- log(dens$prior)
  
  print("removing prior")
  dens2 <- dens %>% 
    #remove the prior, standardize
    mutate(across(starts_with(parameter_type), 
                  ~ remove_prior(densityd = .x,
                                 prior = dens$prior,
                                 floor = floor) ))
  
  dens3 <- dens2
  
  for(p in parameters){
    print(p)  
    dens3 <- sum_over_loci(df = dens3, parameter = p)
  }
  # trying to do the above loop with purrr:map_dfc
  #dens4 <- parameters %>% 
  #        map_dfc(.f = ~ sum_over_loci(df = dens2, parameter = .x))
  return(dens3)
}


# functions for extracting tables and figures from a dataframe created by summarize_posterior()
#
#

parameter_stats <- function(summarized_posterior, quantiles =  c(0.025, 0.25, 0.5, 0.75, 0.975), parameter){
  # this function calculates quantiles, plus mode, mean and sd for a single parameter (e.g. Theta_1)
  # from a data frame that has been summarized for a given set of parameters 
  # with summarize_posterior() from a migrate-n bayesallfile
  require(spatstat.explore)
  require(dplyr)
  require(magrittr)
  p <- summarized_posterior %>% dplyr::select(x,parameter) %>% as.list(bw = 6.33e-05)
  names(p) <- c("x", "y")
  p$bw <- 6.33e-5
  attr(p, "class") <- "density"
  qu <- quantile.density(p,quantiles)
  wmo <- p$x[which(p$y==max(p$y))]
  wme <- weighted.mean(p$x, p$y)
  wsd <- sqrt(weighted.var(p$x, p$y))
  stats <- tibble_row(quantiles=list(qu),mode = wmo, mean = wme, sd = wsd[,1]) %>% unnest_wider(quantiles)
  return(stats)
}

posterior_stats <- function(summarized_posterior, quantiles =  c(0.025, 0.25, 0.5, 0.75, 0.975), 
                            parameter_key, parameter_type = c("Theta","M","Nm","D")){
  # this function uses purrr to loop over all parameters in a dataframe (summarized_posterior)
  # that has been summarized for a given set of parameters with summarize_posterior() from a migrate-n bayesallfile
  # and calculate quantiles, plus mode, mean and sd for all of them using parameter_stats()
  # and rename them using parameter_key that has been created with make_parameter_key()
  require(magrittr)
  require(purrr)
  require(stringr)
  #find the parameter type among the *names* of the parameter key
  parameters <- parameter_key[grep(parameter_type, parameter_key)]  
  
  #setnames swaps the values and names of parameters, map creates a list output, list_rbind binds them into a tibble
  df <- parameters %>%  
        map(\(x) parameter_stats(thetas, quantiles = quantiles, parameter = x)) %>% list_rbind(names_to = "parameter")
  return(df)
}

plot_parameter_posteriors <- function(summarized_posterior, parameter_key,
                                      parameter_type = c("Theta","M","Nm","D"),
                                      x_limits = c(0,1), param_colors = NULL){
  #this function will create a ridge plot of parameters for a given parameter_type
  #set x_limits to where you want the boundaries of the plot to fall on the x-axis
  #param_colors can be a vector of colors the same length as the number of parameters
  #you are plotting. If it is shorter, the function will generate colors from the
  #turbo pallete of the viridis package.
  require(dplyr)
  require(ggplot2)
  require(ggridges)
  require(tidyr)
  
  parameters <- parameter_key[grep(parameter_type, parameter_key)]
  
  if(length(param_colors) < length(parameters)){
    require(viridis)
    param_colors <- turbo(n = length(parameters))
  }
  
  # pivot to a long-format data frame
  sp_long <- summarized_posterior %>% dplyr::select(all_of(c("x","prior",parameters))) %>% 
               pivot_longer(cols = !c(x,prior), 
               names_to = "Parameter",
               values_to = "Density")
  
 ridgeplot <- ggplot(sp_long, aes(x = x, y = Parameter, height = Density, 
                                         fill = Parameter)) +
              geom_density_ridges(stat = "identity") + 
              scale_fill_discrete(type = param_colors) +
              scale_y_discrete(labels = NULL, breaks = NULL) +
              labs(x = parameter_type) + xlim(x_limits) +
              guides(fill = guide_legend(reverse=T))
 return(ridgeplot)
}

# TO DO
#add coda-based functions to test for convergence
#add locus-by-locus plot function
#add single locus capability to harvester function
#figure out why plot function is plotting things out of order


