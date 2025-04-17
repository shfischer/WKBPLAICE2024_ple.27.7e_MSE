### ------------------------------------------------------------------------ ###
### script for running MSE locally (not on HPC) ####
### ------------------------------------------------------------------------ ###
suppressMessages(library(FLCore))
suppressMessages(library(FLasher))
suppressMessages(library(FLBRP))
suppressMessages(library(mse))
suppressMessages(library(FLfse))
suppressMessages(library(GA))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stockassessment))
suppressMessages(library(doParallel))

req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", 
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) 
  suppressMessages(library(package = i, character.only = TRUE))
source("funs.R")
source("funs_GA.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

cl1 <- FALSE

### ------------------------------------------------------------------------ ###
### parallelisation (optional) ####
### ------------------------------------------------------------------------ ###

if (FALSE) {
  req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", 
                 "GA", "doParallel", "doRNG",
                 "tidyr", "dplyr", "stockassessment")
  req_scripts <- c("funs.R", "funs_GA.R", "funs_WKNSMSE.R", "funs_OM.R")
  
  n_workers2 <- 10
  cl2 <- makeCluster(n_workers2)
  registerDoParallel(cl2)
  print(cl2)
  cl2_length <- length(cl2)
  . <- foreach(i = seq(n_workers2)) %dopar% {
    for (i in req_pckgs) 
      suppressPackageStartupMessages(
        library(package = i, character.only = TRUE, warn.conflicts = FALSE, 
                verbose = FALSE, quietly = TRUE))
    for (i in req_scripts) source(i)
  }
}

### ------------------------------------------------------------------------ ###
### list OMs ####
### ------------------------------------------------------------------------ ###
OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")


### ------------------------------------------------------------------------ ###
### all OMs - F=Fmsy & F=0 ####
### ------------------------------------------------------------------------ ###

. <- foreach (OM = OMs) %:% foreach(Ftrgt = list("MSY", 0)) %do% {
  #browser()
  print(paste0("OM=", OM, " - Ftrgt=", Ftrgt))
  
  if (identical(OM, "refset")) {
    n_blocks <- 10; n_workers <- 10; mp_parallel <- TRUE
  } else {
    n_blocks <- 1; n_workers <- 1; mp_parallel <- FALSE
  }
  
  stock_id <- "ple.27.7e_revision"
  scenario <- ""
  MP <- "constF"
  n_yrs <- 100
  ga_search <- FALSE
  save_MP <- TRUE
  
  args_local <- c("a=1")
  source("MP_run.R")
  
}

### ------------------------------------------------------------------------ ###
### SAM - all OMs ####
### ------------------------------------------------------------------------ ###
# args_local <- c("n_blocks=10", "n_workers=10", "mp_parallel=TRUE",
#                 "scenario=''", "MP='ICES_SAM'",
#                 "n_yrs=20", "check_file=FALSE",
#                 "ga_search=FALSE", "OM='baseline'", "save_MP=TRUE", 
#                 "collate=FALSE", "stat_yrs='multiple'"
#                 )
# source("MP_run.R")
# stopCluster(cl); rm(cl); gc()
# rm(args_local)
# ### other OMs
# alt_OMs <- c("Catch_no_disc", "Catch_no_surv", "migr_none", "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", "R_failure", "overcatch", "undercatch", "Idx_higher")
# for (OM in alt_OMs) {
#   print(paste0("OM=", OM))
#   args_local <- c(paste0("OM='", OM, "'"))
#   source("MP_run.R")
#   stopCluster(cl); rm(cl); gc()
# }

### ------------------------------------------------------------------------ ###
### optimised chr rule (MP5) from WKBPLAICE 2024 ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_x_w %>% 
  filter(MP == 5)

OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")

#. <- foreach(x = split(df_x_w, f = seq(nrow(df_x_w)))) %:% 
. <- foreach(x = split(df_x_w, f = seq(nrow(df_x_w)))[6]) %:% 
  foreach(OM = OMs) %do% {
    #browser()
    print(x)
    print(paste0("OM=", OM, ", index = ", x$index))
    
    if (identical(OM, "refset")) {
      n_blocks <- 10; n_workers <- 10; mp_parallel <- TRUE
    } else {
      n_blocks <- 1; n_workers <- 1; mp_parallel <- FALSE
    }
    
    ### stock
    stock_id <- "ple.27.7e_revision"
    
    ### chr params
    idxB_lag <- 1
    idxB_range_3 <- 2
    exp_b <- 1
    comp_b_multiplier <- 3.7
    interval <- 2
    multiplier <- 0.66
    upper_constraint <- 1.2
    lower_constraint <- 0.7
    biomass_index <- "UK-FSP"
    
    ### scenario directory
    scenario <- "multiplier"
    
    args_local <- c("MP='hr'",
                    "n_yrs=20", "check_file=FALSE",
                    "ga_search=TRUE", "MP='hr'", "save_MP=TRUE", 
                    "popSize=1", "maxiter=1",
                    "add_suggestions=FALSE", "collate=FALSE")
    source("MP_run.R")
    
    if (isTRUE(is(cl, "cluster"))) stopCluster(cl)
}

### ------------------------------------------------------------------------ ###
### optimised chr rule (MP5) corrected for data revisions  ####
### ------------------------------------------------------------------------ ###
### MP5 control parameters are defined based on values in specific year(s)
### I_trigger:
###   I_loss = I_2008
###   I_trigger = I_loss * 3.7
###   -> no change because index value did not change in this year
### HR_target:
###   HR (harvest rate) = (dead) catch / index
###   HR_target = HR(2003-2023) * 0.66
###   -> HR values changed because of revisions to catch and index
###   => adapt multiplier of 0.66



OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")

. <- foreach(OM = OMs) %do% {
    #browser()
    print(x)
    print(paste0("OM=", OM, ", index = ", x$index))
    
    if (identical(OM, "refset")) {
      n_blocks <- 10; n_workers <- 10; mp_parallel <- TRUE
    } else {
      n_blocks <- 1; n_workers <- 1; mp_parallel <- FALSE
    }
    
    ### stock
    stock_id <- "ple.27.7e_revision"
    
    ### chr params
    idxB_lag <- 1
    idxB_range_3 <- 2
    exp_b <- 1
    comp_b_multiplier <- 3.7
    interval <- 2
    multiplier <- 0.66/1.06682634070264015235807164572179317474365234375 #0.66
    upper_constraint <- 1.2
    lower_constraint <- 0.7
    biomass_index <- "UK-FSP"
    
    ### scenario directory
    scenario <- "multiplier"
    
    args_local <- c("MP='hr'",
                    "n_yrs=20", "check_file=FALSE",
                    "ga_search=FALSE", "MP='hr'", "save_MP=TRUE")
    source("MP_run.R")
    
    if (isTRUE(is(cl, "cluster"))) stopCluster(cl)
  }


### ------------------------------------------------------------------------ ###
### refset - MP5 100 years ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
n_blocks <- 7 
n_workers <- 7 
mp_parallel <- TRUE

### stock
stock_id <- "ple.27.7e_revision"
    
### chr params
idxB_lag <- 1
idxB_range_3 <- 2
exp_b <- 1
comp_b_multiplier <- 3.7
interval <- 2
multiplier <- 0.66/1.06682634070264015235807164572179317474365234375 #0.66
upper_constraint <- 1.2
lower_constraint <- 0.7
biomass_index <- "UK-FSP"
    
### scenario directory
scenario <- "multiplier"

args_local <- c("MP='hr'", "OM='refset'",
                "n_yrs=100", "check_file=FALSE",
                "ga_search=FALSE", "MP='hr'", "save_MP=TRUE", 
                "collate=FALSE")
source("MP_run.R")

