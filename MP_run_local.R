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
         "R_no_AC", "R_higher", "R_lower", "R_h_lower",
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
args_local <- c("n_blocks=10", "n_workers=10", "mp_parallel=TRUE",
                "scenario=''", "MP='ICES_SAM'",
                "n_yrs=20", "check_file=FALSE",
                "ga_search=FALSE", "OM='baseline'", "save_MP=TRUE", 
                "collate=FALSE", "stat_yrs='multiple'"
                )
source("MP_run.R")
stopCluster(cl); rm(cl); gc()
rm(args_local)
### other OMs
alt_OMs <- c("Catch_no_disc", "Catch_no_surv", "migr_none", "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", "R_h_lower", "R_failure", "overcatch", "undercatch", "Idx_higher")
for (OM in alt_OMs) {
  print(paste0("OM=", OM))
  args_local <- c(paste0("OM='", OM, "'"))
  source("MP_run.R")
  stopCluster(cl); rm(cl); gc()
}

### ------------------------------------------------------------------------ ###
### rfb (default: multiplier=0.95) - all OMs ####
### ------------------------------------------------------------------------ ###
args_local <- c("scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "ga_search=FALSE", "OM='baseline'", "save_MP=TRUE", 
                "collate=FALSE", "stat_yrs='multiple'"
)
source("MP_run.R")
rm(args_local)
### other OMs
alt_OMs <- c("Catch_no_disc", "Catch_no_surv", "migr_none", "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", "R_failure", "overcatch",
             "undercatch", "Idx_higher")
for (OM in alt_OMs) {
  print(paste0("OM=", OM))
  args_local <- c(paste0("OM='", OM, "'"))
  source("MP_run.R")
}

### ------------------------------------------------------------------------ ###
### refset x & w -> all OMs ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_optima <- readRDS("output/paper/refset_x_w_grid_opt.rds")

OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower",
         "R_h_lower", "R_failure", "overcatch", "undercatch", "Idx_higher")

. <- foreach(x = split(df_optima, f = seq(nrow(df_optima)))) %:% 
  foreach(OM = OMs) %do% {
    #browser()
    print(x)
    print(paste0("OM=", OM, ", index = ", x$index))
    
    if (identical(OM, "refset")) {
      n_blocks <- 10; n_workers <- 10; mp_parallel <- TRUE
    } else {
      n_blocks <- 1; n_workers <- 1; mp_parallel <- FALSE
    }
    
    ### chr params
    idxB_lag <- x$idxB_lag
    idxB_range_3 <- x$idxB_range_3
    exp_b <- x$exp_b
    comp_b_multiplier <- x$comp_b_multiplier
    interval <- x$interval
    multiplier <- x$multiplier
    upper_constraint <- x$upper_constraint
    lower_constraint <- x$lower_constraint
    biomass_index <- x$index
    
    #OM <- OM
    ### scenario directory
    scenario <- ifelse(identical(x$index, "Q1SWBeam"), 
                       "multiplier_Q1SWBeam", "multiplier") 
    
    args_local <- c("MP='hr'",
                    "n_yrs=20", "check_file=FALSE",
                    "ga_search=TRUE", "MP='hr'", "save_MP=TRUE", 
                    "popSize=1", "maxiter=1",
                    "add_suggestions=FALSE", "collate=FALSE")
    source("MP_run.R")
    
    if (isTRUE(is(cl, "cluster"))) stopCluster(cl)
}

### ------------------------------------------------------------------------ ###
### refset - CHR2 100 years ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_optima <- readRDS("output/paper/refset_x_w_grid_opt.rds")
pars <- df_optima[df_optima$MP == 2, ]

n_blocks <- 7 
n_workers <- 7 
mp_parallel <- TRUE
    
### chr params
idxB_lag <- pars$idxB_lag
idxB_range_3 <- pars$idxB_range_3
exp_b <- pars$exp_b
comp_b_multiplier <- pars$comp_b_multiplier
interval <- pars$interval
multiplier <- pars$multiplier
upper_constraint <- pars$upper_constraint
lower_constraint <- pars$lower_constraint
biomass_index <- pars$index
    
#OM <- OM
### scenario directory
scenario <- "multiplier"

args_local <- c("MP='hr'", "OM='refset'",
                "n_yrs=100", "check_file=FALSE",
                "ga_search=TRUE", "MP='hr'", "save_MP=TRUE", 
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")


### ------------------------------------------------------------------------ ###
### sensitivity - index uncertainty ####
### ------------------------------------------------------------------------ ###

df_optima <- readRDS("output/paper/refset_x_w_grid_opt.rds")

OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_h_lower",
         "R_lower", "R_failure", "overcatch", "undercatch", "Idx_higher")

. <- foreach(x = split(df_optima, f = seq(nrow(df_optima)))) %:% 
  foreach(OM = OMs[1]) %do% {
    #browser()
    print(x)
    print(paste0("OM=", OM, ", index = ", x$index))
    
    if (identical(OM, "refset")) {
      n_blocks <- 10; n_workers <- 10; mp_parallel <- TRUE
    } else {
      n_blocks <- 1; n_workers <- 1; mp_parallel <- FALSE
    }
    
    ### chr params
    idxB_lag <- x$idxB_lag
    idxB_range_3 <- x$idxB_range_3
    exp_b <- x$exp_b
    comp_b_multiplier <- x$comp_b_multiplier
    interval <- x$interval
    multiplier <- x$multiplier
    upper_constraint <- x$upper_constraint
    lower_constraint <- x$lower_constraint
    biomass_index <- x$index
    
    OM <- OM
    scenario <- "sensitivity_idx"
    
    idx_unc <- seq(0, 2, 0.1)
    
    args_local <- c("ga_search=FALSE")
    source("MP_run.R")
    
    if (isTRUE(is(cl, "cluster"))) stopCluster(cl)
}

