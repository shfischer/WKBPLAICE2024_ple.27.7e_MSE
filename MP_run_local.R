### ------------------------------------------------------------------------ ###
### script for running MSE locally (not on HPC) ####
### ------------------------------------------------------------------------ ###
suppressMessages(library(FLCore))
suppressMessages(library(FLasher))
suppressMessages(library(FLBRP))
suppressMessages(library(mse))
suppressMessages(library(FLfse))
suppressMessages(library(FLXSA))
suppressMessages(library(GA))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stockassessment))
suppressMessages(library(doParallel))

req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", "FLXSA",
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
### baseline - hr - x & w & v & n0 ####
### ------------------------------------------------------------------------ ###
args_local <- c("n_blocks=1", "n_workers=0", 
                "scenario='mult_comp_b_mult'", "MP='hr'",
                "n_yrs=20", "check_file=FALSE",
                "ga_search=TRUE", "OM='baseline'", "save_MP=TRUE", 
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "idxB_lag=1", "idxB_range_3=1", "exp_b=1", 
                "comp_b_multiplier=0.66", "interval=1", 
                "multiplier=0.64",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")

rm(args_local)
args_local <- "n_blocks=1"
idxB_range_3 <- 1
interval <- 2
comp_b_multiplier <- 0.68
multiplier <- 0.61
source("MP_run.R")

idxB_range_3 <- 2
interval <- 1
comp_b_multiplier <- 0.92
multiplier <- 0.59
source("MP_run.R")

idxB_range_3 <- 2
interval <- 2
comp_b_multiplier <- 0.87
multiplier <- 0.58
source("MP_run.R")



### ------------------------------------------------------------------------ ###
### reference set - hr optimised multiplier ####
### ------------------------------------------------------------------------ ###
res <- readRDS("output/ple.27.7e/refset/1000_20/multiplier/hr/multiplier-upper_constraint1.2-lower_constraint0.7--obj_ICES_res_11-20.rds")
res@solution[, "multiplier"]
### 0.6
### but slightly above 5% risk -> use 0.59
args_local <- c("n_blocks=10", "n_workers=10", "mp_parallel=TRUE",
                "scenario='multiplier'", "MP='hr'",
                "n_yrs=20", "check_file=FALSE",
                "ga_search=TRUE", "OM='refset'", "MP='hr'", "save_MP=TRUE", 
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "idxB_lag=1", "idxB_range_3=1", "exp_b=1", 
                "comp_b_multiplier=1.4", "interval=1", 
                "multiplier=0.59",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")

### with Q1SWBeam index
rm(args_local)
args_local <- "a=1"
biomass_index <- "Q1SWBeam"
scenario <- "multiplier_Q1SWBeam"
multiplier <- 0.75
source("MP_run.R")

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
alt_OMs <- c("Catch_no_disc", "Catch_no_surv", "migr_none", "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", "R_failure", "overcatch", "undercatch", "Idx_higher")
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
### baseline - hr - optimise multiplier - Q1SWBeam ####
### ------------------------------------------------------------------------ ###
args_local <- c("n_blocks=1", "n_workers=15", "ga_parallel=TRUE",
                "scenario='multiplier_Q1SWBeam'", "MP='hr'",
                "biomass_index='Q1SWBeam'",
                "n_yrs=20", "check_file=FALSE",
                "ga_search=TRUE", "OM='baseline'", "MP='hr'", "save_MP=FALSE", 
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "idxB_lag=1", "idxB_range_3=1", "exp_b=1", 
                "comp_b_multiplier=1.4", "interval=1", 
                "multiplier=seq(0,2,0.01)",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")

### ------------------------------------------------------------------------ ###
### refset x & w -> all OMs ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)

OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")

. <- foreach(x = split(df_x_w, f = seq(nrow(df_x_w)))) %:% 
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
### refset - MP5 100 years ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
pars <- df_x_w[df_x_w$MP == 5, ]

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

df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)

OMs <- c("refset", "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")

. <- foreach(x = split(df_x_w, f = seq(nrow(df_x_w)))) %:% 
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

