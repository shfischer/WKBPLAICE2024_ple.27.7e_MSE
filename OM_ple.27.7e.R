### ------------------------------------------------------------------------ ###
### create OM for western English Channel plaice ple.27.7e ####
### ------------------------------------------------------------------------ ###
### base OM on SAM model fit
### follow OM routines developed for ICES WKNSMSE 2018

library(ggplot2)
library(FLCore)
library(FLAssess)
library(FLXSA)
library(FLasher)
library(FLfse)
library(ggplotFL)
library(stockassessment)
library(foreach)
library(dplyr)
library(tidyr)
library(doParallel)
library(mse)
library(patchwork)

source("funs.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

### input data, including discard estimates
stk_data <- readRDS("input/ple.27.7e/preparation/model_input_stk.RDS")
idx_data <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
ALKs <- readRDS("input/ple.27.7e/preparation/ALK_MSE.rds")
refpts <- list(
  ### ICES style EqSim reference points (run with SAM fit)
  EqSim_Btrigger = 3265.991, EqSim_Fmsy = 0.2110553, EqSim_Fpa = 0.2436579, 
  EqSim_Bpa = 3265.991, EqSim_Blim = 2332.851,
  ### real OM MSY values
  Fmsy = 0.222, Bmsy = 5500, Cmsy = 1210, Blim = 2352,
  ### length reference points
  Lc = 25, ### from simulated length data
  Lref = 0.75*25 + 0.25*64.53349 ### Linf: quarterly ALK with 2019-2023 data
)
# round(min(ssbtable(fit)[, "Estimate"]))

### intermediate year catch advice
int_yr_catch <- 1219 ### from 2023 advice sheet

### baseline OM
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.5, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "baseline")

### alternative OMs - discard survival ####
### all discards die
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.0, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "Catch_no_surv")
### all discards survive
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 1, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "Catch_no_disc")




### alternative OMs - natural mortality M ####
### M_high: M +50%
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.5, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, 
          int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_high",
          M_alternative_mult = TRUE, M_alternative = 1.5)

### M_low: M -50%
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.5, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, 
          int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_low",
          M_alternative_mult = TRUE, M_alternative = 0.5)

### M_Gislason: age-dependent M according to Gislason et al. (2010)
M_Gislason <- readRDS("input/ple.27.7e/preparation/M_Gislason.rds")
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.5, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, 
          int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_Gislason",
          M_alternative_mult = FALSE, M_alternative = M_Gislason)

### alternative OMs - Recruitment ####
### R_no_AC: no auto-correlation in recruitment residuals
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 0.5, disc_survival_MP = 0.5,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch, 
          int_yr_catch_split = TRUE,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = FALSE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "R_no_AC")

### R_higher: recruitment +20%
### use files from baseline OM and adapt recruitment residuals
### instead of copying file, create a hard link (on Windows)
path_baseline <- "input/ple.27.7e/baseline/1000_100/"
# list.files(path_baseline)
files_link <- c("ALKs.rds", "catch_res.rds", "idx.rds", "idx_dev.rds", 
                "idx_dev_raw.rds", 
                "proc_res.rds", "SAM_conf.rds", "SAM_fit.rds",
                "SAM_initial.rds", "SAM_uncertainty.rds", 
                "stk.rds", "stk_oem.rds")
path_new <- "input/ple.27.7e/R_higher/1000_100/"
dir.create(path_new, recursive = TRUE)
file.link(from = paste0(path_baseline, files_link), 
          to = paste0(path_new, files_link))
### copy some other files so that they can be changed
files_copy <- c("refpts_mse.rds")
file.copy(from = paste0(path_baseline, files_copy), 
          to = paste0(path_new, files_copy))
### adapt recruitment model
sr_R_higher <- readRDS(paste0(path_baseline, "sr.rds"))
params(sr_R_higher)["a"] <- params(sr_R_higher)["a"] * 1.2
saveRDS(sr_R_higher, file = paste0(path_new, "sr.rds"))


### R_lower: recruitment -20%
path_new <- "input/ple.27.7e/R_lower/1000_100/"
dir.create(path_new, recursive = TRUE)
file.link(from = paste0(path_baseline, files_link), 
          to = paste0(path_new, files_link))
### copy some other files so that they can be changed
files_copy <- c("refpts_mse.rds")
file.copy(from = paste0(path_baseline, files_copy), 
          to = paste0(path_new, files_copy))
### adapt recruitment model
sr_R_lower <- readRDS(paste0(path_baseline, "sr.rds"))
params(sr_R_lower)["a"] <- params(sr_R_lower)["a"] * 0.8
saveRDS(sr_R_lower, file = paste0(path_new, "sr.rds"))


### R_failure: recruitment failure 2025-2029
path_new <- "input/ple.27.7e/R_failure/1000_100/"
dir.create(path_new, recursive = TRUE)
file.link(from = paste0(path_baseline, files_link), 
          to = paste0(path_new, files_link))
### copy some other files so that they can be changed
files_copy <- c("refpts_mse.rds")
file.copy(from = paste0(path_baseline, files_copy), 
          to = paste0(path_new, files_copy))
### adapt recruitment residuals
sr_R_failure <- readRDS(paste0(path_baseline, "sr.rds"))
residuals(sr_R_failure)[, ac(2025:2029)] <- 
  residuals(sr_R_failure)[, ac(2025:2029)] * 0.1
saveRDS(sr_R_failure, file = paste0(path_new, "sr.rds"))


### alternative OMs - migration ####

### migr_none: no migration
### load new stock object 
### - landings are 7.e catches, discards are 7.d catches
### - both assuming 50% discard survival
### -> discard rate of this stock object controls migration element
stk_migr <- readRDS("input/ple.27.7e/preparation/model_input_stk_migration_LD.rds")
### update intermediate year catch: 
### -> from 2022 (biennial) advice sheet, for 2024, for area 7.e: 1104
### -> correct for 50% discard survival (in 7.e)
int_yr_catch_7e <- 1104 * (1 - 0.08771373 * 0.5)
create_OM(stk_data = stk_migr, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2023, 
          disc_survival_OM = 1, disc_survival_MP = 1,
          int_yr_add = TRUE, 
          int_yr_catch = int_yr_catch_7e, 
          int_yr_catch_split = FALSE, ### don't split
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, 
          idxB = "UK-FSP", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs_sample = 2019:2023, 
          length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "migr_none")

### ------------------------------------------------------------------------ ###
### MSY reference points ####
### ------------------------------------------------------------------------ ###
### usually called from OM_MSY.pbs -> OM_MSY.R

if (FALSE) {
  ### set up parallel processing
  req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", "FLXSA",
                 "tidyr", "dplyr", "ggplot2")
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  ### load additional functions
  req_scripts <- c("funs.R", "funs_GA.R", "funs_WKNSMSE.R", "funs_OM.R")
  for (i in req_scripts) source(i)
  ### parallelisation with doFuture
  plan(multisession, workers = 5)
  ### load packages and functions into parallel workers
  . <- foreach(i = seq(5)) %dofuture% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    for (i in req_scripts) source(i)
  }
  
  ### baseline OM
  res <- est_MSY(OM = "baseline", stock_id = "ple.27.7e", yr_start = 2025, 
                 n_blocks = 5, n_iter = 1000)
  res$result[which.max(res$result$catch), ]
  #        Ftrgt     catch       ssb       tsb      rec
  # 15 0.1638845   1702.92  9536.053  11136.77 6542.729 paper
  # 24 0.2334596 1200.3999  5278.768  6147.824 6497.135 WKBPLAICE - MSY
  #  1 0.0000000    0.0000 20982.801 22037.208 7195.251 WKBPLAICE - unfished
  
  plan(sequential)
}

### ------------------------------------------------------------------------ ###
### update MSY reference points for alternative OMs ####
### ------------------------------------------------------------------------ ###

### find ratio of R(SSB=Blim)/R0 -> definition of Blim
stk_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
Blim <- min(iterMedians(ssb(stk_baseline)), na.rm = TRUE)
dimnames(stk_baseline)$year[which.min(iterMedians(ssb(stk_baseline)))]
### Blim is SSB in 2008
sr_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
RR0 <- c(((iterMedians(params(sr_baseline)["a"])*Blim) /
            (iterMedians(params(sr_baseline))["b"] + Blim)) /
           iterMedians(params(sr_baseline))["a"])
RR0
### Blim corresponds to SSB at ~ 79.5% of R0

refpts <- FLPar(refpts, iter = 1000, unit = "")
update_refpts <- function(OM, refpts, RR0) {
  ### get MSY levels 
  refpts_MSY <- readRDS(paste0("input/ple.27.7e/", OM,
                               "/1000_100/MSY_trace.rds"))
  refpts_MSY <- refpts_MSY[[which.max(sapply(refpts_MSY, function(x) x$catch))]]
  ### update
  refpts["Fmsy"] <- refpts_MSY$Ftrgt
  refpts["Bmsy"] <- refpts_MSY$ssb
  refpts["Cmsy"] <- refpts_MSY$catch
  ### load recruitment model and estimate Blim
  sr_mse <- readRDS(paste0("input/ple.27.7e/", OM, "/1000_100/sr.rds"))
  pars <- iterMedians(params(sr_mse))
  refpts["Blim"] <- c(pars["b"])*(RR0/(1 - RR0))
  print(refpts)
  ### save updated values
  saveRDS(refpts, file = paste0("input/ple.27.7e/", OM,
                                "/1000_100/refpts_mse.rds"))
}

### baseline
update_refpts(OM = "baseline", refpts = refpts, RR0 = RR0)
### alternative OMs
update_refpts(OM = "Catch_no_disc", refpts = refpts, RR0 = RR0)
update_refpts(OM = "Catch_no_surv", refpts = refpts, RR0 = RR0)
update_refpts(OM = "migr_none", refpts = refpts, RR0 = RR0)
update_refpts(OM = "M_low", refpts = refpts, RR0 = RR0)
update_refpts(OM = "M_high", refpts = refpts, RR0 = RR0)
update_refpts(OM = "M_Gislason", refpts = refpts, RR0 = RR0)
update_refpts(OM = "R_no_AC", refpts = refpts, RR0 = RR0)
update_refpts(OM = "R_higher", refpts = refpts, RR0 = RR0)
update_refpts(OM = "R_lower", refpts = refpts, RR0 = RR0)




### ------------------------------------------------------------------------ ###
### alternative Blim values ####
### ------------------------------------------------------------------------ ###
stk_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
Blim <- min(iterMedians(ssb(stk_baseline)), na.rm = TRUE)
sr_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
### find SSB at R=x
bh <- function(alpha, beta, B) (alpha * B)/(beta + B)
bh_min <- function(alpha, beta, B, R0, prop) {
  (bh(alpha = alpha, beta = beta, B = B) - prop * R0)^2
}
### 0.7R0
B0.7R0 <- optim(par = c(B = Blim), fn = bh_min, 
                alpha = c(iterMedians(params(sr_baseline)["a"])), 
                beta = c(iterMedians(params(sr_baseline)["b"])),
                prop = 0.7, R0 = c(iterMedians(params(sr_baseline))["a"]),
                method = "Brent", lower = 0, upper = 40000)
B0.7R0$par
### 0.3R0
B0.3R0 <- optim(par = c(B = Blim), fn = bh_min, 
                alpha = c(iterMedians(params(sr_baseline)["a"])), 
                beta = c(iterMedians(params(sr_baseline)["b"])),
                prop = 0.3, R0 = c(iterMedians(params(sr_baseline))["a"]),
                method = "Brent", lower = 0, upper = 40000)
B0.3R0$par
### many values
Blim_RR0 <- data.frame(RR0 = seq(0, 1, 0.001))
Blim_RR0$Blim <- sapply(Blim_RR0$RR0, function(x) {
  tmp <- optim(par = c(B = Blim), fn = bh_min, 
               alpha = c(iterMedians(params(sr_baseline)["a"])), 
               beta = c(iterMedians(params(sr_baseline)["b"])),
               prop = x, R0 = c(iterMedians(params(sr_baseline))["a"]),
               method = "Brent", lower = 0, upper = 40000)
  return(tmp$par)
})
saveRDS(Blim_RR0, "input/ple.27.7e/baseline/1000_100/Blim_RR0.rds")



### ------------------------------------------------------------------------ ###
### for harvest rate: check mean catch length history ####
### ------------------------------------------------------------------------ ###

### load stock
stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")

### indices 
### use observed values - equivalent to simulated plus added uncertainty
idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
idx$`FSP-7e`@index ### 2003-2020

### aggregated biomass index
idxB <- quantSums(idx$`FSP-7e`@index * catch.wt(stk)[ac(2:8), ac(2003:2020)])
plot(idxB) + ylim(c(0, NA))
### corresponding catch
idxC <- catch(stk)[, ac(2003:2020)]

### harvest rate
plot(idxC/idxB) + ylim(c(0, NA))

### get mean catch length from WGCSE 2021
Lc <- 26
LFeM <- 36
lmean <- read.csv("input/ple.27.7e/preparation/lmean.csv")
### always below LFeM

### plot mean length
ggplot() +
  geom_hline(yintercept = LFeM, size = 0.4, colour = "red") +
  geom_line(data = lmean, aes(x = Year, y = Lmean),
            size = 0.3) +
  ylim(c(0, NA)) + xlim(c(2010, 2020)) +
  labs(y = "mean catch length [cm]") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_ple_mean_length.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_mean_length.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)

### plot harvest rate
df_hr <- as.data.frame(idxC/idxB)
ggplot() +
  geom_line(data = df_hr, aes(x = year, y = data),
            size = 0.3) +
  geom_point(data = df_hr %>% filter(year %in% 2014), 
             aes(x = year, y = data),
             size = 0.5, colour = "red") +
  ylim(c(0, NA)) + xlim(c(2010, 2020)) +
  labs(y = "harvest rate (catch/index)") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_ple_mean_length_hr_target.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_mean_length_hr_target.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)


