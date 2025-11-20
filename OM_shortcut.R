### ------------------------------------------------------------------------ ###
### SAM shortcut MP preparation ####
### ------------------------------------------------------------------------ ###

library(ggplot2)
library(FLCore)
library(FLAssess)
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


### input data, including discard estimates
stk_data <- readRDS("input/ple.27.7e/preparation/model_input_stk.rds")
idx_data <- readRDS("input/ple.27.7e/preparation/model_input_idx.rds")

### apply discard survival
disc_survival_OM <- 0.5
discards.n(stk_data)[is.na(discards.n(stk_data))] <- 0
discards.wt(stk_data)[is.na(discards.wt(stk_data))] <- 0
discards.n(stk_data)[] <- discards.n(stk_data) * (1 - disc_survival_OM)
discards(stk_data) <- computeDiscards(stk_data)
catch(stk_data) <- computeCatch(stk_data, slot = "all")

### fit SAM
fit <- FLR_SAM(stk_data, idx_data, conf = NULL, conf_full = FALSE,
               idx_weight = FALSE, NA_rm = TRUE)

### ------------------------------------------------------------------------ ###
### SAM retro - define estimator error ####
### ------------------------------------------------------------------------ ###
### assessment estimate error
### based on retro analysis

retro <- retro(fit = fit, year = 10, ncores = 10)

opar <- par()
par(mar = c(2, 4.5, 0.5, 0.5))
plot(retro)
par(opar)
png(filename = "output/plots/shortcut/preparation/SAM_retro.png", 
    width = 20, height = 12, units = "cm", res = 300)
par(mar = c(2, 4.5, 0.5, 0.5))
plot(retro)
dev.off()

ssb_table <- function(x) {
  tmp <- data.frame(SSB = ssbtable(x)[, "Estimate"])
  tmp$year <- as.numeric(rownames(tmp))
  rownames(tmp) <- NULL
  tmp$assessment <- max(tmp$year)
  return(tmp)
}
retro_SSB <- lapply(retro, ssb_table)
retro_SSB <- do.call(rbind, retro_SSB)

### plot SSB retro
df <- rbind(retro_SSB, ssb_table(fit))
p_ssb <- df %>% filter(assessment < 2023) %>%
  ggplot(aes(x = year, y = SSB/1000, colour = as.factor(assessment))) +
  geom_line(data = df %>% filter(assessment == 2023),
            colour = "black") +
  geom_line(show.legend = FALSE, linewidth = 0.2) +
  theme_bw(base_size = 8) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Year", y = "SSB (1000 t)")
p_ssb
ggsave(filename = "output/plots/shortcut/preparation/SAM_retro_SSB.png", 
       width = 8.5, height = 4, units = "cm", dpi = 600)
p_ssb + xlim(c(2000, NA))
ggsave(filename = "output/plots/shortcut/preparation/SAM_retro_SSB_zoom.png", 
       width = 8.5, height = 4, units = "cm", dpi = 600)

### SSB error
SSB_error <- retro_SSB %>%
  filter(year == assessment) %>%
  mutate(SSB_retro = SSB, 
         SSB = NULL, 
         assessment = NULL) %>%
  left_join(ssb_table(fit) %>% 
              mutate(assessment = NULL)) %>%
  mutate(SSB_ratio = SSB_retro/SSB)
sd(SSB_error$SSB_ratio)
# 0.1400751
### plot
trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  scales::trans_new("from", trans, inv, 
            domain = c(from, Inf))
}
p_SSB_error <- SSB_error %>%
  select(year, SSB_ratio) %>%
  ggplot(aes(x = year, y = SSB_ratio)) +
  geom_col(fill = "darkgrey", colour = "darkgrey") +
  theme_bw(base_size = 8) +
  labs(x = "Year", y = "Terminal SSB / SSB") +
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = trans_from(),
                     limits = c(NA, NA), breaks = c(0.8, 0.9, 1, 1.1, 1.2, 1.3)) +
  scale_x_continuous(breaks = seq(2014, 2022, 2))
p_SSB_error
ggsave(filename = "output/plots/shortcut/preparation/SSB_retro_res.png", 
       width = 8.5, height = 4, units = "cm", dpi = 600)

### auto-correlation
SSB_err_acf <- acf(SSB_error$SSB_ratio)
plot(SSB_err_acf)
c(SSB_err_acf$acf)[2]
# 0.1041058
### plot
p_acf <- data.frame(acf = SSB_err_acf$acf, lag = seq(SSB_err_acf$n.used) - 1) %>%
  filter(lag != 0) %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_col(fill = "darkgrey", colour = "darkgrey") +
  theme_bw(base_size = 8) +
  labs(x = "Lag", y = "Auto-correlation") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = 1:9)
p_acf
ggsave(filename = "output/plots/shortcut/preparation/SSB_retro_acf.png", 
       width = 8.5, height = 4, units = "cm", dpi = 600)

### combine plots
p1 <- p_ssb +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  coord_cartesian(xlim = c(2010, 2024), ylim = c(0, 9.9), expand = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2 <- SSB_error %>%
  select(year, SSB_ratio) %>%
  ggplot(aes(x = year, y = SSB_ratio)) +
  annotate(geom = "rect", xmin = 2000, xmax = 2030, 
           ymin = 1 - sd(SSB_error$SSB_ratio), 
           ymax = 1 + sd(SSB_error$SSB_ratio),
           fill = "red", alpha = 0.1) +
  geom_linerange(aes(ymin = 1, ymax = SSB_ratio)) +
  geom_point(aes(colour = as.factor(year)), shape = 19, show.legend = FALSE) +
  theme_bw(base_size = 8) +
  labs(x = "Year", y = "Terminal SSB / SSB") +
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = trans_from(),
                     breaks = c(0.8, 0.9, 1, 1.1, 1.2, 1.3)) +
  scale_x_continuous(breaks = seq(2012, 2022, 2)) +
  coord_cartesian(xlim = c(2010, 2024), ylim = c(0.74, 1.35), expand = FALSE)

p1/p2
ggsave(filename = "output/plots/shortcut/preparation/retro_smry.png", 
       width = 8.5, height = 6, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### prepare residuals ####
### ------------------------------------------------------------------------ ###

### get stock template
stk <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
n <- dim(stk)[6]

### create auto-correlated SSB error
set.seed(3)
ssb_res <- rlnoise(n, ssb(stk) %=% 0, sd = sd(SSB_error$SSB_ratio), 
                   b = c(SSB_err_acf$acf)[2])
plot(ssb_res)
### replicate SSB error into age structure
n_res <- stock.n(stk) %=% NA_real_
n_res[] <- rep(c(ssb_res), each = dim(stk)[1])

### save for all OMs
OMs <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", "M_low", "M_high", "M_Gislason", "R_no_AC", "R_higher", "R_lower", "R_failure")
for (OM_i in OMs) {
  saveRDS(n_res, file = paste0("input/ple.27.7e/", OM_i, 
                               "/1000_100/shortcut_n_res.rds"))
}


### ------------------------------------------------------------------------ ###
### trial runs ####
### ------------------------------------------------------------------------ ###
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse",
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) 
  suppressMessages(library(package = i, character.only = TRUE))

### load additional functions
req_scripts <- c("funs.R", "funs_GA.R", "funs_WKNSMSE.R", "funs_OM.R")
for (i in req_scripts) source(i, local = TRUE)

input <- input_mp(stock_id = "ple.27.7e", OM = "baseline", n_iter = 1000,
                  n_yrs = 20, yr_start = 2025, n_blocks = 1,
                  MP = "ICES_SAM_shortcut")
refpts <- input_refpts(stock_id = "ple.27.7e", OM = "baseline", n_iter = 1000)

res_mp <- do.call(mp, input)
