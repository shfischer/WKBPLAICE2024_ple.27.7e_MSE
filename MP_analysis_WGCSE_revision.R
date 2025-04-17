### ------------------------------------------------------------------------ ###
### analyse MSE results ####
### ------------------------------------------------------------------------ ###
library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(foreach)
library(doParallel)
library(doFuture)
source("funs.R")
source("funs_GA.R")
source("funs_analysis.R")
source("funs_OM.R")

### ------------------------------------------------------------------------ ###
### define operating model (OM) names, reference and robustness set ####
### ------------------------------------------------------------------------ ###
### list of all operating models
OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")
OMs_label <- c("Reference set\n(combined)", 
               "Baseline", "Catch:\nno discards", "Catch:\n100% discards", 
               "Catch:\nno migration", 
               "M: -50%", "M: +50%", "M: Gislason", 
               "R: no AC", "R: +20%", "R: -20%", 
               "R: failure", "Catch: +10%", "Catch: -10%", 
               "Uncertainty:\nindex +20%")

OMs_group <- c("refset (combined)", rep("refset", 7), rep("robset", 7))

OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")

### ------------------------------------------------------------------------ ###
### compare harvest rate before and after WGCSE 2025 revision ####
### ------------------------------------------------------------------------ ###

path_input_original <- "input/ple.27.7e/baseline/1000_100/"
path_input_revision <- "input/ple.27.7e_revision/baseline/1000_100/"

### load observed stock (50% discard survival)
stk_oem_original <- readRDS(paste0(path_input_original, "stk_oem.rds"))
stk_oem_revision <- readRDS(paste0(path_input_revision, "stk_oem.rds"))
### extract catch
catch_oem_original <- iter(window(catch(stk_oem_original), 
                                  start = 2003, end = 2023), 1)
catch_oem_revision <- iter(window(catch(stk_oem_revision), 
                                  start = 2003, end = 2023), 1)

### load indices
idx_original <- readRDS(paste0(path_input_original, "idx.rds"))
idx_dev_original <- readRDS(paste0(path_input_original, "idx_dev.rds"))
idx_revision <- readRDS(paste0(path_input_revision, "idx.rds"))
idx_dev_revision <- readRDS(paste0(path_input_revision, "idx_dev.rds"))

### derive (observed) biomass index
idxB_original <- quantSums(index(idx_original$`UK-FSP`) *
                             idx_dev_original$`UK-FSP` * 
                             catch.wt(idx_original$`UK-FSP`))
idxB_original <- iter(window(idxB_original, end = 2023), 1)
idxB_revision <- quantSums(index(idx_revision$`UK-FSP`) *
                             idx_dev_revision$`UK-FSP` * 
                             catch.wt(idx_revision$`UK-FSP`))
idxB_revision <- iter(window(idxB_revision, end = 2023), 1)
idxB_revision/idxB_original

### calculate (relative) harvest rate
hr_original <- catch_oem_original/idxB_original
hr_revision <- catch_oem_revision/idxB_revision

hr_ref_original <- mean(hr_original) * 0.66
hr_ref_revision <- mean(hr_revision) * 0.66

hr_ratio <- hr_ref_revision/hr_ref_original
hr_ratio
# 1.066826
# 1.06682634070264015235807164572179317474365234375

### ------------------------------------------------------------------------ ###
### violin plots for final MP with revised OMs ####
### ------------------------------------------------------------------------ ###

### MP parameters
df <- data.frame(
  idxB_lag = 1, 
  idxB_range_3 = 2, 
  exp_b = 1, 
  comp_b_multiplier = 3.7, 
  interval = 2, 
  multiplier = 0.618657390447733, 
  upper_constraint = 1.2, 
  lower_constraint = 0.7,
  idx_unc = 1, 
  index = "UK-FSP",
  MP = 5
)
df <- df %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"), "-1_UK-FSP",
                       ".rds"))

### get stats
stats <- foreach(i = split(df, f = 1), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows)  %:%
  foreach(period = c("long-term", "short-term", "all"),
          period_yrs = list(2035:2044, 2025:2034, 2025:2044),
          .combine = bind_rows) %do% {
  #browser()
  print(paste("OM =", OM, "period = ", period))
  ### get projection
  path_i <- paste0("output/ple.27.7e_revision/", OM, "/1000_20/",
                   ifelse(identical(i$index, "Q1SWBeam"),
                          "multiplier_Q1SWBeam", "multiplier"),
                   "/hr/")
  mp_i <- readRDS(paste0(path_i, i$file))
  stk <- mp_i@om@stock
  
  ### get reference points
  refpts <- input_refpts(stock_id = "ple.27.7e_revision", OM = OM)
  
  ### extract metrics
  yr_min <- min(period_yrs)
  yr_max <- max(period_yrs)
  stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
  stk <- window(stk, start = yr_min, end = yr_max)
  ssb_i <- c(ssb(stk)/refpts["Bmsy"])
  catch_i <- c(catch(stk)/refpts["Cmsy"])
  fbar_i <- c(fbar(stk)/refpts["Fmsy"])
  risk_i <- c(apply(ssb(stk) < rep(c(refpts["Blim"]), 
                                   each = dim(ssb(stk))[2]), 2, mean))
  icv_i <- c(iav(catch(stk_icv), period = i$interval))
  icv_annual_i <- c(iav(catch(stk_icv), period = 1))
  ### combine
  df_i <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = catch_i, metric = "catch"),
                              data.frame(val = fbar_i, metric = "Fbar"),
                              data.frame(val = icv_i, metric = "ICV"),
                              data.frame(val = icv_annual_i, 
                                         metric = "ICV_annual"),
                              data.frame(val = risk_i, metric = "risk")
  ))
  df_i <- df_i %>%
    mutate(n1 = i$idxB_range_3,
           v = i$interval,
           x = i$multiplier,
           w = i$comp_b_multiplier,
           index = i$index,
           #group = i$group,
           OM = OM, OM_group = OM_group,
           #optimum = i$optimum,
           period = period,
           MP = i$MP
           )
  return(df_i)
}
saveRDS(stats, file = "output/revision_MP5_stats.rds")
# stats <- readRDS("output/revision_MP5_stats.rds")

### go through all solutions and plots stats
stats_plot <- stats %>%
  mutate(OM = factor(OM, levels = OMs,
                     labels = OMs_label),
         OM_group = factor(OM_group,
                           levels = OMs_group,
                           labels = c("", 
                                      rep("Reference set", 7), 
                                      rep("Robustness set", 7)))) %>%
  mutate(group = paste0("MP", MP, " - ", index, " - ",
                        case_when(v == 1 ~ "annual",
                                  v == 2 ~ "biennial"),
                        " - ",
                        case_when(w == 1.4 ~ "x",
                                  w != 1.4 ~ "x & w"),
                        " - ",
                        case_when(period == "long-term" ~ "long term",
                                  period == "short-term" ~ "short term",
                                  period == "all" ~ "all years"))) %>%
  mutate(group = paste0(group, " - revised OMs")) %>%
  mutate(group_label = paste0(index, "_",
                              case_when(v == 1 ~ "annual",
                                        v == 2 ~ "biennial"),
                              "_",
                              case_when(w == 1.4 ~ "x",
                                        w != 1.4 ~ "x_w"),
                              "_", 
                              case_when(period == "long-term" ~ "long",
                                        period == "short-term" ~ "short",
                                        period == "all" ~ "all")))


cols <- scales::hue_pal()(15)

. <- foreach(group_i = unique(stats_plot$period)) %do% {
  #browser()
  stats_plot_i <- stats_plot %>%
    filter(period == group_i)
  title_i <- stats_plot_i$group[1]
  file_i <- stats_plot_i$group_label[1]
  risk_max <- stats_plot_i %>%
    filter(metric == "risk") %>%
    filter(val == max(val)) %>%
    select(val) %>% unlist()
  risk_max <- ifelse(risk_max <= 0.3, 0.3, NA)
  p_risk <- stats_plot_i %>%
    filter(metric == "risk") %>%
    ggplot() +
    geom_col(data = . %>%
               group_by(OM, OM_group) %>%
               summarise(val = max(val)),
             aes(x = OM, y = val, fill = OM),
             show.legend = FALSE, width = 0.8, colour = "black", size = 0.2,
             position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(x = OM, y = val),
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    geom_hline(yintercept = 0.05, colour = "red", size = 0.4, 
               linetype = "1111") +
    stat_summary(aes(x = OM, y = val),
                 fun = "mean", geom = "point", shape = 4, size = 1,
                 stroke = 0.25) +
    scale_fill_manual("", values = cols) +
    facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
    labs(y = expression(max.~B[lim]~risk), title = title_i) +
    coord_cartesian(ylim = c(0, risk_max)) +
    theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  #p_risk
  p_catch <- stats_plot_i %>%
    filter(metric == "catch") %>%
    ggplot(aes(x = OM, y = val)) +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    stat_summary(aes(x = OM, y = val),
                 fun = "mean", geom = "point", shape = 4, size = 1,
                 stroke = 0.25) +
    geom_hline(yintercept = 1, colour = "#ebebeb", linewidth = 0.4,
               linetype = "1111") +
    scale_fill_manual("", values = cols) +
    facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
    labs(y = expression(Catch/MSY)) +
    coord_cartesian(ylim = c(0, 2.5)) +
    theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_blank())
  #p_catch
  p_ssb <- stats_plot_i %>%
    filter(metric == "SSB") %>%
    ggplot(aes(x = OM, y = val)) +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    stat_summary(aes(x = OM, y = val),
                 fun = "mean", geom = "point", shape = 4, size = 1,
                 stroke = 0.25) +
    geom_hline(yintercept = 1, colour = "#ebebeb", linewidth = 0.4,
               linetype = "1111") +
    scale_fill_manual("", values = cols) +
    facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
    labs(y = expression(SSB/B[MSY])) +
    coord_cartesian(ylim = c(0, 2.5)) +
    theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_blank())
  #p_ssb
  p_icv <- stats_plot_i %>%
    filter(metric == "ICV") %>%
    ggplot(aes(x = OM, y = val)) +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    stat_summary(aes(x = OM, y = val),
                 fun = "mean", geom = "point", shape = 4, size = 1,
                 stroke = 0.25) +
    scale_fill_manual("", values = cols) +
    facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
    labs(y = "ICV") +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          strip.text.x = element_blank())
  #p_icv
  p <- p_risk / p_catch / p_ssb / p_icv
  #p
  ggsave(filename = paste0("output/plots_revision/MP/refset_stats_MP5_revised_",
                           file_i, ".png"), 
         plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
         type = "cairo", bg = "white")
  ggsave(filename = paste0("output/plots_revision/MP/refset_stats_MP5_revised_",
                           file_i, ".pdf"), 
         plot = p, width = 16, height = 10, units = "cm", bg = "white")
}


### ------------------------------------------------------------------------ ###
### compare trajectories - MP5 - original vs revised OMs ####
### ------------------------------------------------------------------------ ###

### get projection
mp_list <- list(
  WKBPLAICE = readRDS(paste0("output/ple.27.7e/refset/1000_20/multiplier/hr/",
                             "mp_1_2_1_3.7_2_0.66_1.2_0.7.rds")),
  WGCSE2025 = readRDS(paste0("output/ple.27.7e_revision/refset/1000_20/",
                             "multiplier/hr/",
                             "mp_1_2_1_3.7_2_0.618657390447733_1.2_0.7-1_UK-FSP.rds"))
)
stk_list <- lapply(mp_list, function(x) x@om@stock)

### history
input_list <- list(
  WKBPLAICE = input_mp(stock_id = "ple.27.7e", OM = "refset", n_yrs = 20, 
                       MP = "hr"),
  WGCSE2025 = input_mp(stock_id = "ple.27.7e_revision", OM = "refset", n_yrs = 20, 
                       MP = "hr")
)
stk_hist <- lapply(input_list, function(x) x$om@stock)
  input$om@stock
stk_hist <- list(stk_hist, stk_hist)

refpts <- input_refpts(OM = "refset")
refpts[] <- NA

#debugonce(plot_worm_comparison)
p <- plot_worm_comparison(stk = stk_list, stk_hist = stk_hist, 
                          names = c("WKBPLAICE 2024",
                                    "WGCSE 2025\n(revision)"), 
                          refpts = refpts)
p

ggsave(filename = "output/plots_revision/wormplots/refset_comparison_revision.png",
       width = 16, height = 8, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots_revision/wormplots/refset_comparison_revision.pdf",
       width = 16, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### OM reference set trajectory ####
### ------------------------------------------------------------------------ ###

### plot refset
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_refset_label <- c("Baseline", "Catch:\nno discards", 
                      "Catch:\n100% discards", 
                      "Catch:\nno migration", 
                      "M: -50%", "M: +50%", "M: Gislason")

### get projection
path_i <- paste0("output/ple.27.7e_revision/", OMs_refset, "/1000_20/", 
                 "multiplier/hr/")
stk <- lapply(path_i, function(y) {
  readRDS(paste0(y, "mp_1_2_1_3.7_2_0.618657390447733",
                 "_1.2_0.7-1_UK-FSP.rds"))@om@stock
})

### historical stock
stk_hist <- lapply(OMs_refset, function(y) {
  input_mp(stock_id = "ple.27.7e_revision", OM = y, n_yrs = 20, 
           MP = "hr")$om@stock
})

### get reference points
refpts <- lapply(OMs_refset, function(y) {
  input_refpts(OM = y)
})

### plot
p <- plot_worm_distr_mult(stk = stk, stk_hist = stk_hist, refpts = refpts,
                          stk_labels = OMs_refset_label,
                          title = "")

ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM.png",
       plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM.pdf",
       plot = p, width = 16, height = 7.5, units = "cm")
    

### same but for 100-year projection
### get projection
stk <- readRDS(
  paste0("output/ple.27.7e_revision/refset/1000_100/multiplier/hr/",
         "mp_1_2_1_3.7_2_0.618657390447733_1.2_0.7-1_UK-FSP.rds"))@om@stock
stk <- lapply(seq_along(OMs_refset), function(x) {
  iter(stk, seq(from = (x - 1) * 1000 + 1, to = (x - 1) * 1000 + 1000))
})
### historical stock
stk_hist <- lapply(OMs_refset, function(y) {
  input_mp(stock_id = "ple.27.7e_revision", OM = y, n_yrs = 100, 
           MP = "hr")$om@stock
})
### get reference points
refpts <- lapply(OMs_refset, function(y) {
  input_refpts(stock_id = "ple.27.7e", OM = y)
})

### plot
p <- plot_worm_distr_mult(stk = stk, stk_hist = stk_hist, refpts = refpts,
                          stk_labels = OMs_refset_label,
                          title = "",
                          yr_end = 2124, xintercept = c(2024, 2044.5))

ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM_100.png",
       plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM_100.pdf",
       plot = p, width = 16, height = 7.5, units = "cm")

### Blim risk over time
stk_combined <- readRDS(
  paste0("output/ple.27.7e_revision/refset/1000_100/multiplier/hr/",
         "mp_1_2_1_3.7_2_0.618657390447733_1.2_0.7-1_UK-FSP.rds"))@om@stock
refpts_combined <- Reduce(refpts, f = FLCore::combine)
SSBs <- FLCore::window(ssb(stk_combined), start = 2025)
yrs <- dim(SSBs)[2]
its <- dim(SSBs)[6]
### collapse correction - not needed, all above threshold
Blim <- c(refpts_combined["Blim"])
Blim_ts <- SSBs %=% rep(c(Blim), each = dim(SSBs)[2])
risk <- iterMeans(SSBs/Blim_ts < 1)

p <- as.data.frame(risk) %>%
  ggplot(aes(x = year, y = data)) +
  geom_vline(xintercept = c(2034.5, 2044.5), colour = "grey") +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "1111") + 
  labs(x = "Year", y = expression(B[lim]~risk)) + 
  coord_cartesian(xlim = c(2025, 2124), ylim = c(0, 0.15), expand = FALSE) + 
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 1, size = 8))
ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM_100_risk.png",
       plot = p, width = 10, height = 5, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = "output/plots_revision/wormplots/refset_MP5_by_OM_100_risk.pdf",
       plot = p, width = 10, height = 5, units = "cm")
  
