### ------------------------------------------------------------------------ ###
### analyse MSE results ####
### ------------------------------------------------------------------------ ###
library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(stringr)
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
source("funs_WKNSMSE.R")

### ------------------------------------------------------------------------ ###
### define operating model (OM) names, reference and robustness set ####
### ------------------------------------------------------------------------ ###
### list of all operating models
OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", "R_h_lower",
         "R_failure", "overcatch", "undercatch", "Idx_higher")
OMs_label <- c("Reference set\n(combined)", 
               "Baseline", "Catch:\nno discards", "Catch:\n100% discards", 
               "Catch:\nno migration", 
               "M: -50%", "M: +50%", "M: Gislason", 
               "R: no AC", "R: +20%", "R: -20%", "R: h -20%",
               "R: failure", "Catch: +10%", "Catch: -10%", 
               "Uncertainty:\nindex +20%")

OMs_group <- c("refset (combined)", rep("refset", 7), rep("robset", 8))

OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")

### ------------------------------------------------------------------------ ###
### refset - x & w - grid search summary ####
### ------------------------------------------------------------------------ ###

### load all results and combine
df_runs <- foreach(index = c("Q1SWBeam", "UK-FSP"), .combine = bind_rows) %do% {
  #browser()
  path <- paste0("output/ple.27.7e/refset/1000_20/",
                 ifelse(identical(index, "Q1SWBeam"),
                        "x_w_n1_v_Q1SWBeam", "x_w_n1_v"),
                 "/hr/")
  tmp_runs <- readRDS(paste0(path, "runs.rds"))
  tmp_runs <- lapply(tmp_runs, function(x) {
    bind_cols(data.frame(t(x$pars)), data.frame(x$stats))
  })
  tmp_runs <- do.call(bind_rows, tmp_runs)
  tmp_runs$index <- index
  return(tmp_runs)
}
### add fitness value
df_runs$fitness <- df_runs$X11.20_Catch_rel -
  penalty(x = df_runs$X11.20_risk_Blim_max, 
          negative = FALSE, max = 1, 
          inflection = 0.06, 
          steepness = 1000)
### strict 5% risk limit
df_runs$fitness2 <- df_runs$X11.20_Catch_rel -
  ifelse(df_runs$X11.20_risk_Blim_max <= 0.05, 0, 1)
### groups
df_runs <- df_runs %>%
  mutate(group = paste(interval, index)) %>%
  mutate(group = factor(group, 
                        levels = c("1 UK-FSP", "1 Q1SWBeam", 
                                   "2 UK-FSP", "2 Q1SWBeam"),
                        labels = c("UK-FSP (annual)", "Q1SWBeam (annual)",
                                   "UK-FSP (biennial)", "Q1SWBeam (biennial)")))

### save results
saveRDS(df_runs, file = "output/paper/chr_refset_grid.rds")
write.csv(df_runs, file = "output/paper/chr_refset_grid.csv", row.names = FALSE)
# df_runs <- readRDS("output/paper/chr_refset_grid.rds")

### runs by chr version
table(df_runs[, c("index", "interval")])

### find optima
df_optima <- bind_rows(
  # df_runs %>%
  #   group_by(group) %>%
  #   filter(comp_b_multiplier <= 1.5) %>%
  #   #& !(index == "Q1SWBeam" & interval == 1)
  #   filter(X11.20_risk_Blim_max <= 0.05) %>%
  #   filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
  #   mutate(optimum = "local"),
  df_runs %>%
    group_by(group) %>%
    filter(X11.20_risk_Blim_max <= 0.05) %>%
    filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
    mutate(optimum = "global")
)

### number the MPs
df_optima$MP <- c(3, 4, 1, 2)
df_optima <- df_optima %>%
  arrange(MP) %>%
  relocate(MP)

saveRDS(df_optima, file = "output/paper/chr_refset_tuned.rds")
write.csv(df_optima, "output/paper/chr_refset_tuned.csv", row.names = FALSE)

### plot raw data
p_raw <- df_runs %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(y = multiplier, x = comp_b_multiplier, label = catch,
             fill = catch_col)) +
  geom_point(alpha = 0.8, shape = 21, stroke = NA, size = 0.4) +
  #geom_text(aes(label = round(X11.20_Catch_rel, 3)), size = 2.5) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_vline(data = df_optima,
             aes(xintercept = comp_b_multiplier, colour = "optimum"),
             linewidth = 0.2, linetype = "1111") +
  geom_hline(data = df_optima,
             aes(yintercept = multiplier, colour = "optimum"),
             linewidth = 0.2, linetype = "1111") +
  scale_colour_manual("Optimum",
                      values = c(optimum = "red")) +
  labs(y = "Harvest rate multiplier (x)", 
       x = expression(I[trigger]~multiplier~"(w)")) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +
  coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
  theme_bw(base_size = 8)
p_raw
# ggsave(filename = "output/plots/MP/refset_x_w_grid.png", plot = p_raw,
#        width = 16, height = 10, units = "cm", dpi = 600, type = "cairo",
#        bg = "white")
# ggsave(filename = "output/plots/MP/refset_x_w_grid.pdf", plot = p_raw,
#        width = 16, height = 10, units = "cm", 
#        bg = "white")


### interpolate missing cells by group (for plotting only)
df_runs_int <- foreach(group_i = levels(df_runs$group),
                       .combine = bind_rows) %do% {
  data_i <- df_runs %>%
    filter(group == group_i) %>%
    mutate(catch = X11.20_Catch_rel) %>%
    select(w = comp_b_multiplier, x = multiplier, 
           catch = X11.20_Catch_rel, risk = X11.20_risk_Blim_max,
           group)
  ### remove some cells to avoid NAs - no change to optimum
  # if (identical(group_i, "UK-FSP (annual)"))
  #   data_i <- data_i %>% filter(!(x %in% c(0.51, 0.52, 0.53, 0.54) &
  #                                   w %in% c(1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 
  #                                            1.95)))
  
  n_x <- length(seq(min(data_i$x), max(data_i$x), 0.01))
  n_w <- length(seq(min(data_i$w), max(data_i$w), 0.01))
  
  out_catch <- akima::interp(x = data_i$x, y = data_i$w,
                             z = data_i$catch, 
                             nx = n_x, ny = n_w,
                             linear = TRUE)
  out_risk <- akima::interp(x = data_i$x, y = data_i$w,
                            z = data_i$risk, 
                            nx = n_x, ny = n_w,
                            linear = TRUE)
  
  ### format
  df_out <- expand.grid(x = out_catch$x, w = out_catch$y)
  df_out <- data.frame(df_out)
  df_out$catch <- as.vector(out_catch$z)
  df_out$risk <- as.vector(out_risk$z)
  df_out$group <- group_i
  return(df_out)
}
df_runs_int$group <- factor(df_runs_int$group,
                            levels = c("UK-FSP (annual)", 
                                       "Q1SWBeam (annual)",
                                       "UK-FSP (biennial)", 
                                       "Q1SWBeam (biennial)"))

p_int <- df_runs_int %>%
  ### remove NAs for risk at very low x to avoid grey cells at border
  mutate(catch = ifelse(x == 0 & is.na(catch), 0, catch)) %>%
  mutate(risk = ifelse(x == 0 & is.na(risk), 0, risk)) %>%
  mutate(catch_col = ifelse(risk <= 0.05, catch, NA)) %>%
  ggplot(aes(y = x, x = w, fill = catch_col)) +
  geom_raster(alpha = 0.8, interpolate = FALSE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_vline(data = df_optima,
             aes(xintercept = comp_b_multiplier, colour = "optimum"),
             linewidth = 0.3, linetype = "1111") +
  geom_hline(data = df_optima,
             aes(yintercept = multiplier, colour = "optimum"),
             linewidth = 0.3, linetype = "1111") +
  scale_colour_manual("",
                      values = c(optimum = "red")) +
  labs(y = "Harvest rate multiplier (x)", 
       x = expression(I[trigger]~multiplier~"(w)")) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
  theme_bw(base_size = 8)
p_int
ggsave(filename = "output/paper/plots/refset_x_w_grid_int.png", plot = p_int,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/paper/plots/refset_x_w_grid_int.pdf", plot = p_int,
       width = 17, height = 10, units = "cm",
       bg = "white")

### summary table
df_x_w <- readRDS("output/paper/chr_refset_tuned.rds")
df_smry <- df_x_w %>%
  ungroup() %>%
  select(MP, index, 
         n1 = idxB_range_3, v = interval, x = multiplier,
         w = comp_b_multiplier,
         risk = X11.20_risk_Blim_max,
         catch = X11.20_Catch_rel,
         ssb = X11.20_SSB_rel,
         icv = X11.20_ICV) %>%
  arrange(MP)
df_smry
write.csv(df_smry, file = "output/paper/chr_refset_tuned_smry.csv", 
          row.names = FALSE)
saveRDS(df_smry, file = "output/paper/chr_refset_tuned_smry.rds")

### ------------------------------------------------------------------------ ###
### refset - visualise chr for all chr rule versions ####
### ------------------------------------------------------------------------ ###
df_smry <- readRDS("output/paper/chr_refset_tuned_smry.rds")

df_plot <- foreach(x = split(df_smry, df_smry$MP), 
                   .combine = bind_rows) %do% {
  #browser()
  ### points of chr rule
  x_values <- c(0, x$w, 5.9)
  y_values <- c(0, x$x, x$x)
  x <- bind_rows(x, x, x)
  x$x_axis <- x_values
  x$y_axis <- y_values
  return(x)
}

df_plot <- df_plot %>%
  mutate(interval = ifelse(v == 1, "annual", "biennial")) %>%
  mutate(interval = factor(interval, levels = c("annual", "biennial"))) %>%
  mutate(index = factor(index, levels = c("UK-FSP", "Q1SWBeam"))) %>%
  mutate(MP = paste0("CHR", MP)) %>%
  group_by(MP) %>%
  mutate(MP_label_y = max(y_axis))

p <- df_plot %>%
  ggplot(aes(x = x_axis, y = y_axis, linetype = MP)) +
  geom_line(show.legend = FALSE) +
  scale_linetype_manual(values = c("1111", "solid", "3131", "4242")) + 
  # geom_text(aes(label = MP,
  #               x = 4.75, y = MP_label_y + 0.045),
            # show.legend = FALSE, size = 2) +
  annotate(geom = "text", size = 2,
           label = paste0("CHR", df_smry$MP),
           x = 6.225, y = df_smry$x + c(0.02, -0.02, -0.01, 0.01)) + 
  #facet_grid(interval ~ index) +
  labs(y = "Harvest rate multiplier (x)", 
       x = expression(I[trigger]~multiplier~"(w)")) +
  coord_cartesian(xlim = c(0, 5.9), ylim = c(0, 0.99), expand = FALSE,
                  clip = "off") + 
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        plot.margin = unit(c(4, 30, 4, 4), units = "pt"))
p
ggsave(filename = "output/paper/plots/chr_illustration.png", plot = p,
       width = 8, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/paper/plots/chr_illustration.pdf", plot = p,
       width = 8, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### chr rule - violin plots ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x_w <- readRDS("output/paper/chr_refset_tuned.rds")
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds")) %>%
  arrange(MP)


### get stats
# , .combine = bind_rows
stats <- foreach(i = split(df_x_w, f = seq(nrow(df_x_w))), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows)  %:%
  foreach(period = c("long-term", "short-term", "all"),
          period_yrs = list(2035:2044, 2025:2034, 2025:2044),
          .combine = bind_rows) %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/",
                     ifelse(identical(i$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    mp_i <- readRDS(paste0(path_i, i$file))
    stk <- mp_i@om@stock

    ### get reference points
    refpts <- input_refpts(OM = OM)

    ### extract metrics
    yr_min <- min(period_yrs)
    yr_max <- max(period_yrs)
    stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
    stk <- window(stk, start = yr_min, end = yr_max)
    ssb_i <- c(ssb(stk)/refpts["Bmsy"])
    ssb_abs_i <- c(ssb(stk))
    catch_i <- c(catch(stk)/refpts["Cmsy"])
    catch_abs_i <- c(catch(stk))
    fbar_i <- c(fbar(stk)/refpts["Fmsy"])
    risk_i <- c(apply(ssb(stk) < rep(c(refpts["Blim"]), 
                                     each = dim(ssb(stk))[2]), 2, mean))
    icv_i <- c(iav(catch(stk_icv), period = i$interval))
    icv_annual_i <- c(iav(catch(stk_icv), period = 1))
    ### combine
    df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = ssb_abs_i, metric = "SSB_abs"),
                              data.frame(val = catch_i, metric = "catch"),
                              data.frame(val = catch_abs_i, metric = "catch_abs"),
                              data.frame(val = fbar_i, metric = "Fbar"),
                              data.frame(val = icv_i, metric = "ICV"),
                              data.frame(val = icv_annual_i, 
                                         metric = "ICV_annual"),
                              data.frame(val = risk_i, metric = "risk")
    ))
    df <- df %>%
      mutate(n1 = i$idxB_range_3,
             v = i$interval,
             x = i$multiplier,
             w = i$comp_b_multiplier,
             index = i$index,
             group = i$group,
             OM = OM, OM_group = OM_group,
             period = period,
             MP = i$MP) %>%
      relocate(MP)
    return(df)
}
saveRDS(stats, file = "output/paper/chr_refset_stats.rds")
# stats <- readRDS("output/paper/chr_refset_stats.rds")

### go through all solutions and plots stats
stats_plot <- stats %>%
  mutate(OM = factor(OM, levels = OMs,
                     labels = OMs_label),
         OM_group = factor(OM_group,
                           levels = OMs_group,
                           labels = c("", 
                                      rep("Reference set", 7), 
                                      rep("Robustness set", 8)))) %>%
  mutate(group = paste0("CHR", MP, " - ", index, " - ",
                        case_when(v == 1 ~ "annual",
                                  v == 2 ~ "biennial"),
                        " - ",
                        case_when(period == "long-term" ~ "long term",
                                  period == "short-term" ~ "short term",
                                  period == "all" ~ "all years")))


cols <- scales::hue_pal()(16)

. <- foreach(group_i = unique(stats_plot$group)) %:%
  foreach(scale = c("relative", "absolute")) %do% {
  #browser()
  stats_plot_i <- stats_plot %>%
    filter(group == group_i)
  MP_i <- stats_plot_i$MP[1]
  title_i <- stats_plot_i$group[1]
  file_i <- gsub(x = stats_plot_i$group[1], pattern = " - | ", 
                 replacement = "_")
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
             show.legend = FALSE, width = 0.8, colour = "black", linewidth = 0.2,
             position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(x = OM, y = val),
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    geom_hline(yintercept = 0.05, colour = "red", linewidth = 0.4, 
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
  ggsave(filename = paste0("output/paper/plots/MP/stats_",  file_i, 
                           "_refset", ".png"), 
         plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
         type = "cairo", bg = "white")
  ggsave(filename = paste0("output/paper/plots/MP/stats_",  file_i, 
                           "_refset", ".pdf"), 
         plot = p, width = 16, height = 10, units = "cm", bg = "white")
}

### ------------------------------------------------------------------------ ###
### ICES MSY - violin plots ####
### ------------------------------------------------------------------------ ###
df_altMPs <- data.frame(MP = c("ICES_SAM", "ICES_SAM"),
                        MP_label = c("MSY1", "MSY2"),
                        file = c("mp.rds", "mp_0.2_5400.rds"),
                        interval = c(1, 1))

### get stats
stats <- foreach(i = split(df_altMPs, f = seq(nrow(df_altMPs))), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows)  %:%
  foreach(period = c("long-term", "short-term", "all"),
          period_yrs = list(2035:2044, 2025:2034, 2025:2044),
          .combine = bind_rows) %do% {
    #browser()
    MP_i <- i$MP
    MP_label_i <- i$MP_label
    file_i <- i$file
    ### refset OM - combine manually
    if (identical(OM, "refset")) {
      stks <- lapply(OMs_refset, function(OM_i) {
        path_i <- paste0("output/ple.27.7e/", OM_i, "/1000_20/", MP_i, "/")
        mp_i <- readRDS(paste0(path_i, file_i))
        stk_i <- mp_i@om@stock
        return(stk_i)
      })
      stk <- Reduce(FLCore::combine, stks)
      
    } else {
      ### get projection
      path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", MP_i, "/")
      if (!file.exists(paste0(path_i, file_i))) return(NULL)
      mp_i <- readRDS(paste0(path_i, file_i))
      stk <- mp_i@om@stock
    }
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### extract metrics
    yr_min <- min(period_yrs)
    yr_max <- max(period_yrs)
    stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
    stk <- window(stk, start = yr_min, end = yr_max)
    ssb_i <- c(ssb(stk)/refpts["Bmsy"])
    ssb_abs_i <- c(ssb(stk))
    catch_i <- c(catch(stk)/refpts["Cmsy"])
    catch_abs_i <- c(catch(stk))
    fbar_i <- c(fbar(stk)/refpts["Fmsy"])
    risk_i <- c(apply(ssb(stk) < rep(c(refpts["Blim"]), 
                                     each = dim(ssb(stk))[2]), 2, mean))
    icv_i <- c(iav(catch(stk_icv), period = i$interval))
    ### combine
    df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = ssb_abs_i, metric = "SSB_abs"),
                              data.frame(val = catch_i, metric = "catch"),
                              data.frame(val = catch_abs_i, metric = "catch_abs"),
                              data.frame(val = fbar_i, metric = "Fbar"),
                              data.frame(val = icv_i, metric = "ICV"),
                              data.frame(val = risk_i, metric = "risk")
    ))
    df <- df %>%
      mutate(MP = MP_label_i, OM = OM, OM_group = OM_group, period = period)
    return(df)
}
saveRDS(stats, file = "output/paper/ICES_MSY_refset_stats.rds")
# stats <- readRDS("output/paper/ICES_MSY_refset_stats.rds")

### go through all solutions and plots stats
stats_plot <- stats %>%
  mutate(
   OM = factor(OM, levels = OMs,
                     labels = OMs_label),
   OM_group = factor(OM_group,
                     levels = OMs_group,
                     labels = c("", 
                                rep("Reference set", 7), 
                                rep("Robustness set", 8)))) %>%
  mutate(group = paste0(MP,
                        " - ",
                        case_when(period == "long-term" ~ "long term",
                                  period == "short-term" ~ "short term",
                                  period == "all" ~ "all years")))

cols <- scales::hue_pal()(16)

. <- foreach(group_i = unique(stats_plot$group)) %do% {
  #browser()
  stats_plot_i <- stats_plot %>%
    filter(group == group_i)
  title_i <- stats_plot_i$group[[1]]
  file_i <- gsub(x = stats_plot_i$group[1], pattern = " - | ", 
                 replacement = "_")
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
    geom_hline(yintercept = 0.05, colour = "red", linewidth = 0.4,
               linetype = "1111") +
   stat_summary(aes(x = OM, y = val),
                fun = "mean", geom = "point", shape = 4, size = 1,
                stroke = 0.25) +
   scale_fill_manual("", values = cols) +
   facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
   labs(y = expression(max.~B[lim]~risk), 
        title = title_i) +
   coord_cartesian(ylim = c(0, NA)) +
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
  ggsave(filename = paste0("output/paper/plots/MP/stats_", file_i,
                           "_refset.png"), 
        plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
        type = "cairo", bg = "white")
  ggsave(filename = paste0("output/paper/plots/MP/stats_", file_i,
                           "_refset.pdf"), 
        plot = p, width = 16, height = 10, units = "cm", bg = "white")
}

### summary table
stats_smry <- foreach(stats_i = split(stats_plot, f = stats_plot$group), 
                      .combine = bind_rows) %do% {
  #browser()
  stats_y <- foreach(metric_i = unique(stats_i$metric),
                     .combine = bind_rows) %do% {
    #browser()
    if (!identical(metric_i, "risk")) {
      res_y <- stats_i %>% 
        filter(metric == metric_i) %>%
        group_by(metric, MP, OM, period) %>%
        summarise(val = median(val))
    } else {
      res_y <- stats_i %>% 
        filter(metric == metric_i) %>%
        group_by(metric, MP, OM, period) %>%
        summarise(val = max(val))
    }
    return(res_y)
  }
}
stats_smry <- stats_smry %>%
  pivot_wider(names_from = metric, values_from = val)
stats_smry <- stats_smry %>%
  select(MP, period, OM, risk, catch, SSB, Fbar, ICV) %>%
  mutate(MP = ifelse(MP == "ICES_SAM", "ICES MSY", MP),
         period = case_when(period == "long-term" ~ "long term",
                            period == "short-term" ~ "short term",
                            period == "all" ~ "all years"),
         OM = gsub(x = OM, pattern = "\n", replacement = " "),
         OM = gsub(x = OM, pattern = "\\(combined\\)", replacement = ""),
         OM = trimws(OM)) 
saveRDS(stats_smry, file = "output/paper/ICES_MSY_stats_smry.rds")
write.csv(stats_smry, file = "output/paper/ICES_MSY_stats_smry.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### violin plots - CHR2 vs MSY2 - all OMs ####
### ------------------------------------------------------------------------ ###
### optimised chr rule - CHR2
stats_CHR2 <- readRDS("output/paper/chr_refset_stats.rds")
stats_CHR2 <- stats_CHR2 %>%
  filter(MP == 2 & period == "long-term") %>%
  mutate(OM = factor(OM, levels = OMs,
                     labels = OMs_label),
         OM_group = factor(OM_group,
                           levels = OMs_group,
                           labels = c("", 
                                      rep("Reference set", 7), 
                                      rep("Robustness set", 8))))

### ICES MSY rule with reduced Fmsy - MSY2
stats_MSY2 <- readRDS("output/paper/ICES_MSY_refset_stats.rds")
stats_MSY2 <- stats_MSY2 %>%
  filter(MP == "MSY2" & period == "long-term") %>%
  mutate(OM = factor(OM, levels = OMs,
                     labels = OMs_label),
         OM_group = factor(OM_group,
                           levels = OMs_group,
                           labels = c("", 
                                      rep("Reference set", 7), 
                                      rep("Robustness set", 8))))

### plot
cols <- scales::hue_pal()(16)

### risk
risk_max <- 0.315
p_CHR2_risk <- stats_CHR2 %>%
  filter(metric == "risk") %>%
  ggplot() +
  geom_col(data = . %>%
             group_by(OM, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM, y = val, fill = OM),
           show.legend = FALSE, width = 0.8, colour = "black", linewidth = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM, y = val),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  geom_hline(yintercept = 0.05, colour = "red", linewidth = 0.4, 
             linetype = "1111") +
  stat_summary(aes(x = OM, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk),
       title = "(a) Empirical: chr rule (CHR2)") +
  coord_cartesian(ylim = c(0, risk_max)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
#p_CHR2_risk
p_MSY2_risk <- stats_MSY2 %>%
  filter(metric == "risk") %>%
  ggplot() +
  geom_col(data = . %>%
             group_by(OM, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM, y = val, fill = OM),
           show.legend = FALSE, width = 0.8, colour = "black", linewidth = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM, y = val),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  geom_hline(yintercept = 0.05, colour = "red", linewidth = 0.4, 
             linetype = "1111") +
  stat_summary(aes(x = OM, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk),
       title = "(b) Model-based: ICES MSY rule (MSY2)") +
  coord_cartesian(ylim = c(0, risk_max)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        #strip.text.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
#p_MSY2_risk

### catch
p_CHR2_catch <- stats_CHR2 %>%
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
#p_CHR2_catch
p_MSY2_catch <- stats_MSY2 %>%
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
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_blank())
#p_MSY2_catch

### Catch - absolute
p_CHR2_catch_abs <- stats_CHR2 %>%
  filter(metric == "catch_abs") %>%
  ggplot(aes(x = OM, y = val/1000)) +
  geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = OM), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = OM, y = val/1000),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = "Catch\n(1000t)") +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
#p_CHR2_catch_abs
p_MSY2_catch_abs <- stats_MSY2 %>%
  filter(metric == "catch_abs") %>%
  ggplot(aes(x = OM, y = val/1000)) +
  geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = OM), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = OM, y = val/1000),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", linewidth = 0.4,
             linetype = "1111") +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_blank())
#p_MSY2_catch_abs

### SSB
p_CHR2_ssb <- stats_CHR2 %>%
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
#p_CHR2_ssb
p_MSY2_ssb <- stats_MSY2 %>%
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
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank())
#p_MSY2_ssb

### SSB - absolute scale
p_CHR2_ssb_abs <- stats_CHR2 %>%
  filter(metric == "SSB_abs") %>%
  ggplot(aes(x = OM, y = val/1000)) +
  geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = OM), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = OM, y = val/1000),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = "SSB\n(1000t)") +
  coord_cartesian(ylim = c(0, 24)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
#p_CHR2_ssb_abs
p_MSY2_ssb_abs <- stats_MSY2 %>%
  filter(metric == "SSB_abs") %>%
  ggplot(aes(x = OM, y = val/1000)) +
  geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = OM), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = OM, y = val/1000),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 24)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank())
#p_MSY2_ssb_abs

### ICV
p_CHR2_icv <- stats_CHR2 %>%
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
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())
#p_CHR2_icv
### ICV CHR2 biennial
p_CHR2_icv_biennial <- stats_CHR2 %>%
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
  labs(y = "ICV\n(biennial)") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
#p_CHR2_icv_biennial
### ICV CHR2 annual
p_CHR2_icv_annual <- stats_CHR2 %>%
  filter(metric == "ICV_annual") %>%
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
  labs(y = "ICV\n(annual)") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())
#p_CHR2_icv_annual
p_MSY2_icv <- stats_MSY2 %>%
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
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
#p_MSY2_icv

### combine plots
p <- p_CHR2_risk + p_MSY2_risk +
  p_CHR2_ssb + p_MSY2_ssb +
  p_CHR2_catch + p_MSY2_catch +
  p_CHR2_icv + p_MSY2_icv +
  plot_layout(ncol = 2)
p
ggsave(filename = "output/paper/plots/MP/stats_comp_CHR2_MSY2.png", 
       plot = p, width = 18, height = 10, units = "cm", dpi = 600, 
       type = "cairo", bg = "white")
ggsave(filename = "output/paper/plots/MP/stats_comp_CHR2_MSY2.pdf", 
       plot = p, width = 18, height = 10, units = "cm", bg = "white")

### with absolute catch and SSB
### and annual/biennial ICV
p <- p_CHR2_risk + p_MSY2_risk +
  p_CHR2_ssb + p_MSY2_ssb +
  p_CHR2_ssb_abs + p_MSY2_ssb_abs +
  p_CHR2_catch + p_MSY2_catch +
  p_CHR2_catch_abs + p_MSY2_catch_abs +
  p_CHR2_icv_biennial + plot_spacer() +
  p_CHR2_icv_annual + p_MSY2_icv +
  plot_layout(ncol = 2)
p
ggsave(filename = "output/paper/plots/MP/stats_comp_CHR2_MSY2_abs.png", 
       plot = p, width = 18, height = 16, units = "cm", dpi = 600, 
       type = "cairo", bg = "white")
ggsave(filename = "output/paper/plots/MP/stats_comp_CHR2_MSY2_abs.pdf", 
       plot = p, width = 18, height = 16, units = "cm", bg = "white")

### ------------------------------------------------------------------------ ###
### CHR2 - sensitivity to index uncertainty ####
### ------------------------------------------------------------------------ ###

### get runs
### runs summarised in MP_analysis.R
df_idx <- readRDS("output/ple.27.7e/refset/1000_20/sensitivity_idx/hr/runs.rds")

### plot
df_plot <- df_idx %>%
  filter(MP == 5) %>%
  select(MP, idx_unc, risk = `11:20_risk_Blim_max`,
         SSB = `11:20_SSB_rel`, catch = `11:20_Catch_rel`) %>%
  pivot_longer(-1:-2) %>%
  mutate(MP_label = factor(MP, levels = 1:10,
                           labels = paste0("MP", 1:10)),
         name = factor(name, levels = c("risk", "catch", "SSB"),
                       labels = c("B[lim]~risk", "Catch/MSY", 
                                  "SSB/B[MSY]"
                       )))
df_risk <- data.frame(value = 0.05,
                      name = "B[lim]~risk") %>%
  mutate(name = factor(name, levels = c("B[lim]~risk", "Catch/MSY", 
                                        "SSB/B[MSY]"
  )))

p <- df_plot %>%
  ggplot(aes(x = idx_unc, y = value)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "black") +
  geom_hline(data = df_risk,
             aes(yintercept = value),
             colour = "red", linewidth = 0.3) +
  geom_point(size = 0.1) + 
  geom_smooth(span = 0.4, linewidth = 0.3) +
  facet_grid(name ~ "CHR2", scales = "free_y", switch = "y", 
             labeller = label_parsed) + 
  labs(x = "Change to index uncertainty (%)") +
  scale_x_continuous(breaks = seq(0, 2, 0.25),
                     labels = sprintf(fmt = "%+3d", seq(-100, 100, 25))) +
  ylim(c(0, NA)) + 
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8))
p
ggsave(filename = "output/paper/plots/MP/idx_unc_CHR2.png",
       plot = p, width = 8, height = 6, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = "output/paper/plots/MP/idx_unc_CHR2.pdf", 
       plot = p, width = 8, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### compare trajectories - chr/ICES MSY ####
### ------------------------------------------------------------------------ ###

### get optimised solutions for chr rule
df_chr <- readRDS("output/paper/chr_refset_tuned.rds")
df_chr <- df_chr %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds"))


### ICES MSY rule 
stk_MSY2 <- lapply(paste0("output/ple.27.7e/", OMs_refset, 
                              "/1000_20/ICES_SAM/"),
                       function(y) {
                         readRDS(paste0(y, "mp_0.2_5400.rds"))@om@stock
                       })
stk_MSY2 <- Reduce(FLCore::combine, stk_MSY2)

### tuned chr rule CHR2
stk_CHR2 <- readRDS(paste0("output/ple.27.7e/refset/1000_20/multiplier/hr/",
                           df_chr$file[2]))@om@stock


MPs <- c("chr rule (CHR2)", "ICES MSY rule (MSY2)")
names(MPs) <- MPs
stk_list <- list(stk_CHR2, stk_MSY2)
input <- input_mp(OM = "refset", n_yrs = 20, MP = "hr")
stk_hist <- input$om@stock
stk_hist <- list(stk_hist, stk_hist)
refpts <- input_refpts(OM = "refset")
refpts[] <- NA

p <- plot_worm_comparison(stk = stk_list, stk_hist = stk_hist, 
                          names = names(MPs), refpts = refpts) +
  theme(legend.position = "bottom")
p
ggsave(filename = "output/paper/plots/wormplots/refset_comp_chr_MSY.png",
       width = 14, height = 8, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/paper/plots/wormplots/refset_comp_chr_MSY.pdf",
       width = 14, height = 8, units = "cm")

### catch and SSB
p <- plot_worm_comparison(stk = stk_list, stk_hist = stk_hist, 
                          names = names(MPs), refpts = refpts,
                          qnts_show = c("catch", "ssb"),
                          ncol = 1) +
  theme(legend.position = "bottom")
p
ggsave(filename = "output/paper/plots/wormplots/refset_comp_chr_MSY_catch_ssb.png",
       width = 8.5, height = 7, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/paper/plots/wormplots/refset_comp_chr_MSY_catch_ssb.pdf",
       width = 8.5, height = 7, units = "cm")

### ------------------------------------------------------------------------ ###
### compare trajectories - CHR1/2 ####
### ------------------------------------------------------------------------ ###

### get optimised solutions for chr rule
df_chr <- readRDS("output/paper/chr_refset_tuned.rds")
df_chr <- df_chr %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds"))


stk_list <- lapply(df_chr$file[1:2], function(x) {#browser()
  readRDS(paste0("output/ple.27.7e/refset/1000_20/multiplier/hr/",
                 x))@om@stock
})

MPs <- c("CHR1 (UK-FSP, annual)", "CHR2 (UK-FSP, biennial)")
names(MPs) <- MPs
input <- input_mp(OM = "refset", n_yrs = 20, MP = "hr")
stk_hist <- input$om@stock
stk_hist <- list(stk_hist, stk_hist)
refpts <- input_refpts(OM = "refset")
refpts[] <- NA

p <- plot_worm_comparison(stk = stk_list, stk_hist = stk_hist, 
                          names = names(MPs), refpts = refpts) +
  theme(legend.position = "bottom")
p

ggsave(filename = "output/paper/plots/wormplots/refset_comp_CHR1-2.png",
       width = 14, height = 8, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/paper/plots/wormplots/refset_comp_CHR1-2.pdf",
       width = 14, height = 8, units = "cm")


### plot SSB and catch, and Blim risk trajectories
stk_list <- window(FLStocks(stk_list), start = 1980)
stk_list[[1]][, ac(1980:2024)] <- input$om@stock[, ac(1980:2024)]
stk_list[[2]][, ac(1980:2024)] <- input$om@stock[, ac(1980:2024)]
refpts <- input_refpts(OM = "refset")
MPs <- c("CHR1 (UK-FSP, annual)", "CHR2 (UK-FSP, biennial)")
qnts <- lapply(seq_along(stk_list), function(x) {#browser()
  Blim_ts <- ssb(stk_list[[x]]) %=% rep(c(refpts["Blim"]), each = 65)
  qnts_i <- FLQuants(catch = catch(stk_list[[x]])/1000, 
                     ssb = ssb(stk_list[[x]])/1000,
                     risk = apply((ssb(stk_list[[x]])/Blim_ts) < 1, 2, 
                           mean, na.rm = TRUE))
  qnts_i <- lapply(qnts_i, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                   na.rm = TRUE)
  qnts_i_perc <- as.data.frame(FLQuants(qnts_i))
  qnts_i_perc$source <- MPs[x]
  return(qnts_i_perc)
})
qnts <- do.call(rbind, qnts)
qnts_perc <- qnts %>% select(year, iter, data, qname, source) %>%
  filter(!(iter != "50%" & qname == "risk")) %>%
  filter(!(qname == "risk" & year < 2024)) %>%
  #filter(data = ifelse(iter != "50%" & qname == "risk"), NA, data) %>%
  #mutate()
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "ssb", "risk"),
                        labels = c("'Catch (1000t)'", "'SSB (1000t)'",
                                   "B[lim]~risk")))
cols <- c(scales::pal_brewer(palette = "Dark2")(1),
          scales::pal_brewer(palette = "Set1")(4)[4])
p <- qnts_perc %>%
  ggplot(aes(x = year, y = `50%`, colour = source, fill = source, 
             linetype = source)) +
  geom_vline(xintercept = 2024, colour = "grey", size = 0.5) +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
              show.legend = FALSE, linewidth = 0) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE, linewidth = 0) +
  geom_line() +
  geom_hline(data = data.frame(y = 0.05, source = NA,
                               qname = factor("B[lim]~risk")),
             aes(yintercept = y), colour = "red", linewidth = 0.3) +
  scale_colour_manual("", values = cols) + 
  scale_fill_manual("", values = cols) + 
  scale_linetype_manual("", values = c("solid", "1111")) + 
  facet_wrap(~ qname, scales = "free_y", strip.position = "left", 
             labeller = label_parsed, ncol = 1) +
  labs(x = "Year") +
  coord_cartesian(ylim = c(0, NA), xlim = c(2010, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_blank())
p
ggsave(filename = "output/paper/plots/wormplots/refset_comp_CHR1-2_risk.png",
       width = 8, height = 8, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/paper/plots/wormplots/refset_comp_CHR1-2_risk.pdf",
       width = 8, height = 8, units = "cm")



### ------------------------------------------------------------------------ ###
### ICES MSY - full tuning ####
### ------------------------------------------------------------------------ ###
### tuned with full refset

input <- input_mp(OM = "refset", n_iter = 1000, MP = "ICES_SAM")
refpts <- input_refpts(OM = "refset")

file.copy("output/ICES_SAM_tuning_stats.rds",
          "output/paper/ICES_SAM_tuning_stats.rds", overwrite = TRUE)
stats <- readRDS("output/paper/ICES_SAM_tuning_stats.rds")

### find files
res_files <- list.files("output/ple.27.7e/refset/1000_20/ICES_SAM/",
                        pattern = "mp_")
res_files <- res_files[!res_files %in% unique(stats$file)]
res_files <- data.frame(file = res_files) %>%
  mutate(tmp = str_remove_all(file, "mp_|\\.rds")) %>%
  separate_wider_delim(tmp, delim = "_", names = c("Ftrgt", "Btrigger")) %>%
  mutate(Ftrgt = as.numeric(Ftrgt), Btrigger = as.numeric(Btrigger))



### go through files and get summary
stats_add <- foreach(file = res_files$file, Ftrgt = res_files$Ftrgt, 
                 Btrigger = res_files$Btrigger,
                 .combine = bind_rows, .errorhandling = "remove") %:%
  foreach(period = c("long-term", "short-term", "all"),
          period_yrs = list(2035:2044, 2025:2034, 2025:2044),
          .combine = bind_rows, .errorhandling = "remove") %do% {
  #browser()
  
  ### get projection
  mp_i <- readRDS(paste0("output/ple.27.7e/refset/1000_20/ICES_SAM/", file))
  stk <- mp_i@om@stock
  rm(mp_i)
  
  ### refpts
  Bmsy <- c(refpts["Bmsy"])
  Fmsy <- c(refpts["Fmsy"])
  Cmsy <- c(refpts["Cmsy"])
  Blim <- c(refpts["Blim"])
  
  ### extract metrics
  yr_min <- min(period_yrs)
  yr_max <- max(period_yrs)
  stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
  stk <- window(stk, start = yr_min, end = yr_max)
  
  SSBs <- ssb(stk)
  Fs <- fbar(stk)
  Cs <- catch(stk)
  Cs_long <- catch(stk_icv)
  
  ### account for OM/iteration-specific values
  Bmsy_ts <- SSBs %=% rep(c(Bmsy), each = dim(SSBs)[2])
  Fmsy_ts <- Fs %=% rep(c(Fmsy), each = dim(Fs)[2])
  Cmsy_ts <- Cs %=% rep(c(Cmsy), each = dim(Cs)[2])
  Blim_ts <- SSBs %=% rep(c(Blim), each = dim(SSBs)[2])
  
  data.frame(file = file,
             MP = "ICES_SAM",
             Ftrgt = Ftrgt, Btrigger = Btrigger,
             period = period,
             risk = max(apply((SSBs/Blim_ts) < 1, 2, mean, na.rm = TRUE), 
                        na.rm = TRUE),
             SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), 
                                                                na.rm = TRUE),
             Catch = median(c(Cs), na.rm = TRUE),
             Fbar = median(c(Fs), na.rm = TRUE),
             SSB_rel = median(c(SSBs/Bmsy_ts), na.rm = TRUE),
             Catch_rel = median(c(Cs/Cmsy_ts), na.rm = TRUE),
             Fbar_rel = median(c(Fs/Fmsy_ts), na.rm = TRUE),
             ICV = iav(Cs_long, period = 1, summary_all = median)
  )
}
stats <- unique(bind_rows(stats, stats_add))
stats <- stats %>% arrange(Ftrgt, Btrigger)
saveRDS(stats, file = "output/paper/ICES_SAM_tuning_stats.rds")
write.csv(stats, file = "output/paper/ICES_SAM_tuning_stats.csv", 
          row.names = FALSE)
View(stats %>% filter(period == "long-term"))

### duplicate Ftarget=0 values (Btrigger doesn't matter if no fishing)
# stats <- stats <- stats %>%
#   bind_rows(stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 1500),
#             stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 2500),
#             stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 3500),
#             stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 4500),
#             stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 5500),
#             stats %>%
#               filter(Ftrgt == 0 & Btrigger == 0) %>% mutate(Btrigger = 500)
#               )


### find optimum
df_optimum <- stats %>%
  group_by(period) %>%
  filter(risk <= 0.05) %>%
  filter(Catch_rel == max(Catch_rel))
df_refpts <- df_optimum %>%
  mutate(type = "Optimum") %>%
  bind_rows(data.frame(type = "ICES",
                       Ftrgt = 0.21106, Btrigger = 3265.99)) %>%
  mutate(type = factor(type, levels = c("Optimum", "ICES"),
                       labels = c("Tuned (MSY2)", "ICES (MSY1)"))) %>%
  select(type, Ftrgt, Btrigger)

### plot raw data
p_raw <- stats %>%
  filter(period == "long-term") %>%
  mutate(catch = Catch_rel,
         catch_col = ifelse(risk <= 0.05, 
                            Catch_rel, NA)) %>%
  ggplot(aes(x = Btrigger, y = Ftrgt, label = catch,
             fill = catch_col)) +
  geom_point(alpha = 0.8, shape = 21, stroke = NA, size = 2) +
  geom_tile(alpha = 0.8) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_refpts,
             aes(yintercept = Ftrgt, colour = type, linetype = type),
             linewidth = 0.2) +
  geom_vline(data = df_refpts,
             aes(xintercept = Btrigger, colour = type, linetype = type),
             linewidth = 0.2) +
  scale_colour_manual("", values = c("Tuned (MSY2)" = "red", 
                                     "ICES (MSY1)" = "black")) +
  scale_linetype_manual("", values = c("Tuned (MSY2)" = "1111", 
                                       "ICES (MSY1)" = "solid")) +
  labs(x = expression(B[trigger]), y = expression(F[target])) +
  coord_cartesian(#expand = TRUE, 
    xlim = c(0, NA), ylim = c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p_raw
ggsave(filename = "output/paper/plots/MSY_grid_raw.png", plot = p_raw,
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/paper/plots/MSY_grid_raw.pdf", plot = p_raw,
       width = 8.5, height = 5, units = "cm",
       bg = "white")


### interpolation to find area where to focus on

### manual linear interpolation to get surface
### format into wide data.frame
x <- seq(0, 7000, 100)
y <- seq(0, 0.3, 0.01)
stats_catch <- stats %>%
  filter(period == "long-term") %>%
  select(Btrigger, Ftrgt, Catch_rel) %>%
  full_join(expand.grid(Btrigger = x,
                        Ftrgt = y)) %>%
  arrange(Btrigger, Ftrgt) %>%
  pivot_wider(names_from = Btrigger, values_from = Catch_rel) %>%
  #arrange(rev(Ftrgt)) %>%
  tibble::column_to_rownames("Ftrgt") %>%
  as.matrix()
stats_risk <- stats %>%
  filter(period == "long-term") %>%
  select(Btrigger, Ftrgt, risk) %>%
  full_join(expand.grid(Btrigger = x,
                        Ftrgt = y)) %>%
  arrange(Btrigger, Ftrgt) %>%
  pivot_wider(names_from = Btrigger, values_from = risk) %>%
  #arrange(rev(Ftrgt)) %>%
  tibble::column_to_rownames("Ftrgt") %>%
  as.matrix()
### 1st - outside edges
stats_catch[y == 0] <- approx(x = x, y = stats_catch[y == 0], xout = x, na.rm = TRUE)$y
stats_catch[y == 0.3] <- approx(x = x, y = stats_catch[y == 0.3], xout = x, na.rm = TRUE)$y
stats_catch[, x == 0] <- approx(x = y, y = stats_catch[, x == 0], xout = y, na.rm = TRUE)$y
stats_catch[, x == 7000] <- approx(x = y, y = stats_catch[, x == 7000], xout = y, na.rm = TRUE)$y
stats_risk[y == 0] <- approx(x = x, y = stats_risk[y == 0], xout = x, na.rm = TRUE)$y
stats_risk[y == 0.3] <- approx(x = x, y = stats_risk[y == 0.3], xout = x, na.rm = TRUE)$y
stats_risk[, x == 0] <- approx(x = y, y = stats_risk[, x == 0], xout = y, na.rm = TRUE)$y
stats_risk[, x == 7000] <- approx(x = y, y = stats_risk[, x == 7000], xout = y, na.rm = TRUE)$y
### inner horizontal lines
for (i in y) {
  if (isTRUE(i %in% c(0, 0.3))) next()
  if (isTRUE(sum(!is.na(stats_catch[y == i])) > 2)) {
    stats_catch[y == i] <- approx(x = x, y = stats_catch[y == i], xout = x, na.rm = TRUE)$y
    stats_risk[y == i] <- approx(x = x, y = stats_risk[y == i], xout = x, na.rm = TRUE)$y
  }
}
### vertical lines
for (i in x) {
  if (isTRUE(i %in% c(0, 7000))) next()
  if (isTRUE(sum(!is.na(stats_catch[, x == i])) > 2)) {
    stats_catch[, x == i] <- approx(x = y, y = stats_catch[, x == i], xout = y, na.rm = TRUE)$y
    stats_risk[, x == i] <- approx(x = y, y = stats_risk[, x == i], xout = y, na.rm = TRUE)$y
  }
}
### try plotting
image(t(stats_catch))
image(t(stats_risk))
### convert back into data.frame
stats_int <- full_join(as.data.frame(stats_catch) %>% 
            tibble::rownames_to_column("Ftrgt") %>%
            pivot_longer(-Ftrgt, names_to = "Btrigger", values_to = "Catch_rel"),
          as.data.frame(stats_risk) %>% 
            tibble::rownames_to_column("Ftrgt") %>%
            pivot_longer(-Ftrgt, names_to = "Btrigger", values_to = "risk")) %>%
  mutate(Ftrgt = as.numeric(Ftrgt),
         Btrigger = as.numeric(Btrigger))
saveRDS(stats_int, file = "output/paper/grid_int_cells.rds")
### plot
p_int <- stats_int %>%
  mutate(catch = Catch_rel,
         catch_col = ifelse(risk <= 0.05, 
                            Catch_rel, NA)) %>%
  ggplot(aes(x = Btrigger, y = Ftrgt, label = catch,
             fill = catch_col)) +
  geom_raster(alpha = 0.8, interpolate = FALSE) +
  # geom_raster(#data = . %>% filter(risk <= 0.05), 
  #             alpha = 0.8, interpolate = TRUE) +
  # geom_tile(data = . %>% filter(risk > 0.05), alpha = 0.8) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
  geom_hline(data = df_refpts,
             aes(yintercept = Ftrgt, colour = type, linetype = type),
             linewidth = 0.2) +
  geom_vline(data = df_refpts,
             aes(xintercept = Btrigger, colour = type, linetype = type),
             linewidth = 0.2) +
  scale_colour_manual("", values = c("Tuned (MSY2)" = "red", 
                                     "ICES (MSY1)" = "black")) +
  scale_linetype_manual("", values = c("Tuned (MSY2)" = "1111", 
                                       "ICES (MSY1)" = "solid")) +
  labs(x = expression(B[trigger]^MP), y = expression(F[target]^MP)) +
  coord_cartesian(#expand = TRUE, 
    xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p_int
ggsave(filename = "output/paper/plots/MSY_grid_int.png", plot = p_int,
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/paper/plots/MSY_grid_int.pdf", plot = p_int,
       width = 8.5, height = 5, units = "cm",
       bg = "white")


### interpolate with akima - doesn't work well...
# df_int <- stats %>%
#   filter(period == "long-term") %>%
#   dplyr::select(Btrigger, Ftrgt, risk, Catch_rel)
# x_val <- seq(min(df_int$Btrigger), max(df_int$Btrigger), 100)
# y_val <- seq(min(df_int$Ftrgt), max(df_int$Ftrgt), 0.01)
# n_x <- length(x_val)
# n_y <- length(y_val)
# 
# out_catch <- akima::interp(x = df_int$Btrigger/1000, y = df_int$Ftrgt,
#                            z = df_int$Catch_rel, 
#                            xo = x_val/1000, yo = y_val,
#                            #nx = n_x, ny = n_w,
#                            linear = TRUE, extrap = TRUE)
# out_risk <- akima::interp(x = df_int$Btrigger/1000, y = df_int$Ftrgt,
#                           z = df_int$risk, 
#                           xo = x_val/1000, yo = y_val,,
#                           #nx = n_x, ny = n_w,
#                           linear = TRUE, extrap = TRUE)
# 
# ### format
# df_int <- expand.grid(Btrigger = out_catch$x * 1000, Ftrgt = out_catch$y)
# df_int <- data.frame(df_int)
# df_int$catch <- as.vector((out_catch$z))
# df_int$risk <- as.vector((out_risk$z))
# 
# df_int %>%
#   mutate(catch = catch,
#          catch_col = ifelse(risk <= 0.05, 
#                             catch, NA)) %>%
#   ggplot(aes(x = Btrigger, y = Ftrgt, label = catch,
#              fill = catch_col)) +
#   geom_point(alpha = 0.8, shape = 21, stroke = NA, size = 2) +
#   geom_raster(alpha = 0.8) +
#   scale_fill_gradientn(paste0("Catch/MSY"),
#                        colours = hcl.colors(10),
#                        values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
#                                   1), 
#                        breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   labs(x = expression(B[trigger]), y = expression(F[trgt])) +
#   coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
#   theme_bw(base_size = 8)
# df_int %>%
#   mutate(catch = catch,
#          catch_col = catch) %>%
#   ggplot(aes(x = Btrigger, y = Ftrgt, label = catch,
#              fill = catch_col)) +
#   geom_point(alpha = 0.8, shape = 21, stroke = NA, size = 2) +
#   geom_raster(alpha = 0.8) +
#   scale_fill_gradientn(paste0("Catch/MSY"),
#                        colours = hcl.colors(10),
#                        values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
#                                   1), 
#                        breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   labs(x = expression(B[trigger]), y = expression(F[trgt])) +
#   coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
#   theme_bw(base_size = 8)
# df_int %>%
#   filter(risk <= 0.05) %>%
#   filter(catch == max(catch, na.rm = TRUE))

### ------------------------------------------------------------------------ ###
### wormplots - chr and ICES MSY rule ####
### ------------------------------------------------------------------------ ###
smry_chr <- readRDS("output/paper/chr_refset_tuned.rds")

### chr by OM
. <- foreach(x = split(smry_chr, seq(nrow(smry_chr)))) %:%
  foreach(OM = OMs[-1], OM_label = OMs_label[-1])  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", 
                     ifelse(identical(x$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    file_i <- paste("mp", x$idxB_lag, x$idxB_range_3, x$exp_b,
                    x$comp_b_multiplier, x$interval, x$multiplier,
                    x$upper_constraint, x$lower_constraint, sep = "_", 
                    collapse = "")
    mp_i <- readRDS(paste0(path_i, file_i, ".rds"))
    stk <- mp_i@om@stock
    
    ### historical stock
    input <- input_mp(OM = OM, n_yrs = 20, MP = "hr")
    stk_hist <- input$om@stock
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### plot
    p <- plot_worm_distr(stk = stk, stk_hist = stk_hist, refpts = refpts,
                         title = paste0("CHR", x$MP, " - ", x$group, " - ",
                                        OM_label))
    ggsave(filename = paste0("output/paper/plots/wormplots/all/CHR", x$MP,
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/paper/plots/wormplots/all/CHR", x$MP,
                             "_", OM, ".pdf"),
           plot = p, width = 16, height = 7.5, units = "cm")
}
### ICES MSY by OM
. <- foreach(MP = c("MSY1", "MSY2"), file = c("mp.rds", "mp_0.2_5400.rds")) %:%
  foreach(OM = OMs[-1], OM_label = OMs_label[-1])  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/ICES_SAM/")
    file_i <- file
    mp_i <- readRDS(paste0(path_i, file_i))
    stk <- mp_i@om@stock
    
    ### historical stock
    input <- input_mp(OM = OM, n_yrs = 20, MP = "ICES_SAM")
    stk_hist <- input$om@stock
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### plot
    p <- plot_worm_distr(stk = stk, stk_hist = stk_hist, refpts = refpts,
                         title = paste0(MP, " - ", OM_label))
    ggsave(filename = paste0("output/paper/plots/wormplots/all/", MP,
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/paper/plots/wormplots/all/", MP,
                             "_", OM, ".pdf"),
           plot = p, width = 16, height = 7.5, units = "cm")
}

### plot refset
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_refset_label <- c("Baseline", "Catch:\nno discards", 
                      "Catch:\n100% discards", 
                      "Catch:\nno migration", 
                      "M: -50%", "M: +50%", "M: Gislason")
### chr refset
. <- foreach(x = split(smry_chr, seq(nrow(smry_chr)))) %:%
  foreach(OM = "refset", OM_label = "Reference set (combined)")  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OMs_refset, "/1000_20/", 
                     ifelse(identical(x$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    file_i <- paste("mp", x$idxB_lag, x$idxB_range_3, x$exp_b,
                    x$comp_b_multiplier, x$interval, x$multiplier,
                    x$upper_constraint, x$lower_constraint, sep = "_", 
                    collapse = "")
    stk <- lapply(path_i, function(y) {
      readRDS(paste0(y, file_i, ".rds"))@om@stock
    })
    
    ### historical stock
    stk_hist <- lapply(OMs_refset, function(y) {
      input_mp(OM = y, n_yrs = 20, MP = "hr")$om@stock
    })
    
    ### get reference points
    refpts <- lapply(OMs_refset, function(y) {
      input_refpts(OM = y)
    })
    
    ### plot
    p <- plot_worm_distr_mult(stk = stk, stk_hist = stk_hist, refpts = refpts,
                              stk_labels = OMs_refset_label,
                              title = paste0("CHR", x$MP, " - ", x$group, " - ",
                                             OM_label))

    ggsave(filename = paste0("output/paper/plots/wormplots/all/CHR", x$MP,
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/paper/plots/wormplots/all/CHR", x$MP,
                             "_", OM, ".pdf"), 
           plot = p, width = 16, height = 7.5, units = "cm")
}
### ICES MSY refset
. <- foreach(MP = c("MSY1", "MSY2"), file = c("mp.rds", "mp_0.2_5400.rds")) %:%
  foreach(OM = "refset", OM_label = "Reference set (combined)")  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OMs_refset, "/1000_20/ICES_SAM/")
    file_i <- file
    stk <- lapply(path_i, function(y) {
      readRDS(paste0(y, file_i))@om@stock
    })
    
    ### historical stock
    stk_hist <- lapply(OMs_refset, function(y) {
      input_mp(OM = y, n_yrs = 20, MP = "hr")$om@stock
    })
    
    ### get reference points
    refpts <- lapply(OMs_refset, function(y) {
      input_refpts(OM = y)
    })
    
    ### plot
    p <- plot_worm_distr_mult(stk = stk, stk_hist = stk_hist, refpts = refpts,
                              stk_labels = OMs_refset_label,
                              title = paste0(MP, " - ", OM_label))
    ggsave(filename = paste0("output/paper/plots/wormplots/all/", MP,
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/paper/plots/wormplots/all/", MP,
                             "_", OM, ".pdf"), 
           plot = p, width = 16, height = 7.5, units = "cm")
}

