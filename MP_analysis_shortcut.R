### ------------------------------------------------------------------------ ###
### analyse shortcut MSE results ####
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
### collate shortcut runs ####
### ------------------------------------------------------------------------ ###

### get stats
Ftrgt_shortcut <- sort(unique(c(seq(0, 0.5, 0.025), seq(0.225, 0.23, 0.001))))
stats_shortcut <- foreach(MP = c("ICES_SAM_shortcut"), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_label = OMs_label, OM_group = OMs_group, 
          .combine = bind_rows)  %:%
  foreach(Ftrgt = Ftrgt_shortcut) %do% {
    #browser()
    print(paste0("OM=", OM, " Ftrgt=", Ftrgt)); flush.console()

    ### refset OM - combine manually
    if (identical(OM, "refset")) {
      stks <- lapply(OMs_refset, function(OM_i) {
        path_i <- paste0("output/ple.27.7e/", OM_i, "/1000_20/", MP, "/")
        mp_i <- readRDS(paste0(path_i, "mp_", Ftrgt, ".rds"))
        stk_i <- mp_i@om@stock
        return(stk_i)
      })
      stk <- Reduce(FLCore::combine, stks)
      
    } else {
      ### get projection
      path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", MP, "/")
      mp_i <- readRDS(paste0(path_i, "mp_", Ftrgt, ".rds"))
      stk <- mp_i@om@stock
    }
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### extract metrics
    yr_min <- 2035
    yr_max <- 2044
    stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
    stk_i <- window(stk, start = yr_min, end = yr_max)
    
    data.frame(
      ssb = median(ssb(stk_i)),
      catch = median(catch(stk_i)),
      icv = median(iav(catch(stk_icv), period = 1)),
      risk = max(iterMeans(ssb(stk_i) < rep(c(refpts["Blim"]), 
                                            each = dim(ssb(stk_i))[2]))),
      MP = MP, OM_label = OM_label, OM = OM, OM_group = OM_group,
      Ftrgt = Ftrgt)
}
saveRDS(stats_shortcut, file = "output/shortcut/shortcut_stats.rds")

stats_shortcut %>%
  filter(OM == "refset") %>%
  select(catch, risk, Ftrgt) %>%
  pivot_longer(c(catch, risk)) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  theme_bw(base_size = 8)


stats_shortcut %>%
  filter(OM_group %in% c("refset (combined)", "refset")) %>%
  mutate(OM_group = factor(OM_group, 
                              levels = c("refset (combined)", "refset"),
                              labels = c("Reference set", "Individual OM"))) %>%
  select(OM, OM_group, catch, risk, Ftrgt) %>%
  pivot_longer(c(catch, risk)) %>%
  mutate(name = factor(name, levels = c("risk", "catch"),
                       labels = c("B[lim]~risk", "Catch~(t)"))) %>%
  ggplot(aes(x = Ftrgt, y = value, linewidth = OM_group, linetype = OM_group,
             group = OM)) +
  geom_line() +
  scale_linewidth_manual("Operating model",
                         values = c("Reference set" = 0.5, 
                                    "Individual OM" = 0.1)) +
  scale_linetype_manual("Operating model",
                        values = c("Reference set" = "solid", 
                                   "Individual OM" = "1111")) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             labeller = label_parsed) +
  labs(x = expression(F[target])) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        legend.key.height = unit(0.6, "lines"))


### ------------------------------------------------------------------------ ###
### collate SAM runs ####
### ------------------------------------------------------------------------ ###


### get stats
Ftrgt_SAM <- sort(unique(c(seq(0, 0.5, 0.05), seq(0.025, 0.2, 0.025),
                           0.19, 0.198, 0.199, "eqsim")))

stats_SAM <- foreach(MP = c("ICES_SAM"), 
                          .combine = bind_rows) %:%
  foreach(OM = OMs[1:8], OM_label = OMs_label[1:8], OM_group = OMs_group[1:8], 
          .combine = bind_rows)  %:%
  foreach(Ftrgt = Ftrgt_SAM) %do% {
    #browser()
    print(paste0("OM=", OM, " Ftrgt=", Ftrgt)); flush.console()
    
    ### refset OM - combine manually
    if (identical(OM, "refset")) {
      stks <- lapply(OMs_refset, function(OM_i) {
        path_i <- paste0("output/ple.27.7e/", OM_i, "/1000_20/", MP, "/")
        file_i <- ifelse(identical(Ftrgt, "eqsim"),
                         "mp.rds", paste0("mp_", Ftrgt, ".rds"))
        mp_i <- readRDS(paste0(path_i, file_i))
        stk_i <- mp_i@om@stock
        return(stk_i)
      })
      stk <- Reduce(FLCore::combine, stks)
      
    } else {
      ### get projection
      path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", MP, "/")
      file_i <- ifelse(identical(Ftrgt, "eqsim"),
                       "mp.rds", paste0("mp_", Ftrgt, ".rds"))
      mp_i <- readRDS(paste0(path_i, file_i))
      stk <- mp_i@om@stock
    }
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### extract metrics
    yr_min <- 2035
    yr_max <- 2044
    stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
    stk_i <- window(stk, start = yr_min, end = yr_max)
    
    cat("\n"); flush.console()
    data.frame(
      ssb = median(ssb(stk_i)),
      catch = median(catch(stk_i)),
      icv = median(iav(catch(stk_icv), period = 1)),
      risk = max(iterMeans(ssb(stk_i) < rep(c(refpts["Blim"]), 
                                            each = dim(ssb(stk_i))[2]))),
      MP = MP, OM_label = OM_label, OM = OM, OM_group = OM_group,
      Ftrgt = Ftrgt)
}
saveRDS(stats_SAM, file = "output/shortcut/SAM_stats.rds")

stats_SAM %>%
  filter(OM == "refset") %>%
  mutate(Ftrgt = as.numeric(ifelse(Ftrgt == "eqsim", 0.211, Ftrgt))) %>%
  select(catch, risk, Ftrgt) %>%
  pivot_longer(c(catch, risk)) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8)

### ------------------------------------------------------------------------ ###
### plot shortcut and full SAM together ####
### ------------------------------------------------------------------------ ###
stats_shortcut <- readRDS("output/shortcut/shortcut_stats.rds")
stats_SAM <- readRDS("output/shortcut/SAM_stats.rds")

df_plot <- bind_rows(
  stats_shortcut %>% mutate(source = "Shortcut"),
  stats_SAM %>% mutate(source = "Full") %>%
    mutate(Ftrgt = as.numeric(ifelse(Ftrgt == "eqsim", 0.211, Ftrgt)))) %>%
  filter(OM_group %in% c("refset (combined)", "refset")) %>%
  mutate(OM_group = factor(OM_group, 
                           levels = c("refset (combined)", "refset"),
                           labels = c("Reference set", "Individual OM"))) %>%
  mutate(source = factor(source, levels = c("Full", "Shortcut"))) %>%
  select(OM, OM_group, source, catch, risk, Ftrgt) %>%
  pivot_longer(c(catch, risk)) %>%
  mutate(name = factor(name, levels = c("risk", "catch"),
                       labels = c("B[lim]~risk", "'Catch (t)'")))

df_plot %>%
  ggplot(aes(x = Ftrgt, y = value, linewidth = OM_group, linetype = OM_group,
             colour = source, group = interaction(OM, source)
             )) +
  geom_line() +
  geom_hline(data = data.frame(y = 0.05,
                               name = factor("B[lim]~risk", 
                                             levels = c("B[lim]~risk", "Catch~(t)"))),
             aes(yintercept = y, alpha = "5% risk limit"), 
             colour = "red", linetype = "2121", linewidth = 0.3) +
  scale_linewidth_manual("Operating model",
                         values = c("Reference set" = 0.5,
                                    "Individual OM" = 0.1)) +
  scale_linetype_manual("Operating model",
                        values = c("Reference set" = "solid",
                                   "Individual OM" = "1111")) +
  scale_colour_manual("MSE type",
                         values = c("Full" = "black",
                                    "Shortcut" = "blue")) +
  scale_alpha_manual("", values = c(1)) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             labeller = label_parsed) +
  labs(x = expression(F[target])) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        legend.key.height = unit(0.6, "lines"))


p1 <- df_plot %>%
  filter(name == "B[lim]~risk") %>%
  ggplot(aes(x = Ftrgt, y = value, linewidth = OM_group, linetype = OM_group,
             colour = source, group = interaction(OM, source)
  )) +
  annotate(geom = "rect", xmin = 0.199, xmax = 1, ymin = -1, ymax = 1e+4,
           fill = "red", alpha = 0.1) +
  geom_line() +
  geom_hline(data = data.frame(y = 0.05,
                               name = factor("B[lim]~risk", 
                                             levels = c("B[lim]~risk", "Catch~(t)"))),
             aes(yintercept = y, alpha = "5% risk limit"), 
             colour = "red", linetype = "2121", linewidth = 0.3) +
  geom_vline(xintercept = 0.229, colour = "blue", linetype = "solid", 
             linewidth = 0.2, alpha = 0.3) +
  geom_vline(xintercept = 0.199, colour = "black", linetype = "solid", 
             linewidth = 0.2, alpha = 0.3) +
  scale_linewidth_manual("Operating model",
                         values = c("Reference set" = 0.5,
                                    "Individual OM" = 0.1)) +
  scale_linetype_manual("Operating model",
                        values = c("Reference set" = "solid",
                                   "Individual OM" = "1212")) +
  scale_colour_manual("MSE type",
                      values = c("Full" = "black",
                                 "Shortcut" = "blue")) +
  scale_alpha_manual("", values = c(1)) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             labeller = label_parsed) +
  labs(x = expression(F[target])) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 1), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        legend.key.height = unit(0.6, "lines"))
p2 <- df_plot %>%
  filter(name == "'Catch (t)'") %>%
  ggplot(aes(x = Ftrgt, y = value, linewidth = OM_group, linetype = OM_group,
             colour = source, group = interaction(OM, source)
  )) +
  annotate(geom = "rect", xmin = 0.199, xmax = 1, ymin = -1, ymax = 1e+4,
           fill = "red", alpha = 0.1) +
  geom_vline(xintercept = 0.229, colour = "blue", linetype = "solid", 
             linewidth = 0.2, alpha = 0.3) +
  geom_vline(xintercept = 0.199, colour = "black", linetype = "solid", 
             linewidth = 0.2, alpha = 0.3) +
  geom_line() +
  scale_linewidth_manual("Operating model",
                         values = c("Reference set" = 0.5,
                                    "Individual OM" = 0.1)) +
  scale_linetype_manual("Operating model",
                        values = c("Reference set" = "solid",
                                   "Individual OM" = "1212")) +
  scale_colour_manual("MSE type",
                      values = c("Full" = "black",
                                 "Shortcut" = "blue")) +
  scale_alpha_manual("", values = c(1)) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             labeller = label_parsed) +
  labs(x = expression(F[target])) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 1650), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 1500, 500)) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        legend.key.height = unit(0.6, "lines"))
p1 + p2 + plot_layout(guides = "collect", ncol = 2)

ggsave(filename = "output/plots/shortcut/SAM_shortcut_stats_comparison.png", 
       width = 16, height = 6, units = "cm", dpi = 600)
ggsave(filename = "output/plots/shortcut/SAM_shortcut_stats_comparison.pdf", 
       width = 16, height = 6, units = "cm")

### vertical plot
p1 + theme(axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           axis.ticks.x = element_blank()) +
  p2 + plot_layout(guides = "collect", ncol = 1)
ggsave(filename = "output/plots/shortcut/SAM_shortcut_stats_comparison_long.png", 
       width = 10, height = 12, units = "cm", dpi = 600)

### compare OMs with high risk at F=0.5
stats_shortcut %>%
  filter(Ftrgt == 0.5 & OM_group %in% c("refset", "refset (combined)")) %>%
  select(risk, OM) %>%
  arrange(-risk)
#        risk            OM
# 1 0.8100000         M_low
# 2 0.6760000    M_Gislason
# 3 0.5210000 Catch_no_surv
# 4 0.4375714        refset
# 5 0.3350000      baseline
# 6 0.2900000        M_high
# 7 0.2560000     migr_none
# 8 0.1720000 Catch_no_disc

stats_SAM %>%
  filter(Ftrgt == 0.5) %>%
  select(risk, OM) %>%
  arrange(-risk)
#        risk            OM
# 1 0.9490000         M_low
# 2 0.4790000 Catch_no_surv
# 3 0.3597143        refset
# 4 0.3550000    M_Gislason
# 5 0.3420000      baseline
# 6 0.2850000     migr_none
# 7 0.2100000 Catch_no_disc
# 8 0.0410000        M_high








### ------------------------------------------------------------------------ ###
### ICES MSY shortcut - full tuning ####
### ------------------------------------------------------------------------ ###
### tuned with full refset
### CPU hours per cell: 00:48:09 ... 00:56:23

refpts <- input_refpts(OM = "refset")

stats <- readRDS("output/shortcut/ICES_SAM_shortcut_tuning_stats.rds")

### find files
res_files <- list.files("output/ple.27.7e/refset/1000_20/ICES_SAM_shortcut/",
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
  mp_i <- readRDS(paste0("output/ple.27.7e/refset/1000_20/ICES_SAM_shortcut/",
                         file))
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
             MP = "ICES_SAM_shortcut",
             Ftrgt = Ftrgt, Btrigger = Btrigger,
             period = period,
             risk = max(apply((SSBs/Blim_ts) < 1, 2, mean, na.rm = TRUE), 
                        na.rm = TRUE),
             SSB = median(c(SSBs), na.rm = TRUE),
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
saveRDS(stats, file = "output/shortcut/ICES_SAM_shortcut_tuning_stats.rds")
write.csv(stats, file = "output/shortcut/ICES_SAM_shortcut_tuning_stats.csv", 
          row.names = FALSE)
View(stats %>% filter(period == "long-term"))

### find optimum
df_optimum <- stats %>%
  group_by(period) %>%
  filter(risk <= 0.05) %>%
  filter(Catch_rel == max(Catch_rel))
df_refpts <- df_optimum %>%
  mutate(type = "Optimum") %>%
  # bind_rows(data.frame(type = "ICES",
  #                      Ftrgt = 0.21106, Btrigger = 3265.99)) %>%
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
ggsave(filename = "output/shortcut/plots/MSY_grid_raw.png", plot = p_raw,
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/shortcut/plots/MSY_grid_raw.pdf", plot = p_raw,
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
saveRDS(stats_int, file = "output/shortcut/grid_int_cells.rds")
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
                       limits = c(0, 1.012)) +
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
ggsave(filename = "output/shortcut/plots/MSY_grid_int.png", plot = p_int,
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/shortcut/plots/MSY_grid_int.pdf", plot = p_int,
       width = 8.5, height = 5, units = "cm",
       bg = "white")


### ------------------------------------------------------------------------ ###
### compare full and shortcut grid ####
### ------------------------------------------------------------------------ ###
stats_int_shortcut <- readRDS("output/shortcut/grid_int_cells.rds")
stats_int_full <- readRDS("output/paper/grid_int_cells.rds")

stats_int_comp <- stats_int_shortcut %>% 
  mutate(type = "Shortcut") %>%
  bind_rows(stats_int_full %>%
              mutate(type = "Full")) %>%
  mutate(type = factor(type, levels = c("Full", "Shortcut"),
                       labels = c("Full MSE", "Shortcut MSE")))

stats_opt <- stats_int_comp %>%
  group_by(type) %>%
  filter(risk <= 0.05) %>%
  filter(Catch_rel == max(Catch_rel))



p_grids <- stats_int_comp %>%
  mutate(catch = Catch_rel,
         catch_col = ifelse(risk <= 0.05, 
                            Catch_rel, NA)) %>%
  ggplot() +
  geom_raster(aes(x = Btrigger, y = Ftrgt, fill = catch_col),
              alpha = 0.8, interpolate = FALSE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1.012)) +
  geom_point(data = stats_opt %>% rename(type2 = type),
             aes(x = Btrigger, y = Ftrgt, colour = type2), 
             shape = 3, stroke = 0.7, size = 5) +
  # geom_hline(data = stats_opt,
  #            aes(yintercept = Ftrgt, colour = type, linetype = type),
  #            linewidth = 0.2) +
  # geom_vline(data = stats_opt,
  #            aes(xintercept = Btrigger, colour = type, linetype = type),
  #            linewidth = 0.2) +
  scale_colour_manual("Tuned", values = c("Full MSE" = "black",
                                     "Shortcut MSE" = "blue")) +
  # scale_linetype_manual("Tuned", values = c("Full MSE" = "solid",
  #                                      "Shortcut MSE" = "1111")) +
  facet_wrap(~ type) +
  labs(x = expression(B[trigger]), y = expression(F[target])) +
  coord_cartesian(#expand = TRUE, 
    xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p_grids

p_grid_full <- stats_int_comp %>%
  filter(type == "Full MSE") %>%
  mutate(catch = Catch_rel,
         catch_col = ifelse(risk <= 0.05, 
                            Catch_rel, NA)) %>%
  ggplot() +
  geom_raster(aes(x = Btrigger, y = Ftrgt, fill = catch_col),
              alpha = 0.8, interpolate = FALSE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1.012)) +
  geom_point(data = stats_opt %>% rename(type2 = type),
             aes(x = Btrigger, y = Ftrgt, colour = type2), 
             shape = 3, stroke = 0.7, size = 5) +
  scale_colour_manual("Tuned", values = c("Full MSE" = "black",
                                          "Shortcut MSE" = "blue")) +
  facet_grid("F[target]" ~ "'Full MSE'", switch = "y", 
             labeller = "label_parsed") +
  labs(x = expression(B[trigger])) +
  coord_cartesian(#expand = TRUE, 
    xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"),
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        axis.title.y = element_blank())

p_grid_shortcut <- stats_int_comp %>%
  filter(type == "Shortcut MSE") %>%
  mutate(catch = Catch_rel,
         catch_col = ifelse(risk <= 0.05, 
                            Catch_rel, NA)) %>%
  ggplot() +
  geom_raster(aes(x = Btrigger, y = Ftrgt, fill = catch_col),
              alpha = 0.8, interpolate = FALSE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1.012)) +
  geom_point(data = stats_opt %>% rename(type2 = type),
             aes(x = Btrigger, y = Ftrgt, colour = type2), 
             shape = 3, stroke = 0.7, size = 5) +
  scale_colour_manual("Tuned", values = c("Full MSE" = "black",
                                          "Shortcut MSE" = "blue")) +
  facet_grid("F[target]" ~ "'Shortcut MSE'", switch = "y", 
             labeller = "label_parsed") +
  labs(x = expression(B[trigger]), y = expression(F[target])) +
  coord_cartesian(#expand = TRUE, 
    xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank())


p_grid_full + p_grid_shortcut

### ------------------------------------------------------------------------ ###
### stats of optimised solutions by OM ####
### ------------------------------------------------------------------------ ###
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
OMs_label2 <- c("Reference set",
                paste0("OM", 1:15))
OMs_group <- c("refset (combined)", rep("refset", 7), rep("robset", 8))




### get stats
stats_comp <- foreach(MSE = c("full", "shortcut"), .combine = bind_rows) %:%
  foreach(optimum = c("full", "shortcut"), .combine = bind_rows) %:% 
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows) %do% {

  #browser()
  file <- paste0("output/ple.27.7e/", OM, "/1000_20/",
                 case_when(MSE == "full" ~ "ICES_SAM",
                           MSE == "shortcut" ~ "ICES_SAM_shortcut"),
                 "/mp_",
                 case_when(optimum == "full" ~ "0.2_5400",
                           optimum == "shortcut" ~ "0.23_5100"),
                 ".rds")
  #if (!isTRUE(file.exists(file))) return(NULL)
  mp <- readRDS(file)
  stk <- mp@om@stock
  
  ### get reference points
  refpts <- input_refpts(OM = OM)
  
  ### extract metrics
  yr_min <- 2035
  yr_max <- 2044
  stk_icv <- window(stk, start = yr_min - 1, end = yr_max)
  stk <- window(stk, start = yr_min, end = yr_max)
  
  catch_i <- median(catch(stk)/refpts["Cmsy"])
  risk_i <- max(iterMeans(ssb(stk) < rep(c(refpts["Blim"]), 
                                         each = dim(ssb(stk))[2])))
  
  # catch_refset <- median(catch_i)
  # catch_OMs <- sapply(1:7, function(x) {
  #   median(iter(catch_i, seq((x - 1)*1000, x*1000)))
  # })
  # 
  # risk_refset <- max(iterMeans(risk_i))
  # risk_OMs <- sapply(1:7, function(x) {
  #   max(iterMeans(iter(risk_i, seq((x - 1)*1000, x*1000))))
  # })
  
  data.frame(MSE = MSE, optimum = optimum,
             OM = OM, OM_group = OM_group,
             catch = catch_i,
             risk = risk_i)

}

df_plot <- stats_comp %>%
  mutate(MSE = factor(MSE, levels = c("full", "shortcut"),
                      labels = c("Full MSE", "Shortcut MSE")),
         optimum = factor(optimum, levels = c("full", "shortcut"),
                          labels = c("Full MSE", "Shortcut MSE")),
         OM = factor(OM, levels = OMs),
         ### OM numbers (ID) from WKBPLAICE paper
         OM_label = factor(OM, levels = OMs,
                           labels = c("Reference\nset",
                                      paste0("OM", c(1:3, 6, 7:9, 10, 12:13,
                                                     14, 11, 4:5, 15)))),
         OM_group = factor(OM_group, 
                           levels = c("refset (combined)", "refset", "robset"),
                           labels = c("' '", "'Reference set'", 
                                      "'Robustness set'"))
         ) %>%
  pivot_longer(catch:risk) %>%
  mutate(name = factor(name, levels = c("catch", "risk"),
                       labels = c("Catch/MSY", "B[lim]~risk")))

df_plot_max <- df_plot %>%
  group_by(MSE, optimum, OM, OM_group, OM_label, name) %>%
  summarise(value = max(value))

df_plot_ref <- df_plot %>% 
  filter(name == "B[lim]~risk") %>% 
  slice_head(n = 1) %>%
  select(MSE, name, value) %>%
  mutate(value = 0.05)

p_stats_full <- df_plot %>%
  filter(optimum == "Full MSE") %>%
  ggplot(aes(x = OM_label, y = value, colour = MSE, fill = MSE)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_blank(data = df_plot %>% filter(optimum == "Shortcut MSE")) +
  geom_hline(data = df_plot_ref,
             aes(yintercept = value), 
             colour = "red", linetype = "2222", alpha = 0.5, 
             show.legend = FALSE) +
  scale_colour_manual("MSE type", values = c("Full MSE" = "black",
                                          "Shortcut MSE" = "blue")) +
  scale_fill_manual("MSE type", values = c("Full MSE" = "black",
                                             "Shortcut MSE" = "blue")) +
  facet_grid(name ~ OM_group, scales = "free", switch = "y",
             labeller = "label_parsed", space = "free_x") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.position = "none",
        strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing.x = unit(0, "lines"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_stats_shortcut <- df_plot %>%
  filter(optimum == "Shortcut MSE") %>%
  ggplot(aes(x = OM_label, y = value, colour = MSE, fill = MSE)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_blank(data = df_plot %>% filter(optimum == "Full MSE")) +
  geom_hline(data = df_plot_ref,
             aes(yintercept = value), 
             colour = "red", linetype = "2222", alpha = 0.5, 
             show.legend = FALSE) +
  scale_colour_manual("MSE type", values = c("Full MSE" = "black",
                                             "Shortcut MSE" = "blue")) +
  scale_fill_manual("MSE type", values = c("Full MSE" = "black",
                                           "Shortcut MSE" = "blue")) +
  facet_grid(name ~ OM_group, scales = "free", switch = "y",
             labeller = "label_parsed", space = "free_x") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing.x = unit(0, "lines"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank())

p_stats_full + p_stats_shortcut

p_full <- p_grid_full + theme(legend.position = "none") + 
  p_grid_shortcut +
  p_stats_full + theme(legend.position = "none") +
  p_stats_shortcut +
  plot_layout(axes = "keep", heights = c(1, 0.8))
p_full
ggsave(filename = "output/shortcut/plots/ple_comparison.png", plot = p_full,
       width = 18, height = 12, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/shortcut/plots/ple_comparison.pdf", plot = p_full,
       width = 18, height = 12, units = "cm",
       bg = "white")

