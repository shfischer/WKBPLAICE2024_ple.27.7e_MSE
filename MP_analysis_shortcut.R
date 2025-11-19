### ------------------------------------------------------------------------ ###
### analyse shortcut MSE results ####
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
  labs(x = "Target F") +
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
Ftrgt_SAM <- sort(unique(c(seq(0, 0.5, 0.1), seq(0.05, 0.25, 0.1), 0.19, 0.198, 0.199, "eqsim")))

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
  labs(x = "Target F") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        legend.key.height = unit(0.6, "lines"))

ggsave(filename = "output/plots/shortcut/SAM_shortcut_stats_comparison.png", 
       width = 16, height = 8, units = "cm", dpi = 600)

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