### ------------------------------------------------------------------------ ###
### Operating model figures ####
### ------------------------------------------------------------------------ ###

library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(tibble)
library(cowplot)
library(patchwork)
library(ggplot2)
library(foreach)
library(stockassessment)
library(FLfse)
source("funs.R")
source("funs_GA.R")
source("funs_analysis.R")
source("funs_OM.R")

### ------------------------------------------------------------------------ ###
### plot OM trajectories vs. SAM assessment - baseline OM ####
### ------------------------------------------------------------------------ ###
### MSY reference points
refpts <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
refpts <- iterMedians(refpts)
### OM stk
stk <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/stk.rds"))
### SAM model fit
fit <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/SAM_fit.rds"))

### OM metrics
qnts <- FLQuants(catch = catch(stk)/1000, rec = rec(stk)/1e+03,
                 ssb = ssb(stk)/1000, fbar = fbar(stk))
qnts <- window(qnts, end = 2023)
### percentiles
qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    na.rm = TRUE)
qnts_perc <- FLQuants(qnts_perc)
qnts_perc <- as.data.frame(qnts_perc)
qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "OM")

### get assessment summary
df_SAM <- bind_rows(list(
  as.data.frame(catchtable(fit)/1000) %>%
    mutate(qname = "catch") %>% 
    rownames_to_column(var = "year"),
  as.data.frame(rectable(fit)/1e+03) %>%
    mutate(qname = "rec") %>% 
    rownames_to_column(var = "year"),
  as.data.frame(ssbtable(fit)/1000) %>%
    mutate(qname = "ssb") %>% 
    rownames_to_column(var = "year"),
  as.data.frame(fbartable(fit)) %>%
    mutate(qname = "fbar") %>% 
    rownames_to_column(var = "year")
)) %>%
  select(year = year, qname = qname, `2.5%` = Low, `50%` = Estimate, 
         `97.5%` = High) %>%
  mutate(source = "SAM",
         year = as.numeric(year))

### combine and format
df <- bind_rows(qnts_perc, df_SAM) %>%
  filter(year <= 2024) %>%
  mutate(source = factor(source, levels = c("OM", "SAM"),
                       labels = c("Operating model", "SAM assessment"))) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "rec", "fbar", "ssb"),
                        labels = c("Catch (1000t)", "Recruitment (1000s)",
                                   "F (ages 3-6)", "SSB (1000t)")))
### MSY levels
df_refs <- data.frame(
  qname = factor(c("catch", "rec", "fbar", "ssb"), 
                 levels = c("catch", "rec", "fbar", "ssb"),
                 labels = c("Catch (1000t)", "Recruitment (1000s)",
                            "F (ages 3-6)", "SSB (1000t)")),
  value = c(c(refpts["Cmsy"])/1000, 6538.489/1000, 
            c(refpts["Fmsy"]), c(refpts["Bmsy"])/1000),
  source = "MSY")


### plot
p <- ggplot() +
  geom_ribbon(data = df %>%
                filter(source == "Operating model"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df %>%
                filter(source == "Operating model"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df %>% filter(source == "SAM assessment"),
            aes(x = year, y = `2.5%`),
            colour = "red", linetype = "2121", linewidth = 0.2, alpha = 0.3) + 
  geom_line(data = df %>% filter(source == "SAM assessment"),
            aes(x = year, y = `97.5%`),
            colour = "red", linetype = "2121", linewidth = 0.2, alpha = 0.3)+ 
  geom_line(data = df,
            aes(x = year, y = `50%`, colour = source, linetype = source)) +
  geom_hline(data = df_refs,
             aes(yintercept = value, colour = source, linetype = source)) +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual("", values = c("Operating model" = "black", 
                                     "SAM assessment" = "red",
                                     "MSY" = "darkgrey")) +
  scale_linetype_manual("", values = c("Operating model" = "solid", 
                                       "SAM assessment" = "2121",
                                       "MSY" = "1111")) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),  
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.1), 
        legend.background = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_vs_SAM.png", plot = p, 
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_vs_SAM.pdf", plot = p, 
       width = 16, height = 8, units = "cm")

### plot again but with some iterations
df_iters <- as.data.frame(qnts) %>% 
  select(year, iter, data, qname) %>%
  filter(year <= 2024) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "rec", "fbar", "ssb"),
                        labels = c("Catch (1000t)", "Recruitment (1000s)",
                                   "F (ages 3-6)", "SSB (1000t)")))
p <- ggplot() +
  geom_ribbon(data = df %>%
                filter(source == "Operating model"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df %>%
                filter(source == "Operating model"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df %>% filter(source == "Operating model"),
            aes(x = year, y = `50%`, colour = "Median"),
            colour = "black") +
  geom_line(data = df_iters %>% filter(iter %in% 1:5),
            aes(x = year, y = data, colour = iter),
            linewidth = 0.1,
            show.legend = FALSE) +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  scale_y_continuous(limits = c(0, NA)) +
  #scale_colour_manual("", values = c("black")) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_vs_SAM_iters.png", plot = p, 
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_vs_SAM_iters.pdf", plot = p, 
       width = 16, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### plot OM trajectories vs. SAM assessment - reference set OMs ####
### ------------------------------------------------------------------------ ###

### reference set OMs
OMs <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason")
OMs_label <- c("Baseline", "Catch: no discards", "Catch: 100% discards", 
               "Catch: no migration", 
               "M: -50%", "M: +50%", "M: Gislason")

### SAM fits
fit_list <- readRDS("input/ple.27.7e/preparation/SAM_all_fits.rds")

res <- foreach(OM_no = seq_along(OMs), OM = OMs, OM_label = OMs_label, 
               .combine = bind_rows) %do% {
  #browser()
  ### OM stk
  stk <- readRDS(paste0("input/ple.27.7e/", OM, "/1000_100/stk.rds"))
  ### SAM model fit
  fit_i <- fit_list[[OM_no]]
  
  ### OM metrics
  qnts <- FLQuants(catch = catch(stk)/1000, rec = rec(stk)/1e+03,
                   ssb = ssb(stk)/1000, fbar = fbar(stk))
  qnts <- window(qnts, end = 2023)
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(source = "OM")
  
  df_SAM <- bind_rows(list(
    as.data.frame(catchtable(fit_i)/1000) %>%
      mutate(qname = "catch") %>% 
      rownames_to_column(var = "year"),
    as.data.frame(rectable(fit_i)/1e+03) %>%
      mutate(qname = "rec") %>% 
      rownames_to_column(var = "year"),
    as.data.frame(ssbtable(fit_i)/1000) %>%
      mutate(qname = "ssb") %>% 
      rownames_to_column(var = "year"),
    as.data.frame(fbartable(fit_i)) %>%
      mutate(qname = "fbar") %>% 
      rownames_to_column(var = "year")
  )) %>%
    select(year = year, qname = qname, `2.5%` = Low, `50%` = Estimate, 
           `97.5%` = High) %>%
    mutate(source = "SAM",
           year = as.numeric(year))
  
  ### combine and format
  df <- bind_rows(qnts_perc, df_SAM) %>%
    filter(year <= 2024) %>%
    mutate(source = factor(source, levels = c("OM", "SAM"),
                           labels = c("Operating model", "SAM assessment"))) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c("Catch (1000t)", "Recruitment (1000s)",
                                     "F (ages 3-6)", "SSB (1000t)")))
  
  df$OM <- OM
  df$OM_label <- OM_label
  
  return(df)
  
}

p_catch <- res %>%
  filter(qname == "Catch (1000t)" & source == "Operating model") %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = res %>%
              filter(qname == "Catch (1000t)" & source == "SAM assessment"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "SAM assessment" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "SAM assessment" = "2121")) +
  labs(y = "Catch (1000t)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text.y = element_blank())
p_R <- res %>%
  filter(qname == "Recruitment (1000s)" & source == "Operating model") %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = res %>%
              filter(qname == "Recruitment (1000s)" & source == "SAM assessment"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "SAM assessment" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "SAM assessment" = "2121")) +
  labs(y = "Recruitment (1000)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position.inside = c(0.5, 0.985),
        legend.position = "inside",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ssb <- res %>%
  filter(qname == "SSB (1000t)" & source == "Operating model") %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = res %>%
              filter(qname == "SSB (1000t)" & source == "SAM assessment"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "SAM assessment" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "SAM assessment" = "2121")) +
  labs(y = "SSB (1000t)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text.y = element_blank())
p_fbar <- res %>%
  filter(qname == "F (ages 3-6)" & source == "Operating model") %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = res %>%
              filter(qname == "F (ages 3-6)" & source == "SAM assessment"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "SAM assessment" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "SAM assessment" = "2121")) +
  labs(y = "Mean F (ages 3-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")

p <- p_catch + p_R + p_ssb + p_fbar + plot_layout(nrow = 1)
p

ggsave(filename = "output/plots/OM/OM_vs_SAM_all_fit.png", plot = p, 
       width = 16, height = 18, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_vs_SAM_all_fit.pdf", plot = p, 
       width = 16, height = 18, units = "cm")


### ------------------------------------------------------------------------ ###
### biological data - baseline OM ####
### ------------------------------------------------------------------------ ###
### catch weights at age
### stock weights at age
### natural mortality
### maturity

### load stk
stk <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/stk.rds"))
stk_median <- iterMedians(stk)

df_biol <- FLQuants(cw = catch.wt(stk_median),
                    lw = landings.wt(stk_median),
                    dw = discards.wt(stk_median),
                    sw = stock.wt(stk_median),
                    m = m(stk_median),
                    mat = mat(stk_median))
df_biol <- as.data.frame(df_biol) %>% 
  filter(year <= 2023)
df_biol <- df_biol %>%
  mutate(qname = factor(qname,
                        levels = c("cw", "lw", "dw", "sw", "m", "mat"),
                        labels = c("Catch weights (kg)",
                                   "Landings weights (kg)",
                                   "Discard weights (kg)",
                                   "Stock weights (kg)",
                                   "Natural mortality",
                                   "Maturity"))) %>%
  mutate(age = factor(age, levels = 10:2,
                      labels = c("10+", 9:2)))
cols <- scales::hue_pal()(length(2:10))
cols <- cols[c(seq(from = 1, to = length(cols), by = 3),
               seq(from = 2, to = length(cols), by = 3),
               seq(from = 3, to = length(cols), by = 3))]

p <- df_biol %>%
  mutate(data = ifelse(data == 0, NA, data)) %>%
  ggplot(aes(x = year, y = data, colour = age)) +
  annotate("rect", xmin = 2018.5, xmax = 2023.5, ymin = -Inf, ymax = Inf,
           alpha = 0.05, fill = "red") + 
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.3) + 
  geom_vline(xintercept = 2018.5, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = 2023.5, colour = "red", linewidth = 0.3) +
  scale_colour_manual("Age (years)", values = cols) + 
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") + 
  coord_cartesian(xlim = c(2000, 2023), ylim = c(0, NA)) +
  labs(x = "Year") + 
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_biol.png", plot = p, 
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_biol.pdf", plot = p, 
       width = 16, height = 8, units = "cm")

### fisheries selectivity
sel <- harvest(stk_median) ### median
sel <- propagate(sel, (dim(stk)[6] + 1))
sel[,,,,, -1] <- harvest(stk) ### iterations

df_sel <- as.data.frame(sel)
df_sel <- df_sel %>%
  group_by(year, iter) %>%
  mutate(data = data/max(data))  %>%
  mutate(source = ifelse(iter == 1, "Median", "Simulation\nreplicates"))

p1 <- df_sel %>%
  filter(year %in% 2014:2018) %>%
  ggplot(aes(x = age, y = data, group = iter, 
             colour = source, linewidth = source, alpha = source)) +
  geom_line(data = . %>% filter(source == "Simulation\nreplicates")) +
  geom_line(data = . %>% filter(source == "Median")) +
  scale_colour_manual("", values = c("Median" = "red",
                                     "Simulation\nreplicates" = "black")) +
  scale_alpha_manual("", values = c("Median" = 1,
                                    "Simulation\nreplicates" = 0.05)) +
  scale_linewidth_manual("", values = c("Median" = 0.5,
                                        "Simulation\nreplicates" = 0.1)) +
  facet_grid("Historical\n" ~ year) +
  ylim(c(0, NA)) +
  labs(x = "Age (years)", y = "Fisheries selectivity") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank())
p2 <- df_sel %>%
  filter(year %in% 2019:2023) %>%
  ggplot(aes(x = age, y = data, group = iter, 
             colour = source, linewidth = source, alpha = source)) +
  geom_line(data = . %>% filter(source == "Simulation\nreplicates")) +
  geom_line(data = . %>% filter(source == "Median")) +
  scale_colour_manual("", values = c("Median" = "red",
                                     "Simulation\nreplicates" = "black")) +
  scale_alpha_manual("", values = c("Median" = 1,
                                    "Simulation\nreplicates" = 0.05)) +
  scale_linewidth_manual("", values = c("Median" = 0.5,
                                        "Simulation\nreplicates" = 0.1)) +
  facet_grid("Historical\n(sampled from)" ~ year) +
  ylim(c(0, NA)) +
  labs(x = "Age (years)", y = "Fisheries selectivity") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p3 <- df_sel %>%
  filter(year %in% 2024:2028) %>%
  ggplot(aes(x = age, y = data, group = iter, 
             colour = source, linewidth = source, alpha = source)) +
  geom_line(data = . %>% filter(source == "Simulation\nreplicates")) +
  geom_line(data = . %>% filter(source == "Median")) +
  scale_colour_manual("", values = c("Median" = "red",
                                     "Simulation\nreplicates" = "black")) +
  scale_alpha_manual("", values = c("Median" = 1,
                                    "Simulation\nreplicates" = 0.05)) +
  scale_linewidth_manual("", values = c("Median" = 0.5,
                                        "Simulation\nreplicates" = 0.1)) +
  facet_grid("Projection\n(re-sampled)" ~ year) +
  ylim(c(0, NA)) +
  labs(x = "Age (years)", y = "Fisheries selectivity") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.position = "none",
        axis.title.y = element_blank())
p <- p1 / p2 / p3
p
ggsave(filename = "output/plots/OM/OM_biol_sel.png", plot = p, 
       width = 16, height = 9.5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_biol_sel.pdf", plot = p, 
       width = 16, height = 9.5, units = "cm")

### discard rates by age
### use OM with 100% discards
stkD <- readRDS(paste0("input/ple.27.7e/Catch_no_surv/1000_100/stk.rds"))
stkD_median <- iterMedians(stkD)

disc <- discards.n(stkD_median)/catch.n(stkD_median)
disc <- propagate(disc, (dim(stkD)[6] + 1))
disc[,,,,, -1] <- discards.n(stkD)/catch.n(stkD) ### iterations
disc <- window(disc, start = 2019, end = 2024)

df_disc <- as.data.frame(disc)
df_disc <- df_disc %>%
  group_by(year, iter) %>%
  mutate(source = ifelse(iter == 1, "Median", "Simulation\nreplicates")) %>%
  mutate(year = as.character(year)) %>%
  mutate(year = ifelse(year <= 2023, year, "Projection"))

p <- df_disc %>%
  filter(source == "Median") %>%
  ggplot(aes(x = age, y = data, group = iter, 
             colour = source, linewidth = source, alpha = source)) +
  geom_line() +
  geom_line(data = . %>% filter(source == "Median")) +
  scale_colour_manual("", values = c("Median" = "red")) +
  scale_alpha_manual("", values = c("Median" = 1)) +
  scale_linewidth_manual("", values = c("Median" = 0.5)) +
  facet_wrap(~ year, nrow = 1) +
  ylim(c(0, NA)) +
  labs(x = "Age (years)", y = "Discard rate") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"))
p
ggsave(filename = "output/plots/OM/OM_biol_disc_rate.png", plot = p, 
       width = 16, height = 4, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_biol_disc_rate.pdf", plot = p, 
       width = 16, height = 4, units = "cm")



### ------------------------------------------------------------------------ ###
### plot OM trajectories vs. ICES assessment - alternative OMs ####
### ------------------------------------------------------------------------ ###

OM_list <- c("baseline",
             "Catch_no_disc", "Catch_no_surv", "migr_none",
             "M_low", "M_high", "M_Gislason",
             "R_no_AC", "R_higher", "R_lower")
OM_labels <- c("Baseline",
               "Catch: no discards", "Catch: 100% discards", "Catch: no migration",
               "M: -50%", "M: +50%", "M: Gislason",
               "R: no AC", "R: +20%", "R: -20%")

### get values from operating models
df_OM <- foreach(OM = OM_list,
                 OM_label = OM_labels,
                 .combine = bind_rows) %do% {
  #browser()
  stk_i <- readRDS(paste0("input/ple.27.7e/", OM, "/1000_100/stk.rds"))
  ### get metrics
  qnts <- FLQuants(ssb = ssb(stk_i)/1000, fbar = fbar(stk_i),
                   catch = catch(stk_i)/1000, rec = rec(stk_i)/1000)
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                       na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(OM = OM, OM_label = OM_label,
           source = "Operating model")
  return(qnts_perc)
}
df_OM <- df_OM %>% filter(year <= 2023)

df_OM_baseline <- df_OM %>% 
  filter(OM == "baseline") %>%
  select(-OM, -OM_label) %>%
  mutate(source = "Baseline")

### load reference points
df_refpts <- foreach(OM = OM_list,
                     OM_label = OM_labels, .combine = bind_rows) %do% {#browser()
  file <- paste0("input/ple.27.7e/", OM, "/1000_100/refpts_mse.rds")
  if (isFALSE(file.exists(file))) return(NULL)
  refpts_i <- readRDS(file)
  refpts_i <- iterMedians(refpts_i)
  data.frame(OM = OM, OM_label = OM_label,
             Fmsy = c(refpts_i["Fmsy"]),
             Bmsy = c(refpts_i["Bmsy"]),
             Cmsy = c(refpts_i["Cmsy"]),
             Blim = c(refpts_i["Blim"]))
}

df_OM <- df_OM %>%
  filter(OM %in% OM_list[1:7])
df_refpts <- df_refpts  %>%
  filter(OM %in% OM_list[1:7])

p_catch <- df_OM %>%
  filter(qname == "catch") %>%
  ggplot(aes(x = year)) +
  geom_hline(data = df_refpts, 
             aes(yintercept = Cmsy/1000), 
             colour = "black", linewidth = 0.3, linetype = "2222",
             show.legend = FALSE) + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = df_OM_baseline %>% filter(qname == "catch"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "Baseline" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "Baseline" = "2121")) +
  labs(y = "Catch (1000t)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text.y = element_blank())
p_R <- df_OM %>%
  filter(qname == "rec") %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = df_OM_baseline %>% filter(qname == "rec"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "Baseline" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "Baseline" = "2121")) +
  labs(y = "Recruitment (1000)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position.inside = c(0.5, 0.985),
        legend.position = "inside",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ssb <- df_OM %>%
  filter(qname == "ssb") %>%
  ggplot(aes(x = year)) +
  geom_hline(data = df_refpts, 
             aes(yintercept = Bmsy/1000), 
             colour = "black", linewidth = 0.3, linetype = "2222",
             show.legend = FALSE) + 
  geom_hline(data = df_refpts, 
             aes(yintercept = Blim/1000), 
             colour = "black", linewidth = 0.3, linetype = "1212",
             show.legend = FALSE) + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = df_OM_baseline %>% filter(qname == "ssb"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "Baseline" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "Baseline" = "2121")) +
  labs(y = "SSB (1000t)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text.y = element_blank())
p_fbar <- df_OM %>%
  filter(qname == "fbar") %>%
  ggplot(aes(x = year)) +
  geom_hline(data = df_refpts, 
             aes(yintercept = Fmsy), 
             colour = "black", linewidth = 0.3, linetype = "2222",
             show.legend = FALSE) + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) + 
  geom_line(data = df_OM_baseline %>% filter(qname == "fbar"),
            aes(x = year, y = `50%`, colour = source, linetype = source),
            linewidth = 0.4) +
  facet_wrap(~ OM_label, ncol = 1, strip.position = "right") + 
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "black",
                                     "Baseline" = "red")) + 
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "Baseline" = "2121")) +
  labs(y = "Mean F (ages 3-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")

p <- p_catch + p_R + p_ssb + p_fbar + plot_layout(nrow = 1)
p

ggsave(filename = "output/plots/OM/OM_vs_SAM_all.png", plot = p, 
       width = 16, height = 18, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_vs_SAM_all.pdf", plot = p, 
       width = 16, height = 18, units = "cm")

### ------------------------------------------------------------------------ ###
### baseline OM - historical SSB and Blim ####
### ------------------------------------------------------------------------ ###
### MSY reference points
refpts <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
refpts <- iterMedians(refpts)
### OM stk
stk <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/stk.rds"))

### SSB
qnt <- ssb(stk)/1000
qnt <- window(qnt, end = 2023)
### percentiles
qnt_perc <- quantile(qnt, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                     na.rm = TRUE)
qnt_perc <- as.data.frame(qnt_perc)
qnt_perc <- qnt_perc %>% select(year, iter, data) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "Operating model (median)")

df_iters <- as.data.frame(qnt) %>% 
  select(year, iter, data)

df_refs <- data.frame(source = c("Bmsy", "Blim"),
                      data = c(refpts$Bmsy/1000, refpts$Blim/1000))
p <- ggplot() +
  geom_ribbon(data = qnt_perc,
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = qnt_perc,
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = qnt_perc,
            aes(x = year, y = `50%`),
            colour = "black", show.legend = FALSE) +
  geom_line(data = df_iters %>% filter(iter %in% 1:5),
            aes(x = year, y = data, colour = iter),
            linewidth = 0.1,
            show.legend = FALSE) +
  geom_hline(yintercept = c(refpts$Bmsy/1000), 
             colour = "red", linetype = "2121",
             show.legend = FALSE) + 
  geom_hline(yintercept = c(refpts$Blim/1000), 
             colour = "red", linetype = "1111",
             show.legend = FALSE) + 
  annotate("text", x = 2020, y = c(refpts$Bmsy/1000), 
           label = expression(B[MSY]), size = 3, colour = "red", 
           hjust = 0, vjust = -0.2) +
  annotate("text", x = 2020, y = c(refpts$Blim/1000), 
           label = expression(B[lim]), size = 3, colour = "red", 
           hjust = 0, vjust = -0.2) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "SSB (1000t)") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p
ggsave(filename = "output/plots/OM/OM_SSB_Blim.png", plot = p, 
       width = 10, height = 4, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_SSB_Blim.pdf", plot = p, 
       width = 10, height = 4, units = "cm")


### ------------------------------------------------------------------------ ###
### baseline OM - MSY search ####
### ------------------------------------------------------------------------ ###
MSY_trace <- readRDS("input/ple.27.7e/baseline/1000_100/MSY_trace.rds")
MSY_trace <- as.data.frame(do.call(rbind, MSY_trace))[-length(MSY_trace), ]
MSY_trace <- as.data.frame(apply(MSY_trace, 2, unlist))
MSY_trace <- unique(MSY_trace)

Fmsy <- MSY_trace$Ftrgt[which.max(MSY_trace$catch)]

df_MSY_trace <- MSY_trace %>%
  mutate(catch = catch/1000, ssb = ssb/1000, tsb = tsb/1000, rec = rec/1000) %>%
  pivot_longer(-Ftrgt) %>%
  mutate(name = factor(name, 
                       levels = c("catch", "rec", "tsb", "ssb"),
                       labels = c("Catch (1000t)", "Recruitment (1000s)",
                                  "TSB (1000t)", "SSB (1000t)")))
### MSY level
df_MSY_trace_MSY <- df_MSY_trace %>%
  filter(Ftrgt == Fmsy)

p <- ggplot(data = df_MSY_trace,
       aes(x = Ftrgt, y = value)) +
  geom_hline(data = df_MSY_trace_MSY,
             aes(yintercept = value, colour = "MSY", linetype = "MSY"),
             linewidth = 0.3,
             show.legend = FALSE) +
  geom_vline(data = df_MSY_trace_MSY,
             aes(xintercept = Fmsy, colour = "MSY", linetype = "MSY"),
             linewidth = 0.3,
             show.legend = FALSE) +
  geom_blank(data = df_MSY_trace %>%
               filter(name == "TSB (1000t)") %>%
               mutate(name = "SSB (1000t)")) +
  stat_smooth(data = df_MSY_trace,
              aes(colour = "Loess smoother", linetype = "Loess smoother"), 
              linewidth = 0.5,
              se = FALSE, span = 0.3, n = 100, show.legend = TRUE) + 
  geom_point(size = 0.5, show.legend = FALSE) +
  scale_linetype_manual("", values = c("Loess smoother" = "solid",
                                       "MSY" = "2121")) +
  scale_colour_manual("", values = c("Loess smoother" = "blue",
                                     "MSY" = "red")) +
  facet_wrap(~name, scales = "free_y", strip.position = "left") +
  ylim(c(0, NA)) +
  labs(x = "F (ages 3-6)") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.7),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.6, "lines"))
p
ggsave(filename = "output/plots/OM/OM_baseline_MSY_search.png", plot = p, 
       width = 10, height = 7, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_baseline_MSY_search.pdf", plot = p, 
       width = 10, height = 7, units = "cm")

### ------------------------------------------------------------------------ ###
### baseline OM - projection wormplot - MSY and F0 ####
### ------------------------------------------------------------------------ ###

stk <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/stk.rds"))

### Fmsy
mp_f <- readRDS(paste0("output/ple.27.7e/baseline/1000_100/OM/constF/mp_MSY.rds"))
stk[, ac(2024:2124)] <- mp_f@om@stock
qnts <- FLQuants(catch = catch(stk)/1000, rec = rec(stk)/1000,
                 ssb = ssb(stk)/1000, fbar = fbar(stk))
### percentiles
qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    na.rm = TRUE)
qnts_perc <- FLQuants(qnts_perc)
qnts_perc <- as.data.frame(qnts_perc)
qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
  pivot_wider(names_from = iter, values_from = data)
### individual iterations
qnts_iter <- as.data.frame(iter(qnts, 1:5))
levels <- c(paste0("Catch (1000t)"), "Recruitment (1000s)",
            "SSB (1000t)", 
            paste0("F (ages ", range(stk)[["minfbar"]], "-",
                   range(stk)[["maxfbar"]], ")"))
qnts_perc$qname <- factor(qnts_perc$qname,
                          levels = c("catch", "rec", "ssb", "fbar"),
                          labels = levels)
qnts_iter$qname <- factor(qnts_iter$qname,
                          levels = c("catch", "rec", "ssb", "fbar"),
                          labels = levels)
p <- ggplot() +
  geom_ribbon(data = qnts_perc,
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(data = qnts_perc,
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_line(data = qnts_iter,
            aes(x = year, y = data, colour = iter), 
            linewidth = 0.1, show.legend = FALSE) + 
  geom_line(data = qnts_perc, mapping = aes(x = year, y = `50%`),
            linewidth = 0.4) +
  geom_vline(xintercept = 2024.5, linewidth = 0.3) + 
  geom_vline(xintercept = 2115, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = 2124, colour = "red", linewidth = 0.3) +
  annotate("rect", xmin = 2115, xmax = 2124, ymin = -Inf, ymax = Inf,
           alpha = .2, fill = "red") + 
  facet_wrap(~ qname, scales = "free_y", ncol = 1, strip.position = "left") +
  ylim(c(0, NA)) + 
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_baseline_MSY_worm.png", plot = p, 
       width = 16, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_baseline_MSY_worm.pdf", plot = p, 
       width = 16, height = 13, units = "cm")

### same for F=0
mp_0 <- readRDS(paste0("output/ple.27.7e/baseline/1000_100/OM/constF/mp_0.rds"))
stk[, ac(2024:2124)] <- mp_0@om@stock
qnts <- FLQuants(catch = catch(stk)/1000, rec = rec(stk)/1000,
                 ssb = ssb(stk)/1000, fbar = fbar(stk))
### percentiles
qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    na.rm = TRUE)
qnts_perc <- FLQuants(qnts_perc)
qnts_perc <- as.data.frame(qnts_perc)
qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
  pivot_wider(names_from = iter, values_from = data)
### individual iterations
qnts_iter <- as.data.frame(iter(qnts, 1:5))
levels <- c(paste0("Catch (1000t)"), "Recruitment (1000s)",
            "SSB (1000t)", 
            paste0("F (ages ", range(stk)[["minfbar"]], "-",
                   range(stk)[["maxfbar"]], ")"))
qnts_perc$qname <- factor(qnts_perc$qname,
                          levels = c("catch", "rec", "ssb", "fbar"),
                          labels = levels)
qnts_iter$qname <- factor(qnts_iter$qname,
                          levels = c("catch", "rec", "ssb", "fbar"),
                          labels = levels)
p <- ggplot() +
  geom_ribbon(data = qnts_perc,
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(data = qnts_perc,
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_line(data = qnts_iter,
            aes(x = year, y = data, colour = iter), 
            linewidth = 0.1, show.legend = FALSE) + 
  geom_line(data = qnts_perc, mapping = aes(x = year, y = `50%`),
            linewidth = 0.4) +
  geom_vline(xintercept = 2024.5, linewidth = 0.3) + 
  geom_vline(xintercept = 2115, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = 2124, colour = "red", linewidth = 0.3) +
  annotate("rect", xmin = 2115, xmax = 2124, ymin = -Inf, ymax = Inf,
           alpha = .2, fill = "red") + 
  facet_wrap(~ qname, scales = "free_y", ncol = 1, strip.position = "left") +
  ylim(c(0, NA)) + 
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_baseline_F0_worm.png", plot = p, 
       width = 16, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_baseline_F0_worm.pdf", plot = p, 
       width = 16, height = 13, units = "cm")

### ------------------------------------------------------------------------ ###
### all OMs - visualisation of OM MSY values ####
### ------------------------------------------------------------------------ ###
OM_list <- c("baseline",
             "Catch_no_disc", "Catch_no_surv", "migr_none",
             "M_low", "M_high", "M_Gislason",
             "R_no_AC", "R_higher", "R_lower", "R_h_lower")
OM_labels <- c("Baseline",
               "Catch: no discards", "Catch: 100% discards", "Catch: no migration",
               "M: -50%", "M: +50%", "M: Gislason",
               "R: no AC", "R: +20%", "R: -20%", "R: h -20%")
MSY_runs <- foreach(OM = OM_list, .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/ple.27.7e/", OM, "/1000_100/MSY_trace.rds")
  if (!file.exists(file_i)) return(NULL)
  runs_i <- readRDS(file_i)
  runs_i <- bind_rows(runs_i)
  runs_i$OM <- OM
  return(runs_i)
}
MSY_runs <- MSY_runs %>%
  group_by(OM) %>%
  mutate(MSY = ifelse(catch == max(catch), TRUE, FALSE))

MSY_runs_plot <- MSY_runs %>%
  pivot_longer(c(catch, ssb)) %>%
  mutate(value = value/1000) %>%
  mutate(label = factor(name, levels = c("catch", "ssb"),
                        labels = c("Catch (1000t)", "SSB (1000t)")),
         OM_label = factor(OM, 
                           levels = OM_list,
                           labels = OM_labels))

### dummy df to get same x-axis limit
df_x1 <- head(MSY_runs_plot, 2) %>%
  mutate(value = c(1.7, 50))
df_x2 <- tail(MSY_runs_plot %>% filter(OM == "M_low"), 2) %>%
  mutate(value = c(1.7, 50))
df_x3 <- tail(MSY_runs_plot, 2) %>%
  mutate(value = c(1.7, 50))

### plot - split into two panels
p1 <- MSY_runs_plot %>%
  filter(OM_label %in% OM_labels[1:4]) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_blank(data = df_x1) +
  geom_point(size = 0.4) +
  geom_smooth(span = 0.4, method = "loess",
              se = FALSE, linewidth = 0.2, colour = "blue") +
  geom_vline(data = . %>%
               filter(MSY == TRUE),
             aes(xintercept = Ftrgt), linewidth = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "mean F (ages 3-6)") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title = element_blank())
p2 <- MSY_runs_plot %>%
  filter(OM_label %in% OM_labels[5:8]) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_blank(data = df_x2) +
  geom_point(size = 0.4) +
  geom_smooth(span = 0.4, method = "loess",
              se = FALSE, linewidth = 0.2, colour = "blue") +
  geom_vline(data = . %>%
               filter(MSY == TRUE),
             aes(xintercept = Ftrgt), linewidth = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "mean F (ages 3-6)") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title = element_blank())
p3 <- MSY_runs_plot %>%
  filter(OM_label %in% OM_labels[9:11]) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_blank(data = df_x3) +
  geom_point(size = 0.4) +
  geom_smooth(span = 0.4, method = "loess",
              se = FALSE, linewidth = 0.2, colour = "blue") +
  geom_vline(data = . %>%
               filter(MSY == TRUE),
             aes(xintercept = Ftrgt), linewidth = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "mean F (ages 3-6)") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank())
design <- "
  11
  22
  3#
"
p <- p1 + p2 + p3 + plot_layout(design = design, widths = c(1, 0.325))
p
ggsave(filename = "output/plots/OM/OM_all_MSY_search.png", plot = p, 
       width = 16, height = 12.5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_all_MSY_search.pdf", plot = p, 
       width = 16, height = 12.5, units = "cm")

### reference points table for working document
### get Blim
Blims <- foreach(OM = OM_list, .combine = bind_rows) %do% {#browser()
    file_i <- paste0("input/ple.27.7e/", OM, "/1000_100/refpts_mse.rds")
    if (!file.exists(file_i)) return(NULL)
    i <- readRDS(file_i)
    return(data.frame(OM = OM, Blim = median(c(i["Blim"]))))
}

### table
### MSY values
refpts_table <- MSY_runs %>%
  filter(MSY == TRUE) %>%
  select(-tsb) %>% unique() %>%
  ### unfished values
  full_join(MSY_runs %>%
              filter(Ftrgt == 0) %>%
              select(ssb, rec, OM) %>%
              rename(B0 = ssb, R0 = rec)) %>%
  ### Blim
  full_join(Blims) %>%
  select(OM, B0, R0, Ftrgt, catch, ssb, rec, Blim) %>%
  rename(FMSY = Ftrgt, MSY = catch, BMSY = ssb, RMSY = rec) 

write.csv(refpts_table, file = "input/OM_refpts.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### all OMs - Blim alternatives ####
### ------------------------------------------------------------------------ ###
### create table with Blim alternatives

### get B0
df_B0 <- MSY_runs %>%
  filter(Ftrgt == 0) %>%
  select(OM, B0 = ssb)

### get Bmsy
df_Bmsy <- MSY_runs %>%
  filter(MSY == TRUE) %>%
  select(OM, Bmsy = ssb)


### get Bloss
df_Bloss <- foreach(OM = OM_list, .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/ple.27.7e/", OM, "/1000_100/stk.rds")
  i <- readRDS(file_i)
  ssb_i <- iterMedians(window(ssb(i), end = 2023))
  which.min(ssb_i)
  return(data.frame(OM = OM, 
                    Bloss_year = dimnames(i)$year[which.min(ssb_i)],
                    Bloss = min(ssb_i)))
}

### Blim based on RR0
df_BlimRR0 <- foreach(OM = OM_list, .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/ple.27.7e/", OM, "/1000_100/refpts_mse.rds")
  i <- readRDS(file_i)
  return(data.frame(OM = OM, 
                    Blim_RR0 = median(c(i["Blim_RR0"]))))
}

### Blim based on BB0
df_BlimBB0 <- foreach(OM = OM_list, .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/ple.27.7e/", OM, "/1000_100/refpts_mse.rds")
  i <- readRDS(file_i)
  return(data.frame(OM = OM, 
                    Blim_BB0 = median(c(i["Blim_BB0"]))))
}


### add relative values
df_Bmsy <- df_Bmsy %>% full_join(df_B0) %>%
  mutate(BmsyB0 = Bmsy/B0) %>%
  select(-B0)
df_Bloss <- df_Bloss %>% full_join(df_B0) %>% full_join(df_Bmsy) %>%
  mutate(BlossB0 = Bloss/B0,
         BlossBmsy = Bloss/Bmsy) %>%
  select(-B0, -Bmsy)
df_BlimRR0 <- df_BlimRR0 %>% full_join(df_B0) %>% full_join(df_Bmsy) %>%
  mutate(Blim_RR0B0 = Blim_RR0/B0,
         Blim_RR0Bmsy = Blim_RR0/Bmsy) %>%
  select(-B0, -Bmsy)
df_BlimBB0 <- df_BlimBB0 %>% full_join(df_B0) %>% full_join(df_Bmsy) %>%
  mutate(Blim_BB0B0 = Blim_BB0/B0,
         Blim_BB0Bmsy = Blim_BB0/Bmsy) %>%
  select(-B0, -Bmsy)

df_Blim <- df_B0 %>%
  full_join(df_Bmsy) %>%
  full_join(df_Bloss) %>%
  full_join(df_BlimRR0) %>%
  full_join(df_BlimBB0)
write.csv(df_Blim, file = "input/OM_refpts_Blim.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### Recruitment model and residual visualisation (baseline OM) ####
### ------------------------------------------------------------------------ ###

### SR model fit
fit <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_fit.rds")
sr_om <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
stk <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
stk <- window(iterMedians(stk), end = 2023)

fit_stk <- SAM2FLStock(object = fit, stk = stk)
sr <- as.FLSR(fit_stk, model = "bevholtSV")
rec(sr) <- rec(sr)/1000
ssb(sr) <- ssb(sr)/1000
sr <- fmle(sr, method = 'L-BFGS-B', fixed = list(), 
           control = list(trace = 0))
sr_params <- abPars("bevholt", s = params(sr)["s"], v = params(sr)["v"], 
                    spr0 = params(sr)["spr0"])
sr_params <- list(a = c(sr_params$a), b = c(sr_params$b))
sr_model <- function(ssb, a, b) {(a*ssb)/(b + ssb)}

### median - data.frame with data
df_sr <- data.frame(year = as.numeric(dimnames(sr)$year),
                    ssb = c(ssb(sr)), 
                    rec = c(rec(sr)), 
                    fitted = c(fitted(sr)),
                    residuals = c(residuals(sr)))
### median - get density
dens <- density(x = df_sr$residuals)
df_dens <- data.frame(x = dens$x, y = dens$y)

### iterations - get density
df_dens_iter <- lapply(seq(dim(sr_om)[6]), function(i) {
  dens <- density(log(c(iter(sr_om@residuals, i))))
  data.frame(x = dens$x, y = dens$y, iter = i)
})
df_dens_iter <- bind_rows(df_dens_iter)

p_res <- df_sr %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_blank(data = data.frame(ssb = 9, rec = 18)) +
  geom_linerange(aes(x = ssb, ymin = rec, ymax = fitted), 
                 colour = "blue", linewidth = 0.2, 
                 linetype = "2121") +
  geom_point(size = 0.5) +
  geom_function(fun = sr_model, args = sr_params, size = 0.4) +
  scale_colour_manual(values = c(observation = "black", model = "black")) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     limits = c(0, NA), expand = c(0, 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     limits = c(0, NA), expand = c(0, 0)) +
  labs(x = "SSB (1000t)", y = "Recruitment (1000s)") +
  theme_bw(base_size = 8)
df_dens_max <- max(df_dens$y)
p_hist <- df_sr %>%
  ggplot(aes(residuals)) +
  geom_histogram(bins = 13, colour = "black", fill = "white", size = 0.4,
                 linewidth = 0.4) +
  geom_line(aes(x = x, y = y * 10 * df_dens_max, colour = "kernel"), 
            data = df_dens, show.legend = TRUE, linewidth = 0.4) + 
  scale_x_continuous(breaks = c(-1, 0, 1),
                     limits = c(-1.5, 1.5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     sec.axis = sec_axis(~./(10 * df_dens_max), 
                                         name = "Density", 
                                         breaks = scales::pretty_breaks()),
                     limits = c(0, 10.5), expand = c(0, 0)) +
  scale_colour_manual("", values = c(kernel = "red"), 
                      labels = c("kernel\ndensity")) +
  labs(y = "Residual count", x = "Log residuals") +
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank())
df_dens_iter_max <- max(df_dens_iter$y)
p_hist_iter <- ggplot() +
  geom_line(data = df_dens_iter, 
            aes(x = x, y = y, group = iter), 
            linewidth = 0.1, alpha = 0.05) + 
  scale_y_continuous(position = "right",
                     breaks = scales::pretty_breaks(),
                     limits = c(0, df_dens_iter_max*1.05), expand = c(0, 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     limits = c(-2.5, 2.5)) +
  labs(y = "Density", x = "Log residuals") +
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90))
### same plot but rasterised to reduce PDF file size
p_hist_iter_raster <- ggplot() +
  ggrastr::rasterise(geom_line(data = df_dens_iter, 
                               aes(x = x, y = y, group = iter), 
                               linewidth = 0.1, alpha = 0.05), dpi = 1200) + 
  scale_y_continuous(position = "right",
                     breaks = scales::pretty_breaks(),
                     limits = c(0, df_dens_iter_max*1.05), expand = c(0, 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     limits = c(-2.5, 2.5)) +
  labs(y = "Density", x = "Log residuals") +
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90))

p <- p_res + p_hist + p_hist_iter + plot_annotation(tag_levels = "a", 
                                 tag_prefix = "(", tag_suffix = ")")  &
  theme(plot.tag = element_text(face = "bold"))
p
ggsave(filename = "output/plots/OM/OM_rec_res.png", plot = p,
       width = 16, height = 5, units = "cm", dpi = 600, type = "cairo")
### rasterised version for PDF
p <- p_res + p_hist + p_hist_iter_raster + 
  plot_annotation(tag_levels = "a", 
                  tag_prefix = "(", tag_suffix = ")")  &
  theme(plot.tag = element_text(face = "bold"))
p
ggsave(filename = "output/plots/OM/OM_rec_res.pdf", plot = p,
       width = 16, height = 5, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### Recruitment models of baseline OM ####
### ------------------------------------------------------------------------ ###

sr_ple <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
df_sr_median <- data.frame(year = as.numeric(dimnames(ssb(sr_ple))$year),
                               ssb = c(iterMedians(ssb(sr_ple))), 
                               rec = c(iterMedians(rec(sr_ple)))) %>%
  filter(year <= 2023) %>%
  mutate(type = "median")
df_sr_iter <- data.frame(year = as.numeric(dimnames(ssb(sr_ple))$year),
                             ssb = c(ssb(sr_ple)), 
                             rec = c(rec(sr_ple)),
                             iter = c(as.numeric(dimnames(sr_ple)$iter))) %>%
  filter(year <= 2023) %>%
  mutate(type = "replicates")
df_sr <- bind_rows(df_sr_median, df_sr_iter) %>%
  mutate(type = factor(type, levels = c("median", "replicates"),
                       labels = c("median (stock-\nrecruit pairs)", 
                                  "replicates (stock-\nrecruit pairs)")))
sr_pars_ple <- abPars("bevholt", s = params(sr_ple)["s"], 
                      v = params(sr_ple)["v"], 
                      spr0 = params(sr_ple)["spr0"])
sr_pars_ple_med <- list(a = c(iterMedians(sr_pars_ple$a)/1000), 
                    b = c(iterMedians(sr_pars_ple$b)/1000))
sr_model_bevholt <- function(ssb, a, b) {(a*ssb)/(b + ssb)}
ssbs_ple <- seq(from = 0, to = 10000, by = 10)
sr_ple_df <- lapply(1:1000, function(x) {
  data.frame(ssb = ssbs_ple,
             rec = sr_model_bevholt(ssb = ssbs_ple, 
                                    a = c(sr_pars_ple$a[, x]),
                                    b = c(sr_pars_ple$b[, x])),
             iter = x)
})
sr_ple_df <- bind_rows(sr_ple_df)
sr_df <- bind_rows(data.frame(ssb = ssbs_ple,
                              rec = sr_model_bevholt(ssb = ssbs_ple, 
                                                     a = c(iterMedians(sr_pars_ple$a)),
                                                     b = c(iterMedians(sr_pars_ple$b))),
                              iter = 0) %>%
                     mutate(type = "median"),
                   sr_ple_df %>%
                     mutate(type = "replicates")) %>%
  mutate(type = factor(type, levels = c("median", "replicates"),
                       labels = c("median (model)", 
                                  "replicates (model)")))

p <- ggplot() +
  ### print points/lines separately so that median is on top
  geom_point(data = df_sr %>% 
               filter(type == "replicates (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type), 
             alpha = 0.01, size = 0.3) +
geom_line(data = sr_df %>% filter(type == "replicates (model)"),
          aes(x = ssb/1000, y = rec/1000, group = iter,
              colour = type),
          linewidth = 0.1, alpha = 0.05) +
  geom_point(data = df_sr %>% filter(type == "median (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type),
             size = 0.4) +
  geom_line(data = sr_df %>% filter(type == "median (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            linewidth = 0.4, alpha = 1) +
  scale_colour_manual("", 
    values = c("median (stock-\nrecruit pairs)" = "red", 
               "replicates (stock-\nrecruit pairs)" = "black",
               "median (model)" = "red", 
               "replicates (model)" = "black")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous("SSB (1000 t)", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Recruitment (1000s)") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 17.5), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))

p
ggsave(filename = "output/plots/OM/OM_rec_baseline.png", plot = p,
       width = 10, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_baseline.pdf", plot = p,
       width = 10, height = 5, units = "cm", dpi = 600)

### same plot but rasterised to reduce PDF file size
p <- ggplot() +
  ### print points/lines separately so that median is on top
  ggrastr::rasterise(geom_point(data = df_sr %>% 
                                  filter(type == "replicates (stock-\nrecruit pairs)"),
                                aes(x = ssb/1000, y = rec/1000, 
                                    colour = type), 
                                alpha = 0.01, size = 0.3), dpi = 600) +
  ggrastr::rasterise(geom_line(data = sr_df %>% filter(type == "replicates (model)"),
                               aes(x = ssb/1000, y = rec/1000, group = iter,
                                   colour = type),
                               size = 0.1, alpha = 0.05), dpi = 600) +
  geom_point(data = df_sr %>% filter(type == "median (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type),
             size = 0.4) +
  geom_line(data = sr_df %>% filter(type == "median (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            size = 0.4, alpha = 1) +
  scale_colour_manual("", 
                      values = c("median (stock-\nrecruit pairs)" = "red", 
                                 "replicates (stock-\nrecruit pairs)" = "black",
                                 "median (model)" = "red", 
                                 "replicates (model)" = "black")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous("SSB (1000 t)", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Recruitment (1000s)") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 17.5), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p
ggsave(filename = "output/plots/OM/OM_rec_baseline.pdf", plot = p,
       width = 10, height = 5, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### Recruitment - SAM estimates vs simulated ####
### ------------------------------------------------------------------------ ###

sr <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")

### auto-correlation from baseline OM
sr_rho <- 0.543945377878967573082036324194632470607757568359375
yrs_res <- 1982:2023

### create simulated residuals for historical period
### same approach as used in create_OM()
res_new <- foreach(iter_i = seq(dim(sr)[6])) %do% {#browser()
  set.seed(iter_i)
  ### get residuals for current iteration
  ### log residuals here - create_OM() used exp() to get multiplicative values
  res_i <- c(FLCore::iter(log(residuals(sr)), iter_i))
  res_i <- res_i[!is.na(res_i)]
  ### calculate kernel density of residuals
  density <- density(x = res_i)
  ### sample residuals
  mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
  ### "smooth", i.e. sample from density distribution
  res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
  ### "add" autocorrelation
  if (TRUE) {
    sr_acf_i <- acf(res_i, lag.max = 1, plot = FALSE, na.action = na.exclude)
    sr_rho_i <- sr_acf_i$acf[2]
    res_ac <- rep(0, length(yrs_res))
    res_ac[1] <- sr_rho * tail(res_i, 1) + sqrt(1 - sr_rho^2) * res_new[1]
    for (r in 2:length(res_ac)) {
      res_ac[r] <- sr_rho * res_ac[r - 1] + sqrt(1 - sr_rho^2) * res_new[r]
    }
  }
  return(list(default = res_new, ac = res_ac))
}
### residuals
res_default <- res_ac <- window(residuals(sr), end = 2023) %=% NA_real_
res_default[] <- exp(unlist(lapply(res_new, "[", "default")))
res_ac[] <- exp(unlist(lapply(res_new, "[", "ac")))
### modelled recruitment (fitted value plus noise)
r_default <- window(fitted(sr), end = 2023) * res_default
r_ac <- window(fitted(sr), end = 2023) * res_ac
### combine data and format
qnts_r <- FLQuants(r_data = window(rec(sr), end = 2023),
                   r_modelled = r_default,
                   r_modelled_ac = r_ac)
df <- as.data.frame(window(ssb(sr), start = 1982, end = 2023)) %>%
  select(year, SSB = data, iter) %>%
  full_join(as.data.frame(qnts_r) %>%
              select(year, R = data, source = qname, iter),
            relationship = "many-to-many")
df <- bind_rows(df %>%
                  filter(iter %in% 1:3),
                df %>%
                  mutate(iter = "all")) %>%
  mutate(iter_group = ifelse(iter == "all", "all", "iter")) %>%
  mutate(iter_group = factor(iter_group, levels = c("iter", "all"))) %>%
  mutate(source = factor(source, 
                         levels = c("r_modelled", "r_modelled_ac", "r_data"),
                         labels = c("Modelled",
                                    "Modelled (with\nauto-correlation)",
                                    "Data")))
# df <- df %>%
#   mutate(source = factor(source, levels = rev(levels(source))))

### plot SR pairs
p_pairs <- df %>%
  ggplot(aes(x = SSB/1000, y = R/1000, colour = source, alpha = iter_group)) +
  geom_point(shape = 19, size = 0.4) + 
  scale_alpha_manual(values = c("iter" = 0.5, "all" = 0.05), 
                     guide = "none") +
  scale_colour_manual("", values = rev(scales::hue_pal()(3))) + 
  facet_wrap(~ iter, nrow = 1) +
  xlim(c(0, NA)) +
  ylim(c(0, NA)) + 
  labs(x = "SSB (1000t)", y = "Recruitment (millions)") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p_pairs
### rasterised version for PDF
p_pairs_raster <- df %>%
  #filter(source == "Data") %>%
  ggplot(aes(x = SSB/1000, y = R/1000, colour = source, alpha = iter_group)) +
  ggrastr::rasterise(geom_point(shape = 19, size = 0.4), dpi = 600) + 
  scale_alpha_manual(values = c("iter" = 0.5, "all" = 0.05), 
                     guide = "none") +
  scale_colour_manual("", values = rev(scales::hue_pal()(3))) + 
  facet_wrap(~ iter, nrow = 1) +
  xlim(c(0, NA)) +
  ylim(c(0, NA)) + 
  labs(x = "SSB (1000t)", y = "Recruitment (millions)") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p_pairs_raster

### plot empirical cumulative distribution
p_ecdf <- df %>%
  ggplot() +
  stat_ecdf(aes(R/1000, colour = source), linewidth = 0.4, alpha = 0.8) +
  scale_colour_manual("", values = rev(scales::hue_pal()(3))) + 
  facet_wrap(~ iter, nrow = 1) +
  xlim(c(0, NA)) +
  ylim(c(0, NA)) + 
  labs(x = "Recruitment (millions)", y = "Empirical cumulative\ndistribution") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p_ecdf  

p <- p_pairs / p_ecdf
p
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model.png", plot = p,
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
p <- p_pairs_raster / p_ecdf
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model.pdf", plot = p,
       width = 16, height = 8, units = "cm", dpi = 600)

### average by SSB blocks
p <- df %>%
  filter(iter == "all") %>%
  mutate(SSB = SSB/1000,
         R = R/1000) %>%
  filter(SSB <= 9) %>%
  mutate(SSB_group = ceiling(SSB)) %>%
  group_by(source, SSB_group) %>%
  summarise(R = median(R)) %>%
  #filter(SSB_group <= 9) %>%
  ggplot(aes(x = SSB_group, y = R, colour = source, linetype = source)) +
  geom_step() +
  scale_colour_manual("", values = rev(scales::hue_pal()(3))) + 
  scale_linetype("") +
  scale_x_continuous(breaks = seq(0, 9, 2), limits = c(0, NA)) +
  labs(x = "SSB (1000t)", y = "Recruitment (1000s)") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model_blocks.png", plot = p,
       width = 10, height = 4, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model_blocks.pdf", plot = p,
       width = 10, height = 4, units = "cm")

### average by year
df %>%
  filter(iter == "all") %>%
  mutate(SSB = SSB/1000,
         R = R/1000) %>%
  group_by(source, year) %>%
  summarise(R = median(R)) %>%
  ggplot(aes(x = year, y = R, colour = source)) +
  geom_step()# +

### boxplots
p <- df %>%
  filter(iter == "all" & source %in% c("Data", "Modelled (with\nauto-correlation)")) %>%
  mutate(SSB = SSB/1000,
         R = R/1000) %>%
  group_by(source, year) %>%
  ggplot(aes(x = year, y = R, colour = source, group = interaction(year, source))) +
  geom_boxplot(outliers = FALSE, linewidth = 0.1) +
  scale_colour_brewer("", palette = "Set1") +
  labs(x = "Year", y = "Recruitment (1000s)") +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model_ts.png", plot = p,
       width = 10, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model_ts.pdf", plot = p,
       width = 10, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### Recruitment models of alternative OMs ####
### ------------------------------------------------------------------------ ###

OM_list <- c("baseline",
             "Catch_no_disc", "Catch_no_surv", "migr_none",
             "M_low", "M_high", "M_Gislason",
             "R_no_AC", "R_failure", "R_higher", "R_lower", "R_h_lower")
OM_labels <- c("Baseline",
               "Catch: no discards", "Catch: 100% discards", "Catch: no migration",
               "M: -50%", "M: +50%", "M: Gislason",
               "R: no AC", "R: failure", "R: +20%", "R: -20%", "R: h -20%")

### load recruitment values and stock-recruit pairs - use medians
df_OM <- foreach(OM = OM_list,
                 OM_label = OM_labels) %do% {
  #browser()
  sr_i <- readRDS(paste0("input/ple.27.7e/", OM,
                         "/1000_100/sr.rds"))
  ### SSB-recruitment pairs used for model fitting
  df_pairs <- data.frame(year = as.numeric(dimnames(ssb(sr_i))$year),
                         ssb = c(iterMedians(ssb(sr_i))), 
                         rec = c(iterMedians(rec(sr_i)))) %>%
    filter(!is.na(ssb)) %>%
    mutate(OM = OM, OM_label = OM_label)
  ### table with modelled recruitment
  sr_i_med <- iter(sr_i, 1)
  params(sr_i_med) <- iterMedians(params(sr_i))
  ssb_max <- max(c(iterMedians(ssb(sr_i))), na.rm = TRUE)
  
  ssb_vals <- seq(from = 0, to = ssb_max * 3, length.out = 1500)
  rec_vals <- c(predict(sr_i_med, ssb = FLQuant(ssb_vals)))

  df_pairs_model <- data.frame(ssb = ssb_vals,
                               rec = rec_vals) %>%
    mutate(OM = OM, OM_label = OM_label, rec_failure = FALSE)
  if (isTRUE(OM == "R_failure")) {
    df_pairs_model <- bind_rows(
      df_pairs_model %>%
        mutate(rec_failure = FALSE),
      df_pairs_model %>%
        mutate(rec_failure = TRUE)) %>%
      mutate(rec = ifelse(rec_failure, rec/10, rec))
  }
  ### recruitment model parameters
  pars_i <- data.frame(a = c(iterMedians(params(sr_i)$a)),
                       b = c(iterMedians(params(sr_i)$b))) %>%
    mutate(OM = OM, OM_label = OM_label)
  return(list(pairs = df_pairs, pairs_add = df_pairs_model, pars = pars_i))
}
df_pairs <- bind_rows(lapply(df_OM, "[[", 1))
df_pairs_add <- bind_rows(lapply(df_OM, "[[", 2))
df_pars <- bind_rows(lapply(df_OM, "[[", 3))

df_pairs <- df_pairs %>%
  mutate(OM_label = factor(OM_label, levels = OM_labels))
df_pairs_add <- df_pairs_add %>%
  mutate(OM_label = factor(OM_label, levels = OM_labels))
df_pars <- df_pars %>%
  mutate(OM_label = factor(OM_label, levels = OM_labels))

p <- ggplot() +
  geom_point(data = df_pairs,
             aes(x = ssb/1000, y = rec/1000),
             size = 0.3) +
  geom_line(data = df_pairs_add,
            aes(x = ssb/1000, y = rec/1000, linetype = rec_failure),
            linewidth = 0.4, show.legend = FALSE) + 
  geom_line(data = df_pairs_add %>%
              filter(OM == "baseline") %>%
              select(-OM, -OM_label),
            aes(x = ssb/1000, y = rec/1000),
            colour = "red", linetype = "2121", linewidth = 0.4) +
  facet_wrap(~ OM_label, nrow = 2) +
  scale_x_continuous("SSB (1000t)",
                     breaks = seq(0, 10, 2)) +
  # scale_y_continuous("Recruitment (1000s)",
  #                    breaks = c(0, 10000, 20000, 30000),
  #                    labels = c(0, 10, 20, 30)) +
  coord_cartesian(xlim = c(0, 11.5), ylim = c(0, 29), expand = FALSE) +
  labs(x = "SSB (1000t)", y = "Recruitment (1000s)") +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/OM/OM_rec_OMs.png", plot = p,
       width = 16, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_OMs.pdf", plot = p,
       width = 16, height = 6, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### Recruitment models - baseline OM - steepness ####
### ------------------------------------------------------------------------ ###
sr <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")

### get steepness
sr_s <- c(params(sr)["s"])
summary(sr_s)
sd(sr_s)
mean(sr_s > 0.95)

p <- data.frame(s = sr_s) %>%
  ggplot(aes(x = s)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)),
                 binwidth = 0.1, fill = "grey", colour = "black", 
                 linewidth = 0.3) +
  labs(x = "Recruitment steepness", y = "Proportion") +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/OM/OM_rec_steepness.png", plot = p,
       width = 8, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_steepness.pdf", plot = p,
       width = 8, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### Recruitment models of OM with lower steepness ####
### ------------------------------------------------------------------------ ###

sr_ple_lowh <- readRDS("input/ple.27.7e/R_h_lower/1000_100/sr.rds")
df_sr_median_lowh <- data.frame(year = as.numeric(dimnames(ssb(sr_ple_lowh))$year),
                           ssb = c(iterMedians(ssb(sr_ple_lowh))), 
                           rec = c(iterMedians(rec(sr_ple_lowh)))) %>%
  filter(year <= 2023) %>%
  mutate(type = "median")
df_sr_iter_lowh <- data.frame(year = as.numeric(dimnames(ssb(sr_ple_lowh))$year),
                         ssb = c(ssb(sr_ple_lowh)), 
                         rec = c(rec(sr_ple_lowh)),
                         iter = c(as.numeric(dimnames(sr_ple_lowh)$iter))) %>%
  filter(year <= 2023) %>%
  mutate(type = "replicates")
df_sr_lowh <- bind_rows(df_sr_median_lowh, df_sr_iter_lowh) %>%
  mutate(type = factor(type, levels = c("median", "replicates"),
                       labels = c("median (stock-\nrecruit pairs)", 
                                  "replicates (stock-\nrecruit pairs)")))
sr_pars_ple_lowh <- abPars("bevholt", s = params(sr_ple_lowh)["s"], 
                      v = params(sr_ple_lowh)["v"], 
                      spr0 = params(sr_ple_lowh)["spr0"])
sr_pars_ple_med_lowh <- list(a = c(iterMedians(sr_pars_ple_lowh$a)/1000), 
                        b = c(iterMedians(sr_pars_ple_lowh$b)/1000))
sr_model_bevholt <- function(ssb, a, b) {(a*ssb)/(b + ssb)}
ssbs_ple <- seq(from = 0, to = 10000, by = 10)
sr_ple_df_lowh <- lapply(1:1000, function(x) {
  data.frame(ssb = ssbs_ple,
             rec = sr_model_bevholt(ssb = ssbs_ple, 
                                    a = c(sr_pars_ple_lowh$a[, x]),
                                    b = c(sr_pars_ple_lowh$b[, x])),
             iter = x)
})
sr_ple_df_lowh <- bind_rows(sr_ple_df_lowh)
sr_df_lowh <- bind_rows(data.frame(ssb = ssbs_ple,
                              rec = sr_model_bevholt(ssb = ssbs_ple, 
                                                     a = c(iterMedians(sr_pars_ple_lowh$a)),
                                                     b = c(iterMedians(sr_pars_ple_lowh$b))),
                              iter = 0) %>%
                     mutate(type = "median"),
                   sr_ple_df_lowh %>%
                     mutate(type = "replicates")) %>%
  mutate(type = factor(type, levels = c("median", "replicates"),
                       labels = c("median (model)", 
                                  "replicates (model)")))

### combine baseline OM rec and lower h rec
sr_df_combined <- bind_rows(sr_df %>% mutate(model = "Baseline"),
                            sr_df_lowh %>% mutate(model = "Lower steepness")) %>%
  mutate(model = factor(model, levels = c("Baseline", "Lower steepness")))
df_sr_combined <- bind_rows(df_sr %>% mutate(model = "Baseline"),
                            df_sr_lowh %>% mutate(model = "Lower steepness")) %>%
  mutate(model = factor(model, levels = c("Baseline", "Lower steepness")))

p_lowh <- ggplot() +
  ### print points/lines separately so that median is on top
  geom_point(data = df_sr_combined %>% 
               filter(type == "replicates (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type), 
             alpha = 0.01, size = 0.3) +
  geom_line(data = sr_df_combined %>% filter(type == "replicates (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            linewidth = 0.1, alpha = 0.05) +
  geom_point(data = df_sr_combined %>% filter(type == "median (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type),
             size = 0.4) +
  geom_line(data = sr_df_combined %>% filter(type == "median (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            linewidth = 0.4, alpha = 1) +
  scale_colour_manual("", 
                      values = c("median (stock-\nrecruit pairs)" = "red", 
                                 "replicates (stock-\nrecruit pairs)" = "black",
                                 "median (model)" = "red", 
                                 "replicates (model)" = "black")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ model) +
  scale_x_continuous("SSB (1000 t)", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Recruitment (1000s)") +
  coord_cartesian(xlim = c(0, 9.9), ylim = c(0, 17.5), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p_lowh
ggsave(filename = "output/plots/OM/OM_rec_lowh.png", plot = p_lowh,
       width = 15.9, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_lowh.pdf", plot = p_lowh,
       width = 15.9, height = 5, units = "cm", dpi = 600)

### same plot but rasterised to reduce PDF file size
p_lowh <- ggplot() +
  ### print points/lines separately so that median is on top
  ggrastr::rasterise(geom_point(data = df_sr_combined %>% 
                                  filter(type == "replicates (stock-\nrecruit pairs)"),
                                aes(x = ssb/1000, y = rec/1000, 
                                    colour = type), 
                                alpha = 0.01, size = 0.3), dpi = 600) +
  ggrastr::rasterise(geom_line(data = sr_df_combined %>% filter(type == "replicates (model)"),
                               aes(x = ssb/1000, y = rec/1000, group = iter,
                                   colour = type),
                               linewidth = 0.1, alpha = 0.05), dpi = 600) +
  geom_point(data = df_sr_combined %>% filter(type == "median (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type),
             size = 0.4) +
  geom_line(data = sr_df_combined %>% filter(type == "median (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            linewidth = 0.4, alpha = 1) +
  scale_colour_manual("", 
                      values = c("median (stock-\nrecruit pairs)" = "red", 
                                 "replicates (stock-\nrecruit pairs)" = "black",
                                 "median (model)" = "red", 
                                 "replicates (model)" = "black")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ model) +
  scale_x_continuous("SSB (1000 t)", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Recruitment (1000s)") +
  coord_cartesian(xlim = c(0, 9.9), ylim = c(0, 17.5), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.7, "lines"))
p_lowh
ggsave(filename = "output/plots/OM/OM_rec_lowh.pdf", plot = p_lowh,
       width = 15.9, height = 5, units = "cm", dpi = 600)




### ------------------------------------------------------------------------ ###
### surveys (catchability and weights at age) - baseline OM ####
### ------------------------------------------------------------------------ ###
### catchability
### index weights at age

### load data
uncertainty <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_uncertainty.rds")
idx <- readRDS("input/ple.27.7e/baseline/1000_100/idx.rds")
idx_dev <- readRDS("input/ple.27.7e/baseline/1000_100/idx_dev.rds")

### catchability
q <- uncertainty$survey_catchability
q <- lapply(q, function(x) propagate(iterMedians(x, 1), dim(x)[6] + 1))
q <- FLQuants(q)
names(q) <- c("UK-FSP", "Q1SWBeam")
### add iterations
q$`UK-FSP`[,,,,, -1] <- uncertainty$survey_catchability[[1]]
q$Q1SWBeam[,,,,, -1] <- uncertainty$survey_catchability[[2]]
### format
df_q <- as.data.frame(q) %>%
  mutate(qname = factor(qname, levels = c("Q1SWBeam", "UK-FSP"))) %>%
  mutate(source = ifelse(iter == 1, "Median", "Simulation\nreplicates"))

p <- df_q %>%
  ggplot(aes(x = age, y = data, group = iter, 
             colour = source, linewidth = source, alpha = source)) +
  geom_line(data = . %>% filter(source == "Simulation\nreplicates")) +
  geom_line(data = . %>% filter(source == "Median")) +
  scale_colour_manual("", values = c("Median" = "red",
                                     "Simulation\nreplicates" = "black")) +
  scale_alpha_manual("", values = c("Median" = 1,
                                    "Simulation\nreplicates" = 0.05)) +
  scale_linewidth_manual("", values = c("Median" = 0.5,
                                        "Simulation\nreplicates" = 0.1)) +
  facet_wrap(~ qname, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  labs(x = "Age (years)", y = "Survey catchability (q)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"))
p
ggsave(filename = "output/plots/OM/OM_idx_q.png", plot = p, 
       width = 12, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_idx_q.pdf", plot = p, 
       width = 12, height = 5, units = "cm")


### index weights at age
df_idx_wts <- FLQuants(Q1SWBeam = FLQuant(iterMedians(catch.wt(idx$Q1SWBeam))),
                       "UK-FSP" = FLQuant(iterMedians(catch.wt(idx$`UK-FSP`))))
df_idx_wts <- as.data.frame(df_idx_wts) %>% 
  filter(year <= 2023) %>%
  mutate(age = factor(age, levels = 9:2))
cols <- scales::hue_pal()(length(2:9))
cols <- cols[c(seq(from = 1, to = length(cols), by = 3),
               seq(from = 2, to = length(cols), by = 3),
               seq(from = 3, to = length(cols), by = 3))]

p <- df_idx_wts %>%
  mutate(data = ifelse(data == 0, NA, data)) %>%
  ggplot(aes(x = year, y = data, colour = age)) +
  annotate("rect", xmin = 2018.5, xmax = 2023.5, ymin = -Inf, ymax = Inf,
           alpha = 0.05, fill = "red") + 
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.3) + 
  geom_vline(xintercept = 2018.5, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = 2023.5, colour = "red", linewidth = 0.3) +
  scale_colour_manual("Age (years)", values = cols) + 
  facet_wrap(~ qname) + 
  coord_cartesian(xlim = c(NA, 2023), ylim = c(0, NA)) +
  labs(x = "Year", y = "Index catch weight (kg)") + 
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p
ggsave(filename = "output/plots/OM/OM_idx_wts.png", plot = p, 
       width = 16, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_idx_wts.pdf", plot = p, 
       width = 16, height = 6, units = "cm")

### survey weights at age in projection
df_idx_wts_sim <- FLQuants(Q1SWBeam = FLQuant(iterMedians(catch.wt(idx$Q1SWBeam))),
                       "UK-FSP" = FLQuant(iterMedians(catch.wt(idx$`UK-FSP`))))
df_idx_wts_sim <- as.data.frame(df_idx_wts_sim) %>% 
  filter(year == 2024)
p_age <- df_idx_wts_sim %>%
  ggplot(aes(x = age, y = data, colour = qname, linetype = qname)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ "Projection") +
  scale_colour_brewer("", palette = "Dark2") + 
  scale_linetype("") +
  labs(x = "Age (years)", y = "Index catch weight (kg)") +
  scale_x_continuous(breaks = seq(2, 9, 2), limits = c(0, NA)) +
  ylim(c(0, 1.587459)) + 
  theme_bw(base_size = 8) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.3, 0.8),
        legend.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key = element_blank())
p_age

p_new <- p + p_age + plot_layout(widths = c(1, 0.5))
p_new
ggsave(filename = "output/plots/OM/OM_idx_wts_age.png", plot = p_new, 
       width = 16, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_idx_wts_age.pdf", plot = p_new, 
       width = 16, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### surveys - simulated historical values ####
### ------------------------------------------------------------------------ ###
idx <- readRDS("input/ple.27.7e/baseline/1000_100/idx.rds")
idx_dev <- readRDS("input/ple.27.7e/baseline/1000_100/idx_dev.rds")
idx_dev_raw <- readRDS("input/ple.27.7e/baseline/1000_100/idx_dev_raw.rds")

### observed numbers and biomass
idx_Q1 <- window(index(idx$Q1SWBeam) * idx_dev$Q1SWBeam, end = 2023)
idx_FSP <- window(index(idx$`UK-FSP`) * idx_dev$`UK-FSP`, end = 2023)
idxB_Q1 <- window(quantSums(index(idx$Q1SWBeam) * catch.wt(idx$Q1SWBeam) *
                              idx_dev$Q1SWBeam), end = 2023)
idxB_FSP <- window(quantSums(index(idx$`UK-FSP`) * catch.wt(idx$`UK-FSP`) *
                               idx_dev$`UK-FSP`), end = 2023)

### simulated values
idx_sim_Q1 <- window(index(idx$Q1SWBeam) * idx_dev_raw$Q1SWBeam, end = 2023)
idx_sim_FSP <- window(index(idx$`UK-FSP`) * idx_dev_raw$`UK-FSP`, end = 2023)
idxB_sim_Q1 <- window(quantSums(index(idx$Q1SWBeam) * catch.wt(idx$Q1SWBeam) *
                                  idx_dev_raw$Q1SWBeam), end = 2023)
idxB_sim_FSP <- window(quantSums(index(idx$`UK-FSP`) * catch.wt(idx$`UK-FSP`) *
                                   idx_dev_raw$`UK-FSP`), end = 2023)

qnts_sim <- FLQuants(Q1SWBeam = idx_sim_Q1,
                     `UK-FSP` = idx_sim_FSP)
qnts_sim_perc <- lapply(qnts_sim, quantile, 
                        probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    na.rm = TRUE)
qnts_sim_perc <- FLQuants(qnts_sim_perc)
df_sim_perc <- as.data.frame(qnts_sim_perc)
df_sim_perc <- df_sim_perc %>% 
  select(year, age, iter, data, survey = qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "simulated")
df_obs_perc <- FLQuants(Q1SWBeam = idx_Q1,
                        `UK-FSP` = idx_FSP) |>
  lapply(quantile, probs = c(0.5)) |>
  as(Class = "FLQuants") |>
  as.data.frame() %>%
  select(year, age, iter, data, survey = qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "observed")
df_perc <- bind_rows(df_sim_perc, df_obs_perc) %>%
  mutate(survey = factor(survey, levels = c("Q1SWBeam", "UK-FSP"),
                         labels = c("Q1SWBeam", "UK-FSP")),
         source = factor(source, levels = c("simulated", "observed")))

p_N_Q1 <- df_perc %>%
  filter(survey == "Q1SWBeam") %>%
  ggplot() +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `50%`, colour = source, linetype = source)) +
  facet_grid(paste0("Age ", age) ~ survey, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual("", values = c("simulated" = "black", 
                                     "observed" = "red")) +
  scale_linetype_manual("", values = c("simulated" = "solid", 
                                       "observed" = "2121")) +
  labs(x = "Year", y = "Index numbers") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.y = element_blank(),
        legend.position = "none")
p_N_FSP <- df_perc %>%
  filter(survey == "UK-FSP") %>%
  ggplot() +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `50%`, colour = source, linetype = source)) +
  facet_grid(paste0("Age ", age) ~ survey, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual("", values = c("simulated" = "black", 
                                     "observed" = "red")) +
  scale_linetype_manual("", values = c("simulated" = "solid", 
                                       "observed" = "2121")) +
  labs(x = "Year", y = "Index numbers") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_N_Q1 + p_N_FSP


### plot biomass
qntsB_sim <- FLQuants(Q1SWBeam = idxB_sim_Q1,
                     `UK-FSP` = idxB_sim_FSP)
qntsB_sim_perc <- lapply(qntsB_sim, quantile, 
                        probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                        na.rm = TRUE)
qntsB_sim_perc <- FLQuants(qntsB_sim_perc)
dfB_sim_perc <- as.data.frame(qntsB_sim_perc)
dfB_sim_perc <- dfB_sim_perc %>% 
  select(year, age, iter, data, survey = qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "simulated")
dfB_obs_perc <- FLQuants(Q1SWBeam = idxB_Q1,
                        `UK-FSP` = idxB_FSP) |>
  lapply(quantile, probs = c(0.5)) |>
  as(Class = "FLQuants") |>
  as.data.frame() %>%
  select(year, age, iter, data, survey = qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = "observed")
dfB_perc <- bind_rows(dfB_sim_perc, dfB_obs_perc) %>%
  mutate(survey = factor(survey, levels = c("Q1SWBeam", "UK-FSP"),
                         labels = c("Q1SWBeam", "UK-FSP")),
         source = factor(source, levels = c("simulated", "observed")))

p_B_Q1 <- dfB_perc %>%
  filter(survey == "Q1SWBeam") %>%
  ggplot() +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `50%`, colour = source, linetype = source)) +
  facet_grid("Biomass" ~ survey, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual("", values = c("simulated" = "black", 
                                     "observed" = "red")) +
  scale_linetype_manual("", values = c("simulated" = "solid", 
                                       "observed" = "2121")) +
  labs(x = "Year", y = "Index biomass") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        legend.position = "none")
p_B_FSP <- dfB_perc %>%
  filter(survey == "UK-FSP") %>%
  ggplot() +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = . %>%
                filter(source == "simulated"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `50%`, colour = source, linetype = source)) +
  facet_grid("Biomass" ~ survey, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual("", values = c("simulated" = "black", 
                                     "observed" = "red")) +
  scale_linetype_manual("", values = c("simulated" = "solid", 
                                       "observed" = "2121")) +
  labs(x = "Year", y = "Index numbers") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")
p <- p_N_Q1 + p_N_FSP + p_B_Q1 + p_B_FSP + 
  plot_layout(ncol = 2, heights = c(1, 0.3))
p
ggsave(filename = "output/plots/OM/OM_idx_sim.png", plot = p, 
       width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_idx_sim.pdf", plot = p, 
       width = 16, height = 16, units = "cm")


### ------------------------------------------------------------------------ ###
### MCMC ####
### ------------------------------------------------------------------------ ###
### use MCMC and compare output with OM
library(tmbstan)

stk_data <- readRDS("input/ple.27.7e/preparation/model_input_stk.RDS")
fit <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_fit.rds")
stk_fit <- SAM2FLStock(object = fit, stk = stk_data)
stk_om <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")

### create template stock for storing results
MCMC_iter <- 1000
MCMC_warmup <- 10000
MCMC_chains <- 10
#stk_MCMC <- propagate(stk_data, MCMC_iter)
stk_MCMC <- propagate(window(stk_om, end = 2023)[,,,,, 1], 10*1000)
harvest(stk_MCMC)[] <- NA
stock.n(stk_MCMC)[] <- NA

### run MCMC
system.time(mcmc <- tmbstan(fit$obj, chains = MCMC_chains, 
                            iter = MCMC_warmup + MCMC_iter,
                            warmup = MCMC_warmup,
                            seed = 1, control = list(max_treedepth = 15)))
saveRDS(mcmc, file = "input/ple.27.7e/baseline/1000_100/MCMC.rds")
mcmc <- readRDS("input/ple.27.7e/baseline/1000_100/MCMC.rds")



methods(class="stanfit")
pairs(mcmc)
traceplot(mcmc, pars = names(fit$obj$par), inc_warmup = TRUE)
traceplot(mcmc, pars = mcmc@model_pars, inc_warmup = TRUE)

traceplot(mcmc, pars = names(fit$obj$par), inc_warmup = TRUE, window = c(1, 2))
#pairs(mcmc, pars=names(obj$par))
#traceplot(mcmc, pars=names(obj$par), inc_warmup=TRUE)



### extract
mc <- extract(mcmc, inc_warmup = FALSE, permuted = TRUE)
idxF <- fit$conf$keyLogFsta[1, ] + 1

### in the MCMC results age and years are mixed within the same row,
### each row represents one iteration
### this needs to be reformatted to be useful...

### fishing mortality
harvest(stk_MCMC)[unique(idxF)] <- 
  aperm(array(data = c(exp(mc$logF)),
              dim = c(dim(stk_MCMC)[6], length(unique(idxF)), dim(stk_MCMC)[2],
                      1, 1, 1)),
        perm = c(2:6, 1))
harvest(stk_MCMC)[] <- harvest(stk_MCMC)[idxF] ### coupled ages
#plot(harvest(stk_MCMC))

### stock numbers at age
stock.n(stk_MCMC)[] <- aperm(array(data = c(exp(mc$logN)),
                                   dim = dim(stock.n(stk_MCMC))[c(6, 1:5)]),
                             perm = c(2:6, 1))
stock(stk_MCMC) <- computeStock(stk_MCMC)
#plot(stock.n(stk_MCMC))
#plot(window(stock.n(stk_om), end = 2023))

### plot MCMC stock
plot(stk_MCMC)
### compare with original SAM fit
plot(window(FLStocks(MCMC = stk_MCMC, original = stk_fit), end = 2023))
### compare with first uncertainty approach
plot(window(FLStocks(MCMC = stk_MCMC, VarCov = stk_om), end = 2023))

### plot SSB and F
### OM metrics
qnts <- FLQuants(OM_ssb = ssb(stk_om)/1000, OM_fbar = fbar(stk_om),
                 MCMC_ssb = ssb(stk_MCMC)/1000, MCMC_fbar = fbar(stk_MCMC))
qnts <- window(qnts, end = 2023)
### percentiles
qnts_perc <- lapply(qnts, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    na.rm = TRUE)
qnts_perc <- FLQuants(qnts_perc)
df_perc <- as.data.frame(qnts_perc) %>%
  select(year, iter, data, qname) %>%
  separate_wider_delim(qname, delim = "_", names = c("source", "quant")) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(source = factor(source, levels = c("OM", "MCMC"),
                         labels = c("Operating model", "MCMC"))) %>%
  mutate(quant = factor(quant,
                        levels = c("ssb", "fbar"),
                        labels = c("SSB (1000t)", "F (ages 3-6)")))
### plot
p <- df_perc %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = source), 
              alpha = 0.15, linewidth = 0,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = source), 
              alpha = 0.15, linewidth = 0,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `2.5%`, colour = source, linetype = source),
            linewidth = 0.2, alpha = 0.15) +
  geom_line(aes(x = year, y = `97.5%`, colour = source, linetype = source),
            linewidth = 0.2, alpha = 0.15) +
  geom_line(aes(x = year, y = `25%`, colour = source, linetype = source),
            linewidth = 0.3, alpha = 0.3) +
  geom_line(aes(x = year, y = `75%`, colour = source, linetype = source),
            linewidth = 0.3, alpha = 0.3) +
  geom_line(aes(x = year, y = `50%`, colour = source, linetype = source)) +
  facet_wrap(~ quant, scales = "free_y", strip.position = "left", ncol = 2) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_colour_manual("", values = c("Operating model" = "red", 
                                     "MCMC" = "blue")) +
  scale_fill_manual("", values = c("Operating model" = "red", 
                                   "MCMC" = "blue")) +
  scale_linetype_manual("", values = c("Operating model" = "solid", 
                                       "MCMC" = "2121")) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),  
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.15), 
        legend.background = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/OM/OM_MCMC_vs_OM.png", plot = p, 
       width = 16, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_MCMC_vs_OM.pdf", plot = p, 
       width = 16, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### M scenarios ####
### ------------------------------------------------------------------------ ###
stk_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
M_baseline <- iterMeans(yearMeans(m(stk_baseline)))
M_Gislason <- readRDS("input/ple.27.7e/preparation/M_Gislason.rds")
M_Gislason <- M_baseline %=% c(yearMeans(M_Gislason)) ### formatting...

qnts <- FLQuants("Baseline" = M_baseline,
                 "M: +50%" = M_baseline * 1.5,
                 "M: -50%" = M_baseline * 0.5,
                 "M: Gislason" = M_Gislason)
df <- as.data.frame(qnts)

p <- df %>%
  ggplot() +
  geom_line(aes(x = age, y = data, colour = qname, linetype = qname),
            linewidth = 0.4) +
  geom_point(aes(x = age, y = data, colour = qname, shape = qname),
             size = 0.4) +
  labs(x = "Age (years)", y = "Natural mortality (M)") +
  scale_colour_discrete("") + 
  scale_linetype_discrete("") +
  scale_shape_discrete("") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"))
p
ggsave(filename = "output/plots/OM/OM_M.png", plot = p, 
       width = 8, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_M.pdf", plot = p, 
       width = 8, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### reference harvest rate (historical data) ####
### ------------------------------------------------------------------------ ###
### use baseline OM to illustrate principle
stk <- readRDS("input/ple.27.7e/baseline/1000_100/stk_oem.rds")
stk_OM <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.rds")

### UK-FSP
idx_FSP <- FLQuants(catch = iterMedians(catch(stk)),
                    index = quantSums(idx$`UK-FSP`@catch.wt * idx$`UK-FSP`@index),
                    rate = quantSums(idx$`UK-FSP`@catch.wt) %=% NA_real_)
idx_FSP <- window(idx_FSP, start = 2003, end = 2023)
idx_FSP$rate <- idx_FSP$catch/idx_FSP$index
### Q1SWBeam
idx_Q1 <- FLQuants(catch = iterMedians(catch(stk)),
                    index = quantSums(idx$Q1SWBeam@catch.wt * idx$Q1SWBeam@index),
                    rate = quantSums(idx$Q1SWBeam@catch.wt) %=% NA_real_)
idx_Q1 <- window(idx_Q1, start = 2006, end = 2023)
idx_Q1$rate <- idx_Q1$catch/idx_Q1$index

p_FSP <- idx_FSP |>
  as.data.frame() %>%
  select(year, data, qname) %>%
  mutate(data = ifelse(qname == "catch", data/1000, data)) %>%
  mutate(qname = factor(qname,
    levels = c("catch", "index", "rate"),
    labels = c("Catch (1000t)",
               "Biomass index\n(kg/hr m beam)",
               "Relative harvest rate\n(catch/biomass index)"))) %>%
  ggplot(aes(x = year, y = data)) +
  geom_line(linewidth = 0.4) +
  geom_hline(data = . %>%
               filter(qname == "Relative harvest rate\n(catch/biomass index)") %>%
               mutate(data = mean(data)),
             aes(yintercept = data), 
             colour = "red", linewidth = 0.3, linetype = "1111") +
  facet_grid(qname ~ "UK-FSP", scales = "free_y", switch = "y", ) +
  ylim(c(0, NA)) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

p_Q1 <- idx_Q1 |>
  as.data.frame() %>%
  select(year, data, qname) %>%
  mutate(data = ifelse(qname == "catch", data/1000, data)) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "index", "rate"),
                        labels = c("Catch (1000t)",
                                   "Biomass index\n(kg/km2)",
                                   "Relative harvest rate\n(catch/biomass index)"))) %>%
  ggplot(aes(x = year, y = data)) +
  geom_line(linewidth = 0.4) +
  geom_hline(data = . %>%
               filter(qname == "Relative harvest rate\n(catch/biomass index)") %>%
               mutate(data = mean(data)),
             aes(yintercept = data), 
             colour = "red", linewidth = 0.3, linetype = "1111") +
  facet_grid(qname ~ "Q1SWBeam", scales = "free_y", switch = "y", ) +
  ylim(c(0, NA)) +
  xlim(c(2003, NA)) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

p <- p_FSP + p_Q1
p
ggsave(filename = "output/plots/OM/OM_idx_hr.png", plot = p, 
       width = 16, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_idx_hr.pdf", plot = p, 
       width = 16, height = 12, units = "cm")

### ------------------------------------------------------------------------ ###
### reference set & robustness set - plot OM time series ####
### ------------------------------------------------------------------------ ###

### plot refset
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_refset_label <- c("Baseline", "Catch: no discards", 
                      "Catch: 100% discards", 
                      "Catch: no migration", 
                      "M: -50%", "M: +50%", "M: Gislason")

title_catch <- "Catch (1000t)"
title_ssb <- "SSB (1000t)"
title_fbar <- paste0("Mean F (ages 3-6)")

### historical stock(s)
stk_hist <- lapply(OMs_refset, function(y) {
  input_mp(OM = y, n_yrs = 20, MP = "hr")$om@stock
})
  
### get reference points
refpts <- lapply(OMs_refset, function(y) {
  input_refpts(OM = y)
})
  
stk_lst <- FLStocks(stk_hist)
qnts_list <- lapply(stk_lst, function(x) {
  qnts_i <- FLQuants(catch = window(catch(x)/1000, end = 2023),
                     fbar = window(fbar(x), end = 2023), 
                     ssb = window(ssb(x)/1000, end = 2024))
  return(qnts_i)
})
df_median <- lapply(seq(length(stk_hist)), function(x) {
  as.data.frame(FLQuants(lapply(qnts_list[[x]], iterMedians))) %>%
    mutate(stk = OMs_refset_label[[x]])
})
df_median <- do.call(bind_rows, df_median) %>%
  select(year, qname, data, stk) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "fbar", "ssb"),
                        labels = c(title_catch, title_fbar, title_ssb))) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))

qnts_combined <- FLQuants(
  catch = window(Reduce(FLCore::combine, lapply(stk_lst, catch))/1000, end = 2023),
  fbar = window(Reduce(FLCore::combine, lapply(stk_lst, fbar)), end = 2023),
  ssb = window(Reduce(FLCore::combine, lapply(stk_lst, ssb))/1000, end = 2024)
)
df_perc <- lapply(qnts_combined, quantile, 
                  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                  na.rm = TRUE)
df_perc <- FLQuants(df_perc)
df_perc <- as.data.frame(df_perc)
df_perc <- df_perc %>% select(year, iter, data, qname) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "fbar", "ssb"),
                        labels = c(title_catch, title_fbar, title_ssb)))

### distribution in last year
### by stk
df_distr_stk <- lapply(seq(length(stk_lst)), function(x) {
  as.data.frame(qnts_list[[x]]) %>%
    filter((year == 2023 & qname %in% c("catch", "fbar")) | 
             (year == 2024 & qname == "ssb")) %>%
    mutate(stk = OMs_refset_label[[x]])
})
df_distr_stk <- do.call(bind_rows, df_distr_stk) %>%
  select(qname, data, stk) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "fbar", "ssb"),
                        labels = c(title_catch, title_fbar, title_ssb))) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))
### combined
df_distr_combined <- as.data.frame(qnts_combined) %>%
  filter((year == 2023 & qname %in% c("catch", "fbar")) | 
           (year == 2024 & qname == "ssb")) %>%
  select(qname, data) %>%
  mutate(qname = factor(qname,
                        levels = c("catch", "fbar", "ssb"),
                        labels = c(title_catch, title_fbar, title_ssb)))

### reference points
val_Cmsy <- sapply(refpts, function(x) {
  c(iterMedians(x["Cmsy"]))/1000
})
df_Cmsy <- data.frame(qname = title_catch,
                      data = val_Cmsy,
                      stk = OMs_refset_label) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))
val_Bmsy <- sapply(refpts, function(x) {
  c(iterMedians(x["Bmsy"]))/1000
})
df_Bmsy <- data.frame(qname = "SSB (1000t)",
                      data = val_Bmsy,
                      stk = OMs_refset_label) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))
val_Blim <- sapply(refpts, function(x) {
  c(iterMedians(x["Blim"]))/1000
})
df_Blim <- data.frame(qname = "SSB (1000t)",
                      data = val_Blim,
                      stk = OMs_refset_label) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))
val_Fmsy <- sapply(refpts, function(x) {
  c(iterMedians(x["Fmsy"]))
})
df_Fmsy <- data.frame(qname = "Mean F (ages 3-6)",
                      data = val_Fmsy,
                      stk = OMs_refset_label) %>%
  mutate(stk = factor(stk, levels = OMs_refset_label))

### max value by qname
df_max <- bind_rows(df_perc %>%
                      group_by(qname) %>%
                      summarise(max = max(`95%`, na.rm = TRUE)) %>%
                      mutate(source = "perc"),
                    df_median %>%
                      group_by(qname) %>%
                      summarise(max = max(data, na.rm = TRUE)) %>%
                      mutate(source = "iter")) %>%
  group_by(qname) %>%
  summarise(max = max(max)) %>%
  mutate(max = max * 1.05)

max_catch <- df_max$max[df_max$qname == title_catch]
max_fbar <- df_max$max[df_max$qname == title_fbar]
max_ssb <- df_max$max[df_max$qname == title_ssb]

yr_start <- 1980
yr_end <- 2024

### catch plot
p_catch_distr <- ggplot() +
  geom_density(data = df_distr_stk %>% 
                 filter(qname == title_catch),
               aes(x = data, y = after_stat(scaled), colour = stk), 
               show.legend = FALSE, size = 0.2,
               outline.type = "full") +
  scale_colour_brewer(palette = "Dark2") +
  geom_density(data = df_distr_combined %>%
                 filter(qname == title_catch),
               aes(x = data, y = after_stat(scaled)), 
               show.legend = FALSE, size = 0.3,
               outline.type = "full") +
  theme_void() +
  coord_flip(xlim = c(0, max_catch), ylim = c(0, 1.05), expand = FALSE)
p_catch <- df_perc %>%
  filter(qname == title_catch) %>%
  ggplot(aes(x = year, y = `50%`)) +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_hline(data = df_Cmsy,
             aes(yintercept = data, colour = stk),
             linewidth = 0.2, linetype = "dashed",
             show.legend = FALSE) +
  geom_line(data = df_median %>%
              filter(qname == title_catch),
            aes(x = year, y = data, colour = stk),
            linewidth = 0.3, show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") + 
  geom_line() +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  labs(x = "Year") +
  coord_cartesian(ylim = c(0, max_catch), 
                  xlim = c(yr_start, yr_end), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(4, 10, 2, 4), "pt")) +
  annotation_custom(grob = ggplotGrob(p_catch_distr),
                    xmin = yr_end, xmax = yr_end + 2,
                    ymin = 0, ymax = max_catch)
#p_catch

### fbar plot
p_fbar_distr <- ggplot() +
  geom_density(data = df_distr_stk %>% 
                 filter(qname == title_fbar),
               aes(x = data, y = after_stat(scaled), colour = stk), 
               show.legend = FALSE, size = 0.2,
               outline.type = "full") +
  scale_colour_brewer(palette = "Dark2") +
  geom_density(data = df_distr_combined %>%
                 filter(qname == title_fbar),
               aes(x = data, y = after_stat(scaled)), 
               show.legend = FALSE, size = 0.3,
               outline.type = "full") +
  theme_void() +
  coord_flip(xlim = c(0, max_fbar), ylim = c(0, 1.05), expand = FALSE)
p_fbar <- df_perc %>%
  filter(qname == title_fbar) %>%
  ggplot(aes(x = year, y = `50%`)) +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_hline(data = df_Fmsy,
             aes(yintercept = data, colour = stk),
             linewidth = 0.2, linetype = "dashed",
             show.legend = FALSE) +
  geom_line(data = df_median %>%
              filter(qname == title_fbar),
            aes(x = year, y = data, colour = stk),
            linewidth = 0.3, show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") + 
  geom_line() +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  labs(x = "Year") +
  coord_cartesian(ylim = c(0, max_fbar), 
                  xlim = c(yr_start, yr_end), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 10, 2, 4), "pt")) +
  annotation_custom(grob = ggplotGrob(p_fbar_distr),
                    xmin = yr_end, xmax = yr_end + 2,
                    ymin = 0, ymax = max_fbar)
p_fbar

### ssb plot
p_ssb_distr <- ggplot() +
  geom_density(data = df_distr_stk %>% 
                 filter(qname == title_ssb),
               aes(x = data, y = after_stat(scaled), colour = stk), 
               show.legend = FALSE, size = 0.2,
               outline.type = "full") +
  scale_colour_brewer(palette = "Dark2") +
  geom_density(data = df_distr_combined %>%
                 filter(qname == title_ssb),
               aes(x = data, y = after_stat(scaled)), 
               show.legend = FALSE, size = 0.3,
               outline.type = "full") +
  theme_void() +
  coord_flip(xlim = c(0, max_ssb), ylim = c(0, 1.05), expand = FALSE)
p_ssb <- df_perc %>%
  filter(qname == title_ssb) %>%
  ggplot(aes(x = year, y = `50%`)) +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_hline(data = df_Bmsy,
             aes(yintercept = data, colour = stk),
             linewidth = 0.2, linetype = "dashed",
             show.legend = FALSE) +
  geom_hline(data = df_Blim,
             aes(yintercept = data, colour = stk),
             linewidth = 0.2, linetype = "dotted",
             show.legend = FALSE) +
  geom_line(data = df_median %>%
              filter(qname == title_ssb),
            aes(x = year, y = data, colour = stk),
            linewidth = 0.3, show.legend = TRUE) +
  scale_colour_brewer("Operating models", palette = "Dark2") + 
  geom_line() +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  labs(x = "Year") +
  coord_cartesian(ylim = c(0, max_ssb), 
                  xlim = c(yr_start, yr_end), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 7),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(2, 10, 4, 4), "pt"),
        legend.key.height = unit(0.4, "lines"),
        legend.key.spacing.y = unit(0.1, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title.position = "top") +
  guides(colour = guide_legend(nrow = 3)) +
  annotation_custom(grob = ggplotGrob(p_ssb_distr),
                    xmin = yr_end, xmax = yr_end + 2,
                    ymin = 0, ymax = max_ssb)
p_ssb

### combine all plots 
p <- p_catch + p_fbar + p_ssb + 
  plot_layout(ncol = 1)
p
ggsave(filename = "output/plots/paper/OM_refset.png", plot = p, 
       width = 8, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper/OM_refset.pdf", plot = p, 
       width = 8, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### recruitment - check +-20% scenarios ####
### ------------------------------------------------------------------------ ###
### SR model fit
fit <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_fit.rds")
sr_om <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
stk <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
stk <- window(iterMedians(stk), end = 2023)

fit_stk <- SAM2FLStock(object = fit, stk = stk)
sr <- as.FLSR(fit_stk, model = "bevholtSV")
rec(sr) <- rec(sr)/1000
ssb(sr) <- ssb(sr)/1000
sr <- fmle(sr, method = 'L-BFGS-B', fixed = list(), 
           control = list(trace = 0))
sr_params <- abPars("bevholt", s = params(sr)["s"], v = params(sr)["v"], 
                    spr0 = params(sr)["spr0"])

### try with 20% higher recruitment
fit_stk20 <- fit_stk
stock.n(fit_stk20)[1, ] <- stock.n(fit_stk20)[1, ] * 1.2
stock.n(fit_stk20)[1, ]/stock.n(fit_stk)[1, ]
sr20 <- as.FLSR(fit_stk20, model = "bevholtSV")
rec(sr20) <- rec(sr20)/1000
ssb(sr20) <- ssb(sr20)/1000
sr20 <- fmle(sr20, method = 'L-BFGS-B', fixed = list(), 
             control = list(trace = 0))
sr_params20 <- abPars("bevholt", s = params(sr20)["s"], v = params(sr20)["v"], 
                      spr0 = params(sr20)["spr0"])

params(sr)
params(sr20)


### ------------------------------------------------------------------------ ###
### TODO ####
### ------------------------------------------------------------------------ ###
### number of iterations

