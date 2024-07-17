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


### ------------------------------------------------------------------------ ###
### __TODO plot OM trajectories vs. ICES assessment - alternative OMs ####
### ------------------------------------------------------------------------ ###
### find alternative OMs
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_alt_OMs <- res_alt %>%
  select(stock, OM_group, OM) %>%
  unique() %>%
  mutate(
    OM_label = factor(OM, 
      levels = c("baseline", "M_low", "M_high", "M_Gislason",
                 "M_dd", "M_no_migration", "no_discards",
                 "rec_higher", "rec_no_AC", "rec_failure"),
      labels = c("baseline", "low", "high", "Gislason",
                 "dens. dep.", "no migration", "no discards",
                 "higher", "no AC", "failure")), .after = "OM") %>%
  mutate(
    OM_label2 = factor(OM_label, 
       levels = c("baseline", "low", "high", "Gislason",
                  "dens. dep.", "no migration", "no discards",
                  "higher", "no AC", "failure"),
       labels = c("baseline", 
                  "M: low", "M: high", "M: Gislason",
                  "M: dens. dep.", 
                  "M: no migration", "Catch: no discards",
                  "Rec: higher", "Rec: no AC", "Rec: failure")),
    .after = "OM_label") %>%
  mutate(OM_group = factor(OM_group, c("baseline", "M", "Catch", "Rec"))) %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))

### load reference points for all stocks and OMs
df_refpts <- foreach(i = seq(nrow(res_alt_OMs)), .combine = bind_rows) %do% {
  OM_i <- res_alt_OMs$OM[i]
  if (identical(OM_i, "rec_failure")) OM_i <- "baseline"
  stock_i <- res_alt_OMs$stock[i]
  refpts_i <- readRDS(paste0("input/", stock_i, "/",
                             OM_i, "/1000_100/refpts_mse.rds"))
  refpts_i <- iterMedians(refpts_i)
  res_alt_OMs[i, ] %>%
    mutate(Fmsy = c(refpts_i["Fmsy"]),
           Bmsy = c(refpts_i["Bmsy"]),
           Cmsy = c(refpts_i["Cmsy"]),
           Blim = c(refpts_i["Blim"])) %>%
    rename(stock_id = stock, stock = stock_label)
}

### get values from operating models
df_OM <- foreach(stock_id = res_alt_OMs$stock,
                 OM = res_alt_OMs$OM,
                 OM_label2 = res_alt_OMs$OM_label2,
                 stock = res_alt_OMs$stock_label,
                 .combine = bind_rows) %do% {
  #browser()
  stk_i <- readRDS(paste0("input/", stock_id, "/", 
                          ifelse(identical(OM, "rec_failure"),
                                 "baseline", OM), 
                          "/1000_100/stk.rds"))
  ### get metrics
  qnts <- FLQuants(ssb = ssb(stk_i)/1000, fbar = fbar(stk_i),
                   catch = catch(stk_i)/1000, rec = rec(stk_i)/1000)
  if (isTRUE(OM_label2 == "Catch: no discards"))
    qnts$catch <- landings(stk_i)/1000
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                       na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2,
           source = "OM")
  return(qnts_perc)
}
df_OM <- df_OM %>% filter(year <= 2020)
### get ICES assessment summary
df_ICES <- foreach(stock = c("Plaice", "Cod", "Herring"),
                   stock_id = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                   .combine = bind_rows) %do% {
  #browser()
  smry <- read.csv(paste0("input/", stock_id, "/preparation/",
                          "ices_assessment_summary.csv"))
  names(smry)[1] <- "year"
  smry <- smry %>% 
    select(year, SSB, F, Catches, Landings, Recruitment) %>%
    mutate(SSB = SSB/1000, Catches = Catches/1000, Landings = Landings/1000,
           Recruitment = Recruitment/1000) %>%
    rename(fbar = F, ssb = SSB, catch = Catches, landings = Landings,
           rec = Recruitment) %>%
    pivot_longer(c(ssb, fbar, catch, landings, rec), names_to = "qname") %>%
    mutate(stock = stock, stock_id = stock_id, source = "ICES")
  return(smry)
}
df_ICES <- df_ICES %>% filter(year <= 2020)
### add OM levels
df_ICES <- df_ICES %>%
  full_join(df_OM %>% 
              select(stock, stock_id, OM, OM_label2) %>% 
              unique())
# ### use landings for plaice scenario
# df_ICES$value[df_ICES$stock == "Plaice" & 
#                 df_ICES$OM_label2 == "Catch: no discards" &
#                 df_ICES$qname == "catch"] <- 
#   df_ICES$value[df_ICES$stock == "Plaice" & 
#                   df_ICES$OM_label2 == "Catch: no discards" &
#                   df_ICES$qname == "landings"]

df_combined <- 
  bind_rows(df_ICES %>%
              select(year, qname, value, source, stock, stock_id,
                     OM, OM_label2),
            df_OM %>% 
              select(year, qname, value = `50%`, stock, stock_id, source, 
                     OM, OM_label2)) %>%
  mutate(source = factor(source, levels = c("OM", "ICES"),
                         labels = c("Operating model",
                                    "ICES assessment")))


p_ple_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 3.2), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "rec"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "rec"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "rec"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 39), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [1000s]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.5, 0.98),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 23), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 0.9), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 3-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())

p_cod_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 610), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "rec"),
              aes(x = year, ymin = `5%`/1000, ymax = `95%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "rec"),
              aes(x = year, ymin = `25%`/1000, ymax = `75%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "rec"),
            aes(x = year, y = value/1000, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 3.4), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [millions]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 290), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 1.3), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 2-4)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())

p_her_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 1450), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "rec"),
              aes(x = year, ymin = `5%`/1000, ymax = `95%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "rec"),
              aes(x = year, ymin = `25%`/1000, ymax = `75%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "rec"),
            aes(x = year, y = value/1000, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 245), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [millions]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 11000), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 1.59), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 2-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())


p_ple <- p_ple_catch + 
  ggtitle(label = "(a) Plaice") + 
  theme(plot.title = element_text(face = "bold")) +
  p_ple_rec + p_ple_ssb + p_ple_fbar + 
  plot_layout(nrow = 1)
ggsave(filename = "output/plots/risk_OM_vs_ICES_ple.png", plot = p_ple, 
       width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES_ple.pdf", plot = p_ple, 
       width = 16, height = 16, units = "cm")

p_cod_her <- p_cod_catch + 
  ggtitle(label = "(b) Cod") + 
  theme(plot.title = element_text(face = "bold")) +
  p_cod_rec + p_cod_ssb + p_cod_fbar + 
  p_her_catch + 
  ggtitle(label = "(c) Herring") + 
  theme(plot.title = element_text(face = "bold")) +
  p_her_rec + p_her_ssb + p_her_fbar + 
  plot_layout(nrow = 2, heights = c(5, 5))
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.png", plot = p_cod_her, 
       width = 16, height = 22, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.pdf", plot = p_cod_her, 
       width = 16, height = 22, units = "cm")


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
### __TODO visualisation of OM MSY values ####
### ------------------------------------------------------------------------ ###

MSY_runs <- foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                    .combine = bind_rows) %:% 
  foreach(OM = c("baseline", "M_low", "M_high", "M_Gislason",
                 "no_discards_not_hidden",
                 "rec_no_AC", "rec_higher", "M_dd", "M_no_migration"),
          .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/", stock, "/", OM, "/1000_100/MSY_trace.rds")
  if (!file.exists(file_i)) return(NULL)
  runs_i <- readRDS(file_i)
  runs_i <- bind_rows(runs_i)
  runs_i$stock <- stock
  runs_i$OM <- OM
  return(runs_i)
}
MSY_runs <- MSY_runs %>%
  group_by(stock, OM) %>%
  mutate(MSY = ifelse(catch == max(catch), TRUE, FALSE))


MSY_runs_plot <- MSY_runs %>%
  pivot_longer(c(catch, ssb)) %>%
  mutate(value = value/1000) %>%
  mutate(label = factor(name, levels = c("catch", "ssb"),
                        labels = c("Catch [1000t]", "SSB [1000t]")),
         OM_label = factor(OM, 
                           levels = c("baseline", 
                                      "rec_no_AC", "rec_higher",
                                      "M_low", "M_high", "M_Gislason",
                                      "M_dd", "M_no_migration",
                                      "no_discards_not_hidden"),
                           labels = c("baseline",
                                      "R: no AC", "R: higher",
                                      "M: low", "M: high", "M: Gislason",
                                      "M: dens. dep.", "M: no migration",
                                      "Catch: no discards")))

p_MSY_ple <- MSY_runs_plot %>%
  filter(stock == "ple.27.7e") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  geom_smooth(span = 0.4,
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "ple.27.7e" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(a) Plaice") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))
p_MSY_cod <- MSY_runs_plot %>%
  filter(stock == "cod.27.47d20") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  stat_smooth(span = 0.2, 
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "cod.27.47d20" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(b) Cod") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))
p_MSY_her <- MSY_runs_plot %>%
  filter(stock == "her.27.3a47d") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  stat_smooth(span = 0.25, 
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "her.27.3a47d" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(c) Herring") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))

p_MSY <- p_MSY_ple + p_MSY_cod + p_MSY_her +
  plot_layout(design = "AA\nB#\nC#\n", widths = c(4, 2))
ggsave(filename = "output/plots/risk_MSY_all.png", plot = p_MSY, 
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_MSY_all.pdf", plot = p_MSY, 
       width = 17, height = 15, units = "cm")

### get Blim
Blims <- foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                 .combine = bind_rows) %:% 
  foreach(OM = c("baseline", "M_low", "M_high", "M_Gislason",
                 "no_discards_not_hidden",
                 "rec_no_AC", "rec_higher", "M_dd", "M_no_migration"),
          .combine = bind_rows) %do% {#browser()
    file_i <- paste0("input/", stock, "/", OM, "/1000_100/refpts_mse.rds")
    if (!file.exists(file_i)) return(NULL)
    i <- readRDS(file_i)
    return(data.frame(stock = stock, OM = OM ,Blim = median(c(i["Blim"]))))
}

### table
### MSY values
MSY_runs %>%
  filter(MSY == TRUE) %>%
  select(-tsb) %>% unique() %>%
  ### unfished values
  full_join(MSY_runs %>%
              filter(Ftrgt == 0) %>%
              select(ssb, rec, stock, OM) %>%
              rename(B0 = ssb, R0 = rec)) %>%
  ### Blim
  full_join(Blims) %>%
  select(stock, OM, B0, R0, Ftrgt, catch, ssb, rec, Blim) %>%
  rename(FMSY = Ftrgt, MSY = catch, BMSY = ssb, RMSY = rec) %>%
  write.csv(file = "input/OM_refpts.csv", row.names = FALSE)



### ------------------------------------------------------------------------ ###
### Recruitment model and residual visualisation (baseline OM) ####
### ------------------------------------------------------------------------ ###

### SR model fit
fit <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_fit.rds")
sr_om <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
stk <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
stk <- window(iterMedians(stk), end = 2023)

fit_stk <- SAM2FLStock(object = fit, stk = stk_data)
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
                 colour = "blue", size = 0.2, 
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

p <- p_res + p_hist + p_hist_iter + plot_annotation(tag_levels = "a", 
                                 tag_prefix = "(", tag_suffix = ")")  &
  theme(plot.tag = element_text(face = "bold"))
p
ggsave(filename = "output/plots/OM/OM_rec_res.png", plot = p,
       width = 16, height = 5, units = "cm", dpi = 600, type = "cairo")
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
  geom_point(data = df_sr %>% filter(type == "replicates (stock-\nrecruit pairs)"),
             aes(x = ssb/1000, y = rec/1000, 
                 colour = type), 
             alpha = 0.01, size = 0.3) +
  geom_line(data = sr_df %>% filter(type == "replicates (model)"),
            aes(x = ssb/1000, y = rec/1000, group = iter,
                colour = type),
            size = 0.1, alpha = 0.05) +
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
ggsave(filename = "output/plots/OM/OM_rec_baseline.png", plot = p,
       width = 10, height = 5, units = "cm", dpi = 600, type = "cairo")
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
res_new <- foreach(iter_i = seq(dim(sr)[6])) %do% {
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
ggsave(filename = "output/plots/OM/OM_rec_data_vs_model.pdf", plot = p,
       width = 16, height = 8, units = "cm", dpi = 600)
  
### ------------------------------------------------------------------------ ###
### Recruitment models of alternative OMs (all stocks) ####
### ------------------------------------------------------------------------ ###

### find alternative OMs
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_alt_OMs <- res_alt %>%
  select(stock, OM_group, OM) %>%
  unique() %>%
  mutate(
    OM_label = factor(OM, 
                      levels = c("baseline", "M_low", "M_high", "M_Gislason",
                                 "M_dd", "M_no_migration", "no_discards",
                                 "rec_higher", "rec_no_AC", "rec_failure"),
                      labels = c("baseline", "low", "high", "Gislason",
                                 "dens. dep.", "no migration", "no discards",
                                 "higher", "no AC", "failure")), .after = "OM") %>%
  mutate(
    OM_label2 = factor(OM_label, 
                       levels = c("baseline", "low", "high", "Gislason",
                                  "dens. dep.", "no migration", "no discards",
                                  "higher", "no AC", "failure"),
                       labels = c("baseline", 
                                  "M: low", "M: high", "M: Gislason",
                                  "M: dens. dep.", 
                                  "M: no migration", "Catch: no discards",
                                  "Rec: higher", "Rec: no AC", "Rec: failure")),
    .after = "OM_label") %>%
  mutate(OM_group = factor(OM_group, c("baseline", "M", "Catch", "Rec"))) %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))


### load recruitment values and stock-recruit pairs - use medians
df_OM <- foreach(stock_id = res_alt_OMs$stock,
                 OM = res_alt_OMs$OM,
                 OM_label2 = res_alt_OMs$OM_label2,
                 stock = res_alt_OMs$stock_label) %do% {
  #browser()
  sr_i <- readRDS(paste0("input/", stock_id, "/", 
                         ifelse(identical(OM, "rec_failure"),
                                "baseline", OM), 
                         "/1000_100/sr.rds"))
  ### SSB-recruitment pairs used for model fitting
  df_pairs <- data.frame(year = as.numeric(dimnames(ssb(sr_i))$year),
                         ssb = c(iterMedians(ssb(sr_i))), 
                         rec = c(iterMedians(rec(sr_i)))) %>%
    filter(!is.na(ssb)) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2)
  ### table with modelled recruitment
  sr_i_med <- iter(sr_i, 1)
  params(sr_i_med) <- iterMedians(params(sr_i))
  ssb_max <- max(c(iterMedians(ssb(sr_i))), na.rm = TRUE)
  
  if (isTRUE(model(FLSR(model = bevholt)) == model(sr_i))) {
    ssb_vals <- seq(from = 0, to = ssb_max * 3, length.out = 1500)
    rec_vals <- c(predict(sr_i_med, ssb = FLQuant(ssb_vals)))
  } else if (isTRUE(model(FLSR(model = segreg)) == model(sr_i))) {
    ssb_vals <- c(0, c(params(sr_i_med)$b), ssb_max * 3)
    rec_vals <- c(0, c(params(sr_i_med)$a * params(sr_i_med)$b),
                  c(params(sr_i_med)$a * params(sr_i_med)$b))
  } else {
    stop("unknown model")
  }
  df_pairs_model <- data.frame(ssb = ssb_vals,
                               rec = rec_vals) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2, rec_failure = FALSE)
  if (isTRUE(OM == "rec_failure")) {
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
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2)
  return(list(pairs = df_pairs, pairs_add = df_pairs_model, pars = pars_i))
}
df_pairs <- bind_rows(lapply(df_OM, "[[", 1))
df_pairs_add <- bind_rows(lapply(df_OM, "[[", 2))
df_pars <- bind_rows(lapply(df_OM, "[[", 3))


p_ple <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Plaice"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Plaice"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [1000t]",
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous("Recruitment [1000s]",
                     breaks = c(0, 10000, 20000, 30000),
                     labels = c(0, 10, 20, 30)) +
  coord_cartesian(xlim = c(0, 9500), ylim = c(0, 35000), expand = FALSE) +
  labs(title = "(a) Plaice") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p_cod <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Cod"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Cod"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [1000t]",
                     breaks = c(0, 25000, 50000, 75000, 100000),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous("Recruitment [millions]",
                     breaks = c(0, 500000, 1000000, 1500000),
                     labels = c(0, "0.5", "1.0", "1.5")) +
  coord_cartesian(xlim = c(0, 125000), ylim = c(0, 1500000), expand = FALSE) +
  labs(title = "(b) Cod") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p_her <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Herring"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Herring"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [million t]",
                     breaks = c(0, 2000000, 4000000, 6000000),
                     labels = c(0, 2, 4, 6)) +
  scale_y_continuous("Recruitment [millions]",
                     breaks = c(0, 50000000, 100000000, 150000000),
                     labels = c(0, 50, 100, 150)) +
  coord_cartesian(xlim = c(0, 6200000), ylim = c(0, 180000000), expand = FALSE) +
  labs(title = "(c) Herring") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p <- (p_ple + p_cod + p_her + plot_layout(heights = c(2.4, 1, 1)))
p
ggsave(filename = "output/plots/OM/OM_rec_OMs.png", plot = p,
       width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_OMs.pdf", plot = p,
       width = 16, height = 16, units = "cm", dpi = 600)

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
### TODO ####
### ------------------------------------------------------------------------ ###
### number of iterations

