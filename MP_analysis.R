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
### check Blim risk and projection length - baseline ####
### ------------------------------------------------------------------------ ###
### baseline OM
### use optimised multiplier
res_mult <- readRDS("output/ple.27.7e/baseline/1000_20/multiplier/hr/mp_1_1_1_1.4_1_0.63_1.2_0.7.rds")
refpts_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
Blim <- c(refpts_baseline["Blim"])[1] ### identical for all iterations

### start parallel workers
plan(multisession, workers = 15)

### number of iterations

### risk3: max risk (of last 10 years) -> year 2035
### cumulative risk over iterations
risk3_cum <- foreach(its = 1:1000, .combine = bind_rows) %dofuture% {
  data.frame(iter = its,
             risk3 = c(iterMeans(iter(ssb(res_mult)[, ac(2035)], seq(its)) < 
                                   Blim)))
}
risk3_cum %>%
  ggplot(aes(x = iter, y = risk3)) +
  geom_line() +
  theme_bw(base_size = 8)

### risk1: average risk
### cumulative risk over iterations
risk1_cum <- foreach(its = 1:1000, .combine = bind_rows) %dofuture% {
  data.frame(iter = its,
             risk1 = mean(iter(ssb(res_mult)[, ac(2035:2044)], seq(its)) < 
                            Blim))
}
risk1_cum %>%
  ggplot(aes(x = iter, y = risk1)) +
  geom_line() +
  theme_bw(base_size = 8)

df_risk_cum <- full_join(risk1_cum, risk3_cum)
df_risk_cum %>%
  pivot_longer(-1, names_to = "type", values_to = "risk") %>%
  mutate(type = factor(type, levels = c("risk1", "risk3"),
                       labels = c("Risk 1 (mean)", "Risk 3 (max)"))) %>%
  ggplot(aes(x = iter, y = risk, linetype = type)) +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "1111",
             linewidth = 0.3) +
  scale_linetype("") +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
                     limits = c(0, 0.052)) +
  labs(x = "# simulation replicates", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"))

### same for 10k iterations
res_mult_10k <- readRDS("output/ple.27.7e/baseline/10000_20/multiplier/hr/mp_1_1_1_1.4_1_0.63_1.2_0.7.rds")
### use first 1000 iterations from baseline run above
stk10k <- stock(om(res_mult_10k))
iter(stk10k, 1:1000) <- stock(om(res_mult))

### risk3: max risk (of last 10 years) -> year 2035
### cumulative risk over iterations
risk3_cum_10k <- foreach(its = 1:10000, .combine = bind_rows) %dofuture% {
  data.frame(iter = its,
             risk3 = c(iterMeans(iter(ssb(stk10k)[, ac(2035)], seq(its)) < 
                                   Blim)))
}
risk3_cum_10k %>%
  ggplot(aes(x = iter, y = risk3)) +
  geom_line() +
  theme_bw(base_size = 8)

### risk1: average risk
### cumulative risk over iterations
risk1_cum_10k <- foreach(its = 1:10000, .combine = bind_rows) %dofuture% {
  data.frame(iter = its,
             risk1 = mean(iter(ssb(stk10k)[, ac(2035:2044)], seq(its)) < 
                            Blim))
}

### plot
df_risk_cum_10k <- full_join(risk1_cum_10k, risk3_cum_10k)
p <- df_risk_cum_10k %>%
  pivot_longer(-1, names_to = "type", values_to = "risk") %>%
  mutate(type = factor(type, levels = c("risk1", "risk3"),
                       labels = c("Risk 1 (mean)", "Risk 3 (max)"))) %>%
  ggplot(aes(x = iter, y = risk, linetype = type)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "1111",
             linewidth = 0.3) +
  scale_linetype_manual("", values = c("solid", "1111")) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
                     limits = c(0, 0.052)) +
  labs(x = "# simulation replicates", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(4, 8, 4, 4), "pt")) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0, 10200), expand = FALSE) 
p
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_iter_10k.png", plot = p,
       width = 10, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_iter_10k.pdf", plot = p,
       width = 10, height = 6, units = "cm")
p <- p + coord_cartesian(xlim = c(0, 1000), expand = FALSE)
p
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_iter_1k.png", plot = p,
       width = 10, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_iter_1k.pdf", plot = p,
       width = 10, height = 6, units = "cm")

### projection period
### risk over time
iterMeans(ssb(res_mult) < c(refpts_baseline["Blim"]))

### get 100-year projection
res_mult_100 <- readRDS("output/ple.27.7e/baseline/1000_100/multiplier/hr/mp_1_1_1_1.4_1_0.63_1.2_0.7.rds")
### add fishing history
stk_hist <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
stk <- stk_hist
stk[, ac(2024:2124)] <- stock(om(res_mult_100))

### SSB trajectory
### get metrics
df_ssb_qnt <- quantile(ssb(stk)/1000, 
                       c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
df_ssb_qnt <- as.data.frame(df_ssb_qnt) %>%
  select(year, iter, data) %>%
  pivot_wider(names_from = iter, values_from = data)
df_ssb_iter <- as.data.frame(iter(ssb(stk)/1000, 1:5))
Bmsy <- c(median(refpts_baseline["Bmsy"], na.rm = TRUE))
Blim <- c(median(refpts_baseline["Blim"], na.rm = TRUE))
p_ssb <- ggplot() +
  geom_vline(xintercept = 2024.5, colour = "grey", size = 0.5) +
  geom_ribbon(data = df_ssb_qnt ,
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(data = df_ssb_qnt,
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_line(data = df_ssb_iter,
            aes(x = year, y = data, colour = iter),
            linewidth = 0.1, show.legend = FALSE) +
  geom_line(data = df_ssb_qnt,
            aes(x = year, y = `50%`), linewidth = 0.4) +
  # geom_hline(yintercept = Bmsy/1000,
  #            colour = "black", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = Blim/1000,
             colour = "black", linewidth = 0.5, linetype = "1111") +
  labs(x = "Year", y = "SSB (1000t)") +
  coord_cartesian(xlim = c(2020, 2124), ylim = c(0, 15), expand = FALSE) + 
  theme_bw(base_size = 8)
p_ssb
p_ssb_20 <- p_ssb +
  coord_cartesian(xlim = c(2020, 2044), ylim = c(0, 15), expand = FALSE)
p_ssb_20

### risk trajectory
risk <- iterMeans(ssb(stk) < Blim)
df_risk <- as.data.frame(risk)
risk_max <- max(risk[, ac(2035:2044)])
risk_year <- 2035
p_risk <- df_risk %>%
  filter(year >= 2024) %>%
  ggplot(aes(x = year, y = data)) +
  geom_vline(xintercept = 2024.5, colour = "grey", size = 0.5) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed") + 
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  annotate("rect", xmin = 2034.5, xmax = 2044.5, ymin = -Inf, ymax = Inf,
           alpha = 0.05, fill = "red") + 
  geom_hline(yintercept = risk_max, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = risk_year, colour = "red", linewidth = 0.3) +
  coord_cartesian(xlim = c(2020, 2124), ylim = c(0, 0.055), expand = FALSE) 
p_risk

p_risk_20 <- p_risk +
  coord_cartesian(xlim = c(2020, 2044), ylim = c(0, 0.055), expand = FALSE) 
p_risk_20

p1 <- p_ssb_20 + 
  annotate(geom = "text", x = 2021, y = Blim/1000*0.95, label = "B[lim]",
           parse = TRUE,
           colour = "black", size = 2.5, vjust = 1, hjust = 0) + 
  facet_wrap(~ "20-year projection") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2 <- p_ssb + facet_wrap(~ "100-year projection") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p3 <- p_risk_20 + 
  annotate(geom = "text", x = 2021, y = 0.05*0.98, label = "5% risk threshold",
           colour = "red", size = 2.5, vjust = 1, hjust = 0)
p4 <- p_risk + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
p
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_yrs.png", plot = p,
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr0.63_Blim_yrs.pdf", plot = p,
       width = 16, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### check Blim risk - reference set ####
### ------------------------------------------------------------------------ ###
### same for reference set OM
res <- readRDS("output/ple.27.7e/refset/1000_20/multiplier/hr/multiplier-upper_constraint1.2-lower_constraint0.7--obj_ICES_res_11-20.rds")
res@solution[, "multiplier"]
mp_mult_refset <- readRDS("output/ple.27.7e/refset/1000_20/multiplier/hr/mp_1_1_1_1.4_1_0.61_1.2_0.7.rds")
stk_refset <- stock(om(mp_mult_refset))
refpts_refset <- input_refpts(OM = "refset")
Blim <- c(refpts_refset["Blim"]) ### differ by iterations/OM
n_iter <- 1000
n_OMs <- dim(stk_refset)[6]/n_iter
risk_refset <- iterMeans(ssb(stk_refset) < 
                           rep(Blim, each = dim(ssb(stk_refset))[2]))
which.max(risk_refset[, ac(2035:2044)])
risk_max_year <- 2035

### risk3: max risk (of last 10 years)
### cumulative risk over iterations
risk3_cum <- foreach(its = seq(n_iter), .combine = bind_rows) %dofuture% {#browser()
  ### group OMs, e.g. iteration 1 for all OMs
  its_i <- sort(unlist(lapply(seq(its), function(x) {x + (seq(n_OMs) - 1) * n_iter})))
  data.frame(iter = its,
             risk3 = c(iterMeans(iter(ssb(stk_refset)[, ac(2035)], seq(its_i)) < 
                                   Blim[its_i])))
}
risk3_cum %>%
  ggplot(aes(x = iter, y = risk3)) +
  geom_line() +
  theme_bw(base_size = 8)

### risk1: average risk
### cumulative risk over iterations
risk1_cum <- foreach(its = seq(n_iter), .combine = bind_rows) %dofuture% {
  its_i <- sort(unlist(lapply(seq(its), function(x) {x + (seq(n_OMs) - 1) * n_iter})))
  data.frame(iter = its,
             risk1 = mean(iter(ssb(stk_refset)[, ac(2035:2044)], seq(its_i)) < 
                            rep(Blim[its_i], each = 10)))
}

### plot
df_risk_cum <- full_join(risk1_cum, risk3_cum)
p <- df_risk_cum %>%
  pivot_longer(-1, names_to = "type", values_to = "risk") %>%
  mutate(type = factor(type, levels = c("risk1", "risk3"),
                       labels = c("Risk 1 (mean)", "Risk 3 (max)"))) %>%
  ggplot(aes(x = iter, y = risk, linetype = type)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "1111",
             linewidth = 0.3) +
  scale_linetype("") +
  # scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
  #                    limits = c(0, 0.052)) +
  labs(x = "# simulation replicates", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(4, 8, 4, 4), "pt")) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0, 1000), expand = FALSE) 
p
ggsave(filename = "output/plots/MP/refset_hr0.63_Blim_iter.png", plot = p,
       width = 10, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/refset_hr0.63_Blim_iter.pdf", plot = p,
       width = 10, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### baseline - multiplier ####
### ------------------------------------------------------------------------ ###

### load GA results
res_mult <- readRDS("output/ple.27.7e/baseline/1000_20/multiplier/hr/multiplier-upper_constraint1.2-lower_constraint0.7--obj_ICES_res_11-20.rds")
runs_mult <- readRDS("output/ple.27.7e/baseline/1000_20/multiplier/hr/runs.rds")
### optimum multiplier
as.data.frame(as.list(res_mult@solution[1, ]))

### combine performance statistics of all runs
df_mult <- lapply(runs_mult, function(x) {
  bind_cols(data.frame(t(x$pars)), data.frame(x$stats))
})
df_mult <- do.call(bind_rows, df_mult)
### add fitness value
df_mult$fitness <- df_mult$X11.20_Catch_rel -
  penalty(x = df_mult$X11.20_risk_Blim_max, 
          negative = FALSE, max = 1, 
          inflection = 0.06, 
          steepness = 1000)

df_mult %>%
  select(multiplier, X11.20_risk_Blim_max, X11.20_Catch_rel, fitness) %>%
  View()

### plot
df_plot <- df_mult %>%
  select(multiplier, Blim_risk = X11.20_risk_Blim_max, 
         Catch_rel = X11.20_Catch_rel, SSB_rel = X11.20_SSB_rel) %>%
  pivot_longer(-1) %>%
  mutate(name = factor(name, 
                       levels = c("SSB_rel", "Catch_rel", "Blim_risk"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
df_Blim <- data.frame(name = "B[lim]~risk",
                      value = 0.05) %>%
  mutate(name = factor(name, 
                       levels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
df_label <- data.frame(name = c("B[lim]~risk", "Catch/MSY"),
                       x = c(1.2, 0.65), y = c(0.1, 0.5),
                       label = c("5% risk threshold", "maximum catch")) %>%
  mutate(name = factor(name, 
                       levels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
p <- df_plot %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_hline(data = df_Blim,
             aes(yintercept = value), colour = "red", linewidth = 0.3) +
  geom_point(size = 0.1, shape = 16) +
  geom_smooth(linewidth = 0.4, n = 100, span = 0.2, se = FALSE, 
              colour = "black") +
  geom_vline(xintercept = 0.63, colour = "red", linewidth = 0.3,
             linetype = "1111") +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             ncol = 1, labeller = label_parsed) +
  geom_text(data = df_label,
            aes(x = x, y = y, label = label),
            colour = "red", size = 2.5, hjust = 0) +
  labs(x = "Multiplier") +
  theme_bw(base_size = 8) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/MP/baseline_hr_multiplier.png", plot = p,
       width = 8, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_multiplier.pdf", plot = p,
       width = 8, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### refset - multiplier ####
### ------------------------------------------------------------------------ ###

### load results
df_runs <- foreach(index = c("UK-FSP", "Q1SWBeam"), .combine = bind_rows) %do% {
  path <- paste0("output/ple.27.7e/refset/1000_20/",
                 ifelse(identical(index, "Q1SWBeam"),
                        "multiplier_Q1SWBeam", "multiplier"),
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
  mutate(group = index) %>%
  mutate(group = factor(group, 
                        levels = c("UK-FSP", "Q1SWBeam")))

# df_runs %>%
#   select(group, multiplier, X11.20_risk_Blim_max, X11.20_Catch_rel) %>%
#   View()

### find optima
df_optima <- df_runs %>%
  group_by(group) %>%
  filter(X11.20_risk_Blim_max <= 0.05) %>%
  filter(X11.20_Catch_rel == max(X11.20_Catch_rel))

### save results
saveRDS(df_runs, file = "output/refset_x_runs.rds")
saveRDS(df_optima, file = "output/refset_x_runs_opt.rds")
write.csv(df_optima, "output/refset_x_runs_opt.csv", row.names = FALSE)

### plot
df_runs_plot <- df_runs %>%
  select(group, multiplier, Blim_risk = X11.20_risk_Blim_max, 
         Catch_rel = X11.20_Catch_rel, SSB_rel = X11.20_SSB_rel) %>%
  pivot_longer(-1:-2) %>%
  mutate(name = factor(name, 
                       levels = c("SSB_rel", "Catch_rel", "Blim_risk"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
df_Blim <- data.frame(name = "B[lim]~risk",
                      value = 0.05) %>%
  mutate(name = factor(name,
                       levels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))

df_runs_plot_catch <- bind_rows(
  df_runs %>%
    select(group, multiplier, risk = X11.20_risk_Blim_max,
           catch = X11.20_Catch_rel) %>%
    group_by(group) %>%
    filter(catch == max(catch)) %>%
    mutate(label = "absolute"),
  df_runs %>%
    select(group, multiplier, risk = X11.20_risk_Blim_max,
           catch = X11.20_Catch_rel) %>%
    group_by(group) %>%
    filter(risk <= 0.05) %>%
    filter(catch == max(catch)) %>%
    mutate(label = "within risk limit")
)

p_x <- df_runs_plot %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_hline(data = df_Blim,
             aes(yintercept = value), colour = "red", linewidth = 0.3) +
  geom_point(size = 0.1, shape = 16) +
  geom_smooth(linewidth = 0.4, n = 100, span = 0.2, se = FALSE, 
              colour = "black") +
  geom_vline(data = df_runs_plot_catch,
             aes(xintercept = multiplier, colour = label),
             linewidth = 0.4, linetype = "1111") +
  scale_colour_discrete("Catch maximum") +
  facet_grid(name ~ group, scales = "free_y", labeller = label_parsed, 
             switch = "y") +
  labs(x = "Multiplier (x)") +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.key = element_blank())
p_x
ggsave(filename = "output/plots/MP/refset_x.png", plot = p_x,
       width = 16, height = 7, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/refset_x.pdf", plot = p_x,
       width = 16, height = 7, units = "cm")

### ------------------------------------------------------------------------ ###
### baseline explorations - x, w, n1, v ####
### ------------------------------------------------------------------------ ###
### x - multiplier 
### w - Itrigger
### n1 - number of years in index
### v - advice interval

### load all runs
if (FALSE) {
  runs2 <- readRDS("output/ple.27.7e/baseline/1000_20/mult_comp_b_mult_int2/hr/runs.rds")
  runs1 <- readRDS("output/ple.27.7e/baseline/1000_20/mult_comp_b_mult/hr/runs.rds")
  runs <- c(runs1, runs2)
  length(runs)
  length(unique(runs))
  saveRDS(runs, "output/ple.27.7e/baseline/1000_20/mult_comp_b_mult/hr/runs.rds")
  runs <- readRDS("output/ple.27.7e/baseline/1000_20/mult_comp_b_mult/hr/runs.rds")
  
  ### combine performance statistics of all runs
  df_runs <- lapply(runs, function(x) {
    bind_cols(data.frame(t(x$pars)), data.frame(x$stats))
  })
  df_runs <- do.call(bind_rows, df_runs)
  
  ### add fitness value
  df_runs$fitness <- df_runs$X11.20_Catch_rel -
    penalty(x = df_runs$X11.20_risk_Blim_max, 
            negative = FALSE, max = 1, 
            inflection = 0.06, 
            steepness = 1000)
  max(df_runs$fitness)
  df_runs[which.max(df_runs$fitness), ]
  ### strict 5% risk limit - identical results
  df_runs$fitness2 <- df_runs$X11.20_Catch_rel -
    ifelse(df_runs$X11.20_risk_Blim_max <= 0.05, 0, 1)
  max(df_runs$fitness2)
  df_runs[which.max(df_runs$fitness2), ]
  
  saveRDS(df_runs, file = "output/baseline_results_x_w_v_n1.rds") 
}
df_runs <- readRDS("output/baseline_results_x_w_v_n1.rds")

# df_runs %>%
#   select(multiplier, X11.20_risk_Blim_max, X11.20_Catch_rel, fitness) %>%
#   View()

### multiplier
df_mult <- df_runs %>%
  filter(idxB_range_3 == 1 & 
           interval == 1 &
           comp_b_multiplier == 1.4)
df_mult[which.max(df_mult$fitness), ]

### multiplier & Itrigger
df_mult_w <- df_runs %>%
  filter(idxB_range_3 == 1 & 
           interval == 1)
df_mult_w_best <- df_mult_w[which.max(df_mult_w$fitness), ]

### multiplier & Itrigger & interval
df_mult_w_v <- df_runs %>%
  filter(idxB_range_3 == 1)
df_mult_w_v[which.max(df_mult_w_v$fitness), ]

### plot multiplier
df_plot <- df_mult %>%
  select(multiplier, Blim_risk = X11.20_risk_Blim_max, 
         Catch_rel = X11.20_Catch_rel, SSB_rel = X11.20_SSB_rel) %>%
  pivot_longer(-1) %>%
  mutate(name = factor(name, 
                       levels = c("SSB_rel", "Catch_rel", "Blim_risk"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
df_Blim <- data.frame(name = "B[lim]~risk",
                      value = 0.05) %>%
  mutate(name = factor(name, 
                       levels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
df_label <- data.frame(name = c("B[lim]~risk", "Catch/MSY"),
                       x = c(1.2, 0.65), y = c(0.1, 0.5),
                       label = c("5% risk limit", "maximum catch")) %>%
  mutate(name = factor(name, 
                       levels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")))
p <- df_plot %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_hline(data = df_Blim,
             aes(yintercept = value), colour = "red", linewidth = 0.3) +
  geom_point(size = 0.1, shape = 16) +
  geom_smooth(linewidth = 0.4, n = 100, span = 0.2, se = FALSE, 
              colour = "black") +
  geom_vline(xintercept = 0.63, colour = "red", linewidth = 0.3,
             linetype = "1111") +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", 
             ncol = 1, labeller = label_parsed) +
  geom_text(data = df_label,
            aes(x = x, y = y, label = label),
            colour = "red", size = 2.5, hjust = 0) +
  labs(x = "Multiplier") +
  theme_bw(base_size = 8) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        axis.title.y = element_blank())
p
ggsave(filename = "output/plots/MP/baseline_hr_multiplier.png", plot = p,
       width = 8, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_multiplier.pdf", plot = p,
       width = 8, height = 6, units = "cm")

### plot - multiplier & Itrigger
### first: full area
p <- df_mult_w %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(x = multiplier, y = comp_b_multiplier, label = X11.20_Catch_rel,
             fill = catch)) +
  geom_raster(alpha = 0.8) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       #colours = rainbow(10),
                       #colours = topo.colors(10),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Multiplier (x)", y = expression(I[trigger]~multiplier~"(w)")) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger.png", plot = p,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger.pdf", plot = p,
       width = 16, height = 10, units = "cm")

### add optimum
p <- p +
  geom_hline(data = df_mult_w_best,
             aes(yintercept = comp_b_multiplier),
             linewidth = 0.3, linetype = "1111") +
  geom_vline(data = df_mult_w_best,
             aes(xintercept = multiplier),
             linewidth = 0.3, linetype = "1111")
p
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_opt.png", plot = p,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_opt.pdf", plot = p,
       width = 16, height = 10, units = "cm")
### remove cells with risk > 5%
p <- df_mult_w %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(x = multiplier, y = comp_b_multiplier, label = X11.20_Catch_rel,
             fill = catch_col)) +
  geom_raster(alpha = 0.8) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_mult_w_best,
             aes(yintercept = comp_b_multiplier),
             linewidth = 0.3, linetype = "1111") +
  geom_vline(data = df_mult_w_best,
             aes(xintercept = multiplier),
             linewidth = 0.3, linetype = "1111") +
  labs(x = "Multiplier (x)", y = expression(I[trigger]~multiplier~"(w)")) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_risk.png", plot = p,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_risk.pdf", plot = p,
       width = 16, height = 10, units = "cm")

### plot - multiplier & Itrigger & interval & index years
df_runs %>%
  group_by(idxB_range_3, interval) %>%
  filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
  select(idxB_range_3, interval, comp_b_multiplier, multiplier, 
         X11.20_Catch_rel, X11.20_SSB_rel, X11.20_ICV, X11.20_risk_Blim_max) %>%
  View()
### -> several identical solutions
### filter down by SSB (max), risk (min)
df_runs_opt <- df_runs %>%
  group_by(idxB_range_3, interval) %>%
  filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
  filter(X11.20_SSB_rel == max(X11.20_SSB_rel)) %>%
  filter(X11.20_risk_Blim_max == min(X11.20_risk_Blim_max)) %>%
  select(idxB_range_3, interval, comp_b_multiplier, multiplier, 
         X11.20_Catch_rel, X11.20_SSB_rel, X11.20_ICV, X11.20_risk_Blim_max)
View(df_runs_opt)
p <- df_runs %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(x = multiplier, y = comp_b_multiplier, label = X11.20_Catch_rel,
             fill = catch_col)) +
  geom_raster(alpha = 0.8) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_runs_opt,
             aes(yintercept = comp_b_multiplier),
             linewidth = 0.3, linetype = "1111") +
  geom_vline(data = df_runs_opt,
             aes(xintercept = multiplier),
             linewidth = 0.3, linetype = "1111") +
  facet_grid(paste0("Index years (n0): ", idxB_range_3) ~ 
               paste0("Advice interval (v): ", interval)) + 
  labs(x = "Multiplier (x)", y = expression(I[trigger]~multiplier~"(w)")) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(panel.spacing = unit(0.7, "lines"))
p
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_v_n0.png", plot = p,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/baseline_hr_mult_trigger_v_n0.pdf", plot = p,
       width = 16, height = 10, units = "cm")

### summarise results
df_runs %>%
  group_by(idxB_range_3, interval) %>%
  filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
  filter(X11.20_SSB_rel == max(X11.20_SSB_rel)) %>%
  filter(X11.20_risk_Blim_max == min(X11.20_risk_Blim_max)) %>%
  filter(comp_b_multiplier == max(comp_b_multiplier)) %>%
  select(idxB_range_3, interval, comp_b_multiplier, multiplier, 
         X11.20_Catch_rel, X11.20_SSB_rel, X11.20_ICV, X11.20_risk_Blim_max,
         fitness) %>%
  View()

df_baseline <- df_runs %>%
  filter(comp_b_multiplier == 1.4, interval == 1, idxB_range_3 == 1) %>%
  filter(fitness == max(fitness)) %>%
  bind_rows(
    df_runs %>%
      group_by(idxB_range_3, interval) %>%
      filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
      filter(X11.20_SSB_rel == max(X11.20_SSB_rel)) %>%
      filter(X11.20_risk_Blim_max == min(X11.20_risk_Blim_max)) %>%
      filter(comp_b_multiplier == max(comp_b_multiplier)) %>%
      ungroup()
  ) %>%
  mutate(optimisation = c("v=1_n1=1_w=1.4_x", "v=1_n1=1_x_w", "v=2_n1=1_x_w", 
                          "v=2_n1=2_x_w", "v=1_n1=2_x_w")) %>%
  mutate(path = "output/ple.27.7e/baseline/1000_20/x_w_n1_v/hr/") %>%
  unite("file", 1:8, remove = FALSE) %>%
  mutate(file = paste0("mp_", file, ".rds"))
saveRDS(df_baseline, "output/baseline_smry_v_n1_w_x.rds")
write.csv(df_baseline, "output/baseline_smry_v_n1_w_x.csv", row.names = FALSE)

### wormplots
### go through all runs
for (i in split(df_baseline, f = seq(nrow(df_baseline)))) {#browser()
  ### find OM stock - slot depends on version of mse package
  res_mp <- readRDS(paste0(i$path, i$file))
  . <- try(stk_res <- res_mp@om@stock, silent = TRUE)
  if (is(., "try-error")) stk_res <- res_mp@stock
  path_OM <- paste0("input/ple.27.7e/baseline/1000_100/")
  stk_hist <- readRDS(paste0(path_OM, "stk.rds"))
  refpts <- readRDS(paste0(path_OM, "refpts_mse.rds"))
  p <- plot_worm(stk = stk_res, stk_hist = stk_hist, refpts = refpts)
  ggsave(filename = paste0("output/plots/wormplots/baseline_", i$optimisation, 
                           ".png"), plot = p, 
         width = 16, height = 7, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/wormplots/baseline_", i$optimisation, 
                           ".pdf"), plot = p, 
         width = 16, height = 7, units = "cm")
}


### summary statistics
### get stats
baseline_stats <- foreach(i = split(df_baseline, f = seq(nrow(df_baseline))), 
                          .combine = bind_rows) %do% {#browser()
  ### get projections and reference points
  mp_i <- readRDS(paste0(i$path, i$file))
  path_OM <- paste0("input/ple.27.7e/baseline/1000_100/")
  refpts <- readRDS(paste0(path_OM, "refpts_mse.rds"))
  stk <- mp_i@om@stock
  ### extract metrics
  stk_icv <- window(stk, start = 2034, end = 2044)
  stk <- window(stk, start = 2035, end = 2044)
  ssb_i <- c(ssb(stk)/refpts["Bmsy"])
  ssb20_i <- c(ssb(stk)[, ac(2044)]/refpts["Bmsy"])
  catch_i <- c(catch(stk)/refpts["Cmsy"])
  fbar_i <- c(fbar(stk)/refpts["Fmsy"])
  risk_i <- c(apply(ssb(stk) < c(refpts["Blim"]), 2, mean))
  icv_i <- c(iav(catch(stk_icv), period = i$interval))
  ### combine
  df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                           data.frame(val = ssb20_i, metric = "SSB20"),
                           data.frame(val = catch_i, metric = "catch"),
                           data.frame(val = fbar_i, metric = "Fbar"),
                           data.frame(val = icv_i, metric = "ICV"),
                           data.frame(val = risk_i, metric = "risk")
  ))
  return(bind_cols(i, df))
}
baseline_stats <- baseline_stats %>%
  mutate(
    MP_label = factor(optimisation,
                      levels = c("v=1_n1=1_w=1.4_x",
                                 "v=1_n1=1_x_w",
                                 "v=2_n1=1_x_w",
                                 "v=2_n1=2_x_w",   
                                 "v=1_n1=2_x_w"),
                      labels = c("x\n(v=1, n1=1, w=1.4)", 
                                 "x & w\n(v=1, n1=1)",
                                 "x & w\n(v=2, n1=1)",
                                 "x & w\n(v=2, n1=2)",
                                 "x & w\n(v=1, n1=2)")))
saveRDS(baseline_stats, file = "output/baseline_stats.rds")
# baseline_stats <- readRDS("output/baseline_stats.rds")

### plot
p_risk <- baseline_stats %>%
  filter(metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.05, colour = "red") +
  geom_col(data = baseline_stats %>%
             filter(metric == "risk") %>%
             group_by(MP_label) %>%
             summarise(val = max(val)),
           aes(x = MP_label, y = val, fill = MP_label), 
           show.legend = FALSE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = MP_label, y = val, group = MP_label),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1) +
  scale_fill_brewer(name = "", palette = "Dark2") +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p_catch <- baseline_stats %>%
  filter(metric == "catch") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ssb <- baseline_stats %>%
  filter(metric == "SSB") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_icv <- baseline_stats %>%
  filter(metric == "ICV") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  labs(y = "ICV") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank())
p <- p_risk / p_catch / p_ssb / p_icv
p
ggsave(filename = "output/plots/MP/baseline_stats_x_w_v_n0.png", plot = p,
       width = 16, height = 15, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/MP/baseline_stats_x_w_v_n0.pdf", plot = p,
       width = 16, height = 15, units = "cm", 
       bg = "white")



### ------------------------------------------------------------------------ ###
### refset - x & w - summary ####
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

### find optima
df_optima <- df_runs %>%
  group_by(group) %>%
  filter(X11.20_risk_Blim_max <= 0.05) %>%
  filter(X11.20_Catch_rel == max(X11.20_Catch_rel))

### save results
saveRDS(df_runs, file = "output/refset_x_w_grid.rds")
saveRDS(df_optima, file = "output/refset_x_w_grid_opt.rds")
write.csv(df_optima, "output/refset_x_w_grid_opt.csv", row.names = FALSE)

### plot raw data
p_raw <- df_runs %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(x = multiplier, y = comp_b_multiplier, label = catch,
             fill = catch_col)) +
  geom_tile(alpha = 0.8) +
  #geom_text(aes(label = round(X11.20_Catch_rel, 3)), size = 2.5) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_optima,
             aes(yintercept = comp_b_multiplier),
             linewidth = 0.3, linetype = "1111", colour = "red") +
  geom_vline(data = df_optima,
             aes(xintercept = multiplier),
             linewidth = 0.3, linetype = "1111", colour = "red") +
  labs(x = "Multiplier (x)", y = expression(I[trigger]~multiplier~"(w)")) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +
  coord_cartesian(expand = FALSE, xlim = c(0, 1)) +
  theme_bw(base_size = 8)
p_raw
ggsave(filename = "output/plots/MP/refset_x_w_grid.png", plot = p_raw,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/MP/refset_x_w_grid.pdf", plot = p_raw,
       width = 16, height = 10, units = "cm", 
       bg = "white")


### interpolate missing cells by group (for plotting only)
df_runs_int <- foreach(group_i = levels(df_runs$group),
                       .combine = bind_rows) %do% {
  data_i <- df_runs %>%
    filter(group == group_i) %>%
    mutate(catch = X11.20_Catch_rel) %>%
    select(w = comp_b_multiplier, x = multiplier, 
           catch = X11.20_Catch_rel, risk = X11.20_risk_Blim_max,
           group)
  n_x <- length(seq(min(data_i$x), max(data_i$x), 0.01))
  n_w <- length(seq(min(data_i$w), max(data_i$w), 0.01))
  
  out_catch <- akima::interp(x = data_i$x, y = data_i$w,
                             z = data_i$catch, 
                             nx = n_x, ny = n_w)
  out_risk <- akima::interp(x = data_i$x, y = data_i$w,
                            z = data_i$risk, 
                            nx = n_x, ny = n_w,
                            linear = FALSE)
  
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
  mutate(catch_col = ifelse(risk <= 0.05, catch, NA)) %>%
  ggplot(aes(x = x, y = w, fill = catch_col)) +
  geom_raster(alpha = 0.8, interpolate = TRUE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_optima,
             aes(yintercept = comp_b_multiplier),
             linewidth = 0.3, linetype = "1111", colour = "red") +
  geom_vline(data = df_optima,
             aes(xintercept = multiplier),
             linewidth = 0.3, linetype = "1111", colour = "red") +
  labs(x = "Multiplier (x)", y = expression(I[trigger]~multiplier~"(w)")) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +
  coord_cartesian(expand = FALSE, xlim = c(0, 1)) +
  theme_bw(base_size = 8)
p_int
ggsave(filename = "output/plots/MP/refset_x_w_grid_int.png", plot = p_int,
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/MP/refset_x_w_grid_int.pdf", plot = p_int,
       width = 16, height = 10, units = "cm", 
       bg = "white")

### summary table
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(df_x, df_x_w)
df_smry <- df_x_w %>%
  ungroup() %>%
  select(index, n1 = idxB_range_3, v = interval, x = multiplier,
         w = comp_b_multiplier,
         risk = X11.20_risk_Blim_max,
         catch = X11.20_Catch_rel,
         ssb = X11.20_SSB_rel)
df_smry
write.csv(df_smry, file = "output/refset_x_w_smry.csv", row.names = FALSE)


### ------------------------------------------------------------------------ ###
### refset - x & w - violin plots ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(df_x, df_x_w)
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds"))
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
               "R: failure", "Catch: +10%", "Catch: -20%", 
               "Uncertainty:\nindex +20%")
OMs_group <- c("refset (combined)", rep("refset", 7), rep("robset", 7))

### get stats
# , .combine = bind_rows
stats <- foreach(i = split(df_x_w, f = seq(nrow(df_x_w))), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows) %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/",
                     ifelse(identical(i$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    mp_i <- readRDS(paste0(path_i, i$file))

    ### get reference points
    refpts <- input_refpts(OM = OM)

    ### extract metrics
    stk <- mp_i@om@stock
    stk_icv <- window(stk, start = 2034, end = 2044)
    stk <- window(stk, start = 2035, end = 2044)
    ssb_i <- c(ssb(stk)/refpts["Bmsy"])
    catch_i <- c(catch(stk)/refpts["Cmsy"])
    fbar_i <- c(fbar(stk)/refpts["Fmsy"])
    risk_i <- c(apply(ssb(stk) < rep(c(refpts["Blim"]), 
                                     each = dim(ssb(stk))[2]), 2, mean))
    icv_i <- c(iav(catch(stk_icv), period = i$interval))
    ### combine
    df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = catch_i, metric = "catch"),
                              data.frame(val = fbar_i, metric = "Fbar"),
                              data.frame(val = icv_i, metric = "ICV"),
                              data.frame(val = risk_i, metric = "risk")
    ))
    df <- df %>%
      mutate(n1 = i$idxB_range_3,
             v = i$interval,
             x = i$multiplier,
             w = i$comp_b_multiplier,
             index = i$index,
             group = i$group,
             OM = OM, OM_group = OM_group)
    return(df)
}
saveRDS(stats, file = "output/refset_stats.rds")
# stats <- readRDS("output/refset_stats.rds")

### go through all solutions and plots stats
stats_plot <- stats %>%
  mutate(OM = factor(OM, levels = OMs,
                     labels = OMs_label),
         OM_group = factor(OM_group,
                           levels = OMs_group,
                           labels = c("", 
                                      rep("Reference set", 7), 
                                      rep("Robustness set", 7))),
         group_plot = factor(group,
            levels = c("UK-FSP", "Q1SWBeam",
                       "UK-FSP (annual)", "UK-FSP (biennial)",
                       "Q1SWBeam (annual)", "Q1SWBeam (biennial)"),
            labels = c("UK-FSP - x", "Q1SWBeam - x",
                       "UK-FSP - x & w - annual", "UK-FSP - x & w - biennial",
                       "Q1SWBeam - x & w - annual", 
                       "Q1SWBeam - x & w - biennial")),
         group_file = factor(group,
           levels = c("UK-FSP", "Q1SWBeam",
                      "UK-FSP (annual)", "UK-FSP (biennial)",
                      "Q1SWBeam (annual)", "Q1SWBeam (biennial)"),
           labels = c("UK-FSP_x", "Q1SWBeam_x",
                      "UK-FSP_x_w_annual", "UK-FSP_x_w_biennial",
                      "Q1SWBeam_x_w_annual", 
                      "Q1SWBeam_x_w_biennial")))

. <- foreach(group_i = levels(stats_plot$group),
             group_plot_i = levels(stats_plot$group_plot),
             group_file_i = levels(stats_plot$group_file)) %do% {
  #browser()
  p_risk <- stats_plot %>%
    filter(metric == "risk" & group == group_i) %>%
    ggplot() +
    geom_hline(yintercept = 0.05, colour = "red") +
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
    stat_summary(aes(x = OM, y = val),
                 fun = "mean", geom = "point", shape = 4, size = 1) +
    #scale_fill_brewer(name = "", palette = "Dark2") +
    facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
    labs(y = expression(max.~B[lim]~risk), title = group_plot_i) +
    coord_cartesian(ylim = c(0, 0.3)) +
    theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  #p_risk
  p_catch <- stats_plot %>%
    filter(metric == "catch" & group == group_i) %>%
    ggplot(aes(x = OM, y = val)) +
    geom_hline(yintercept = 1, colour = "grey") +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    #scale_fill_brewer(name = "", palette = "Dark2") +
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
  p_ssb <- stats_plot %>%
    filter(metric == "SSB" & group == group_i) %>%
    ggplot(aes(x = OM, y = val)) +
    geom_hline(yintercept = 1, colour = "grey") +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    #scale_fill_brewer(name = "", palette = "Dark2") +
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
  p_icv <- stats_plot %>%
    filter(metric == "ICV" & group == group_i) %>%
    ggplot(aes(x = OM, y = val)) +
    geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
                position = position_dodge(width = 0.8), scale = "width") +
    geom_boxplot(aes(group = OM), 
                 position = position_dodge(width = 0.8),
                 fill = "white", width = 0.1, size = 0.2,
                 outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                 outlier.fill = "transparent") +
    #scale_fill_brewer(name = "", palette = "Dark2") +
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
  ggsave(filename = paste0("output/plots/MP/refset_stats_",
                           group_file_i, ".png"), 
         plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
         type = "cairo", bg = "white")
  ggsave(filename = paste0("output/plots/MP/refset_stats_",
                           group_file_i, ".pdf"), 
         plot = p, width = 16, height = 10, units = "cm", bg = "white")
}






### ------------------------------------------------------------------------ ###
### OLD CODE from here - not updated ####
### ------------------------------------------------------------------------ ###



### ------------------------------------------------------------------------ ###
### collate results - baseline MPs ####
### ------------------------------------------------------------------------ ###

### baseline OM runs for all MPs (including optimisations)
res <- foreach(MP = c("rfb", "2over3", "2over3_XSA", "hr", "ICES_SAM"), 
               .combine = bind_rows) %:%
  foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"), 
          .combine = bind_rows) %:%
  foreach(OM = c("baseline"), .combine = bind_rows) %:%
  foreach(optimised = c("default", "multiplier", "all"), 
          .combine = bind_rows) %:%
  foreach(period = c("1-20", "11-20"), .combine = bind_rows) %do% {
    #browser()
    res_i <- data.frame(stock = stock, OM = OM, MP = MP, optimised = optimised, 
                        period = period)
    if (MP %in% c("rfb", "hr")) {
      ga_path <- paste0("output/", stock, "/", OM, "/1000_20/multiplier/", MP, 
                        "/")
      if (isTRUE(MP %in% c("rfb", "hr")) & 
          optimised %in% c("default", "multiplier")) {
        ga_prefix <- "multiplier-upper_constraint1.2-lower_constraint0.7"
      } else if (identical(MP, "rfb") & identical(optimised, "all")) {
        ga_prefix <- paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                            "interval-multiplier-upper_constraint1.2-lower_",
                            "constraint0.7")
      } else if (identical(MP, "hr") & identical(optimised, "all")) {
        ga_prefix <- paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval-",
                            "multiplier-upper_constraint1.2-lower_constraint0.7")
      }
      ga_res_file <- paste0(ga_prefix, "--obj_ICES_res_", period, ".rds")
      ga_runs_file <- paste0(ga_prefix, "--obj_ICES_runs.rds")
      
      ### load file
      if (!file.exists(paste0(ga_path, ga_res_file))) return(NULL)
      ga_res <- readRDS(paste0(ga_path, ga_res_file))
      ga_solution <- as.data.frame(as.list(ga_res@solution[1, ]))
      ### round parameters to significance used in optimisation
      if (identical(MP, "rfb")) {
        ga_solution[c(1:4, 8)] <- round(ga_solution[c(1:4, 8)], 0)
        ga_solution[c(5:7)] <- round(ga_solution[c(5:7)], 1)
        ga_solution[c(9:11)] <- round(ga_solution[c(9:11)], 2)
      } else if (identical(MP, "hr")) {
        ga_solution[c(1:2, 5)] <- round(ga_solution[c(1:2, 5)], 0)
        ga_solution[c(3:4)] <- round(ga_solution[c(3:4)], 1)
        ga_solution[c(6:8)] <- round(ga_solution[c(6:8)], 2)
      }
      if (identical(optimised, "default") & identical(MP, "rfb")) {
        if (stock %in% c("ple.27.7e", "cod.27.47d20")) {
          ga_solution$multiplier <- 0.95
        } else if (stock %in% c("her.27.3a47d")) {
          ga_solution$multiplier <- 0.9
        }
      } else if (identical(optimised, "default") & identical(MP, "hr")) {
        ga_solution$multiplier <- 1
      }
      
      ### load stats
      ga_runs <- readRDS(paste0(ga_path, ga_runs_file))
      ga_run_i <- ga_runs[[paste0(ga_solution, collapse = "_")]]
      ga_stats_i <- bind_rows(ga_run_i$stats)
      
      res_i$generations <- ga_res@iter
      res_i$fitness <- ga_res@fitnessValue
      if (identical(optimised, "default")) {
        res_i$generations <- NA
        res_i$fitness <- NA
      }
      ### combine definition and stats
      res_i <- bind_cols(res_i, ga_solution, ga_stats_i)
      ### check if MP results exist
      path_i <- paste0("output/", stock, "/", OM, "/1000_20/", MP, "/")
      file_i <- paste0("mp_", paste0(ga_solution, collapse = "_"), ".rds")
      file_i <- paste0(path_i, file_i)
      res_i$file <- ifelse(file.exists(file_i), file_i, NA)
    ### other MPs - only default exists
    } else if (identical(optimised, "default")) {
      ### load stats
      path_i <- paste0("output/", stock, "/", OM, "/1000_20/", MP, "/")
      file_i <- paste0(path_i, "stats.rds")
      if (!file.exists(file_i)) return(NULL)
      stats_i <- readRDS(file_i)
      stats_i <- bind_rows(stats_i)
      ### check if MP results exist
      path_mp_i <- paste0(path_i, "mp.rds")
      stats_i$file <- ifelse(file.exists(path_mp_i), path_mp_i, NA)
      res_i <- bind_cols(res_i, stats_i)
      ### catch interval
      res_i$interval <- case_when(
        MP %in% c("2over3") ~ 2,
        MP %in% c("2over3_XSA", "ICES_SAM") ~ 1
      )
    } else {
      return(NULL)
    }
    ### calculate fitness
    if (isTRUE(is.na(res_i$fitness)) | is.null(res_i$fitness)) {
      tmp_catch <- unlist(res_i[, paste0(gsub(period, pattern = "-",
                                              replacement = ":"), 
                                         "_Catch_rel") ])
      tmp_risk <- unlist(res_i[, paste0(gsub(period, pattern = "-",
                                             replacement = ":"),
                                        "_risk_Blim_max") ])
      res_i$fitness <- tmp_catch - penalty(x = tmp_risk, 
                                           negative = FALSE, max = 1, 
                                           inflection = 0.06, 
                                           steepness = 1000)
    }
    return(res_i)
}

write.csv(res, file = "output/MPs_baseline.csv", row.names = FALSE)
saveRDS(res, file = "output/MPs_baseline.rds")
res <- readRDS("output/MPs_baseline.rds")

### ------------------------------------------------------------------------ ###
### summary table for baseline OM: optimisation results ####
### ------------------------------------------------------------------------ ###
res <- readRDS("output/MPs_baseline.rds")
res %>% 
  filter(period == "11-20" & MP %in% c("rfb", "hr")) %>%
  select(`stock`:lower_constraint, idxB_lag, idxB_range_3, comp_b_multiplier, 
         -OM, -period) %>%
  group_by(stock, MP) %>%
  ### relative growth -> used because some fitness values might be negative
  mutate(fitness_improvement = ((fitness - min(fitness))/abs(min(fitness)))*100,
         .after = fitness) %>%
  mutate(fitness_improvement = round(fitness_improvement)) %>%
  write.csv("tmp.csv", row.names = FALSE)


### ------------------------------------------------------------------------ ###
### collate results - alternative OMs ####
### ------------------------------------------------------------------------ ###

### include alternative OMs
OMs_ple <- c("baseline", "M_low", "M_high", "M_Gislason", 
             "no_discards", "rec_no_AC", "rec_failure")
OMs_cod <- c("baseline", "rec_higher", "M_dd", "M_no_migration", "rec_failure")
OMs_her <- c("baseline", "rec_higher", "rec_failure", "M_high", "M_low")
res_alt <- res %>%
  filter(period == "11-20") %>%
  select(stock:period, fitness, lag_idx:lower_constraint, file) %>%
  mutate(id = seq(n()))
res_alt <- foreach(OM = unique(c(OMs_ple, OMs_cod, OMs_her)),
                   .combine = bind_rows) %:% 
  foreach(i = split(res_alt, f = seq(nrow(res_alt))), 
          .combine = bind_rows) %do% {
    #browser()
    i$OM <- OM
    if (identical(OM, "baseline")) return(i)
    if ((identical(i$stock, "ple.27.7e") & !OM %in% OMs_ple) |
        (identical(i$stock, "cod.27.47d20") & !OM %in% OMs_cod) |
        (identical(i$stock, "her.27.3a47d") & !OM %in% OMs_her)) 
      return(NULL)
    file_i <- i$file
    if (identical(OM, "rec_failure")) {
      file_i <- gsub(x = file_i, pattern = "baseline/1000_20", 
                     replacement = "baseline/1000_20/rec_failure")
    } else {
      file_i <- gsub(x = file_i, pattern = "baseline", replacement = OM)
    }
    i$file <- ifelse(file.exists(file_i), file_i, NA)
    i$fitness <- NA
    i$OM_group = case_when(
      OM == "baseline" ~ "baseline",
      OM %in% c("M_low", "M_high", "M_Gislason", "M_dd", 
                        "M_no_migration") ~ "M",
      OM %in% c("rec_higher", "rec_no_AC", "rec_failure") ~ "Rec",
      OM %in% c("no_discards") ~ "Catch")
    return(i)
}
res_alt$OM_group[res_alt$OM == "baseline"] <- "baseline"
#sum(table(res_alt$OM_group))
#nrow(res_alt)
#table(res_alt$id)
#any(is.na(res_alt$file))
#View(res_alt)

saveRDS(res_alt, file = "output/MPs_alternative_OMs.rds")
# res_alt <- readRDS("output/MPs_alternative_OMs.rds")
write.csv(res_alt, file = "output/MPs_alternative_OMs.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate results - alternative OMs - stats ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
stats_alt <- foreach(i = split(res_alt, f = seq(nrow(res_alt))), 
        .combine = bind_rows) %do% {#browser()
  ### get projections and reference points
  mp_i <- readRDS(i$file)
  path_OM <- paste0("input/", i$stock, "/", 
                    ifelse(identical(i$OM, "rec_failure"), "baseline", i$OM), 
                    "/1000_100/")
  refpts_i <- iterMedians(readRDS(paste0(path_OM, "refpts_mse.rds")))
  ### find OM stock - slot depends on version of mse package
  . <- try(stk <- mp_i@om@stock, silent = TRUE)
  if (is(., "try-error")) stk <- mp_i@stock
  ### extract metrics
  stk_icv <- window(stk, start = 2030, end = 2040)
  stk <- window(stk, start = 2031, end = 2040)
  ssb_i <- c(ssb(stk)/refpts_i["Bmsy"])
  ssb20_i <- c(ssb(stk)[, ac(2040)]/refpts_i["Bmsy"])
  catch_i <- c(catch(stk)/refpts_i["Cmsy"])
  fbar_i <- c(fbar(stk)/refpts_i["Fmsy"])
  risk_i <- c(apply(ssb(stk) < c(refpts_i["Blim"]), 2, mean))
  icv_i <- c(iav(catch(stk_icv), period = i$interval))
  if (identical(i$OM, "no_discards")) 
    catch_i <- c(landings(stk)/refpts_i["Cmsy"])
  ### combine
  df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                            data.frame(val = ssb20_i, metric = "SSB20"),
                            data.frame(val = catch_i, metric = "catch"),
                            data.frame(val = fbar_i, metric = "Fbar"),
                            data.frame(val = icv_i, metric = "ICV"),
                            data.frame(val = risk_i, metric = "risk")
  ))
  return(bind_cols(i, df))
}
stats_alt <- stats_alt %>%
  mutate(
    OM_label = factor(OM, 
      levels = c("baseline", "M_low", "M_high", "M_Gislason",
                 "M_dd", "M_no_migration", "no_discards",
                 "rec_higher", "rec_no_AC", "rec_failure"),
      labels = c("baseline", "low", "high", "Gislason",
                 "dens. dep.", "no migration", "no discards",
                 "higher", "no AC", "failure")), .after = "OM") %>%
  mutate(OM_group = factor(OM_group, c("baseline", "M", "Catch", "Rec"))) %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
      levels = c("2over3_default", "2over3_XSA_default", 
                 "rfb_default", "rfb_multiplier", "rfb_all", 
                 "hr_default", "hr_multiplier", "hr_all", 
                 "ICES_SAM_default"),
      labels = c("2 over 3", "2 over 3 (XSA)", 
                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                 "hr (generic)", "hr (multiplier)", "hr (all)",
                 "ICES MSY")), .after = "MP") %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))

saveRDS(stats_alt, file = "output/MPs_alternative_OMs_stats.rds")
# stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")

### ------------------------------------------------------------------------ ###
### wormplots for all MPs/OMs ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")

### go through all runs
for (i in split(res_alt, f = seq(nrow(res_alt)))) {#browser()
  ### find OM stock - slot depends on version of mse package
  res_mp <- readRDS(i$file)
  . <- try(stk_res <- res_mp@om@stock, silent = TRUE)
  if (is(., "try-error")) stk_res <- res_mp@stock
  path_OM <- paste0("input/", i$stock, "/", 
                    ifelse(identical(i$OM, "rec_failure"), "baseline", i$OM), 
                    "/1000_100/")
  stk_hist <- readRDS(paste0(path_OM, "stk.rds"))
  refpts <- readRDS(paste0(path_OM, "refpts_mse.rds"))
  p <- plot_worm(stk = stk_res, stk_hist = stk_hist, refpts = refpts)
  ggsave(filename = paste0("output/plots/wormplots/", i$stock, "_", i$OM, "_",
                           i$MP, "_", i$optimised, ".png"), plot = p, 
         width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/wormplots/", i$stock, "_", i$OM, "_",
                           i$MP, "_", i$optimised, ".pdf"), plot = p, 
         width = 17, height = 8, units = "cm")
}

### ------------------------------------------------------------------------ ###
### violin plots - all stocks/OMs/MPs - grouped by OM ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")

col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
p_ple_catch <- stats_alt %>%
  filter(stock == "ple.27.7e" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_SSB <- stats_alt %>%
  filter(stock == "ple.27.7e" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_risk <- stats_alt %>%
  filter(stock == "ple.27.7e" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "ple.27.7e" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())
p_cod_catch <- stats_alt %>%
  filter(stock == "cod.27.47d20" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_SSB <- stats_alt %>%
  filter(stock == "cod.27.47d20" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_risk <- stats_alt %>%
  filter(stock == "cod.27.47d20" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "cod.27.47d20" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())
p_her_catch <- stats_alt %>%
  filter(stock == "her.27.3a47d" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_SSB <- stats_alt %>%
  filter(stock == "her.27.3a47d" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_risk <- stats_alt %>%
  filter(stock == "her.27.3a47d" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "her.27.3a47d" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

### combine catch, SSB & risk plots for each stock
p_ple <- plot_grid(p_ple_catch, p_ple_SSB, 
                   p_ple_risk + theme(legend.position = "none"),
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.6))
p_cod <- plot_grid(p_cod_catch, p_cod_SSB, 
                   p_cod_risk + theme(legend.position = "none"),
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.6))
p_her <- plot_grid(p_her_catch, p_her_SSB, 
                   p_her_risk + theme(legend.position = "none"),
                   NULL,
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.4, 0.2))

p <- plot_grid(NULL,
               plot_grid(p_ple, get_legend(p_cod_risk), 
                         nrow = 1, rel_widths = c(1, 0.22)),
               plot_grid(NULL, NULL, nrow = 1, rel_widths = c(1, 0.6),
                         labels = c("(b) Cod", "(c) Herring"),
                         label_size = 9, label_x = c(-0.025)),
               plot_grid(p_cod, p_her, align = "vh", axis = "tb",
                         nrow = 1, rel_widths = c(1, 0.6)),
               labels = c("(a) Plaice", "", "", ""),
               label_size = 9, label_x = c(-0.025),
               ncol = 1, rel_heights = c(0.07, 1, 0.07, 1), align = "v")
p
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin.png", plot = p,
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin.pdf", plot = p,
       width = 17, height = 15, units = "cm", 
       bg = "white")

### ------------------------------------------------------------------------ ###
### violin plots - all stocks/OMs/MPs - grouped by MP ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")
### create OM labels
stats_alt_MP <- stats_alt %>%
  mutate(OM_label_axis = paste0(OM_group, ": ", OM_label)) %>%
  mutate(OM_label_axis = factor(OM_label_axis, 
    levels = c("baseline: baseline", "M: low", "M: high", "M: Gislason",
               "M: dens. dep.", "M: no migration",
               "Catch: no discards",
               "Rec: no AC", "Rec: higher", "Rec: failure"),
    labels = c("baseline", "M: low", "M:high", "M: Gislason",
               "M: dens. dep.", "M: no migr.",
               "Catch: no disc.",
               "R: no AC", "R: higher", "R: failure")))

col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")


p_ple_risk <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "ple.27.7e" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(a) Plaice") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_ple_catch <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_SSB <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

p_cod_risk <- stats_alt_MP %>%
  filter(stock == "cod.27.47d20" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "cod.27.47d20" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(b) Cod") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_cod_catch <- stats_alt_MP %>%
  filter(stock == "cod.27.47d20" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_SSB <- stats_alt_MP %>%
  filter(stock == "cod.27.47d20" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

p_her_risk <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "her.27.3a47d" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(c) Herring") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_her_catch <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_SSB <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

layout <- "AA\nBB\nCC\n#D\n#E\n#F\n#G\n#H\n#I\n"
p <- p_ple_risk + p_ple_catch + p_ple_SSB +
  p_cod_risk + p_cod_catch + p_cod_SSB +
  p_her_risk + p_her_catch + p_her_SSB +
  plot_layout(design = layout, widths = c(0.195, 1))
p
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin_MP.png", plot = p,
       width = 18, height = 23, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin_MP.pdf", plot = p,
       width = 18, height = 23, units = "cm", 
       bg = "white")


### ------------------------------------------------------------------------ ###
### violin plots - baseline OM & all stocks & MPs  ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")
stats_baseline <- readRDS("output/MPs_baseline.rds")
stats_baseline <- stats_baseline %>%
  filter(period == "11-20")
col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
p_catch <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "catch" & OM == "baseline") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              scale = "width") +
  geom_boxplot(fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())
p_ssb <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "SSB" & OM == "baseline") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              scale = "width") +
  geom_boxplot(fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")
p_risk <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "risk" & OM == "baseline") %>%
  ggplot() +
  geom_hline(yintercept = 0.0525, colour = "red", size = 0.75) +
  geom_col(data = stats_alt %>%
             filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") &
                      metric == "risk" & OM == "baseline") %>%
             group_by(stock_label, MP_label) %>%
             summarise(val = max(val)),
           aes(x = MP_label, y = val, fill = MP_label),
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2) +
  geom_boxplot(aes(x = MP_label, y = val),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p_fitness <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") &
           OM == "baseline") %>%
  select(stock_label, MP_label, fitness) %>%
  unique() %>%
  ggplot() +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_col(aes(x = MP_label, y = fitness, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2) +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  labs(y = expression("fitness "*italic(phi))) +
  scale_y_continuous(#limits = c(-0.27, 1.05), 
                     breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p <- p_risk + p_catch + p_ssb + p_fitness + plot_layout(ncol = 1)

p
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.png", plot = p,
       width = 18, height = 16, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.pdf", plot = p,
       width = 18, height = 16, units = "cm",
       bg = "white")

### ------------------------------------------------------------------------ ###
### projections - baseline OM with all MPs and stocks ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_lookup <- res_alt %>%
  filter(OM == "baseline") %>%
  bind_rows(data.frame(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                       MP = "history"))
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)
refpts_ple <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
refpts_ple <- iterMedians(refpts_ple)
refpts_her <- readRDS("input/her.27.3a47d/baseline/1000_100/refpts_mse.rds")
refpts_her <- iterMedians(refpts_her)


proj <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  if (identical(i$MP, "history")) {
    stk_i <- readRDS(paste0("input/", i$stock, "/baseline/1000_100/stk.rds"))
    stk_i <- window(stk_i, end = 2020)
  } else {
    ### find OM stock - slot depends on version of mse package
    res_mp <- readRDS(i$file)
    . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
    if (is(., "try-error")) stk_i <- res_mp@stock
  }
  ### get metrics
  qnts <- FLQuants(catch = catch(stk_i)/1000, rec = rec(stk_i)/1000,
                   ssb = ssb(stk_i)/1000, fbar = fbar(stk_i))
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  return(bind_cols(i, qnts_perc))
}
proj <- proj %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("2over3_default", "2over3_XSA_default", 
                                 "rfb_default", "rfb_multiplier", "rfb_all", 
                                 "hr_default", "hr_multiplier", "hr_all", 
                                 "ICES_SAM_default"),
                      labels = c("2 over 3", "2 over 3 (XSA)", 
                                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                                 "hr (generic)", "hr (multiplier)", "hr (all)",
                                 "ICES MSY")), .after = "MP") %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))
proj_distr <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  if (identical(i$MP, "history")) {
    return(NULL)
  } else {
    ### find OM stock - slot depends on version of mse package
    res_mp <- readRDS(i$file)
    . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
    if (is(., "try-error")) stk_i <- res_mp@stock
  }
  ### get metrics & distribution in last year
  qnts <- FLQuants(catch = catch(stk_i)[, ac(2040)]/1000, 
                   rec = rec(stk_i)[, ac(2040)]/1000,
                   ssb = ssb(stk_i)[, ac(2040)]/1000, 
                   fbar = fbar(stk_i)[, ac(2040)])
  qnts <- as.data.frame(qnts)
  qnts <- qnts %>% select(year, iter, data, qname)
  return(bind_cols(i, qnts))
}
proj_distr <- proj_distr %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("2over3_default", "2over3_XSA_default", 
                                 "rfb_default", "rfb_multiplier", "rfb_all", 
                                 "hr_default", "hr_multiplier", "hr_all", 
                                 "ICES_SAM_default"),
                      labels = c("2 over 3", "2 over 3 (XSA)", 
                                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                                 "hr (generic)", "hr (multiplier)", "hr (all)",
                                 "ICES MSY")), .after = "MP")


col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("2 over 3" = "solid", 
              "2 over 3 (XSA)" = "2121", 
              "rfb (generic)" = "solid", 
              "rfb (multiplier)" = "3131", 
              "rfb (all)" = "1111", 
              "hr (generic)" = "solid", 
              "hr (multiplier)" = "3131", 
              "hr (all)" = "1111", 
              "ICES MSY" = "solid")

p_ple_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_ple["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") + 
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 4), expand = FALSE) +
  facet_wrap(~ "Plaice") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_catch_distr <- proj_distr %>% 
  filter(stock == "ple.27.7e" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, linetype = MP_label, colour = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 4), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_ple_catch <- p_ple_catch +
  annotation_custom(grob = ggplotGrob(p_ple_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 4)
p_cod_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_cod["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 150), expand = FALSE) +
  facet_wrap(~ "Cod") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_catch_distr <- proj_distr %>% 
  filter(stock == "cod.27.47d20" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 150), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_cod_catch <- p_cod_catch +
  annotation_custom(grob = ggplotGrob(p_cod_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 150)
p_her_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_her["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") + 
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 700), expand = FALSE) +
  facet_wrap(~ "Herring") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_catch_distr <- proj_distr %>% 
  filter(stock == "her.27.3a47d" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, linetype = MP_label, colour = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 700), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_her_catch <- p_her_catch +
  annotation_custom(grob = ggplotGrob(p_her_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 700)

p_ple_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_ple["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_ple["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 33), expand = FALSE) +
  #facet_wrap(~ "Plaice") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_ple_ssb_distr <- proj_distr %>% 
  filter(stock == "ple.27.7e" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 33), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_ple_ssb <- p_ple_ssb +
  annotation_custom(grob = ggplotGrob(p_ple_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 33)
p_cod_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_cod["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_cod["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 620), expand = FALSE) +
  #facet_wrap(~ "Cod") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.3, 0.65),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.key.width = unit(0.7, "lines"),
        axis.title.y = element_blank(),
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_cod_ssb_distr <- proj_distr %>% 
  filter(stock == "cod.27.47d20" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 620), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_cod_ssb <- p_cod_ssb +
  annotation_custom(grob = ggplotGrob(p_cod_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 620)
p_her_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_her["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_her["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 3500), expand = FALSE) +
  #facet_wrap(~ "Cod") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_her_ssb_distr <- proj_distr %>% 
  filter(stock == "her.27.3a47d" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 3500), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_her_ssb <- p_her_ssb +
  annotation_custom(grob = ggplotGrob(p_her_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 3500)


# p <- plot_grid(p_ple_catch, p_cod_catch, p_her_catch,
#                p_ple_ssb, p_cod_ssb, p_her_ssb,
#           ncol = 3, align = "v", axis = "l")
p <- p_ple_catch + p_cod_catch + p_her_catch + 
  p_ple_ssb + p_cod_ssb + p_her_ssb + plot_layout(ncol = 3, widths = 1)
p
ggsave(filename = "output/plots/risk_baseline_MPs_projection.png", plot = p, 
       width = 18, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_baseline_MPs_projection.pdf", plot = p, 
       width = 18, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### projections - rec_failure OM for cod ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_lookup <- res_alt %>%
  filter(OM == "rec_failure" & stock == "cod.27.47d20" &
         ((MP == "rfb" & optimised == "multiplier") | 
            (MP == "hr" & optimised == "multiplier") |
           MP == "ICES_SAM"))
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)

### load projections
proj <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  ### find OM stock - slot depends on version of mse package
  res_mp <- readRDS(i$file)
  . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
  if (is(., "try-error")) stk_i <- res_mp@stock
  ### get metrics
  qnts <- FLQuants(catch = catch(stk_i)/1000, rec = rec(stk_i)/1000,
                   ssb = ssb(stk_i)/1000, fbar = fbar(stk_i),
                   risk = ssb(stk_i) %=% NA_real_)
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  ### annual risks
  qnts_perc$risk[,,,,, "50%"] <- apply(qnts$ssb < c(refpts_cod["Blim"]/1000), 2,
                                       mean, na.rm = TRUE)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  return(bind_cols(i, qnts_perc))
}
proj <- proj %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("rfb_multiplier", "hr_multiplier",
                                 "ICES_SAM_default"),
                      labels = c("rfb (multiplier)", "hr (multiplier)", 
                                 "ICES MSY")), 
    .after = "MP") %>%
  mutate(stock_label = "cod.27.47d20") %>%
  mutate(qname_label = factor(qname,
                              levels = c("rec", "ssb", "risk", "catch", "fbar"),
                              labels = c("'Recruitment\n  [millions]'", 
                                         "'SSB [1000 t]'", "B[lim]~'risk'",
                                         "'Catch [1000 t]'",
                                         "'F (ages 2-4)'")))
proj_plot <- proj %>% 
  filter(MP_label %in% c("rfb (multiplier)", "hr (multiplier)", 
                         "ICES MSY") & 
           qname %in% c("catch", "ssb", "rec", "risk"))
col_vals <- c("rfb (multiplier)" = "#6baed6", 
              "hr (multiplier)" = "#fb6a4a", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("rfb (multiplier)" = "3131", 
              "hr (multiplier)" = "3131", 
              "ICES MSY" = "solid")

df_area <- data.frame(
  xmin = rep(2021, 4), 
  xmax = rep(2025, 4), 
  ymin = rep(0, 4), 
  ymax = c(300, 390, 1, 47),
  qname_label = factor(c("'Recruitment\n  [millions]'", 
                         "'SSB [1000 t]'", "B[lim]~'risk'",
                         "'Catch [1000 t]'"),
                       levels = c("'Recruitment\n  [millions]'",
                                  "'SSB [1000 t]'", "B[lim]~'risk'",
                                  "'Catch [1000 t]'")))
df_text <- data.frame(
  x = rep(2023, 4), y = c(200, 200, 1, 20),
  text = c("recruitment\nfailure:\n2021-2025", "", "", ""),
  qname_label = factor(c("'Recruitment\n  [millions]'", 
                         "'SSB [1000 t]'", "B[lim]~'risk'",
                         "'Catch [1000 t]'"),
                       levels = c("'Recruitment\n  [millions]'",
                                  "'SSB [1000 t]'", "B[lim]~'risk'",
                                  "'Catch [1000 t]'")))
p <- ggplot() +
  geom_rect(data = df_area,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = 0.2, linetype = 0) +
  geom_text(data = df_text, aes(x = x, y = y, label = text),
            size = 7 * 0.35) +
  geom_line(data = proj_plot,
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label)) +
  # facet_wrap(~ qname_label, scales = "free_y", strip.position = "left",
  #            ncol = 1) +
  facet_grid(qname_label ~ "Cod", scales = "free_y", switch = "y", 
             labeller = label_parsed) +
  scale_colour_manual("", values = col_vals) +
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2020, 2040.5), ylim = c(0, NA), expand = FALSE) +
  #ylim(c(0, NA)) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.85),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.width = unit(0.7, "lines"),
        strip.text.y = element_text(margin = unit(c(0, 0, 0, 8), "pt"))
        )
        #strip.text.y = unit(c(4, 4, 4, 20), "pt"))
p
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.png", plot = p, 
       width = 8, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.pdf", plot = p, 
       width = 8, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### Plaice: explore Blim options ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")

### use ICES MSY rule as example
res_tmp <- res_alt %>% 
  filter(stock == "ple.27.7e" & OM == "baseline" & MP == "ICES_SAM" &
           optimised == "default")
mp <- readRDS(res_tmp$file)
SSBs <- ssb(mp@om@stock)[, ac(2031:2040)]
refpts <- readRDS(paste0("input/ple.27.7e/baseline/1000_100/refpts_mse.rds"))

### find Blim (SSB) where risk is 5%
risk_ecdf <- ecdf(sort(c(SSBs)))
risk_ecdf <- data.frame(SSB = sort(c(SSBs)),
                        risk = risk_ecdf(sort(c(SSBs))))
tail(which(risk_ecdf$risk < 0.05), 1)
# 499
risk_ecdf[498:500, ]
risk_ecdf$risk[500]
Blim0.05 <- risk_ecdf$SSB[500]
Blim0.05
# 60.39531

### get various Blim where SSB(R=xR0)
Blim_RR0 <- readRDS("input/ple.27.7e/baseline/1000_100/Blim_RR0.rds")

#ecdf(sort(c(SSBs)))(100)
#mean(c(SSBs) < Blim_RR0$Blim[701])
#mean(c(SSBs) < Blim_RR0$Blim[301])

p_full <- as.data.frame(window(SSBs/1000)) %>%
  ggplot(aes(x = data)) +
  geom_histogram(aes(y = (stat(count)/sum(count))),
                 binwidth = 0.5, colour = "grey", fill = "white",
                 show.legend = TRUE, size = 0.1) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  stat_ecdf(aes(y = ..y..),
           pad = FALSE, geom = "step", size = 0.5) +
  scale_y_continuous(labels = scales::percent, name = expression(cumulative~B[lim]~risk),
  sec.axis = sec_axis(trans = ~ .,
                      name = "SSB distribution\n ",
                      labels = scales::percent),
                      expand = expansion(mult = c(0, 0.05), add = 0)) +
  geom_vline(xintercept = median(c(refpts["Blim"])/1000),
             colour = "black", linetype = "dashed",
             size = 0.4) +
  labs(x = expression(italic(B)[lim]~"[1000t]"),
       title = "(a)") + 
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90),
        plot.title = element_text(face = "bold"))
  
p_zoom <- as.data.frame(window(SSBs/1000)) %>%
  ggplot(aes(x = data)) +
  geom_histogram(aes(y = (stat(count)/sum(count))),
                 binwidth = 0.1, colour = "grey", fill = "white",
                 show.legend = TRUE, size = 0.1) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  stat_ecdf(aes(y = ..y..),
            pad = FALSE, geom = "step", size = 0.5) +
  scale_y_continuous(labels = scales::percent, name = expression(cumulative~B[lim]~risk),
                     sec.axis = sec_axis(trans = ~ .,
                                         name = "SSB distribution",
                                         labels = scales::percent),
                     expand = expansion(mult = c(0, 0.05), add = 0)) +
  geom_vline(xintercept = median(c(refpts["Blim"])/1000),
             colour = "black", linetype = "dashed",
             size = 0.4) +
  annotate(geom = "label",
           x = median(c(refpts["Blim"])/1000) + 0.05, y = 0.90, 
           hjust = 0, size = 2.5, check_overlap = TRUE,
           label = expression("default:"~italic(B)[lim]==2119*"t")) +
  geom_vline(xintercept = Blim_RR0$Blim[701]/1000,
             colour = "black", linetype = "1212",
             size = 0.4) +
  annotate(geom = "label",
           x = Blim_RR0$Blim[701]/1000 + 0.05, y = 0.80, 
           hjust = 0, size = 2.5, check_overlap = TRUE,
           label = expression("SSB("*R == 0.7*R[0]*")")) +
  geom_vline(xintercept = Blim_RR0$Blim[301]/1000,
             colour = "black", linetype = "1212",
             size = 0.4) +
  annotate(geom = "label",
           x = Blim_RR0$Blim[301]/1000 + 0.05, y = 0.70, 
           hjust = 0, size = 2.5, check_overlap = TRUE,
           label = expression("SSB("*R == 0.3*R[0]*")")) +
  geom_vline(xintercept = Blim0.05/1000,
             colour = "black", linetype = "1212",
             size = 0.4) +
  annotate(geom = "label",
           x = Blim0.05/1000 + 0.05, y = 0.60, 
           hjust = 0, size = 2.5, check_overlap = TRUE,
           label = expression("SSB("*P[B[lim]] == 0.05*")")) +
  labs(x = expression(italic(B)[lim]~"[1000t]"),
       title = "(b)") + 
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90),
        plot.title = element_text(face = "bold")) +
  coord_cartesian(xlim = c(0, 4))

p <- p_full + p_zoom 
ggsave(filename = "output/plots/risk_ple_Blim_SAM.png", plot = p, 
       width = 16, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_ple_Blim_SAM.pdf", plot = p, 
       width = 16, height = 8, units = "cm")
