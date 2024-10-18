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
### number and sort
df_optima$MP <- c(1, 6)
df_optima <- df_optima %>%
  arrange(MP) %>%
  relocate(MP)

### save results
saveRDS(df_runs, file = "output/refset_x_runs.rds")
saveRDS(df_optima, file = "output/refset_x_runs_opt.rds")
write.csv(df_optima, "output/refset_x_runs_opt.csv", row.names = FALSE)
# df_optima <- readRDS("output/refset_x_runs_opt.rds")

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
  facet_grid(paste0("Index years (n1): ", idxB_range_3) ~ 
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
  p <- plot_worm_distr(stk = stk_res, stk_hist = stk_hist, refpts = refpts)
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
                                 "v=1_n1=2_x_w",
                                 "v=2_n1=1_x_w",
                                 "v=2_n1=2_x_w"),
                      labels = c("x\n(v=1, n1=1, w=1.4)", 
                                 "x & w\n(v=1, n1=1)",
                                 "x & w\n(v=1, n1=2)",
                                 "x & w\n(v=2, n1=1)",
                                 "x & w\n(v=2, n1=2)")))
saveRDS(baseline_stats, file = "output/baseline_stats.rds")
# baseline_stats <- readRDS("output/baseline_stats.rds")

### plot
p_risk <- baseline_stats %>%
  filter(metric == "risk") %>%
  ggplot() +
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
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4, linetype = "1111") +
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
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", size = 0.4, 
             linetype = "dotted") +
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
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", size = 0.4, 
             linetype = "dotted") +
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
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", size = 0.4, 
             linetype = "dotted") +
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

### find optima
df_optima <- bind_rows(
  df_runs %>%
    group_by(group) %>%
    filter(comp_b_multiplier <= 1.5) %>%
             #& !(index == "Q1SWBeam" & interval == 1)
    filter(X11.20_risk_Blim_max <= 0.05) %>%
    filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
    mutate(optimum = "local"),
  df_runs %>%
    group_by(group) %>%
    filter(X11.20_risk_Blim_max <= 0.05) %>%
    filter(X11.20_Catch_rel == max(X11.20_Catch_rel)) %>%
    mutate(optimum = "global")
) %>%
  mutate(optimum = factor(optimum, levels = c("local", "global")))

### number the MPs
df_optima$MP <- c(7, 9, 2, 4, 8, 10, 3, 5)
df_optima <- df_optima %>%
  arrange(MP) %>%
  relocate(MP)


### save results
saveRDS(df_runs, file = "output/refset_x_w_grid.rds")
# df_runs <- readRDS("output/refset_x_w_grid.rds")
saveRDS(df_optima, file = "output/refset_x_w_grid_opt.rds")
write.csv(df_optima, "output/refset_x_w_grid_opt.csv", row.names = FALSE)

### plot raw data
p_raw <- df_runs %>%
  mutate(catch = X11.20_Catch_rel,
         catch_col = ifelse(X11.20_risk_Blim_max <= 0.05, 
                            X11.20_Catch_rel, NA)) %>%
  ggplot(aes(x = multiplier, y = comp_b_multiplier, label = catch,
             fill = catch_col)) +
  geom_point(alpha = 0.8, shape = 21, stroke = NA, size = 0.4) +
  #geom_text(aes(label = round(X11.20_Catch_rel, 3)), size = 2.5) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_optima,
             aes(yintercept = comp_b_multiplier, colour = optimum),
             linewidth = 0.2, linetype = "1111") +
  geom_vline(data = df_optima,
             aes(xintercept = multiplier, colour = optimum),
             linewidth = 0.2, linetype = "1111") +
  scale_colour_manual("Optimum",
                      values = c(local = "blue", global = "red")) +
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
  mutate(catch_col = ifelse(risk <= 0.05, catch, NA)) %>%
  ggplot(aes(x = x, y = w, fill = catch_col)) +
  geom_raster(alpha = 0.8, interpolate = TRUE) +
  scale_fill_gradientn(paste0("Catch/MSY"),
                       colours = hcl.colors(10),
                       values = c(0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975,
                                  1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(data = df_optima,
             aes(yintercept = comp_b_multiplier, colour = optimum),
             linewidth = 0.3, linetype = "1111") +
  geom_vline(data = df_optima,
             aes(xintercept = multiplier, colour = optimum),
             linewidth = 0.3, linetype = "1111") +
  scale_colour_manual("Optimum",
                      values = c(local = "blue", global = "red")) +
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
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_smry <- df_x_w %>%
  ungroup() %>%
  select(MP, index, optimum, 
         n1 = idxB_range_3, v = interval, x = multiplier,
         w = comp_b_multiplier,
         risk = X11.20_risk_Blim_max,
         catch = X11.20_Catch_rel,
         ssb = X11.20_SSB_rel,
         icv = X11.20_ICV) %>%
  arrange(MP)
df_smry
write.csv(df_smry, file = "output/refset_x_w_smry.csv", row.names = FALSE)
saveRDS(df_smry, file = "output/refset_x_w_smry.rds")

### ------------------------------------------------------------------------ ###
### refset - visualise chr for all chr rule versions ####
### ------------------------------------------------------------------------ ###
df_smry <- readRDS("output/refset_x_w_smry.rds")

df_plot <- df_smry %>%
  mutate(opt_pars = ifelse(MP %in% c(1, 6), "x", "x & w"),
         interval = ifelse(v == 1, "annual", "biennial"))
View(df_plot)
df_plot <- foreach(x = split(df_plot, df_plot$MP), 
                   .combine = bind_rows) %do% {
  #browser()
  ### points of chr rule
  x_values <- c(0, x$w, 10)
  y_values <- c(0, x$x, x$x)
  x <- bind_rows(x, x, x)
  x$x_axis <- x_values
  x$y_axis <- y_values
  return(x)
}

p <- df_plot %>%
  mutate(index = factor(index, levels = c("UK-FSP", "Q1SWBeam")),
         optimum = factor(optimum, levels = c("local", "global"))) %>%
  group_by(MP) %>%
  mutate(y_axis_label = max(y_axis)) %>%
  ggplot(aes(x = x_axis, y = y_axis, group = MP, 
             linetype = optimum, colour = optimum)) +
  geom_line() +
  geom_text(aes(label = paste0("MP", MP),
                x = 4.75, y = y_axis_label + 0.045),
            show.legend = FALSE, size = 2) + 
  scale_colour_manual("Optimum",
                        values = c("local" = "blue", "global" = "red")) + 
  scale_linetype_manual("Optimum",
                        values = c("local" = "solid", "global" = "1111")) + 
  facet_grid(paste0(interval, " - ", opt_pars) ~ index) +
  labs(x = "w", y = "x") + 
  coord_cartesian(xlim = c(0, 5.9), ylim = c(0, 0.99), expand = FALSE) + 
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))
p
ggsave(filename = "output/plots/MP/refset_chr_illustration.png", plot = p,
       width = 8, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/MP/refset_chr_illustration.pdf", plot = p,
       width = 8, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### refset - x & w - violin plots - by MP ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds")) %>%
  arrange(MP)
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
    catch_i <- c(catch(stk)/refpts["Cmsy"])
    fbar_i <- c(fbar(stk)/refpts["Fmsy"])
    risk_i <- c(apply(ssb(stk) < rep(c(refpts["Blim"]), 
                                     each = dim(ssb(stk))[2]), 2, mean))
    icv_i <- c(iav(catch(stk_icv), period = i$interval))
    icv_annual_i <- c(iav(catch(stk_icv), period = 1))
    ### combine
    df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = catch_i, metric = "catch"),
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
             optimum = i$optimum,
             period = period,
             MP = i$MP) %>%
      relocate(MP)
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
                                      rep("Robustness set", 7)))) %>%
  mutate(group = paste0("MP", MP, " - ", index, " - ",
                        case_when(v == 1 ~ "annual",
                                  v == 2 ~ "biennial"),
                        " - ",
                        case_when(w == 1.4 ~ "x",
                                  w != 1.4 ~ "x & w"),
                        case_when(optimum == "local" ~ " (local optimum)",
                                  optimum == "global" ~ " (global optimum)"),
                        " - ",
                        case_when(period == "long-term" ~ "long term",
                                  period == "short-term" ~ "short term",
                                  period == "all" ~ "all years"))) %>%
  mutate(group_label = paste0(index, "_",
                              case_when(v == 1 ~ "annual",
                                        v == 2 ~ "biennial"),
                              "_",
                              case_when(w == 1.4 ~ "x",
                                        w != 1.4 ~ "x_w"),
                              "_", optimum, "_",
                              case_when(period == "long-term" ~ "long",
                                        period == "short-term" ~ "short",
                                        period == "all" ~ "all")))


#scales::show_col(scales::hue_pal()(20))
#cols <- scales::hue_pal()(15)[c(1, 2, 3:5, 8:10, 11:14, 6:7, 15)]
cols <- scales::hue_pal()(15)
#cols <- scales::hue_pal()(20)[c(1, 3, 5:7, 11:13, 14:17, 8:9, 19)]
# cols <- scales::hue_pal()(30)[c(1,
#                                 29,
#                                 16:18, # catch blue
#                                 4:6, # M brown
#                                 24:27, # R purple
#                                 19:20, # catch blue
#                                 11
#                                )]

. <- foreach(group_i = unique(stats_plot$group)) %do% {
  #browser()
  stats_plot_i <- stats_plot %>%
    filter(group == group_i)
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
  ggsave(filename = paste0("output/plots/MP/refset_stats_",
                           file_i, ".png"), 
         plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
         type = "cairo", bg = "white")
  ggsave(filename = paste0("output/plots/MP/refset_stats_",
                           file_i, ".pdf"), 
         plot = p, width = 16, height = 10, units = "cm", bg = "white")
}

### ------------------------------------------------------------------------ ###
### refset - x & w - violin plots - compare MPs ####
### ------------------------------------------------------------------------ ###

stats_plot_MP <- stats_plot %>%
  mutate(MP_label = paste0("MP", MP, " - ", index, " - ",
                        case_when(v == 1 ~ "annual",
                                  v == 2 ~ "biennial"),
                        " - ",
                        case_when(w == 1.4 ~ "x",
                                  w != 1.4 ~ "x & w"),
                        "\n", 
                        case_when(optimum == "local" ~ "(local optimum)",
                                  optimum == "global" ~ "(global optimum)"))) %>%
  mutate(MP_label = factor(MP_label,
    levels = sort(unique(MP_label))[c(1, 3:10, 2)])) %>%
  mutate(period_label = factor(period, 
                               levels = c("long-term", "short-term", "all"),
                               labels = c("long term", "short term",
                                          "all years")))

p_risk <- stats_plot_MP %>%
  filter(metric == "risk" & OM == "Reference set\n(combined)") %>%
  ggplot() +
  geom_col(data = . %>%
             group_by(MP_label, period_label) %>%
             summarise(val = max(val)),
           aes(x = MP_label, y = val), fill = "#F8766D",
           show.legend = FALSE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = MP_label, y = val),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  geom_hline(yintercept = 0.05, colour = "red", linewidth = 0.4,
             linetype = "1111") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  scale_fill_manual("", values = cols) +
  facet_wrap(~ period_label) +
  labs(y = expression(max.~B[lim]~risk)) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
#p_risk
p_catch <- stats_plot_MP %>%
  filter(metric == "catch" & OM == "Reference set\n(combined)") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_violin(fill = "#F8766D",, size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", linewidth = 0.4,
             linetype = "1111") +
  facet_wrap(~ period_label) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
#p_catch
p_ssb <- stats_plot_MP %>%
  filter(metric == "SSB" & OM == "Reference set\n(combined)") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_violin(fill = "#F8766D", size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  geom_hline(yintercept = 1, colour = "#ebebeb", linewidth = 0.4,
             linetype = "1111") +
  facet_wrap(~ period_label) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
#p_ssb
p_icv <- stats_plot_MP %>%
  filter(metric == "ICV" & OM == "Reference set\n(combined)") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_violin(fill = "#F8766D", size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  facet_wrap(~ period_label) +
  labs(y = "ICV") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())
#p_icv
p <- p_risk / p_catch / p_ssb / p_icv
p
ggsave(filename = paste0("output/plots/MP/refset_stats_comparison.png"), 
       plot = p, width = 16, height = 13, units = "cm", dpi = 600, 
       type = "cairo", bg = "white")
ggsave(filename = paste0("output/plots/MP/refset_stats_comparison.pdf"), 
       plot = p, width = 16, height = 13, units = "cm", bg = "white")


### ICV only - default (annual/biennial) and annual - long term
p_icv <- stats_plot_MP %>%
  filter(metric == "ICV" & OM == "Reference set\n(combined)" &
           period == "long-term") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_violin(fill = "#F8766D", size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  facet_wrap(~ period_label) +
  labs(y = "ICV") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.text.x = element_blank())
p_icv_annual <- stats_plot_MP %>%
  filter(metric == "ICV_annual" & OM == "Reference set\n(combined)" &
           period == "long-term") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_violin(fill = "#F8766D", size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = MP_label),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  stat_summary(aes(x = MP_label, y = val),
               fun = "mean", geom = "point", shape = 4, size = 1,
               stroke = 0.25) +
  facet_wrap(~ period_label) +
  labs(y = "ICV (annual)") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p <- p_icv_annual / p_icv
p
ggsave(filename = paste0("output/plots/MP/refset_stats_comparison_ICV.png"), 
       plot = p, width = 8, height = 12, units = "cm", dpi = 600, 
       type = "cairo", bg = "white")
ggsave(filename = paste0("output/plots/MP/refset_stats_comparison_ICV.pdf"), 
       plot = p, width = 8, height = 12, units = "cm", bg = "white")


### ------------------------------------------------------------------------ ###
### refset - x & w - wormplots ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds"))
df_x_w <- df_x_w %>%
  mutate(group = paste0(index, " - ",
                        case_when(interval == 1 ~ "annual",
                                  interval == 2 ~ "biennial"),
                        " - ",
                        case_when(comp_b_multiplier == 1.4 ~ "x",
                                  comp_b_multiplier != 1.4 ~ "x & w"),
                        case_when(optimum == "local" ~ " (local optimum)",
                                  optimum == "global" ~ " (global optimum)"))) %>%
  mutate(group_label = paste0(index, "_",
                              case_when(interval == 1 ~ "annual",
                                        interval == 2 ~ "biennial"),
                              "_",
                              case_when(comp_b_multiplier == 1.4 ~ "x",
                                        comp_b_multiplier != 1.4 ~ "x_w"),
                              "_", optimum))


OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")
OMs_label <- c("Reference set (combined)", 
               "Baseline", "Catch: no discards", "Catch: 100% discards", 
               "Catch: no migration", 
               "M: -50%", "M: +50%", "M: Gislason", 
               "R: no AC", "R: +20%", "R: -20%", 
               "R: failure", "Catch: +10%", "Catch: -10%", 
               "Uncertainty: index +20%")

. <- foreach(x = split(df_x_w, seq(nrow(df_x_w)))) %:%
  foreach(OM = OMs[-1], OM_label = OMs_label[-1])  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", 
                     ifelse(identical(x$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    mp_i <- readRDS(paste0(path_i, x$file))
    stk <- mp_i@om@stock
    
    ### historical stock
    input <- input_mp(OM = OM, n_yrs = 20, MP = "hr")
    stk_hist <- input$om@stock
    
    ### get reference points
    refpts <- input_refpts(OM = OM)
    
    ### plot
    # p <- plot_worm(stk = stk, stk_hist = stk_hist, refpts = refpts,
    #                title = paste0(MP_label, " - ", OM_label))
    p <- plot_worm_distr(stk = stk, stk_hist = stk_hist, refpts = refpts,
                         title = paste0("MP", x$MP, " - ", x$group, " - ",
                                        OM_label))
    
    ggsave(filename = paste0("output/plots/wormplots/hr_", x$group_label, 
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/plots/wormplots/hr_", x$group_label, 
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
. <- foreach(x = split(df_x_w, seq(nrow(df_x_w)))) %:%
  foreach(OM = OMs[1], OM_label = OMs_label[1])  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OMs_refset, "/1000_20/", 
                     ifelse(identical(x$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    stk <- lapply(path_i, function(y) {
      readRDS(paste0(y, x$file))@om@stock
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
                              title = paste0("MP", x$MP, " - ", x$group, " - ",
                                             OM_label))
    
    ggsave(filename = paste0("output/plots/wormplots/hr_", x$group_label, 
                             "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/plots/wormplots/hr_", x$group_label, 
                             "_", OM, ".pdf"), 
           plot = p, width = 16, height = 7.5, units = "cm")
    
}

### ------------------------------------------------------------------------ ###
### refset - x & w - proportion below Itrigger ####
### ------------------------------------------------------------------------ ###
### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste("mp", idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_"),
                       ".rds"))
df_x_w <- df_x_w %>%
  mutate(group = paste0(index, " - ",
                        case_when(interval == 1 ~ "annual",
                                  interval == 2 ~ "biennial"),
                        " - ",
                        case_when(comp_b_multiplier == 1.4 ~ "x",
                                  comp_b_multiplier != 1.4 ~ "x & w"),
                        case_when(optimum == "local" ~ " (local optimum)",
                                  optimum == "global" ~ " (global optimum)"))) %>%
  mutate(group_label = paste0(index, "_",
                              case_when(interval == 1 ~ "annual",
                                        interval == 2 ~ "biennial"),
                              "_",
                              case_when(comp_b_multiplier == 1.4 ~ "x",
                                        comp_b_multiplier != 1.4 ~ "x_w"),
                              "_", optimum))

res_b <- foreach(x = split(df_x_w, seq(nrow(df_x_w))),
                 .combine = bind_rows) %:%
  foreach(OM = "refset", OM_label = "Reference set (combined)",
          .combine = bind_rows)  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", 
                     ifelse(identical(x$index, "Q1SWBeam"),
                            "multiplier_Q1SWBeam", "multiplier"),
                     "/hr/")
    tracking <- readRDS(paste0(path_i, x$file))@tracking[[1]]
    
    b_tmp <- tracking["comp_b", ]
    b_tmp <- iterMeans(b_tmp)
    df_tmp <- as.data.frame(b_tmp) %>%
      select(year, b = data) %>%
      mutate(MP = x$MP)
    
    return(df_tmp)
    
}

p <- res_b %>%
  mutate(prop = 1 - b) %>%
  mutate(MP_label = paste0("MP", MP)) %>%
  mutate(MP_label = factor(MP_label,
                           levels = paste0("MP", 1:10))) %>%
  ggplot(aes(x = year, y = prop)) +
  geom_line() +
  facet_wrap(~ MP_label, nrow = 2) +
  labs(x = "Year", y = expression("Proportion below "*I[trigger])) +
  coord_cartesian(xlim = c(2024.5, NA), ylim = c(-0.01, 0.5), expand = FALSE) +
  theme_bw(base_size = 8)
p

ggsave(filename = paste0("output/plots/MP/refset_prop_b.png"),
       plot = p, width = 16, height = 6, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = paste0("output/plots/MP/refset_prop_b.pdf"), 
       plot = p, width = 16, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### rfb & SAM - violin plots - by OM ####
### ------------------------------------------------------------------------ ###

### list of all operating models
OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_label <- c("Reference set\n(combined)", 
               "Baseline", "Catch:\nno discards", "Catch:\n100% discards", 
               "Catch:\nno migration", 
               "M: -50%", "M: +50%", "M: Gislason", 
               "R: no AC", "R: +20%", "R: -20%", 
               "R: failure", "Catch: +10%", "Catch: -10%", 
               "Uncertainty:\nindex +20%")
OMs_group <- c("refset (combined)", rep("refset", 7), rep("robset", 7))

### get stats
# , .combine = bind_rows
stats <- foreach(MP = c("rfb", "ICES_SAM"), 
                 .combine = bind_rows) %:%
  foreach(OM = OMs, OM_group = OMs_group, .combine = bind_rows)  %:%
  foreach(period = c("long-term", "short-term", "all"),
          period_yrs = list(2035:2044, 2025:2034, 2025:2044),
          .combine = bind_rows) %do% {
    #browser()
    ### refset OM - combine manually
    if (identical(OM, "refset")) {
      stks <- lapply(OMs_refset, function(OM_i) {
        path_i <- paste0("output/ple.27.7e/", OM_i, "/1000_20/", MP, "/")
        mp_i <- readRDS(paste0(path_i, "mp.rds"))
        stk_i <- mp_i@om@stock
        return(stk_i)
      })
      stk <- Reduce(FLCore::combine, stks)
      
    } else {
      ### get projection
      path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", MP, "/")
      #if (!file.exists(paste0(path_i, "mp.rds"))) return(NULL)
      mp_i <- readRDS(paste0(path_i, "mp.rds"))
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
      mutate(MP = MP, OM = OM, OM_group = OM_group, period = period)
    return(df)
}
saveRDS(stats, file = "output/altMPs_stats.rds")
# stats <- readRDS("output/altMPs_stats.rds")

### go through all solutions and plots stats
stats_plot <- stats %>%
  mutate(
   OM = factor(OM, levels = OMs,
                     labels = OMs_label),
   OM_group = factor(OM_group,
                     levels = OMs_group,
                     labels = c("", 
                                rep("Reference set", 7), 
                                rep("Robustness set", 7))),
   MP_title = case_when(MP == "rfb" ~ "rfb (default)",
                        MP == "ICES_SAM" ~ "ICES MSY (with SAM)"),
   period_label = factor(period, levels = c("long-term", "short-term", "all"),
                         labels = c("long term", "short term",
                                    "all years")))

cols <- scales::hue_pal()(15)

. <- foreach(MP_i = unique(stats_plot$MP),
             MP_title_i = unique(stats_plot$MP_title)) %:%
  foreach(period_i = unique(stats_plot$period),
          period_label_i = unique(stats_plot$period_label)) %do% {
  #browser()
  p_risk <- stats_plot %>%
   filter(metric == "risk" & MP == MP_i & period == period_i) %>%
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
                fun = "mean", geom = "point", shape = 4, size = 1) +
   scale_fill_manual("", values = cols) +
   facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
   labs(y = expression(max.~B[lim]~risk), 
        title = paste0(MP_title_i, " - ", period_label_i)) +
   coord_cartesian(ylim = c(0, NA)) +
   theme_bw(base_size = 8) +
   theme(panel.spacing.x = unit(0, "lines"),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(hjust = 0.5))
  #p_risk
  p_catch <- stats_plot %>%
   filter(metric == "catch" & MP == MP_i & period == period_i) %>%
   ggplot(aes(x = OM, y = val)) +
   geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
               position = position_dodge(width = 0.8), scale = "width") +
   geom_boxplot(aes(group = OM), 
                position = position_dodge(width = 0.8),
                fill = "white", width = 0.1, size = 0.2,
                outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                outlier.fill = "transparent") +
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
  p_ssb <- stats_plot %>%
   filter(metric == "SSB" & MP == MP_i & period == period_i) %>%
   ggplot(aes(x = OM, y = val)) +
   geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
               position = position_dodge(width = 0.8), scale = "width") +
   geom_boxplot(aes(group = OM), 
                position = position_dodge(width = 0.8),
                fill = "white", width = 0.1, size = 0.2,
                outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                outlier.fill = "transparent") +
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
  p_icv <- stats_plot %>%
   filter(metric == "ICV" & MP == MP_i & period == period_i) %>%
   ggplot(aes(x = OM, y = val)) +
   geom_violin(aes(fill = OM), size = 0.2, show.legend = FALSE,
               position = position_dodge(width = 0.8), scale = "width") +
   geom_boxplot(aes(group = OM), 
                position = position_dodge(width = 0.8),
                fill = "white", width = 0.1, size = 0.2,
                outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
                outlier.fill = "transparent") +
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
  ggsave(filename = paste0("output/plots/MP/", MP_i, "_stats_",
                           period_i, ".png"), 
        plot = p, width = 16, height = 10, units = "cm", dpi = 600, 
        type = "cairo", bg = "white")
  ggsave(filename = paste0("output/plots/MP/", MP_i, "_stats_",
                           period_i, ".pdf"), 
        plot = p, width = 16, height = 10, units = "cm", bg = "white")
}

### summary table
stats_smry <- foreach(stats_i = split(stats_plot, f = stats_plot$MP), 
                      .combine = bind_rows) %:% 
  foreach(period_i = unique(stats_plot$period),
          period_label_i = unique(stats_plot$period_label)) %do% {
  #browser()
  stats_y <- foreach(metric_i = unique(stats_i$metric),
                     .combine = bind_rows) %do% {
    #browser()
    if (!identical(metric_i, "risk")) {
      res_y <- stats_i %>% 
        filter(metric == metric_i & period == period_i) %>%
        group_by(metric, MP, OM, period) %>%
        summarise(val = median(val))
    } else {
      res_y <- stats_i %>% 
        filter(metric == metric_i & period == period_i) %>%
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
         OM = gsub(x = OM, pattern = "\\(combined\\)", replacement = "")) 
saveRDS(stats_smry, file = "output/altMPs_stats_smry.rds")
write.csv(stats_smry, file = "output/altMPs_stats_smry.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### rfb & SAM - wormplots ####
### ------------------------------------------------------------------------ ###
OMs <- c("refset", 
         "baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
         "M_low", "M_high", "M_Gislason", 
         "R_no_AC", "R_higher", "R_lower", 
         "R_failure", "overcatch", "undercatch", "Idx_higher")
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_label <- c("Reference set (combined)", 
               "Baseline", "Catch: no discards", "Catch: 100% discards", 
               "Catch: no migration", 
               "M: -50%", "M: +50%", "M: Gislason", 
               "R: no AC", "R: +20%", "R: -20%", 
               "R: failure", "Catch: +10%", "Catch: -10%", 
               "Uncertainty: index +20%")

. <- foreach(MP = c("rfb", "ICES_SAM"), 
             MP_label = c("rfb", "ICES MSY (with SAM)")) %:%
  foreach(OM = OMs[-1], OM_label = OMs_label[-1])  %do% {
  #browser()
  ### get projection
  path_i <- paste0("output/ple.27.7e/", OM, "/1000_20/", MP, "/")
  mp_i <- readRDS(paste0(path_i, "mp.rds"))
  stk <- mp_i@om@stock
  
  ### historical stock
  input <- input_mp(OM = OM, n_yrs = 20, MP = "hr")
  stk_hist <- input$om@stock
  
  ### get reference points
  refpts <- input_refpts(OM = OM)
  
  ### plot
  # p <- plot_worm(stk = stk, stk_hist = stk_hist, refpts = refpts,
  #                title = paste0(MP_label, " - ", OM_label))
  p <- plot_worm_distr(stk = stk, stk_hist = stk_hist, refpts = refpts,
                       title = paste0(MP_label, " - ", OM_label))

  ggsave(filename = paste0("output/plots/wormplots/", MP, "_", OM, ".png"),
         plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
         type = "cairo")
  ggsave(filename = paste0("output/plots/wormplots/", MP, "_", OM, ".pdf"), 
         plot = p, width = 16, height = 7.5, units = "cm")

}

### plot refset
OMs_refset <- c("baseline", "Catch_no_disc", "Catch_no_surv", "migr_none", 
                "M_low", "M_high", "M_Gislason")
OMs_refset_label <- c("Baseline", "Catch:\nno discards", 
                      "Catch:\n100% discards", 
                      "Catch:\nno migration", 
                      "M: -50%", "M: +50%", "M: Gislason")
. <- foreach(MP = c("rfb", "ICES_SAM"), 
             MP_label = c("rfb", "ICES MSY (with SAM)")) %:%
  foreach(OM = OMs[1], OM_label = OMs_label[1])  %do% {
    #browser()
    ### get projection
    path_i <- paste0("output/ple.27.7e/", OMs_refset, "/1000_20/", MP, "/")
    stk <- lapply(path_i, function(y) {
      readRDS(paste0(y, "mp.rds"))@om@stock
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
                              title = paste0(MP_label, " - ", OM_label))
    
    ggsave(filename = paste0("output/plots/wormplots/", MP, "_", OM, ".png"),
           plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
           type = "cairo")
    ggsave(filename = paste0("output/plots/wormplots/", MP, "_", OM, ".pdf"), 
           plot = p, width = 16, height = 7.5, units = "cm")
    
}


### ------------------------------------------------------------------------ ###
### refset - x & w - sensitivity to index uncertainty ####
### ------------------------------------------------------------------------ ###

### get optimised solutions
### get optimised solutions
df_x <- readRDS("output/refset_x_runs_opt.rds")
df_x_w <- readRDS("output/refset_x_w_grid_opt.rds")
df_x_w <- bind_rows(
  df_x %>% mutate(optimum = "global"), 
  df_x_w)
df_x_w <- df_x_w %>%
  mutate(file = paste0(paste(idxB_lag, idxB_range_3, exp_b, 
                             comp_b_multiplier, interval, multiplier, 
                             upper_constraint, lower_constraint, 
                             sep = "_")))
runs_idx <- foreach(df_i = split(df_x_w, f = seq(nrow(df_x_w))),
                    .combine = bind_rows) %do% {
  # browser()
  path <- "output/ple.27.7e/refset/1000_20/sensitivity_idx/hr/"
  runs_i <- readRDS(paste0(path, "runs_", df_i$file, "_0-2_",
                           df_i$index, ".rds"))
  runs_i <- runs_i %>%
    mutate(MP = df_i$MP) %>%
    relocate(MP, idx_unc)
  
  return(runs_i)
  
}
saveRDS(runs_idx, 
        file = "output/ple.27.7e/refset/1000_20/sensitivity_idx/hr/runs.rds")

### plot
df_plot <- runs_idx %>%
  select(MP, idx_unc, risk = `11:20_risk_Blim_max`,
         SSB = `11:20_SSB_rel`, catch = `11:20_Catch_rel`) %>%
  pivot_longer(-1:-2) %>%
  mutate(MP_label = factor(MP, levels = 1:10,
                           labels = paste0("MP", 1:10)),
         name = factor(name, levels = c("risk", "catch", "SSB"),
                       labels = c("max.~B[lim]~risk", "Catch/MSY", 
                                  "SSB/B[MSY]"
                                  )))
df_risk <- data.frame(value = 0.05,
                      name = "max.~B[lim]~risk") %>%
  mutate(name = factor(name, levels = c("max.~B[lim]~risk", "Catch/MSY", 
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
  facet_grid(name ~ MP_label, scales = "free_y", switch = "y", 
             labeller = label_parsed) + 
  labs(x = "Index uncertainty multiplier") +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  ylim(c(0, NA)) + 
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8))
p
ggsave(filename = paste0("output/plots/MP/refset_idx_unc.png"),
       plot = p, width = 16, height = 7.5, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = paste0("output/plots/MP/refset_idx_unc.pdf"), 
       plot = p, width = 16, height = 7.5, units = "cm")


### ------------------------------------------------------------------------ ###
### HR hockey-stick principle visualisation ####
### ------------------------------------------------------------------------ ###

data.frame(x = c(0, 1, 2),
           y = c(0, 1, 1)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  scale_x_continuous("Biomass index I", expand = c(0, 0),
                     breaks = c(0, 1), 
                     labels = c(0, expression(italic(I)[trigger]))) +
  scale_y_continuous("Harvest rate", limits = c(0, 1.2), expand = c(0, 0),
                     breaks = c(0, 1), 
                     labels = c(0, expression(italic(H)))) +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = 1,
           linetype = "dotted") +
  annotate(geom = "segment", x = 0, xend = 1, y = 1, yend = 1,
           linetype = "dotted") +
  theme_classic()
ggsave(filename = "output/plots/HR_principle.png",
       width = 8.5, height = 5, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/HR_principle.pdf",
       width = 8.5, height = 5, units = "cm")


### ------------------------------------------------------------------------ ###
### Exceptional circumstances - biomass index range ####
### ------------------------------------------------------------------------ ###
### use MP4

### load mp results
mp <- readRDS("output/ple.27.7e/refset/1000_20/multiplier/hr/mp_1_2_1_1.02_2_0.58_1.2_0.7.rds")
### input data
input <- input_mp(OM = "refset", n_yrs = 20)

### historical index data
### @index slot does not include weight at age
idxB_hist <- quantSums(input$oem@observations$idx$`UK-FSP`@index *
  input$oem@observations$idx$`UK-FSP`@catch.wt *
  input$oem@deviances$idx$`UK-FSP`)

### projected index
### @index slot includes weight at age
idxB_proj <- quantSums(mp@oem@observations$idx$`UK-FSP`@index *
                         window(input$oem@deviances$idx$`UK-FSP`, start = 2024))

### combine history and projection
idxB <- idxB_hist
idxB[, ac(2024:2044)] <- idxB_proj

#plot(idxB) + ylim(c(0, NA))

### get percentiles
idxB_qnt <- quantile(idxB, 
                     probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                     na.rm = TRUE)
df <- as.data.frame(idxB_qnt) %>%
  select(year, iter, data) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(period = ifelse(year < 2024, "Data", "Projection"))

### plot
p <- df %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_line(aes(x = year, y = `50%`), linewidth = 0.4) +
  facet_grid(1 ~ period, shrink = TRUE, space = "free_x", scales = "free_x") +
  coord_cartesian(ylim = c(0, 3.5), expand = FALSE) + 
  labs(x = "Year", y = "UK-FSP biomass index (kg/hr m beam)") + 
  theme_bw(base_size = 8) +
  theme(strip.text.y = element_blank())
ggsave(filename = "output/plots/EC/MP4_idx_hist_proj.png",
       width = 16, height = 8, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/EC/MP4_idx_hist_proj.pdf",
       width = 16, height = 8, units = "cm")
