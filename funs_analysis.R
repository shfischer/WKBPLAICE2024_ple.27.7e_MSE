### ------------------------------------------------------------------------ ###
### wormplot of projection ####
### ------------------------------------------------------------------------ ###

plot_worm_refset <- function(stk, stk_hist, refpts, 
                      n_OMs = 7,
                      scale_ssb = 1/1000, scale_catch = 1/1000, 
                      scale_rec = 1/1000,
                      history = TRUE, 
                      its = 1:5,
                      n_years = 20, yr_end = 2044, yr_start = 2010,
                      xintercept = 2024, 
                      title_rec = "Recruitment (1000s)",
                      title_catch = "Catch (1000t)",
                      title_ssb = "SSB (1000t)",
                      title_fbar = paste0("Mean F (ages ", 
                                          range(stk)[["minfbar"]], "-", 
                                          range(stk)[["maxfbar"]], ")")
) {
  #browser()
  ### load projection
  stk_plot <- stk_hist
  yrs_res <- dimnames(stk)$year
  stk_plot[, ac(yrs_res)] <- stk
  stk <- stk_plot
  stk <- window(stk, max(yr_start - 10, dims(stk)$minyear), end = yr_end)
  ### load reference points
  refpts_plot <- refpts
  #refpts <- iterMedians(refpts)
  ### get metrics
  qnts <- FLQuants(catch = catch(stk), rec = rec(stk),
                   ssb = ssb(stk), fbar = fbar(stk))
  ### scale values
  qnts$catch <- qnts$catch * scale_catch
  qnts$rec <- qnts$rec * scale_rec
  qnts$ssb <- qnts$ssb * scale_ssb
  refpts_plot["Cmsy"] <- refpts_plot["Cmsy"] * scale_catch
  refpts_plot["Bmsy"] <- refpts_plot["Bmsy"] * scale_ssb

  ### percentiles - all OMs
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  
  ### median by OM
  n_iter_total <- dim(stk)[6]
  n_iter_OM <- n_iter_total/n_OMs
  qnts_OMs <- foreach(OM = seq(n_OMs), 
                      OM_iters = split(seq(n_iter_total), 
                            rep(seq(n_OMs), each = n_iter_OM)),
                      .combine = bind_rows) %do% {
    qnts_i <- FLQuants(catch = iter(catch(stk), OM_iters), 
                       rec = iter(rec(stk), OM_iters), 
                       ssb = iter(ssb(stk), OM_iters), 
                       fbar = iter(fbar(stk), OM_iters))
    qnts_med_i <- lapply(qnts_i, iterMedians)
    qnts_med_i <- FLQuants(qnts_med_i)
    qnts_med_i <- as.data.frame(qnts_med_i)
    qnts_med_i <- qnts_med_i %>% 
      select(year, value = data, qname) %>%
      mutate(group = OM)
    return(qnts_med_i)
}
  
  
  df_refpts <- bind_rows(
    data.frame(qname = "catch", 
               value = unique(c(refpts_plot["Cmsy"])),
               group = seq_along(unique(c(refpts_plot["Cmsy"])))),
    data.frame(qname = "rec", 
               value = NA,
               group = NA),
    data.frame(qname = "ssb", 
               value = unique(c(refpts_plot["Bmsy"])),
               group = seq_along(unique(c(refpts_plot["Bmsy"])))),
    data.frame(qname = "fbar", 
               value = unique(c(refpts_plot["Fmsy"])),
               group = seq_along(unique(c(refpts_plot["Fmsy"]))))
  ) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "ssb", "fbar")))

  p <- qnts_perc %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_hline(data = df_refpts, 
               aes(yintercept = value, colour = as.factor(group)),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line() +
    geom_line(data = qnts_OMs,
              aes(x = year, y = value, colour = as.factor(group))) +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, NA), xlim = c(yr_start, NA), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank())
  p
  if (identical(OM, "refset")) {
    p <- p + 
      scale_colour_brewer(palette = "Dark2")
  } else {
    p <- p + geom_hline(data = df_refpts, 
                        aes(yintercept = value),
                        colour = "black", linewidth = 0.5, linetype = "dashed")
  }
  
  ### individual iterations
  qnts_iter <- as.data.frame(iter(qnts, its))
  ### plot
  p_catch <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "catch"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "catch"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "catch"),
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "catch"),
              aes(x = year, y = `50%`), linewidth = 0.4) +
    geom_hline(yintercept = unique(c(refpts_plot["Cmsy"])),
               colour = "black", linewidth = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_catch),
                    expand = FALSE) + 
    labs(y = title_catch) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_rec <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "rec"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "rec"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "rec"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "rec"),
              aes(x = year, y = `50%`), size = 0.4) +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_rec)) + 
    labs(y = title_rec) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_ssb <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "ssb"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "ssb"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "ssb"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "ssb"),
              aes(x = year, y = `50%`), size = 0.4) +
    geom_hline(yintercept = c(median(refpts["Bmsy"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dashed") +
    geom_hline(yintercept = c(median(refpts["Blim"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dotted") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_ssb)) + 
    labs(y = title_ssb, x = "Year") +
    theme_bw(base_size = 8)
  p_fbar <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "fbar"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "fbar"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "fbar"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "fbar"),
              aes(x = year, y = `50%`), size = 0.4) +
    geom_hline(yintercept = c(median(refpts["Fmsy"], na.rm = TRUE)),
               colour = "black", size = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_fbar)) + 
    labs(y = title_fbar, x = "Year") +
    theme_bw(base_size = 8)
  p <-  plot_grid(p_catch, p_rec, p_fbar, p_ssb, align = "v", 
                  rel_heights = c(1, 1.1))
  return(p)
}

### ------------------------------------------------------------------------ ###
### wormplot - individual OM ####
### ------------------------------------------------------------------------ ###

plot_worm <- function(stk, stk_hist, refpts, 
                      title = NULL,
                      scale_ssb = 1/1000, scale_catch = 1/1000, 
                      scale_rec = 1/1000,
                      history = TRUE, 
                      its = 1:5,
                      n_years = 20, yr_end = 2044, yr_start = 2010,
                      xintercept = 2024, 
                      title_rec = "Recruitment (1000s)",
                      title_catch = "Catch (1000t)",
                      title_ssb = "SSB (1000t)",
                      title_fbar = paste0("Mean F (ages ", 
                                          range(stk)[["minfbar"]], "-", 
                                          range(stk)[["maxfbar"]], ")")
) {
  #browser()
  ### load projection
  stk_plot <- stk_hist
  yrs_res <- dimnames(stk)$year
  stk_plot[, ac(yrs_res)] <- stk
  stk <- stk_plot
  stk <- window(stk, max(yr_start - 10, dims(stk)$minyear), end = yr_end)
  ### load reference points
  refpts <- iterMedians(refpts)
  ### get metrics
  qnts <- FLQuants(catch = catch(stk), rec = rec(stk),
                   fbar = fbar(stk), ssb = ssb(stk))
  ### scale values
  qnts$catch <- qnts$catch * scale_catch
  qnts$rec <- qnts$rec * scale_rec
  qnts$ssb <- qnts$ssb * scale_ssb
  refpts["Cmsy"] <- refpts["Cmsy"] * scale_catch
  refpts["Bmsy"] <- refpts["Bmsy"] * scale_ssb
  refpts["Blim"] <- refpts["Blim"] * scale_ssb
  
  ### percentiles - all OMs
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  ### individual iterations
  qnts_iter <- as.data.frame(iter(qnts, its)) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  df_MSY <- data.frame(
    qname = c("catch", "fbar", "ssb"),
    value = c(refpts["Cmsy"], refpts["Fmsy"], refpts["Bmsy"])
  ) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  df_Blim <- data.frame(
    qname = c("ssb"),
    value = c(refpts["Blim"])
  ) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  p <- qnts_perc %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_hline(data = df_MSY, 
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim, 
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line(data = qnts_iter,
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line() +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, NA), xlim = c(yr_start, NA), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank())
  if (!is.null(title)) {
    p <- p +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  #p

  return(p)
}

### ------------------------------------------------------------------------ ###
### wormplot - individual OM - with distribution ####
### ------------------------------------------------------------------------ ###

plot_worm_distr <- function(stk, stk_hist, refpts, 
                      title = NULL,
                      scale_ssb = 1/1000, scale_catch = 1/1000, 
                      scale_rec = 1/1000,
                      history = TRUE, 
                      its = 1:5,
                      n_years = 20, yr_end = 2044, yr_start = 2010,
                      xintercept = 2024, 
                      title_rec = "Recruitment (1000s)",
                      title_catch = "Catch (1000t)",
                      title_ssb = "SSB (1000t)",
                      title_fbar = paste0("Mean F (ages ", 
                                          range(stk)[["minfbar"]], "-", 
                                          range(stk)[["maxfbar"]], ")")
) {
  #browser()
  ### load projection
  stk_plot <- stk_hist
  yrs_res <- dimnames(stk)$year
  stk_plot[, ac(yrs_res)] <- stk
  stk <- stk_plot
  stk <- window(stk, max(yr_start - 10, dims(stk)$minyear), end = yr_end)
  ### load reference points
  refpts <- iterMedians(refpts)
  ### get metrics
  qnts <- FLQuants(catch = catch(stk), rec = rec(stk),
                   fbar = fbar(stk), ssb = ssb(stk))
  ### scale values
  qnts$catch <- qnts$catch * scale_catch
  qnts$rec <- qnts$rec * scale_rec
  qnts$ssb <- qnts$ssb * scale_ssb
  refpts["Cmsy"] <- refpts["Cmsy"] * scale_catch
  refpts["Bmsy"] <- refpts["Bmsy"] * scale_ssb
  refpts["Blim"] <- refpts["Blim"] * scale_ssb
  
  ### percentiles - all OMs
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  ### individual iterations
  qnts_iter <- as.data.frame(iter(qnts, its)) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  df_MSY <- data.frame(
    qname = c("catch", "fbar", "ssb"),
    value = c(refpts["Cmsy"], refpts["Fmsy"], refpts["Bmsy"])
  ) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  df_Blim <- data.frame(
    qname = c("ssb"),
    value = c(refpts["Blim"])
  ) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  ### max value by qname
  df_max <- bind_rows(qnts_perc %>%
                        group_by(qname) %>%
                        filter(year >= 2010) %>%
                        summarise(max = max(`95%`)) %>%
                        mutate(source = "perc"),
                      qnts_iter %>%
                        group_by(qname) %>%
                        filter(year >= 2010) %>%
                        summarise(max = max(data)) %>%
                        mutate(source = "iter")) %>%
    group_by(qname) %>%
    summarise(max = max(max)) %>%
    mutate(max = max * 1.05)
  
  max_catch <- df_max$max[df_max$qname == title_catch]
  max_rec <- df_max$max[df_max$qname == title_rec]
  max_fbar <- df_max$max[df_max$qname == title_fbar]
  max_ssb <- df_max$max[df_max$qname == title_ssb]
  
  ### catch plot
  p_catch_distr <- as.data.frame(qnts$catch[, ac(yr_end)]) %>%
    select(data) %>%
    ggplot(aes(x = data)) +
    geom_density(aes(y = after_stat(scaled)), show.legend = FALSE, size = 0.2,
                 outline.type = "full") +
    theme_void() +
    coord_flip(xlim = c(0, max_catch), ylim = c(0, 1.05), expand = FALSE)
  p_catch <- qnts_perc %>%
    filter(qname == title_catch) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_hline(data = df_MSY %>%
                 filter(qname == title_catch),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim %>%
                 filter(qname == title_catch),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line(data = qnts_iter %>%
                filter(qname == title_catch),
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line() +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, max_catch), 
                    xlim = c(yr_start, yr_end), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(4, 10, 4, 4), "pt")) +
    annotation_custom(grob = ggplotGrob(p_catch_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_catch)
  #p_catch
  
  ### rec plot
  p_rec_distr <- as.data.frame(qnts$rec[, ac(yr_end)]) %>%
    select(data) %>%
    ggplot(aes(x = data)) +
    geom_density(aes(y = after_stat(scaled)), show.legend = FALSE, size = 0.2,
                 outline.type = "full") +
    theme_void() +
    coord_flip(xlim = c(0, max_rec), ylim = c(0, 1.05), expand = FALSE)
  p_rec <- qnts_perc %>%
    filter(qname == title_rec) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_hline(data = df_MSY %>%
                 filter(qname == title_rec),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim %>%
                 filter(qname == title_rec),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line(data = qnts_iter %>%
                filter(qname == title_rec),
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line() +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, max_rec), 
                    xlim = c(yr_start, yr_end), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(4, 10, 4, 4), "pt")) +
    annotation_custom(grob = ggplotGrob(p_rec_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_rec)
  #p_rec
  
  ### fbar plot
  p_fbar_distr <- as.data.frame(qnts$fbar[, ac(yr_end)]) %>%
    select(data) %>%
    ggplot(aes(x = data)) +
    geom_density(aes(y = after_stat(scaled)), show.legend = FALSE, size = 0.2,
                 outline.type = "full") +
    theme_void() +
    coord_flip(xlim = c(0, max_fbar), ylim = c(0, 1.05), expand = FALSE)
  p_fbar <- qnts_perc %>%
    filter(qname == title_fbar) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_hline(data = df_MSY %>%
                 filter(qname == title_fbar),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim %>%
                 filter(qname == title_fbar),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line(data = qnts_iter %>%
                filter(qname == title_fbar),
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line() +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, max_fbar), 
                    xlim = c(yr_start, yr_end), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(4, 10, 4, 4), "pt")) +
    annotation_custom(grob = ggplotGrob(p_fbar_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_fbar)
  #p_fbar
  
  ### ssb plot
  p_ssb_distr <- as.data.frame(qnts$ssb[, ac(yr_end)]) %>%
    select(data) %>%
    ggplot(aes(x = data)) +
    geom_density(aes(y = after_stat(scaled)), show.legend = FALSE, size = 0.2,
                 outline.type = "full") +
    theme_void() +
    coord_flip(xlim = c(0, max_ssb), ylim = c(0, 1.05), expand = FALSE)
  p_ssb <- qnts_perc %>%
    filter(qname == title_ssb) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_hline(data = df_MSY %>%
                 filter(qname == title_ssb),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim %>%
                 filter(qname == title_ssb),
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line(data = qnts_iter %>%
                filter(qname == title_ssb),
              aes(x = year, y = data, colour = iter), 
              linewidth = 0.1, show.legend = FALSE) + 
    geom_line() +
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, max_ssb), 
                    xlim = c(yr_start, yr_end), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(4, 10, 4, 4), "pt")) +
    annotation_custom(grob = ggplotGrob(p_ssb_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_ssb)
  #p_ssb
  
  ### combine all plots 
  p <- p_catch + p_rec + p_fbar + p_ssb + 
    plot_layout(ncol = 2) +
    plot_annotation(title = title, 
                    theme = theme(
                      plot.title = element_text(hjust = 0.5,
                                                size = 10)))
  
  return(p)
}
