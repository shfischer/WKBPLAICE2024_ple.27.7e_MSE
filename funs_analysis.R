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

### ------------------------------------------------------------------------ ###
### wormplot - several OMs - with distribution ####
### ------------------------------------------------------------------------ ###

plot_worm_distr_mult <- function(stk, stk_hist, refpts, stk_labels,
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
                                                range(stk[[1]])[["minfbar"]], "-", 
                                                range(stk[[1]])[["maxfbar"]], ")")
) {
  #browser()
  ### load projection
  #stk_lst <- FLStocks(stk)
  stk_lst <- lapply(seq(length(stk)), function(x) {
    stk_i <- stk_hist[[x]]
    stk_i[, ac(dimnames(stk[[1]])$year)] <- stk[[x]]
    return(stk_i)
  })
  stk_lst <- FLStocks(stk_lst)
  stk_lst <- window(stk_lst, max(yr_start - 10, dims(stk_lst[[1]])$minyear), 
                    end = yr_end)
  
  qnts_list <- lapply(stk_lst, function(x) {
    qnts_i <- FLQuants(catch = catch(x) * scale_catch, 
                       rec = rec(x) * scale_rec,
                       fbar = fbar(x), 
                       ssb = ssb(x) * scale_ssb)
    return(qnts_i)
  })
  df_median <- lapply(seq(length(stk)), function(x) {
    as.data.frame(FLQuants(lapply(qnts_list[[x]], iterMedians))) %>%
      mutate(stk = stk_labels[[x]])
  })
  df_median <- do.call(bind_rows, df_median) %>%
    select(year, qname, data, stk) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb))) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  
  qnts_combined <- FLQuants(
    catch = Reduce(FLCore::combine, lapply(stk_lst, catch)) * scale_catch,
    rec = Reduce(FLCore::combine, lapply(stk_lst, rec)) * scale_rec,
    fbar = Reduce(FLCore::combine, lapply(stk_lst, fbar)),
    ssb = Reduce(FLCore::combine, lapply(stk_lst, ssb)) * scale_ssb
  )
  df_perc <- lapply(qnts_combined, quantile, 
                    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                    na.rm = TRUE)
  df_perc <- FLQuants(df_perc)
  df_perc <- as.data.frame(df_perc)
  df_perc <- df_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))
  
  ### distribution in last year
  ### by stk
  df_distr_stk <- lapply(seq(length(stk)), function(x) {
    as.data.frame(qnts_list[[x]]) %>%
      filter(year == yr_end) %>%
      mutate(stk = stk_labels[[x]])
  })
  df_distr_stk <- do.call(bind_rows, df_distr_stk) %>%
    select(qname, data, stk) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb))) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  ### combined
  df_distr_combined <- as.data.frame(qnts_combined) %>%
    filter(year == yr_end) %>%
    select(qname, data) %>%
    mutate(qname = factor(qname,
                          levels = c("catch", "rec", "fbar", "ssb"),
                          labels = c(title_catch, title_rec, 
                                     title_fbar, title_ssb)))

  ### load reference points
  val_Cmsy <- sapply(refpts, function(x) {
    c(iterMedians(x["Cmsy"])) * scale_catch
  })
  df_Cmsy <- data.frame(qname = title_catch,
                        data = val_Cmsy,
                        stk = stk_labels) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  val_Bmsy <- sapply(refpts, function(x) {
    c(iterMedians(x["Bmsy"])) * scale_ssb
  })
  df_Bmsy <- data.frame(qname = title_ssb,
                        data = val_Bmsy,
                        stk = stk_labels) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  val_Blim <- sapply(refpts, function(x) {
    c(iterMedians(x["Blim"])) * scale_ssb
  })
  df_Blim <- data.frame(qname = title_ssb,
                        data = val_Blim,
                        stk = stk_labels) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  val_Fmsy <- sapply(refpts, function(x) {
    c(iterMedians(x["Fmsy"]))
  })
  df_Fmsy <- data.frame(qname = title_fbar,
                        data = val_Fmsy,
                        stk = stk_labels) %>%
    mutate(stk = factor(stk, levels = stk_labels))
  

  ### max value by qname
  df_max <- bind_rows(df_perc %>%
                        group_by(qname) %>%
                        filter(year >= 2010) %>%
                        summarise(max = max(`95%`)) %>%
                        mutate(source = "perc"),
                      df_median %>%
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
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
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
  p_rec_distr <- ggplot() +
    geom_density(data = df_distr_stk %>% 
                   filter(qname == title_rec),
                 aes(x = data, y = after_stat(scaled), colour = stk), 
                 show.legend = FALSE, size = 0.2,
                 outline.type = "full") +
    scale_colour_brewer(palette = "Dark2") +
    geom_density(data = df_distr_combined %>%
                   filter(qname == title_rec),
                 aes(x = data, y = after_stat(scaled)), 
                 show.legend = FALSE, size = 0.3,
                 outline.type = "full") +
    theme_void() +
    coord_flip(xlim = c(0, max_rec), ylim = c(0, 1.05), expand = FALSE)
  p_rec <- df_perc %>%
    filter(qname == title_rec) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = df_median %>%
                filter(qname == title_rec),
              aes(x = year, y = data, colour = stk),
              linewidth = 0.3, show.legend = TRUE) +
    scale_colour_brewer("Operating model", palette = "Dark2") + 
    #guides(color = guide_legend(ncol = 3)) +
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
          plot.margin = unit(c(4, 10, 4, 4), "pt"),
          legend.key.height = unit(0.5, "lines"),
          legend.key.width = unit(0.6, "lines"),
          legend.background = element_blank()) +
    annotation_custom(grob = ggplotGrob(p_rec_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_rec)
  #p_rec
  
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
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
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
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(4, 10, 4, 4), "pt")) +
    annotation_custom(grob = ggplotGrob(p_fbar_distr),
                      xmin = yr_end, xmax = yr_end + 1,
                      ymin = 0, ymax = max_fbar)
  #p_fbar
  
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
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
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
              linewidth = 0.3, show.legend = FALSE) +
    scale_colour_brewer(palette = "Dark2") + 
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
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = title, 
                    theme = theme(
                      plot.title = element_text(hjust = 0.5,
                                                size = 10)))
  
  return(p)
}

### ------------------------------------------------------------------------ ###
### wormplot - compare MPs ####
### ------------------------------------------------------------------------ ###

plot_worm_comparison <- function(stk, stk_hist, refpts, 
                                 names = NA,
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
                                          range(stk[[1]])[["minfbar"]], "-", 
                                          range(stk[[1]])[["maxfbar"]], ")")
) {
  #browser()
  
  stk <- FLStocks(stk_hist)
  stk_list <- FLStocks(stk_list)
  yrs_res <- dimnames(stk_list[[1]])$year
  
  for (i in seq_along(stk_list))
    stk[[i]][, ac(yrs_res)] <- stk_list[[i]]
  
  stk <- window(stk, max(yr_start - 10, dims(stk[[1]])$minyear), 
                end = yr_end)
  
  ### load reference points
  refpts <- iterMedians(refpts)
  
  ### get metrics
  
  qnts <- lapply(seq_along(stk), function(x) {
    qnts_i <- FLQuants(catch = catch(stk[[x]]), rec = rec(stk[[x]]),
                       fbar = fbar(stk[[x]]), ssb = ssb(stk[[x]]))
    ### scale values
    qnts_i$catch <- qnts_i$catch * scale_catch
    qnts_i$rec <- qnts_i$rec * scale_rec
    qnts_i$ssb <- qnts_i$ssb * scale_ssb
    ### percentiles
    qnts_i <- lapply(qnts_i, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                     na.rm = TRUE)
    qnts_i_perc <- FLQuants(qnts_i)
    qnts_i_perc <- as.data.frame(qnts_i_perc)
    qnts_i_perc$source <- names[x]
    return(qnts_i_perc)
  })
  qnts <- do.call(rbind, qnts)
  
  refpts["Cmsy"] <- refpts["Cmsy"] * scale_catch
  refpts["Bmsy"] <- refpts["Bmsy"] * scale_ssb
  refpts["Blim"] <- refpts["Blim"] * scale_ssb
  
  qnts_perc <- qnts %>% select(year, iter, data, qname, source) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
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
    ggplot(aes(x = year, y = `50%`, colour = source, fill = source)) +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE, linewidth = 0) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE, linewidth = 0) +
    geom_hline(data = df_MSY, 
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dashed",
               show.legend = FALSE) +
    geom_hline(data = df_Blim, 
               aes(yintercept = value),
               linewidth = 0.5, linetype = "dotted",
               show.legend = FALSE) +
    geom_line() +
    scale_colour_brewer("", palette = "Dark2") + 
    scale_fill_brewer("", palette = "Dark2") + 
    facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
    labs(x = "Year") +
    coord_cartesian(ylim = c(0, NA), xlim = c(yr_start, NA), expand = FALSE) +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          legend.key.height = unit(0.6, "lines"))
  if (!is.null(title)) {
    p <- p +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  #p
  
  return(p)
}

