### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
if (exists(x = "args_local")) args <- append(args, args_local)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### set default arguments
  ### parallelisation
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("mp_parallel")) mp_parallel <- FALSE
  if (!exists("ga_parallel")) ga_parallel <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  ### simulation dimensions
  if (!exists("n_iter")) n_iter <- 1000
  if (!exists("n_yrs")) n_yrs <- 20
  if (!exists("yr_start")) yr_start <- 2025
  ### OM
  if (!exists("stock_id")) stock_id <- "ple.27.7e"
  if (!exists("OM")) OM <- "baseline"
  ### MP
  if (!exists("MP")) MP <- "hr"
  ### scenario definition
  if (!exists("scenario")) scenario <- "multiplier"
  if (!exists("Ftrgt")) Ftrgt <- "MSY" # only for constF MP
  if (!exists("biomass_index")) biomass_index <- "UK-FSP" ### which biomass index?
  if (!exists("idx_unc")) idx_unc <- 1 ### increase survey uncertainty?
  # if (!exists("rec_failure")) rec_failure <- FALSE
  # if (!exists("overcatch")) overcatch <- FALSE
  # if (!exists("oem_catch_bias")) oem_catch_bias <- FALSE
  ### GA search
  if (!exists("ga_search")) ga_search <- TRUE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- FALSE
    ### objective function elements
    if (!exists("obj_fun")) obj_fun <- "ICES"
    if (!exists("obj_yrs")) obj_yrs <- "11:20"
    ### penalty function
    if (!exists("pen_neg")) pen_neg <- FALSE
    if (!exists("pen_max")) pen_max <- 1
    if (!exists("pen_infl")) pen_infl <- 0.06
    if (!exists("pen_steep")) pen_steep <- 1000
    ### GA
    if (!exists("add_suggestions")) add_suggestions <- TRUE
    ### local search after GA
    if (!exists("ga_optim")) ga_optim <- FALSE
  }
  if (!exists("stat_yrs")) stat_yrs <- "multiple"
  if (!exists("save_MP")) save_MP <- FALSE
  if (!exists("check_file")) check_file <- TRUE
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", "FLXSA",
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) 
  suppressMessages(library(package = i, character.only = TRUE))

### load additional functions
req_scripts <- c("funs.R", "funs_GA.R", "funs_WKNSMSE.R", "funs_OM.R")
for (i in req_scripts) source(i, local = TRUE)

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

### MPI
if (isTRUE(use_MPI) & isTRUE(ga_parallel)) {
  ### 1st: doMPI cluster with 1 worker per node
  message("starting doMPI")
  library(doMPI)
  cl <- startMPIcluster()
  message("startMPIcluster() succeeded")
  print(cl)
  registerDoMPI(cl)
  cl_length <- cl$workerCount
  cl_length
  ### set up MPI workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
  }
  message("MPI package loading succeeded")
  . <- foreach(i = seq(cl_length)) %dopar% {
    source("funs.R")
    source("funs_GA.R")
    source("funs_WKNSMSE.R")
  }
  message("MPI script loading succeeded")
    
  ### 2nd: doParallel workers inside doMPI workers
  if (isTRUE(mp_parallel)){  
    ### start doFuture inside MPI processes
    # if (isTRUE(n_workers > 1)) {
    #   . <- foreach(i = seq(cl_length_1)) %dopar% {
    #     cl2 <- makeCluster(n_workers)
    #     registerDoParallel(cl2)
    #     cl_length_2 <- length(cl2)
    #     ### load packages and functions into parallel workers
    #     . <- foreach(i = seq(cl_length_2)) %dopar% {
    #       for (i in req_pckgs) library(package = i, character.only = TRUE,
    #                                    warn.conflicts = FALSE, verbose = FALSE,
    #                                    quietly = TRUE)
    #       source("funs.R", echo = FALSE)
    #       source("funs_GA.R", echo = FALSE)
    #       source("funs_WKNSMSE.R", echo = FALSE)
    #     }
    #   }
    # }
    #message("setting up doParallel inside MPI succeeded")
  }
} else if (isTRUE(n_workers > 1)) {
  ### parallelisation within compute node

  ### mp() parallelisation - %dofuture%
  ### ga() parallelisation - %dopar%
  ### mse::mp() uses %dofuture% -> use fork with %dopar%
  ### https://github.com/shfischer/mse/tree/dopar

  if (isTRUE(mp_parallel) & isTRUE(ga_parallel))
    stop("Cannot parallelise mp() and ga() at the same time without MPI!")
  
  ### set up doFuture parallelisation - not used
  #plan(multisession, workers = n_workers)

  ### load packages and functions into parallel workers
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  print(cl)
  cl_length <- length(cl)
  . <- foreach(i = seq(n_workers)) %dopar% {
    for (i in req_pckgs) 
      suppressPackageStartupMessages(
        library(package = i, character.only = TRUE, warn.conflicts = FALSE, 
                verbose = FALSE, quietly = TRUE))
    for (i in req_scripts) source(i)
  }

} else {
  mp_parallel <- ga_parallel <- FALSE
  message("Not using any parallelisation")
}

### pass cluster to ga_parallel for use in ga()
if (isTRUE(ga_parallel)) ga_parallel <- cl

### ------------------------------------------------------------------------ ###
### load OM and create input for MSE ####
### ------------------------------------------------------------------------ ###

input <- input_mp(stock_id = stock_id, OM = OM, n_iter = n_iter,
                  n_yrs = n_yrs, yr_start = yr_start, n_blocks = n_blocks,
                  MP = MP, biomass_index = biomass_index)
refpts <- input_refpts(stock_id = stock_id, OM = OM, n_iter = n_iter)


### ------------------------------------------------------------------------ ###
### genetic algorithm (GA) search ####
### ------------------------------------------------------------------------ ###
if (isTRUE(MP %in% c("rfb", "hr")) & isTRUE(ga_search)) {
  
  ### GA arguments
  if (identical(MP, "rfb")) {
    ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                  "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                  "upper_constraint", "lower_constraint")
    ga_suggestions <- rbind(
      c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1.2, 0.7), ### most current data
      c(0, 1, 1, 1, 1, 1, 1, 2, 1, 1.2, 0.7),
      c(0, 1, 1, 1, 0, 0, 0, 2, 1, Inf, 0), ### constant catch
      c(0, 1, 1, 1, 0, 0, 0, 2, 1, 1.2, 0.7), ### constant catch (dummy cap)
      ### default, annual/biennial, turning off elements
      expand.grid(1, 2, 3, 1, 0:1, 0:1, 0:1, 1:2, 1, 1.2, 0.7),
      c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0), ### default without cap
      c(1, 2, 3, 1, 1, 1, 1, 2, 0, Inf, 0), ### zero catch
      c(1, 2, 3, 1, 1, 1, 1, 2, 0, 1.2, 0.7) ### zero catch but capped
    )
    ga_default <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, 1.2, 0.7)
    ga_lower <- c(0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0)
    ga_upper <- c(1, 5, 5, 1, 2, 2, 2, 5, 2, 5, 1)
    ### ga() samples uniform real (double) values from lower ga_lower to 
    ### ga_upper and these are then rounded to the significant digits
    ### -> adjust ga_lower/upper so that minima/maxima have same probability
    ga_step  <-    c(1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 0.01, 0.01, 0.01)
    ga_rounding <- c(0, 0, 0, 0,   1,   1,   1, 0,    2,    2,    2)
    ga_lower <- ga_lower - (ga_step/2 - .Machine$double.eps)
    ga_upper <- ga_upper + (ga_step/2 - .Machine$double.eps)
  } else if (identical(MP, "hr")) {
    ga_names <- c("idxB_lag", "idxB_range_3", "exp_b", "comp_b_multiplier",
                  "interval", "multiplier",
                  "upper_constraint", "lower_constraint")
    ga_suggestions <- rbind(
      c(1, 1, 1, 1.4, 1, 0, Inf, 0), ### zero catch
      expand.grid(1, 1, 1, c(0, 1, 1.4, 2), 1:2, c(0, 0.5, 1), 
                  c(1.2, Inf), c(0, 0.7))
    )
    ga_default <- c(1, 1, 1, 1.4, 1, 1, 1.2, 0.7)
    ga_lower <-   c(1, 1, 0, 0,   1, 0, 1,   0)
    ga_upper <-   c(1, 3, 2, 5,   2, 5, 5,   1)
    ### ga() samples uniform real (double) values from lower ga_lower to 
    ### ga_upper and these are then rounded to the significant digits
    ### -> adjust ga_lower/upper so that minima/maxima have same probability
    ga_step  <-    c(1, 1, 0.1, 0.01, 1, 0.01, 0.01, 0.01)
    ga_rounding <- c(0, 0,   1,    2, 0,    2,    2,    2)
    ga_lower <- ga_lower - (ga_step/2 - .Machine$double.eps)
    ga_upper <- ga_upper + (ga_step/2 - .Machine$double.eps)
  }
  ### turn of parameters not requested, i.e. limit to default value
  pos_default <- which(sapply(mget(ga_names, ifnotfound = FALSE), isFALSE))
  ga_lower[pos_default] <- ga_default[pos_default]
  ga_upper[pos_default] <- ga_default[pos_default]
  ### remove not requested parameters from suggestions
  ga_suggestions[, pos_default] <- rep(ga_default[pos_default], 
                                       each = nrow(ga_suggestions))
  ga_suggestions <- unique(ga_suggestions)
  ### fix parameters?
  pos_fixed <- which(sapply(mget(ga_names, ifnotfound = FALSE), is.numeric))
  if (isTRUE(length(pos_fixed) > 0)) {
    par_fixed <- names(pos_fixed)
    val_fixed <- mget(ga_names, ifnotfound = FALSE)[pos_fixed]
    ga_suggestions[, pos_fixed] <- NA
    ga_suggestions <- unique(ga_suggestions)
    ### parameters fixed to single value
    par_fixed_single <- names(which(sapply(val_fixed, length) == 1))
    val_fixed_single <- val_fixed[par_fixed_single]
    if (isTRUE(length(par_fixed_single) > 0)) {
      pos_par_fixed_single <- match(par_fixed_single, ga_names)
      for (pos in seq_along(par_fixed_single)) {
        ga_suggestions[, pos_par_fixed_single[pos]] <- val_fixed_single[[pos]]
        ### ensure that fixed parameters are not changed in GA
        ga_lower[pos_par_fixed_single[pos]] <- val_fixed_single[[pos]]
        ga_upper[pos_par_fixed_single[pos]] <- val_fixed_single[[pos]]
      }
      ga_suggestions <- unique(ga_suggestions)
    }
    ### parameters fixed to multiple values
    ### -> get all combinations and run all
    par_fixed_multiple <- names(which(sapply(val_fixed, length) > 1))
    if (isTRUE(length(par_fixed_multiple) > 0)) {
      val_fixed_multiple <- val_fixed[par_fixed_multiple]
      pos_par_fixed_multiple <- match(par_fixed_multiple, ga_names)
      ga_suggestions[, pos_par_fixed_multiple] <- NA
      ga_suggestions <- unique(ga_suggestions)
      ga_suggestions <- lapply(as.list(ga_suggestions), unique)
      for (pos in seq_along(par_fixed_multiple)) {
        ga_suggestions[[pos_par_fixed_multiple[pos]]] <-
          val_fixed[[par_fixed_multiple[pos]]]
      }
      ### get all combinations
      ga_suggestions <- unique(expand.grid(ga_suggestions))
    }
  } else {
    par_fixed_single <- par_fixed_multiple <- NULL
    val_fixed_single <- val_fixed_multiple <- NULL
  }
  ### finalise
  ga_suggestions <- unique(ga_suggestions)
  names(ga_suggestions) <- ga_names
  
  ### if only one parameter modified & fixed, 
  ### or parameters fixed to multiple values:
  ### run all supplied values instead of GA search
  if (exists("par_fixed")) {
    if (isTRUE(length(par_fixed_multiple) > 0 ) |
        (length(pos_default) + length(pos_fixed)) == length(ga_names)) {
      ### adapt GA dimensions
      maxiter <- run <- 1
      popSize <- nrow(ga_suggestions)
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  
  ### output path
  ### set name depending on which GA parameters are used
  scn_pars <- ga_names[setdiff(seq_along(ga_names), pos_default)]
  if (isTRUE(length(pos_fixed) > 0)) {
    if (isTRUE(length(val_fixed_single) > 0)) {
      scn_pars[which(scn_pars %in% par_fixed_single)] <- paste0(
        scn_pars[which(scn_pars %in% par_fixed_single)], val_fixed_single)
    }
  }
  scn_pars_c <- paste0(scn_pars, collapse = "-")
  
  path_out <- paste0("output/", stock_id, "/", paste0(OM, collapse = "-"), "/",
                     n_iter, "_", n_yrs, "/", scenario, "/", MP, "/")
  dir.create(path_out, recursive = TRUE)
  
  ### objective function elements
  obj_desc <- paste0("obj_", paste0(obj_fun, collapse = "_"), collapse = "")
  if (isTRUE(obj_yrs == "all")) obj_yrs <- paste0("1:", n_yrs)
  
  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()
  
  ### ------------------------------------------------------------------------ ###
  ### check if previous solutions can be used as suggestions ####
  ### ------------------------------------------------------------------------ ###
  
  ### years for summary statistics
  file_ext <- ifelse(obj_yrs == "all", "_res",
                     paste0("_res_", 
                            gsub(x = obj_yrs, pattern = ":", replacement = "-")))
  ### suffix if different risk limit used
  file_ext <- ifelse(isTRUE(!identical(pen_infl, 0.06) &
                              any(c("MSYPA", "ICES") %in% obj_fun)),
                     paste0(file_ext, "_", pen_infl),
                     file_ext)
  file_ext <- paste0(file_ext, ".rds")

  if (isTRUE(add_suggestions)) {
    ### find files
    avail <- list.files(path_out, pattern = paste0("--", obj_desc, file_ext))
    avail <- gsub(x = avail, pattern = paste0("--", obj_desc, file_ext),
                  replacement = "")
    avail <- strsplit(x = avail, split = "-")
    ### need to have fewer parameters
    avail <- avail[which(sapply(avail, length) < length(scn_pars))]
    ### skip parameters not used
    if (isTRUE(length(avail) > 0)) {
      avail <- avail[which(sapply(avail, function(x) all(x %in% scn_pars)))]
      ### if some parameters fixed, remove suggestions without them
      if (isTRUE(length(par_fixed) > 0)) {
        if (isTRUE(length(par_fixed_single > 0))) {
          avail <- avail[which(sapply(avail, function(x)
            all(paste0(par_fixed_single, val_fixed_single) %in% x)))]
        }
      }
      if (isTRUE(length(avail) > 0)) {
        ### load results
        res_add <- lapply(avail, function(x) {
          tmp <- readRDS(file =
            paste0(path_out, paste0(x, collapse = "-"), "--", obj_desc, 
                   file_ext))
          tmp <- tmp@solution[1, ]
          if (is.na(tmp[which("upper_constraint" == names(tmp))])) {
            tmp[which("upper_constraint" == names(tmp))] <- Inf
          }
          return(tmp)
        })
        res_add <- do.call(rbind, res_add)
        res_add <- data.frame(res_add, stringsAsFactors = FALSE)
        cat("adding GA suggestions:\n")
        print(res_add)
        ### add to GA suggestions
        ga_suggestions <- rbind(ga_suggestions, res_add)
        ga_suggestions <- unique(ga_suggestions)
      }
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### run MSE with GA ####
  ### ---------------------------------------------------------------------- ###
  
  ### set random seed for reproducibility
  registerDoRNG(123)
  set.seed(1)
  
  ### run GA
  system.time({
    res <- ga(type = "real-valued", fitness = mp_fitness, inp_file = inp_file,
              obj_fun = obj_fun, obj_yrs = obj_yrs, stat_yrs = stat_yrs, 
              pen_neg = pen_neg, pen_max = pen_max,
              pen_infl = pen_infl, pen_steep = pen_steep,
              path = path_out, check_file = check_file, save_MP = save_MP,
              scenario = scenario, MP = MP, refpts = refpts,
              suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
              names = ga_names, ga_names = ga_names, ga_rounding = ga_rounding,
              maxiter = maxiter, popSize = popSize, run = run,
              monitor = TRUE, keepBest = TRUE, parallel = ga_parallel, seed = 1,
              summarise_runs = TRUE, postFitness = mp_postFitness, 
              optim = ga_optim)
  })
  
  ### save result
  saveRDS(object = res, file = paste0(path_out, scn_pars_c, 
                                      "--", obj_desc, file_ext))
  
  ### ---------------------------------------------------------------------- ###
  ### collate runs ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(collate)) {
    files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                        full.names = FALSE)
    files <- files[grep(x = files, pattern = "--", invert = TRUE)]
    names(files) <- sapply(files, function(x) {
      sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
    })
    scns <- lapply(files, function(x) {#browser()
      pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", 
                              fixed = TRUE),
                          split = "_")[[1]])
      names(pars) <- ga_names
      ### only keep scenarios where requested parameters are changed
      if (!all(ga_default[pos_default] == pars[pos_default])) return(NULL)
      if (isTRUE(length(pos_fixed) > 0)) {
        if (isTRUE(length(val_fixed_single) > 0)) {
          if (!all(val_fixed_single == pars[pos_par_fixed_single])) return(NULL)
        }
      }
      stats <- readRDS(paste0(path_out, x))
      list(pars = pars, stats = stats)
    })
    scns[sapply(scns, is.null)] <- NULL
    saveRDS(scns, 
            file = paste0(path_out, scn_pars_c, "--", obj_desc, "_runs", 
                          ".rds"))
  }
  
### ------------------------------------------------------------------------ ###
### run MP without GA ####
### ------------------------------------------------------------------------ ###
} else {
  
  ### output path
  path_out <- paste0("output/", stock_id, "/", OM, "/",
                     n_iter, "_", n_yrs, "/", scenario, "/", MP, "/")
  dir.create(path_out, recursive = TRUE)
  
  if (identical(MP, "hr")) {
    ### MP parameters
    pars_names <- c("idxB_lag", "idxB_range_3", "exp_b", "comp_b_multiplier",
                    "interval", "multiplier",
                    "upper_constraint", "lower_constraint")
    pars <- mget(pars_names, 
                 ifnotfound = c(1, 1, 1, 1.4, 1, 1, 1.2, 0.7))
    pars <- expand.grid(pars)
    ### additional parameters
    pars_more_names <- c("idx_unc")
    pars_more <- expand.grid(idx_unc = idx_unc,
                             biomass_index = biomass_index,
                             stringsAsFactors = FALSE)
    
    res_stats <- foreach(mp_pars_i = split(pars, f = seq(nrow(pars))), 
                         .combine = bind_rows) %:%
      foreach(pars_more_i = split(pars_more, f = seq(nrow(pars_more))),
              .combine = bind_rows) %do% {
        #browser()
        ### create input
        input_i <- input_mp(stock_id = stock_id, OM = OM, n_iter = n_iter,
                            n_yrs = n_yrs, yr_start = yr_start, 
                            n_blocks = n_blocks, MP = MP, 
                            biomass_index = pars_more_i$biomass_index,
                            idx_unc = pars_more_i$idx_unc)
        ### insert MP parameters
        input_i = mp_pars(params = unlist(mp_pars_i), 
                          params_names = pars_names, 
                          input = input_i, MP = MP)
        ### run MP
        res_i <- do.call(mp, input_i)
        if (isTRUE(save_MP))
          saveRDS(res_i, paste0(path_out, "mp_", 
                                paste0(mp_pars_i, collapse = "_"), "-",
                                paste0(pars_more_i, collapse = "_"),
                                ".rds"))
        ### calculate stats
        stats <- mp_stats(input = input_i, res_mp = res_i, stat_yrs = stat_yrs,
                          refpts = refpts,
                          collapse_correction = TRUE)
        saveRDS(stats, paste0(path_out, 
                              paste0(mp_pars_i, collapse = "_"), "-",
                              paste0(pars_more_i, collapse = "_"),
                              ".rds"))
        ### return stats
        return(bind_cols(mp_pars_i, pars_more_i, stats))
    }
    
    file <- paste0(sapply(as.list(bind_cols(pars, pars_more)), function(x) {
      ifelse(identical(length(unique(x)), 1L),
             unique(x), paste0(min(x), "-", max(x)))
    }), collapse = "_")
    saveRDS(res_stats, file = paste0(path_out, "runs_", file, ".rds"))
    
    
  } else {
  
    if (identical(MP, "constF")) {
      if (identical(Ftrgt, "MSY")) {
        input$ctrl$hcr@args$ftrg <- median(c(refpts["Fmsy"]))
      } else {
        input$ctrl$hcr@args$ftrg <- Ftrgt
      }
    }
    
    res_mp <- do.call(mp, input)
    file_name <- paste0("mp_", Ftrgt)
    if (isTRUE(save_MP))
      saveRDS(res_mp, paste0(path_out, file_name, ".rds"))
    stats <- mp_stats(input = input, res_mp = res_mp, refpts = refpts, 
                      stat_yrs = stat_yrs)
    saveRDS(stats, paste0(path_out, "stats", ".rds"))
    
  }
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

if (!exists(x = "args_local")) {
  quit(save = "no")
}

