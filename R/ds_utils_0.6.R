# This file contains utils for distance sampling. Dave Fifield 2014.

# run a single ds model with specification given in dataframe mod
# 
#
# This Function is used by do.ds() to run many candidate models and by do.det.fcn 
# to re-run final model if needed.
# 
# Params:
# mod - model specification with 4 cols: label, key, adj, form
# data - the ovservation data
# folder - folder to store results. If null, then just return results
# logFileConn - destination for logging messsage. Will either be "" (i.e. stdout)
#   or an open file connection. Note, in the latter case it will not work when 
#   using parallel processing - logging output is lost.
# verbose - provide extra debugging info
# min_data_size_limit - minimum number of observations needed in order to attempt
#   fitting a ddf. If nrow(data) < min_data_size_limit model fitting is not attempted
#   and a failed model dummy file is created.
run.ddf.model <-
  function(mod,
           data,
           folder = "",
           logfileConn = "",
           verbose,
           min_data_size_limit = 20,
           ...) {
    
  if (nrow(mod) > 1) stop("run.ddf.model: given more than 1 model to run!")

  cat(paste("Running", mod$label, "..."), file = logfileConn)
  
  # set up for adjustments
  adj = NULL
  if (mod$adj != "")
    adj <- mod$adj
  
  # Make sure we haven't specified both covars and adjustment terms. Need to 
  # deal with case where formula is either a formula or a string.
  if (!is.null(adj) && mod$form != as.formula("~1") &&
      str_replace_all(mod$form, " ", "") != "~1") {
    mess <- sprintf("run.ddf.model: model (%s) has both covars and adjustment terms!",
                    mod$label)
    stop(mess)
  }
  
  monotonicity <- "none"
  if (mod$form == "~1")
    monotonicity <- "strict"

  model <- NULL
  res <- NULL
  if (nrow(data) >= min_data_size_limit) {
    # fit model - sometimes fitting fails
    res <- try({
      # create a new environment with only the things we want in attempt to keep
      # .rda file small when we eventually use save() which stores everything in the
      # environments in which any formula was created
      env <- new.env(parent = globalenv())
      env$data <- data
      env$mod <- mod
      env$adj <- adj
      env$monotonicity <- monotonicity
      env$dots <- list(...)
      
      if (logfileConn == "")
        model <- with(env, {
          model <- do.call('ds',
                           c(
                             list(
                               data = data,
                               key = mod$key,
                               formula = as.formula(mod$form),
                               monotonicity = monotonicity,
                               adjustment = adj
                             ),
                             dots
                           ))
          return(model) # This just returns outside this block
        })
      else
        capture.output(model <- with(env, {
          model <- do.call('ds',
                           c(
                             list(
                               data = data,
                               key = mod$key,
                               formula = as.formula(mod$form),
                               monotonicity = monotonicity,
                               adjustment = adj
                             ),
                             dots
                           ))
          return(model)
        }),
        file = logfileConn)
    })
  }
  
  # model fitting failed. Create dummy output file to stop model from being
  # attempted again.
  if (inherits(res, "try-error") || is.null(model)) {
    file.create(file.path(folder, paste0("failed_model_", mod$label, ".txt")))
    return(NA)
  }

  cat("finished!\n", file = logfileConn)

  # send summary output to file named "AIC_XXX.XXX_`mod$label.txt`"
  if (folder != "") {
    filename = paste(folder, paste0("AIC_", round(model$ddf$criterion, 3), "_", 
                                    mod$label, ".txt"), sep = "/")
    
    if (verbose)
      cat(paste0("Opening model output file '", filename, "' \n"), file = logfileConn)
    modOutFileConn <- file(filename, open = "wt")
    
    # Create model output file
    if (exists("model")) {
      try.res <- try({
        # first spit out all the things that we want in the model summary
        ch <- ddf.gof(model$ddf, qq = FALSE)$chisquare$chi1
        summ <- summary(model)
        cat(paste0("Date: ", date(), "\n"), file = modOutFileConn)
        cat(paste0("Model name: ", mod$label, "\n"), file = modOutFileConn)
        cat(paste0("Key: ", summ$ds$key, "\n"), file = modOutFileConn)
        cat(paste0("Formula: ", model$ddf$ds$aux$ddfobj$scale$formula),
            "\n",
            file = modOutFileConn)
        cat(paste0("AIC: ",  round(model$ddf$criterion, 3), "\n"), file = modOutFileConn)
        #ddf.gof(model$ddf, qq=FALSE)$dsgof$CvM$p,
        cat(paste0("Chi_scores: ", paste(round((ch$observed - ch$expected) ^
                                                 2 / ch$expected, 3
        ), collapse = ", "), "\n"), file = modOutFileConn)
        cat(paste0("Chisquare p: ",  ch$p, "\n"), file = modOutFileConn)
        cat(paste0("Det prob: ", round(summ$ds$average.p, 3), "\n"), file = modOutFileConn)
        # May not exist for unif key model
        if ("average.p.se" %in% names(summ$ds)) {
          cat(paste0("SE(p): ", round(summ$ds$average.p.se, 5), "\n"), file = modOutFileConn)
          cat(paste0(
            "CV(p): ",
            round(summ$ds$average.p.se / summ$ds$average.p, 3),
            "\n"
          ), file = modOutFileConn)
        }
        cat("\n", file = modOutFileConn)
        capture.output(summary(model), file = modOutFileConn)
      })  
      
      if (inherits(try.res, "try-error")) {
        cat("Creating model summary file failed.\n", file = logfileConn)
      }
      
      # Now save the model object as an .RData file
      save(model, file = sub(".txt", ".RData", filename))
      
    } else {
      cat("Model failed to fit.\n", file = logfileConn)
    }
    close(modOutFileConn)
  }
  
  if (verbose) cat("Returning model object\n", file = logfileConn)
  model
}


# For a given key and vector of covars, generate all possible combinations of n
# items from covars.
gen.form.N <- function(n, key, covars){
  require(utils)

  # if n == 0 then no covars, just key function.
  if (n == 0) {
    return(data.frame(
      form = "~1",
      label = key,
      key = key,
      adj = "",
      stringsAsFactors = F
    ))
  }

  # get all combinations of chooseing n items from covars
  covarList <- combn(covars, n, simplify = F)

  plyr::ldply(covarList, function(item) {
    data.frame(form = paste("~", paste(item, collapse = " + ")), key = key,
               label = paste(key, paste(item, collapse = "."), sep = "."), adj = "", stringsAsFactors = F)
  }
  )
}

# run a series of ds models using all combinations of covariates (including none).
#
# inputs:
#   data - same as ds
#   key - key function(s) to use
#   covars - a character vector of covariate names
#   runModels - for debugging. Should models actually be run or just return the dataframe with a list of what models would be run.
#       Does not do logging.
#   parallel - use parallel processing via doSNOW library to speed things up?
#   nCores - number of cores to use in the parallel cluster
#   folder - the name of a folder where summaries of each model will be written (one per file named
#       modDat$label). This is especially useful for monitoring progress when parallel==TRUE,
#       since the progress bar does not work. If summary files already exist in folder for some models,
#       then these will not be re-run, unless rerun == TRUE.
#   logfile - If provided all logging information will be sent to this file. Otherwise it goes to stdout, which will be lost if parallel == T.
#        Note that all message directly from do.ds() (but not run.ddf.model()) are also sent to stdout.
#   rerun - should models that have already been run and output saved in folder be run again.
#     Only applies when folder == T
#   verbose - enable extra output
#   models - dataframe indicating the models to be run. Possibly created from a
#     previous invocation of do.ds with runModels == F.
#   incl.adj - If TRUE(default), include key+adjustment only models in model set.
## Works by generating a dataframe of models to be run first and then running them.
do.ds <-
  function(data,
           key = c("unif", "hn", "hr"),
           covars = NULL,
           runModels = TRUE,
           parallel = FALSE,
           nCores = 2,
           folder = "",
           logfile = "",
           rerun = FALSE,
           verbose = FALSE,
           models = NULL,
           incl.adj = TRUE,
           cleanFolder = FALSE,
           ...) {
    
  # Models pre-supplied?
  if (!is.null(dim(models)[1]))
    modDat <- models
  else
    # generate models but don't add covars to "unif" key
    modDat <- map_dfr(key, function(k, covars) {
      if (k == "unif")
        gen.form.N(0, k, NULL)
      else
        map_dfr(0:length(covars), gen.form.N, k, covars)
    }, covars = covars)
  

  # append the adjustment only models
  if (incl.adj)
    modDat <- bind_rows(adj.models, modDat)

  # setup logfile connection. If it's "" (the default) then logfileConn will be
  # "", which corresponds to stdout.
  if (logfile != "") {
    logfileConn <- file(logfile, "at") # open for append.
    cat(paste0("Logging info to ", logfile, "\n"))
  } else
    logfileConn <- logfile

  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

  cat(paste0("\n\ndoDS: ", date(), "\n"), file = logfileConn)

  if (parallel & logfile != "" & runModels)
    cat("WARNING: No logging to file for individual sub-processes when parallel == TRUE.\n", file = logfileConn)

  # only run models that currently have no results in folder. Useful when a
  # previous run has failed partway through and we want to avoid re-running a
  # bunch of models for which we already have results.
  if (folder != "") {
    if (runModels) cat(paste0("Saving output files in: ", folder, "\n"), file = logfileConn)

    # Clean the results folder?
    if (cleanFolder) {
      fl <- list.files(path = folder, full.names = TRUE)
      message(sprintf(
        "Cleaning results folder '%s'. Removing %d files.",
        folder,
        length(fl)
      ))
      file.remove(fl)
    }
    
      
    # Get list of models that already have result files
    if (!rerun) {
      fl <- sub(".txt", "", 
                list.files(path = folder, 
                           pattern = ".*\\.txt", 
                           full.names = T), fixed = T)
      existMods <- unlist(lapply(strsplit(basename(fl), "_"), function(x) x[3]))
      cat(paste0(sum(modDat$label %in% existMods), " of ", nrow(modDat), " models have already been run.\n"), file = logfileConn)
      modDat <- modDat[!(modDat$label %in% existMods),]
    }
  }

  if (runModels) {

    cat(paste0("Running ", nrow(modDat), " model(s).\n"), file = logfileConn)
    # Avoids setting up parallel cluster if not needed
    if(nrow(modDat) > 0) {
      if (parallel) {
        cat(paste0("Using parallel processing with ", nCores, " cores.\n"), file = logfileConn)
        cl <- makeCluster(nCores, type = "SOCK")
        registerDoSNOW(cl)
        
        models <- plyr::dlply(
          modDat,
          plyr::.(label),
          .parallel = T,
          .paropts = list(
            .packages = "Distance",
            .verbose = TRUE,
            .errorhandling = "pass" # Don't think this makes any difference
          ),
          run.ddf.model,
          data = data,
          folder = folder,
          logfileConn = "",
          verbose = verbose,
          ...
        )
        stopCluster(cl)
        
        # Tried but had problems with some things not defined...
        # future::plan(multisession, workers = nCores)
        # models <- modDat %>%
        #   split(1:nrow(.)) %>%
        #   furrr::future_map(
        #     run.ddf.model,
        #     data = data,
        #     folder = folder,
        #     logfileConn = logfileConn,
        #     verbose = verbose,
        #     ...
        #   )
      } else { # parallel == FALSE
        models <- modDat %>%
          split(1:nrow(.)) %>%
          map(
            run.ddf.model,
            data = data,
            folder = folder,
            logfileConn = logfileConn,
            verbose = verbose,
            .progress = "Detection Function Progress",
            ...
          )
        
        # models <- plyr::dlply(modDat, .(label), .progress = progress_win(title = "Detection Function Progress"), run.ddf.model, data = data,
        #                       folder = folder,logfileConn = logfileConn, verbose = verbose, ...)
      }
    } else { # nrow(modDat) > 0
      models <- list()
    }

    cat(paste0("do.ds: finishing at ", date(), "\n"), file = logfileConn)

    if (logfile != "")
      close(logfileConn)

    return(models)
  } else { # runmodels == FALSE
    if (logfile != "")
      close(logfileConn)

    return(modDat)
  }
}


# stole this from https://github.com/gsk3/taRifx/blob/master/R/Rfunctions.R#L1161. I'm using it to
# turn NA's into "" so as.numeric wont issue warning about NA's introduced by coercion
destring <- function(x,keep="0-9.-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep = ""),"",x)) )
}

# Check a ddf. Returns original distdata augmented with det. probs and adjSize
check.det.fcn <-
  function(model = NULL,
           species = NULL,
           modname = NULL) {
    
  message("Running check.det.fcn...\n")
  print(summary(model))

  mod.dat <- df.get.data(model) # returns a list
  distdata <- mod.dat$distdata
  form <- mod.dat$form
  key <- mod.dat$key
  aic <- mod.dat$aic
  
  if (nrow(distdata) > 0) {
    # Note that fitted is in the right order since df.get.data uses:
    # fitted(dfobject$ddf)[as.character(distdata$object)] so we don't need to 
    # do that here.
    distdata$detProb <- mod.dat$fitted
    distdata$adjSize <- distdata$size / distdata$detProb
    
    # make sure nothing went wrong
    if (any(is.na(distdata$detProb)))
      warning("Some detection probabilities are NA!", immediate. = T)
    
    cat("\n")
    cat(sprintf("Detection prob range: %s\n", paste(round(
      range(distdata$detProb), 4
    ), collapse = " - ")))
    cat(sprintf("Number of detection probs < 0.15: %d\n", nrow(filter(
      distdata, detProb < .15
    ))))
    cat(sprintf("range of size: %s\n", paste(range(distdata$size), collapse = " - ")))
    cat(sprintf("range of adjusted size: %s\n", paste(round(
      range(distdata$adjSize), 2
    ), collapse = " - ")))
    
    # plot hist of det probs
    p <- ggplot(data = distdata, aes(x = detProb))
    p <- p + geom_histogram(binwidth = 0.1)
    p <- p + scale_x_continuous(breaks = seq(0, 1, .1))
    print(p)
  }
  
  # GOF testing
  message("GOF testing")
  if (!("fake_ddf" %in% class(model))){
    print(ddf.gof(model$ddf, asp = 1))
    # plot model
    plot(model, main = sprintf("%s, %s, %s, %s AIC: %.3f - overall detection function",
                               species, modname, key, form, aic))
  } else
    message("Dummy ddf - no GOF or plot!")

  message("End of check.det.fcn\n")
  distdata
}

# XXXX Might be useful but needs to be updated to use current ddftype column
# implementation instead of "det.fcn.type" left over from subArctic analysis
# 
# find segments with observations with more than one df type
get.trouble <- function(distdata, segdata) {
  distdata %>%
    group_by(Sample.Label) %>%
    summarize(ndet = n_distinct(det.fcn.type)) %>%
    filter(ndet > 1) %>%
    left_join(segdata, by = "Sample.Label") %>%
    left_join(distdata, by = "Sample.Label") %>%
    select(Sample.Label, ndet, det.fcn.type, df.type, distance, size, TransectType, Observer, object) %>%
    as.data.frame()
}


# Extract distance sampling data from a detection function model
df.get.data <- function(dfobject) {
  
  
  if ("dsmodel" %in% class(dfobject)) {
    distdata <- dfobject$ddf$data
    list(
      distdata = distdata,
      form = as.character(dfobject$ddf$call$dsmodel[[2]]$formula)[2],
      key = dfobject$ddf$call$dsmodel[[2]]$key,
      aic = dfobject$ddf$criterion,
      fitted = fitted(dfobject$ddf)[as.character(distdata$object)]
    )
  } else if ("fake_ddf" %in% class(dfobject)) {
    list(
      distdata = dfobject$data,
      form = "~1",
      key = "strip_transect",
      aic = NA,
      fitted = 1 # fitted det prob
    )
  } else if ("list" %in% class(dfobject)) {
    # Not sure what this case is for. Put stop() here to see if it's ever called
    dfobject$distdata
    stop("df.get.data: Don't know how to extract form, key, or aic from this object")
  }
  
  
  # distdata <- df.get.data(model)
  # form <- as.character(model$ddf$call$dsmodel[[2]]$formula)[2]
  # key <- model$ddf$call$dsmodel[[2]]$key
  # aic <- model$ddf$criterion
  
}

# XXXX Could be useful but needs some arg changes - behav probably needs removing??
# and args to check.det.fcn need to be updated.
# 
# check top model for given species, behav, and model det fcn type (e.g.
# "normal"). If list.only is TRUE, just print out the filename of the top model.
check.top.model <- function(fold, species, behav, list.only = FALSE) {
  files <- list.files(fold, ".+\\.RData$", full.names = TRUE)
  
  if (list.only) {
    message(sprintf(
      "Top model for %s %s %s is %s",
      basename(fold),
      species,
      behav,
      basename(files[1])
    ))
    return()
  }
  
  if (files[1] == "") {
    message(sprintf(
      "There were no model results for %s %s %s to load",
      basename(fold),
      species,
      behav
    ))
    return()
  }
  
  message(sprintf("Loading model results from %s", files[1]))
  load(files[1])
  try(check.det.fcn(model = model, paste(species, behav)))
}

# Check all top models for given species and behav
# XXXX Could be useful but needs some arg changes - 
# folder naming and layout has changed....
check.top.models <- function(species, behav, ...){
  folder <- here("R", species, "DF Summaries", behav)

  if (!dir.exists(folder)) {
    message(sprintf("Folder %s does not exist!", folder))
    return()
  }

  dirs <- list.dirs(folder, recursive = FALSE)
  walk(dirs, check.top.model, species, behav, ...)
}


# function to summarize a dsm model
summarize.dsm <- function(model){

  # Check if model had error when it ran  
  if (inherits(model, "try-error")) {
    data.frame(
      response = "Error",
      terms    = NA,
      AIC      = NA,
      REML     = NA,
      OverDisp = NA,
      "Deviance_explained" = NA
    )
  } else {
    summ <- summary(model)
    
    data.frame(
      response = model$family$family,
      terms    = paste(c(
        rownames(summ$pTerms.table), rownames(summ$s.table)
      ), collapse = ", "),
      AIC      = AIC(model),
      REML     = model$gcv.ubre,
      OverDisp = sum(resid(model, type = "pearson")) / model$df.residual,
      "Deviance_explained" = paste0(round(summ$dev.expl * 100, 2), "%")
    )
  }
}

# break a predgrid into chunks and apply func in parallel
par.predict.or.var <- function(this.dsm, func, df, nchunks, off.set = 1, show.progress = FALSE, debug = FALSE){
  browser(expr = debug)

  force(this.dsm)
  force(off.set)
  if (!("multisession" %in% class(plan()))) {
    message("Not in multisession mode! Is this what you wanted?")
  }
  # split data into chunks. base on number of cores?
  if (nchunks > 1) {
    dat.split <- split(df, cut(1:nrow(df), nchunks, FALSE))
  } else{
    dat.split <- list(df)
  }

  #think about what data split should be when func == dsm.var.gam
  # if func is dsm.var.gam, then split each chunk in dat.split into sublists with 1 cell per element
  # should probably move this split into a function that is called in parallel too atht does the split and the dsm.var.gam
  if (identical(func, dsm.var.gam)) {
    dat.split <- map(dat.split, ~ split(.x, 1:nrow(.x)))
  }

  # call function func (either dsm.var.gam or predict) and return either the
  # predictions (fro predict) or a list containing the pred.var and pred
  # (dsm.var.gam). Is predict even necessary if I'm getting the prediction from
  # dsm.var.gam anyway?
  apply.my.func <- function(func, this.dsm, dat, off.set){
    res <- func(this.dsm, dat, off.set)

    if (identical(func, dsm.var.gam)) {
      list(pred.var = res$pred.var, pred = unlist(res$pred))
    } else {
      unlist(res, use.names = F)
    }
  }

  # apply the function to the splits
  print(system.time(res <- future_map(dat.split,
                                      ~ apply.my.func(func, this.dsm, .x, off.set = off.set),
                                      .progress = show.progress,
                                      .options = future_options(globals = c("this.dsm"),
                                                                packages = c("Distance", "dsm")))))
  if (identical(func, predict)) {
    res <- unlist(res)
  }
  res
}


# XXX Not currently used
# Given a DSM and a dataframe to predict to, produce a density estimate with
# measures of uncertainty. Previously Used dsm.var.gam to get an ABUNDANCE
# estimate and then turn that into density estimate. Now just gets density
# estimate directly bu supplying 1 as the off.set
get.dens.est <- function(dsm_final, predgrid) {

  # use dsm.var.gam to get estimated abundance params
  densEst <- summary(dsm.var.gam(dsm_final, predgrid, off.set = predgridCellArea))

  # The estimates in dsm.var.gam are summed for the entire predgrid study area
  # so we need to divide by the number of predgrid cells * cellarea.
  densEst %<>%
    map_at(c("pred.est", "se"), ~ .x/(predgridCellArea * nrow(predgrid)))

  #calculate a lognormal CI for the density est.
  cv.square <- densEst$cv^2
  asymp.ci.c.term <- exp(1.96*sqrt(log(1 + cv.square))) # stolen from print.summary.dsm.var.R. Also in buckland et al 2001.
  list(pred.est = densEst$pred.est, CV = densEst$cv, SE = densEst$se, CI = c("%5" = densEst$pred.est/asymp.ci.c.term,
                                                                             "Mean" = densEst$pred.est, "95%" = densEst$pred.est * asymp.ci.c.term))
}

# XXX not currently used
# pretty-print the output from get.dens.est
print.dens.est <- function(densEst){
  cat("Density estimate:\n\n")

  cat("Approximate asymptotic confidence interval:\n")
  print(densEst$CI)

  cat("\n")
  cat("Point estimate                 :", densEst$pred.est,"\n")
  cat("Standard error                 :", densEst$SE,"\n")
  cat("Coefficient of Variation       :", densEst$CV,"\n")
}


# Render one of the Generic_x_xxx.Rmd files passing species as a param
# and save resulting knitted output.
# 
# Eg. do.generic.render("ATPU", "Generic_2_dsm.Rmd")
# 
# If parallel == TRUE then we are being called from in parallel with other ongoing
# renders() which would cause problems overwriting intermediate files, so we 
# create a separate folder for each render.
do.generic.render <- function(species, file, parallel = FALSE){
  
  suffix <- str_replace(file, "^Generic", "") %>%
    str_replace("Rmd$", "html")
  out.file <- file.path(ResultsDir, species, paste0(species, suffix))
  
  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)

  if (parallel == TRUE) {
    tf <- tempfile()
    dir.create(tf)
    
    message(
      sprintf(
        "Rendering generic file %s for %s to %s using intermediate dir %s",
        file,
        species,
        out.file,
        tf
      )
    )
    
    render(
      file.path(RDir, file),
      params = list(species = species),
      output_file = out.file,
      intermediates_dir = tf
    )
    
    unlink(tf)
  } else {
    message(sprintf("Rendering generic file %s for %s to %s", file, species, out.file))
    render(
      file.path(RDir, file),
      params = list(species = species),
      output_file = out.file
    )
  }
}


# Render the extrapolation analysis Rmd for one dataset
do.extrapolation <- function(spill, dataset, debug = FALSE){
  browser(expr = debug)

  out.file <- file.path(here(), paste(spill, dataset, "0_extrapolation.html", sep = "_"))
  message(sprintf("Rendering generic extrapolation  for %s to %s", dataset, out.file))
  render(file.path(here(), "Generic_0_extrapolation.Rmd"),
                  params = list(spill = spill, dataset = dataset),
                  output_file = out.file)
}


# fitted dsm model checking
check.dsm <- function(dsm_final, modname, segdata) {

  # accessing columns with segdata[, termlab] below doesn't work with tbls or
  # sf objects.
  segdata <- as.data.frame(segdata)
  
  # Several dsm-package functions do partial matching so turn off and re-enable
  # at end.
  options(warnPartialMatchDollar = FALSE)

  message("\n===========================================================================\n\nDSM model summary for model ", modname)
  print(summary(dsm_final))

  # Smoother plots
  message("Smoother plots")
  par(mfrow = c(1,1))
  # XXXX Residual plotting is not working b/c of some mess-up with the name of the 
  # offset column in dsm. need to debug some more
  # gratia::draw(dsm_final, wrap = FALSE, residuals = dsm.options$do.dsm.residuals) %>% 
  gratia::draw(dsm_final, wrap = FALSE) %>% 
    print
  
  if (any(grepl("s(x.sc, y.sc", as.character(dsm_final$formula), fixed = T))) {
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), main = "s(x.sc,y.sc) (response scale)",
            type = "response", asp = 1, plot.type = "contour")
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 0, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 60, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = -60, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
  }
  
  if (any(grepl("s(x, y", as.character(dsm_final$formula), fixed = T))) {
    vis.gam(dsm_final,  view = c("x","y"), main = "s(x.sc,y.sc) (response scale)",
            type = "response", asp = 1, plot.type = "contour")
    vis.gam(dsm_final,  view = c("x","y"), theta = 0, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x","y"), theta = 60, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x","y"), theta = -60, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
  }
  
  
  # Sometimes whines about S3 methods.
  message("Gratia checks")
  suppressWarnings(gratia::appraise(dsm_final))
  
  # Remove "dsm" class to make DHARMa happy
  message("DHARMa checks")
  simmod <- dsm_final
  class(simmod) <- class(simmod)[-1] # Simulate residuals doesn't like dsm class
  sims <- simulateResiduals(fittedModel = simmod)
  plot(sims)
  testResiduals(sims)
  testZeroInflation(sims)
  # May fail if x,y locations are not unique
  res <- try(testSpatialAutocorrelation(sims, segdata$x, segdata$y))
  
  # If more than one resid at a given location.
  # Note still use try() since this may fail to allocate enough memory if 
  # size of data is too big.
  if(inherits(res, "try-error")) {
    message("testSpatialAUtocorrelation failed: aggregating resids spatially")
    
    # This always fails trying to allocate more than 256GB (the RAM I have so
    # let's not bother)
    # segdata$loc <- paste(as.character(segdata$x), as.character(segdata$y), sep = "_")
    # recal <- recalculateResiduals(sims, group = segdata$loc)
    # locs <- segdata %>% 
    #   distinct(loc, .keep_all = TRUE)
    # try(testSpatialAutocorrelation(recal, locs$x, locs$y))
    # 
    # # May fail due to not enough memory for big jobs
    # if(inherits(res, "try-error")) {
    #   message("Spatial autocorrecation test failed for aggregated data. Garbage collecting...")
    #   gc()
    # }
  }

  # plot resids vs each term in model
  attr(terms(simmod), "term.labels") %>%
    walk(function(termlab, sims, segdata) {
      plotResiduals(sims, segdata[, termlab], xlab = termlab)
    }, sims = sims, segdata = segdata)

  # Concurvity
  message("Concurvity checks\nEach term with whole of rest of model")
  try(print(concurvity(dsm_final) %>% round(digits = 3)))
  message("Concurvity of pairwise terms ('estimate' measure presented)\nEach row shows how terms in columns depend on the term in that row.")
  try(print(concurvity(dsm_final, full = FALSE)[["estimate"]] %>% round(digits = 3)))
  message("Plot is non-symmetric, showing how terms on y-axis depend on terms on the x-axis")
  try(vis_concurvity(dsm_final))


  # Gam checks from MGCV
  par(mfrow = c(1,1))
  message("MGCV checks")
  try(my.gam.check(dsm_final))
  rqgam_check(dsm_final)

  par(mfrow = c(1,1))
  # check observed vs expected. See Miller et al 2021 pg 11 (of 18) for 
  # reccommendation to use the "platform" variable (ddftype here) to aggregate
  # by.
  message("Observed vs expected plot")
  oe.dens(dsm_final, covar = "ddftype", plotit = T)
  oe.dens(dsm_final, covar = "depth", plotit = T)
  # oe.dens(dsm_final, covar = "depth.g", plotit = T)
  # oe.dens(dsm_final, covar = "sst", plotit = T)
  # oe.dens(dsm_final, covar = "sst.g", plotit = T)
  # # oe.dens(dsm_final, covar = "year", plotit = T)
  # # oe.dens(dsm_final, covar = "yday", plotit = T)

  # check zero infl - not sure this is right
  (OD_dsm_final <- sum(resid(dsm_final, type = "pearson"))/dsm_final$df.res)

  # check autocorellogram
  # create transect label as cruiseid & date & flyswim in order to avoid
  # having the duplicated segments (one for fly one for swim) together in the 
  # same "Transect".
  # create segment label as %H:%M:%S
  dsm_final$data <- dsm_final$data %>% 
    mutate(
      tr.lab = paste(dsm_final$data$CruiseID,
                     dsm_final$data$Date,
                     dsm_final$data$FlySwim,
                     sep = "_"),
      seg.lab = format(as_datetime(dsm_final$data$StartTime), "%H:%M:%S")
    )
  par(mfrow = c(1,1))
  dsm_cor(
    dsm_final,
    Transect.Label = "tr.lab",
    Segment.Label = "seg.lab",
    max.lag = 20
  )

  # check variogram
  #
  mydata <-
    data.frame(
      E = resid(dsm_final, type = "pearson"),
      x = segdata$x / 1000,
      y = segdata$y / 1000
    )
  coordinates(mydata) <- ~ x + y
  
  print(bubble(
    mydata,
    "E",
    col = c("black", "red"),
    main = "Residuals",
    xlab = "X-coords",
    ylab = "y-coords"
  ))
  
  V <- (gstat::variogram(E ~ 1, mydata))
  plot(
    x = V$dist,
    y = V$gamma,
    xlab = "Distance (km)",
    ylab = "Semi-variance",
    pch = 16,
    cex = 2 * V$np / max(V$np)
  )
  
  V <- (gstat::variogram(E ~ 1, mydata, cutoff = 100))
  plot(
    x = V$dist,
    y = V$gamma,
    xlab = "Distance (km)",
    ylab = "Semi-variance",
    pch = 16,
    cex = 2 * V$np / max(V$np)
  )
  
  
  V <- (gstat::variogram(E ~ 1, mydata, cutoff = 10))
  plot(
    x = V$dist,
    y = V$gamma,
    xlab = "Distance (km)",
    ylab = "Semi-variance",
    pch = 16,
    cex = 2 * V$np / max(V$np)
  )
  
  V <- (gstat::variogram(E ~ 1, mydata, cutoff = 5))
  plot(
    x = V$dist,
    y = V$gamma,
    xlab = "Distance (km)",
    ylab = "Semi-variance",
    pch = 16,
    cex = 2 * V$np / max(V$np)
  )
  
  V <- (gstat::variogram(E ~ 1, mydata, cutoff = 2))
  plot(
    x = V$dist,
    y = V$gamma,
    xlab = "Distance (km)",
    ylab = "Semi-variance",
    pch = 16,
    cex = 2 * V$np / max(V$np)
  )
  
  # Turn partial match warnings back on
  options(warnPartialMatchDollar = TRUE)
}



# Augment segdata with summed counts of observations and create zeros for
# segments where there were no observations of the species of interest.
#
# Also add fields for rawCount, density-corrected estAbund, and estDens from a
# given set of observation distdata. These are used to check if my computed
# response is equal to that from dsm, after adjusting for one-sided transects
# with convert.units. This is used in run.dsm.model() to ensure that we are
# building segdata properly, and that convert.units is doing what we think it
# is.
#
# Note that if you want to restrict the analysis to a specific species/group,
# then distdata should already have been filtered for the species of interest
# before calling this function so that this code will fill in the zeros
# properly. This happens in Generic_2_dsm.Rmd in the usual case.
augment.segdata <- function(segdata, distdata) {
  newsegdata <- distdata %>%
    group_by(Sample.Label) %>%
    summarize(estAbund = sum(adjSize), rawCount = sum(size)) %>%
    right_join(segdata, by = "Sample.Label") %>% 
    st_as_sf
  newsegdata$estAbund[is.na(newsegdata$estAbund)] <- 0
  newsegdata$rawCount[is.na(newsegdata$rawCount)] <- 0
  newsegdata$estDens <- newsegdata$estAbund/newsegdata$segment.area
  newsegdata
}


# add 'k =  xxx' modifier to smooth terms for year or yday if there are less
# than k (default 10) unique values of that variable in segdata.
#
# if the resulting k would be < 2 then just remove the term
#
# XXX This should be made generic to go through all model terms. For each term
# just convert the formula to character and use string functions to grep and
# replace the term as needed instead of trying to use update().
#
adjust.time.covars <- function(form, segdata, k = 10) {
  # check if form contains year or yday
  # check unique number of values for each
  # if either is less than default 10 then

  # check for year
  new.k <- length(unique(segdata$year))
  if (grepl("s(year", as.character(form)[3], fixed = TRUE) &&
      (new.k < k)) {
    message("Removing year from formula for lack of data")
    form <- update(form, . ~ . - s(year, bs = "ts"))

    # enough data to have term at all?
    if (new.k >= 2) {
      message("\t... and replacing it with smooth where k = ", new.k)
      form %<>% as.character
      form <- sprintf('%s ~ %s + s(year, bs = "ts", k = %d)', form[2],
                        form[3], new.k) %>%
        as.formula
    }
  }

  # check for yday
  new.k <- length(unique(segdata$yday))
  if (grepl("s(year", as.character(form)[3], fixed = TRUE) &&
      (new.k < k)) {
    message("Removing yday from formula for lack of data")
    form <- update(form, . ~ . - s(yday, bs = "ts"))

    # enough data to have term at all?
    if (new.k >= 2) {
      message("\t... and replacing it with smooth where k = ", new.k)
      form %<>% as.character
      form <- sprintf('%s ~ %s + s(yday, bs = "ts", k = %d)', form[2],
                        form[3], new.k) %>%
        as.formula
    }
  }

  form
}

# Run a dsm model.
#
# mod.def = dataframe with elements modname, formula, and family (the latter is a list columns).
# folder = path to folder to save results in.
# ddf.obj = is the list of ddfs in the correct order
# rerun.dsms = If TRUE (default) then rerun dsm() for mod.def, otherwise
# if .Rdata file exists from previous run of this model then load results from 
# that file.
run.dsm.model <- function(mod.def,
                          ddf.obj,
                          segment.data,
                          observation.data,
                          method = "REML",
                          convert.units = 1,
                          control = list(keepData = T),
                          folder,
                          rerun.dsms = TRUE,
                          ...) {
  
  # needed if being called from future_map() on a worker process. If not,
  # it messes up trying to access any sf object (segment.data) when using 
  # the s2 spherical geometry package.
  sf_use_s2(FALSE) 
  
  filename <- file.path(folder, paste0(mod.def$modname, ".Rdata"))
  
  # Rerun dsm?
  if (rerun.dsms == FALSE && file.exists(filename)) {
    # Load saved model result from file
    message(sprintf(
      "Loading previously saved model results for %s.",
      mod.def$modname
    ))
    load(filename)
  } else {
    # Rerun the dsm
    message("Running dsm model ", mod.def$modname)
    
    # Adjust formula
    form <- adjust.time.covars(mod.def$formula[[1]], segment.data)

    # Call dsm(). Note the list indexing for formula and family since these are
    # list columns in mod.def. Temporarily disable warnings about partial
    # matches of list/dataframe names since mgcv is full of such things.
    options(warnPartialMatchDollar = FALSE)
    print(system.time(model <- try(dsm(
      formula = form,
      ddf.obj = ddf.obj,
      segment.data = segment.data,
      observation.data = observation.data,
      family = mod.def$family[[1]],
      method = method,
      control = control,
      convert.units = convert.units,
      segment.area = segment.data$segment.area,
      ...
    ),
    outFile = stdout())))
    options(warnPartialMatchDollar = TRUE)
    
    
    if(!inherits(model, "try-error")){
      # if we used "bam" then the data is not kept even if keepData == TRUE, so
      # add it back in
      kd <- control[["keepData"]]
      if (!is.null(kd) & isTRUE(kd) & is.na(model$data)){
        model$data <-
          dsm:::make.data(
            response = as.character(mod.def$formula[[1]])[2],
            segdata = segment.data,
            obsdata = observation.data,
            ddfobject = ddf.obj,
            family = mod.def$family[[1]],
            group = FALSE,
            convert.units = 1, 
            availability = 1,
            segment.area = segment.data$segment.area
          )
      }
  
      # check if my computed response is equal to that from dsm, just to make sure
      # we understand how response is being computed. Note the "as.data.frame"
      # is used to remove the atribs from segment.data b/c it's an sf object.
      stopifnot(
        all.equal(
          model$data %>% arrange(Sample.Label) %>% magrittr::extract("estAbund"),
          segment.data %>% arrange(Sample.Label) %>% as.data.frame %>% magrittr::extract("estAbund")
        )
      )
    } else { # Model failed to fit. Give message and then continue on to save.
      message("run.dsm.model: dsm() failed: ", model)
    } 
    
    message(sprintf("Saving model result to %s", filename))
    save(model, file = filename)
  } # End rerun dsm
  
  model
}


# Calculate predictions for a dsm model
#
# modname = name of model to extract from mod.res and predict for
# mod.res = list of dsm model objects
# species = name of species this is for
# predgrid = prediction grid containing one row per cell (polygons) in each of 4 
#     seasons. Must contain all predictors used in the model indicated by modname 
#     (excluding offset which is calculated internally by predict.)
#
# Initially predgrid contains 400240 rows, with one copy of the cells for 
# each combination of season and flyswim (The spatial extent of the study area 
# contains 50030 cells):
# 
#
#               F     W
#     Fall   50030 50030
#     Spring 50030 50030
#     Summer 50030 50030
#     Winter 50030 50030
#
# Value: predgrid augmented with NHat and Density
# 
# At the end ret will contain:
#   - predictions for fly (200120 rows - one for each spatial cell in each season)
#   - predictions for swim (200120 rows - one for each spatial cell in each season)
#   - predictions for combined (swim+fly) (200120 rows - one for each spatial cell in each season)
#   
# 
# NOTE NOTE NOTE if you do not supply newdata arg (ie. predgrid) but rather just
# predict to the same data that was used to fit the model then the "...order of
# the results will not necessarily be the same as the segdata (segment data)
# data.frame that was supplied (it will be sorted by the Segment.Label field)."
# (From predict.dsm manual)
dsm.pred <-
  function(modname,
           species,
           mod.res,
           predgrid) {
    
  message(sprintf("Predicting %s for dsm model %s", species, 
                  modname))

  # Extract model and mark subset inpredgrid with FlySwim   
  model <- mod.res[[modname]]
  ret <- predgrid %>% 
    mutate(subset = FlySwim)
  
  # Predict from model for all seasons and fly or swim
  ret <- ret %>% 
    mutate(NHat = predict(model, newdata=ret, off.set=ret$area),
           Dens = NHat/area)
  
  # Add in rows for combined fly+swim 
  ret <- ret %>% 
    head(nrow(.) / 2) %>% # Peel off first half for template
    mutate(subset = "Combined",
           NHat = head(ret, nrow(predgrid) / 2)$NHat + # flying +
                  tail(ret, nrow(predgrid) / 2)$NHat,  # swimming
           Dens = NHat/area,
           FlySwim = NA) %>% 
    rbind(ret)   # tack on original fly and swim rows
  
  # if (do.plots) {
  #   # Plot combined result
  #   combined_plot <- ggplot() +
  #     geom_sf(
  #       data = filter(ret, subset == "Combined"),
  #       mapping = aes(colour = Dens, fill = Dens)
  #     ) +
  #     labs(x = "", y = "", fill = "Dens") +
  #     ggtitle("Fly+Swim Combined") +
  #     theme_minimal() +
  #     scale_colour_viridis_c(option = "E") +
  #     scale_fill_viridis_c(option = "E")
  #   print(combined_plot)
  #   
  #   # Plot individual fly and swim results
  #   ind_plot <- ggplot() +
  #     geom_sf(
  #       data = filter(ret, subset != "Combined"),
  #       mapping = aes(colour = Dens, fill = Dens)
  #     ) +
  #     labs(x = "", y = "", fill = "Dens") +
  #     facet_wrap(vars(FlySwim)) +
  #     theme_minimal() +
  #     scale_colour_viridis_c(option = "E") +
  #     scale_fill_viridis_c(option = "E")
  #   print(ind_plot)
  # }  
  
  # Save as a shapefile
  ret %>% 
    filter(subset == "Combined") %>% 
    st_write(
      dsn = ShapeDir,
      layer = paste(species, modname, "predictions", sep = "_"),
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
  

    
  # create 4-lyr seasonal raster of fly+swim combined values for plotting and
  # saving
  pp_raster <- season.names %>%
    map(make.season.raster,
        obj = filter(ret, subset == "Combined"),
        variable = "Dens") %>% 
    rast
  names(pp_raster) <- season.names
  
  # produce 4 panel plot. Not necessary since this is done at end of Generic_3_prediction.rmd
  # plot(pp_raster, main = paste(season.names, paste(species, modname, sep = "_")))
  
  # Save as a raster
  writeRaster(pp_raster,
              file.path(
                predDir,
                sprintf(
                  "%s.%s.%s.%d_sqkm.img",
                  species,
                  season.names,
                  modname,
                  predgridCellArea
                )
              ),
              datatype = "FLT4S",
              overwrite = T)
  ret
}



# render ddf 
do.det.fcn.render <- function(species,
                              do.final,
                              rerun,
                              parallel,
                              nCores,
                              do.eda,
                              cleanFolder) {
  
  
  if (do.final)
    suffix <- "01_final_ddf.html"
  else
    suffix <- "01_candidate_ddf.html"
  
  out.file <- file.path(ResultsDir, species, paste(species, suffix, sep = "_"))
  message(sprintf("Rendering generic ddf fitting for %s to %s", species, out.file))
  
  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)
  
  render(
    file.path(here("R"), "Generic_1_ddf_fitting.rmd"),
    params = list(
      species = species,
      do.final = do.final,
      rerun = rerun,
      parallel = parallel,
      nCores = nCores,
      do.eda = do.eda,
      cleanFolder = cleanFolder
    ),
    output_file = out.file
  )
}

# produce a summary table of candidate det. fcns in one folder
candidate.detfcn.summary <- function(df.spec.name, species) {
  
  folder <- file.path(ResultsDir, species, paste("DF Summaries", df.spec.name, sep = "_"))
  do.ds.det.fcn.checks(folder, species = species, dsetname = df.spec.name)
  
  get.ds.res(folder) %>%
    mutate(dataset = df.spec.name, species = species) %>%
    relocate(species, dataset, Model, Key, Formula, DetProb, AIC, deltaAIC, )
}


# Check all .txt files in folder for DS summary and return a dataframe with
# model name and summary items.
get.ds.res <- function(folder) {
  fl <- list.files(path = folder, pattern = "AIC.*\\.txt", full.names = T)
  
  # If all ddf failed there will be no results so we can't use "order"
  # to sort by aic, so just return the empty dataframe
  if (length(fl) == 0) {
    data.frame(
      Model = NA,
      Key  = NA,
      Formula = NA,
      AIC  = NA,
      ChiScores  = NA,
      ChisquareP  = NA,
      DetProb = NA,
      DetSE = NA,
      DetCV = NA,
      NObs = NA,
      NCov = NA,
      NCovSE = NA,
      NCovCV = NA,
      deltaAIC = NA
    )
  } else {
    # XXXX replace with dplyr
    res <- plyr::ldply(fl, get.stats.file)
    res$deltaAIC <- res$AIC - min(res$AIC, na.rm = T)
    res[order(res$AIC), ]
  }
}

#Load each detfcn model in folder in turn and plot it
do.ds.det.fcn.checks <- function(folder, species, dsetname) {
  fl <- list.files(path = folder, pattern = "AIC.*\\.RData", full.names = T)
  
  walk(fl, function(filename) {
    message(sprintf("Checking detection function: %s", filename))
    load(filename)
    check.det.fcn(model, 
                  species = species, 
                  str_replace(basename(folder), fixed("DF Summaries_"), ""))
  })
}


# Get various values from a file containing a dsmodel object summary.
#folder <- "C:/Users/fifieldd/Documents/Offline/R/DS Utils/Test Output"
#path <- paste(folder, "hn.size 19-Jan-2016_145149.txt", sep = "/")
get.stats.file <- function(path) {
  message("Getting stats file ", path)
  
  lines <- readLines(path)
  
  getStat <- function(lines, string){
    if (length(line <- grep(string, lines)) == 0) {
      warning(path, ": No line containing ", string, "!", immediate. = T)
      return("Not available")
    }
    strsplit(lines[line], string)[[1]][2]
  }
  
  res <- data.frame(Model = getStat(lines, "^Model name: "),
                    Key  = getStat(lines, "Key: "),
                    Formula  = getStat(lines, "Formula: "),
                    AIC  = as.numeric(getStat(lines, "AIC: ")),
                    ChiScores  = getStat(lines, "Chi_scores: "),
                    ChisquareP  = as.numeric(destring(getStat(lines, "Chisquare p: "))),
                    DetProb = as.numeric(destring(getStat(lines, "Det prob: "))),
                    DetSE = as.numeric(destring(getStat(lines, "^SE\\(p\\): "))),
                    DetCV = as.numeric(destring(getStat(lines, "CV\\(p\\): "))),
                    NObs = as.numeric(destring(getStat(lines, "Number of observations :"))),
                    NCov = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][1])),
                    NCovSE = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][2])),
                    NCovCV = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][3])),
                    stringsAsFactors = F)
  res
}


# Takes a ddf spec and adds to $fitted.distdata the following:
# 1) ddftype as a suffix to the Sample.Labels in fitted.distdata.
# 2) ddfobj.orig as an integer version of the ddftype factor.
# 
# Used by create.dsm.data to make Sample.Labels unique in each segdata copy, and
# to assign the initial ddfobj needed by dsm(). Note the final ddfobj value must
# be sequential with no missing values and will be assigned in create.dsm.data()
# once we know which ddfs are included in the dsm (ie. ones with no obs are not
# included).
augment.distdata <- function(df.mod.spec) {
  df.mod.spec$fitted.distdata <- df.mod.spec$fitted.distdata %>%
    mutate(Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
           ddfobj.orig = as.integer(df.mod.spec$ddftype))
  df.mod.spec
}  


# Create a copy of the appropriate segdata (according to surveytype -
# either aerial or ship) for a df.mod.spec. Set the Sample.Label suffix and 
# ddftype accordingly.
create.segdata.copy <- function(df.mod.spec, init_segdata) {
  segdata <- switch(
    substr(df.mod.spec$ddftype, 1, 1),
    S = filter(init_segdata, SurveyType == "Ship"),
    A = filter(init_segdata, SurveyType == "Aerial")
  ) %>%
    mutate(
      Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
      ddftype = df.mod.spec$ddftype,
    )
  segdata
}

# Nov 15 2024
#   - ddgobj may have less than the full number of values (currently 6) if some
# ddfs ("platforms" in Miller et al. (2021) parlance) have no obs
#   - nonetheless, there will be segdata copies for each unique "platform", 6
# in this case as below. so need to creata a factor variable in segdata to hold
# this info.
#   - move the basic set of ddfs to analysis_settings.R
#   - get rid of dummy epithet
#   - rename behav to something else??
#   - Sample.Label will have a W (water) , F (fly), or S (strip) suffix XXXX OR 
#   maybe it needs six different levels - one for each segdata copy????   AW, AF,
#   AS, SW, SF, SS (Aerial water, aerial fly, aerial strip, etc)? Go with this for now.
#
#
# 1. Grab all distdata, set ddf.obj, remove empty distdatas, assign Sample.Label suffix
# 2. setup segdata copies with dfd.obj and ddftype (or whatever it will be called), and
#    Sample.Label suffix. How much code will break if I make it two characters?
#
#
#
# Function to create distdata, segdata, and ddfs for a given species to pass to
# dsm().
#
# Segdata needs to contain (even if there are no observations in some of these
# classes):
#
# 1. ship segdata for water birds with perp distances - ddf.obj == 1 
# 2. ship segdata for flying birds with perp distances - ddf.obj == 2 
# 3. aerial segdata for water birds with perp distances - ddf.obj == 3 
# 4. aerial segdata for flying birds with perp distances - ddf.obj == 4 
# 5. ship segdata for birds with NO perp distances - ddf.obj == 5 
# 6. aerial segdata for birds with NO perp distances - ddf.obj == 6
#
# Instead of hardcoding the number and nature of these, I should probably do it
# dynamically based on df.mod.specs.
#
# We need to explicitly calculate the segment.area for each segment b/c aerial
# segs are (usually) two sided and ship ones are not. Normally, this is handled
# by the convert_units= arg to dsm, but it only takes a single scalar value
# whereas we need a different one for each type of segment (aerial vs ship).
# With this setup there is no need for a convert_units= arg to dsm anymore. Note
# that under normal circumstances it would be fine to just divide the ship-based
# effort by 2 and let dsm:::make.data() do the calc of segment.area but it isn't
# that simple when one of the four "platforms" has no observations. In that
# case, there is currently no way to make a dummy_ddf with no obs and so you
# have to arbitrarily assign a ddf.obj to the segments in that platform. If the
# transect geometry differs between these two then make.data() cannot properly
# calculate the segment.area. So I do it here by hand and pass the segment.area
# to dsm(). Note that this calculation can't be done in Create_seg_data.Rmd b/c
# we don't have the ddf info at that point.
#
# Note that perhaps the "right" solution would be for dsm to work properly with
# to dummy_ddfs that have no observations.
#
# Segment area is:
#
# For two-sided transects: Effort * (width - left) * 2 
# - for one-sided transects: Effort * (width - left)
#
# Width is taken from the meta-data for the ddf.
#
# Arguments: 
# species - char string. Needed to chose the correct season for each
#  segment b/c seasons can be species specific. 
# df.mod.specs - list of detection function specifications - 
# init_segdata - the initial segdata computed by Create_seg_data.Rmd. This is
#  subdivided into aerial and ship.
#
# Value: 
#  Returns a list of 3 elements: 
#  
#  list(distdata, segdata, ddfs) 
#  
# distdata will be augmented with ddfobj, and have appropriate suffix added to
#   Sample.Label (see segdata below)
#
# segdata will consist of multiple copies of the appropriate (aerial or ship) 
#   segdata (one for each "platform" aka observation process) augmented with 
#   variables for:
#   segment.area, 
#   ddfobj - tells dsm:::make_data() which det prob to use to modify the response, 
#   platform.factor - ,
#   Season
# Sample.Label will have a W (water) , F (fly), or S (strip) suffix XXXX OR 
#   maybe it needs six different levels - one for each segdata copy????   
#
# ddfs - a list containing up to n detection functions (eg, Ship_fly, ship_swim,
# aerial_fly, aerial_swim, ship_strip, aerial_strip) excluding any for which
# there were no obs.
create.dsm.data <- function(species, df.mod.specs, init_segdata) {
  
  #### Create distdata ---------------------
    
  # Create the observation data that will be passed to dsm(). Note that some
  # of these distdata components may have no obs and will thus not be included
  distdata <- df.mod.specs %>% 
    map(augment.distdata) %>% 
    # Extract list of fitted.distdatas
    map("fitted.distdata") %>% 
    # Remove datasets with no obs - note use of base R Filter()
    Filter(function(x) nrow(x) > 0, .) %>%
    map_dfr( ~ .) %>%  # Convert list to single dataframe
    mutate(
      # Need to renumber the ddfobj values sequentially from 1 to the number
      # of ddfs, whereas they may currently have gaps in the numeric sequence.
      # For example, if there were no obs in aerial_F or aerial_S then the
      # ddfobj.orig would be numbered 1, 2, 3, 5 but should be 1:4.
      ddfobj = dense_rank(ddfobj.orig)
    )

  # Make sure we still have some obs in distdata and all rows  are assigned a
  # ddfobj.
  stopifnot(nrow(distdata) > 0)
  stopifnot(sum(table(distdata$ddfobj)) == nrow(distdata))  

  #### Create segdata --------------------
  
  # Cylce through the list of ddf model specs adding a copy of segdata for each.
  segdata <- df.mod.specs %>% 
    # Insert a copy of the right segdata (ship vs aerial) for each ddf spec.
    # regardless of whether there were any obs in that ddf. If there weren't
    # then all segments for this ddf spec will end up with zero obs in the 
    # response variable created by dsm:::make_data()
    map(create.segdata.copy, init_segdata = init_segdata) %>% 
    map_dfr(~ .) %>% # Convert to one large dataframe
    assign.season(seasons[[species]], datefield = "Date") %>% 
    select(
      SurveyType,
      Sample.Label,
      CruiseID,
      Program,
      ObserverName,
      DistMeth,
      CalcDurMin,
      ObsHeight,
      TransectID,
      TransectSides,
      Dataset,
      year,
      yday,
      Date,
      StartTime,
      EndTime,
      LatStart,
      LatEnd,
      LongStart,
      LongEnd,
      x,
      y,
      depth,
      depth.g,
      sst,
      sst.g,
      Effort,
      segment.area,
      ddftype, # based on platform (ship/aerial) and behav (F, W, S(trip))
      Season,
      ends_with(".sc")
    ) %>% 
    # Remove rows with any NA values in columns used as dsm covars
    filter(!(rowSums(is.na(select(., all_of(dsm.covars.ns)))) > 0)) %>% 
    
    # reset the ddfobj for any segments which were originally associated with a
    # ddf that had no observations. Otherwise, dsm will get upset. Just set it
    # to any other ddfobj that actually exists.The abundance.est in these
    # segments (as computed by dsm:::make.data()) will still be 0 since none of
    # the obs in distdata will match that sample.label in segdata (XXXX But
    # doesn't this bias the gam????). Note it might be better just to drop these
    # rows from segdata altogether, but that would mean that the number of
    # init.segdata copies in the final segdata would vary dynamically and this
    # would need to be kept track of so that downstream code (ie. in prediction
    # step) would know how to parse out the rows properly.
    mutate(ddfobj = case_when(
      as.integer(ddftype) %in% unique(distdata$ddfobj) ~ as.integer(ddftype),
      TRUE ~ unique(distdata$ddfobj)[1]
    )) %>% 
    # dsm() uses observations in segdata and associated detection function to
    # calculate the gam response internally. Calculate it here so I can check it's
    # done right and I understand what's going on internally in dsm(). 
    augment.segdata(distdata)
    
  
  ### Create ddfs ----------------------------------------
  
  # Extract ddfs from df.mod.specs and drop any which don't have any observations
  # or else dsm() will get upset. Note we use the original ddfobj numbering
  # before it was recoded since this extracts the correct ddfs.
  ddfs <- map(df.mod.specs, "fitted.model") %>% 
    magrittr::extract(unique(distdata$ddfobj.orig))

  list(distdata = distdata, segdata = segdata, ddfs = ddfs)
}


# Get a dataframe listing all ddf models that currently have no results files
# in folder. List of models is constructed from all candidate keys, adjustment 
# terms and covars. Useful for figuring out what ddf model fitting went off into
# la-la land...
list.unfinished.ddfs <- function(folder, covars = dfCovars){
  do.ds(
    data = NULL,
    folder = folder,
    runModels = FALSE,
    covars = covars
  )
}

# Creates a failed ddf model file in the given folder for model specified by mod.
# Useful when a given ddf never finishes running and needs to be marked as failed
# by hand.
create.failed.ddf.file <- function(mod, folder) {
    file.create(file.path(folder, paste0("failed_model_", mod$label, ".txt")))
}
