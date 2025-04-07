# Make sure both distance and distbegin/distend columns are not both in data
# otherwise ds() (as of Distance 1.0.9) will get mad. If distbegin and distend
# are both all NA (as for ECSAS ship data) then just remove them. If they're not
# NA then (e.g., in the case of aerial data) then get rid of distance after
# making sure distance == (distbegin + distend)/2
check.distdata.cols <- function(data) {
  if (!is.null(data$distance) && !is.null(data$distbegin) && !is.null(data$distend)) {
    if (all(is.na(data$distbegin)) && all(is.na(data$distend))) {
      message("check.distdata.cols: removing distbegin and distend columns that are all NA from data")
      data <- select(data, -distbegin, -distend)
    } else if (all.equal(data$distance, (data$distbegin + data$distend) / 2)) {
      message(
        "check.distdata.cols: removing distance column because data has distbegin and distend columns"
      )
      data <- select(data, -distance)
    } else {
      stop(
        paste0(
          "run.ddf.model: data has both distance and non-NA distbegin/distend ",
          "columns but data$distance != (data$distbegin + data$distend)/2"
        )
      )
    }
  }
  data
}


# rel.folder = path (relative to here()) to place the output file in
check.problem.data = function(alldat,
                              rel.folder,
                              leave = NULL,
                              fixed_db = NULL,
                              quiet = TRUE,
                              filename = "suspicious_posns.csv",
                              ...) {

  # Look for watches with problematic watches
  probs <- ECSAS.find.suspicious.posn(alldat, leave = leave, ...) %>%
    mutate(prev_fixed = WatchID %in% fixed_db)

  if (nrow(probs) > 1) {
    message (sprintf("%d problem watches found. %d of these were flagged as previously fixed in db",
                     nrow(probs), sum(probs$prev_fixed)))

    if (!quiet) {
      message("Problem watchIDs:")
      print(probs$WatchID)

      message("Problem watchIDs that were flagged as previously fixed in db:")
      print(filter(probs, prev_fixed == TRUE)$WatchID)

      hist(probs$dist_diff_km)
      hist(probs$pct_diff)
    }

    probs %>%
      select(CruiseID, WatchID, ObserverName, PlatformName, Date, StartTime, EndTime,
             LatStart, LongStart, LatEnd, LongEnd, WatchLenKm, PlatformSpeed, CalcDurMin,
             prev_fixed, dist_dr_km, dist_geo_km, dist_diff_km, pct_diff) %T>%
      write_csv(file = here(rel.folder, filename))
  } else {
    probs <- NULL

    # None found - clean up old files
    message("No problem watches found.")

    # create empty file or truncate if it exists
    file.create(file.path(here(rel.folder, filename)))

    #remove shapefile
    # file.remove(list.files(ShapeDir, pattern = paste0(filename, "\\..*"), full.names = T))
  }

  probs
}



# Creates data structure needed for distance sampling and dsm analysis.
#
# Takes a set of raw data and creates separate obs, and watches clipped to the
# study area and stores them along with the distdata objects in a list called
# "the.data" and returns it. file.prefix is a string that is used to form the
# filenames for various shapefiles. Typically either "ECSAS.ship" etc.
#
# dataset = character string identifying the dataset. Good for multi dataset
# analyses.
# inproj, outproj - projections of input and output datasets respectively
# file.prefix - filename prefix for obs and watches shapefiles that gets
#     prepended to "_obs.shp" and "_watches.shp" to form the output shapefile
#     names.
# saveshp - whether to save shapefiles or not.
# intransect.only - only include observations where InTransect was TRUE?
#
# Note that obs with no count, no distmeth, non-birds, FlySwim != "W" or "F" or
#   association == 18 are filtered out.
create.survey.data <- function(raw.dat = NULL,
                               dataset = NULL,
                               file.prefix = NULL,
                               inproj = 4326,
                               outproj = segProj,
                               saveshp = TRUE,
                               intransect.only = TRUE) {


  coll = makeAssertCollection()
  assert_data_frame(raw.dat, add = coll)
  assert(
    check_class(study.area, "sf"),
    add = coll
  )
  assert(
    check_string(file.prefix),
    add = coll
  )
  reportAssertions(coll)


  # create watches
  message("Creating watches...")
  watches <- raw.dat %>%
    select(
      SurveyType,
      TransectID,
      Program,
      CruiseID,
      WatchID,
      ObserverName,
      TransFarEdge,
      DistMeth,
      Date,
      StartTime,
      EndTime,
      LatStart,
      LongStart,
      LatEnd,
      LongEnd,
      WatchLenKm,
      ObsHeight,
      CalcDurMin,
      TotalWidthKm
    ) %>%
    mutate(Sample.Label = WatchID,
           TransectSides = 1) %>%
    distinct() %>%
    arrange(CruiseID, ObserverName, Date, StartTime)

  # SOMEC data has multiple observers in the same watch (which I guess is ok)
  if (dataset == "SOMEC") {
    watches <- watches %>%
      select(-ObserverName) %>%
      distinct()
  }

  # make sure all data for a watch is consistent. Find rows with duplicate watchIDs
  # and remove these watches
  dups <- watches %>%
    group_by(WatchID) %>%
    summarise(nrows = n()) %>%
    filter(nrows > 1) %>%
    pull("WatchID")

  # Remove watches with inconsistent watch info
  if (length(dups) > 0) {
    warning(paste0(sprintf("Removing %d watches due to inconsistent watch info between rows: ", length(dups)),
                   paste(dups, collapse = ", ")), immediate. = TRUE)
    watches <- filter(watches, !(WatchID %in% dups))
  }

  # clip to study area
  watches <- watches %>%
    st_as_sf(coords = c("LongStart", "LatStart"), crs = st_crs(inproj)) %>%
    st_transform(st_crs(4326)) %>% # for ms_clip below
    select(WatchID) %>% # just keep WatchID
    ms_clip(study.area %>% st_transform(st_crs(4326))) %>%   # do the clipping -
    left_join(watches, by = "WatchID") %>%  # add other cols back in
    st_transform(outproj) %>%
    mutate(StartTime = as.character(StartTime),
           EndTime = as.character(EndTime),)

  # Save as shapefile
  if (saveshp) {
    # Will whine about discarded datum and abbreviated field names until these
    # warnings are removed from new rgdal.
    layer.name <- paste0(file.prefix, "_watches.shp")
    message(sprintf(
      "Saving watches to shapefile '%s'",
      file.path(ShapeDir, layer.name)
    ))

    suppressWarnings(
      st_write(
        watches,
        dsn = ShapeDir,
        layer = layer.name,
        driver = "ESRI Shapefile",
        delete_layer = TRUE
      )
    )
  }

  # Create Obs: remove non-birds, convert distances to km,  windforce is
  # converted to windspeed if necessary, convert InTransect to T/F, remove ship
  # followers, remove zeros (ie. watches where Count == NA),
  # only keep Flyswim == W or F (not L or S). Add in sample.Label,
  # rename columns and select ones of interest
  message("Creating observations...")

  obs <- raw.dat %>%
    mutate(
      Region.Label = 1,
      InTransect = case_when(InTransect == -1 ~ TRUE,
                             InTransect == 0 ~ FALSE,
                             TRUE ~ NA),
      Distance = Distance / 1000
    ) %>%
    filter(
      if (isTRUE(intransect.only))
        InTransect == TRUE
      else
        InTransect %in% c(TRUE, FALSE),
      WatchID %in% watches$WatchID,
      !is.na(Count),
      !is.na(DistMeth),
      # !is.na(Distance), # Now using these with dummy_ddf so don't filter out
      FlySwim %in% c("W", "F"),
      Class == "Bird",
      is.na(Association) | Association != 18
    ) %>%
    # Assign DistType (Need to do after filtering FlySwim for W or F),
    # and add FlockID if there isn't one
    mutate(DistType = assign.dist.type(.),
           FlockID = case_when(all(is.na(FlockID)) ~ 1:nrow(.),
                               TRUE ~ FlockID)) %>%
    left_join(watches[, c("WatchID", "Sample.Label")], by = "WatchID") %>%
    rename(object = FlockID,
           size = Count,
           distance = Distance) %>%
    select(
      SurveyType,
      object,
      Sample.Label,
      Region.Label,
      ObserverName,
      DistMeth,
      DistType,
      Program,
      WatchID,
      Date,
      Year,
      Alpha,
      English,
      size,
      FlySwim,
      InTransect,
      distance,
      Visibility,
      Windspeed,
      Windforce,
      Glare,
      SeaState,
      TransFarEdge,
      Swell,
      LatStart,
      LongStart,
      ObsLat,
      ObsLong,
      ObsTime,
      ObsHeight,
      Association,
      Behaviour,
      Age,
      Plumage,
      Sex,
      DistanceCode,
      Seabird,
      Waterbird,
      spec.grp
    ) %>%
    mutate(
      Windspeed = case_when(
        is.na(Windspeed) &
          !is.na(Windforce) ~ left_join(., beaufort_conversion,
                                        by = c("Windforce" = "beaufort"))$speed.kts,
        TRUE ~ Windspeed
      ),
      weights = case_when(
        DistanceCode %in% c("A", "B") ~ 2,
        DistanceCode %in% c("C", "D") ~ 1,
        TRUE ~ 0
      ),
      Visibility = case_when(Visibility > 20 ~ 20,
                             TRUE ~ Visibility),
      FlySwim = as.factor(FlySwim),
      distbegin = NA_integer_,
      distend = NA_integer_
    ) %>%
    droplevels

  # Make sure there's still some left
  if (nrow(obs) == 0)
    warning("No observations left after filtering!", immediate. = TRUE)

  # Clip to study area
  obs <- obs %>%
    st_as_sf(coords = c("LongStart", "LatStart"), crs = st_crs(inproj)) %>%
    st_transform(st_crs(4326)) %>% # for ms_clip below
    select(object) %>% # just keep WatchID
    ms_clip(study.area %>% st_transform(st_crs(4326))) %>%   # do the clipping -
    left_join(obs, by = "object") %>%  # add other cols back in
    st_transform(outproj)

  # Save as shapefile
  if (saveshp) {
    # Will whine about discarded datum and abbreviated field names until these
    # warnings are removed from new rgdal.
    layer.name <- paste0(file.prefix, "_obs.shp")
    message(sprintf(
      "Saving obs to shapefile '%s'",
      file.path(ShapeDir, layer.name)
    ))

    suppressWarnings(
      st_write(
        obs,
        dsn = ShapeDir,
        layer = layer.name,
        driver = "ESRI Shapefile",
        delete_layer = TRUE
      )
    )
  }

  ##### After all that, now just create a single distdata for use in distance
  #sampling.This doesn't actually get used (unless we do analysis of all
  #seabirds). Rather, specific datasets for a given species/group are created by
  #create.dsm.data( in order to properly generate the samples where there were 0
  #observations of a given species.
  message("Creating distdata....")
  distdata <- obs %>%
    st_drop_geometry() %>%
    droplevels

  # Set dataset attribute
  distdata$Dataset <- watches$Dataset <- dataset

  the.data <-
    list(
      distdata = distdata,
      watches = st_drop_geometry(watches)
    )

  message("\nDone\n")
  the.data
}

# Convert watches object to spatialinesdataframe
watches.to.lines <- function(watches) {

  lins <- watches %>%
    sp::split(.$WatchID) %>%
    lapply(function(x) {
      sp::Lines(list(sp::Line(
        matrix(c(x$LongStart, x$LatStart, x$LongEnd, x$LatEnd),
               nrow =  2, byrow = T))),
        x$WatchID[1L])
      }) %>%
    sp::SpatialLines()

  watches  <- as.data.frame(watches)
  rownames(watches) <- watches$WatchID
  l <- sp::SpatialLinesDataFrame(lins, watches)
  sp::proj4string(l) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  l
}


my.gam.check <- function(b, old.style = FALSE, type = c("deviance", "pearson",
  "response"), k.sample = 5000, k.rep = 200, rep = 0, level = 0.9,
  rl.col = 2, rep.col = "gray80", ...)
{
  type <- match.arg(type)
  resid <- residuals(b, type = type)
  linpred <- if (is.matrix(b$linear.predictors) && !is.matrix(resid))
    napredict(b$na.action, b$linear.predictors[, 1])
  else napredict(b$na.action, b$linear.predictors)

  if (old.style)
    qqnorm(resid, ...)
  else qq.gam(b, rep = rep, level = level, type = type, rl.col = rl.col,
    rep.col = rep.col, ...)
  plot(linpred, resid, main = "Resids vs. linear pred.", xlab = "linear predictor",
    ylab = "residuals", ...)
  hist(resid, xlab = "Residuals", main = "Histogram of residuals",
    ...)
  fv <- if (inherits(b$family, "extended.family"))
    predict(b, type = "response")
  else fitted(b)
  if (is.matrix(fv) && !is.matrix(b$y))
    fv <- fv[, 1]
  plot(fv, napredict(b$na.action, b$y), xlab = "Fitted Values",
    ylab = "Response", main = "Response vs. Fitted Values",
    ...)
  gamm <- !(b$method %in% c("GCV", "GACV", "UBRE", "REML",
    "ML", "P-ML", "P-REML", "fREML"))
  if (gamm) {
    cat("\n'gamm' based fit - care required with interpretation.")
    cat("\nChecks based on working residuals may be misleading.")
  }
  else {
    cat("\nMethod:", b$method, "  Optimizer:", b$optimizer)
    if (!is.null(b$outer.info)) {
      if (b$optimizer[2] %in% c("newton", "bfgs")) {
        boi <- b$outer.info
        cat("\n", boi$conv, " after ", boi$iter, " iteration",
          sep = "")
        if (boi$iter == 1)
          cat(".")
        else cat("s.")
        cat("\nGradient range [", min(boi$grad), ",",
          max(boi$grad), "]", sep = "")
        cat("\n(score ", b$gcv.ubre, " & scale ", b$sig2,
          ").", sep = "")
        ev <- eigen(boi$hess)$values
        if (min(ev) > 0)
          cat("\nHessian positive definite, ")
        else cat("\n")
        cat("eigenvalue range [", min(ev), ",", max(ev),
          "].\n", sep = "")
      }
      else {
        cat("\n")
        print(b$outer.info)
      }
    }
    else {
      if (length(b$sp) == 0)
        cat("\nModel required no smoothing parameter selection")
      else {
        cat("\nSmoothing parameter selection converged after",
          b$mgcv.conv$iter, "iteration")
        if (b$mgcv.conv$iter > 1)
          cat("s")
        if (!b$mgcv.conv$fully.converged)
          cat(" by steepest\ndescent step failure.\n")
        else cat(".\n")
        cat("The RMS", b$method, "score gradient at convergence was",
          b$mgcv.conv$rms.grad, ".\n")
        if (b$mgcv.conv$hess.pos.def)
          cat("The Hessian was positive definite.\n")
        else cat("The Hessian was not positive definite.\n")
      }
    }
    if (!is.null(b$rank)) {
      cat("Model rank = ", b$rank, "/", length(b$coefficients),
        "\n")
    }
  }
  cat("\n")
  kchck <- k.check(b, subsample = k.sample, n.rep = k.rep)
  if (!is.null(kchck)) {
    cat("Basis dimension (k) checking results. Low p-value (k-index<1) may\n")
    cat("indicate that k is too low, especially if edf is close to k'.\n\n")
    printCoefmat(kchck, digits = 3)
  }
}

# Do observed/expected for density response or abundance response
oe.dens <- function(model, covar, cut = NULL, plotit = FALSE, debug = FALSE, ...)
{
    if (debug) browser()
    #get data
    oe <- model$data

    # add in predictions as abundances
    resp <- as.character(model$formula[[2]])
    if (resp %in% c("D", "density", "Dhat", "density.est")) {
      oe$N <- predict(model) * oe$segment.area
    } else {
      oe$N <- predict(model)
    }
    if (!is.null(cut)) {
        oe[[covar]] <- cut(oe[[covar]], breaks = cut)
    }
    oe <- plyr::ddply(oe, covar, function(x) {
        data.frame(Observed = sum(x[resp]), Expected = sum(x$N), n = nrow(x))
    })
    cn <- oe[, 1]
    oe <- t(oe[, 2:3])
    colnames(oe) <- cn
    if (plotit) {
      maxlim <- max(oe[1,], oe[2,])
      plot(oe[1,], oe[2,], xlab = "obs", ylab = "exp",
           main = paste0("Obs vs exp at specific values of ", covar),
           xlim = c(0, maxlim), ylim = c(0, maxlim), ...)
    }
    return(oe)
}


# Render a full preliminary dsm analysis Rmd for one species.
# XXX Used?
do.full.prelim <- function(species){

  if(!dir.exists(file.path(ResultsDir, species))) {
    dir.create(file.path(ResultsDir, species))
  }

  out.file <-
    file.path(ResultsDir,
              species,
              paste(species, "00.04_full_prelim.html", sep = "_"))
  message(sprintf("Rendering full prelim %s to %s",
                  species,
                  out.file))
  render(
    file.path(RDir, "Generic_full_test.Rmd"),
    params = list(species = species),
    output_file = out.file
  )
}


# Set default convert.units, cutpoints, etc values when setting up ddf model
# list.
set.def.df.spec.values <- function(ddf.def, nm) {

  if(length(nm) != 1)
    stop("set.def.df.spec.values: length of ddf name is not 1.")


  # Note that the convert_units won't actually get used because we
  # only use ds() to compute the detection function and not
  # abundance.
  if(str_detect(nm, "Ship")) {
    # ECSAS Ship surveys have set cutpoints whereas aerial uses distbegin
    # and distend columns computed in 00.01_Extract_data.Rmd.
    ddf.def$convert_units = ecsas.ship.convert.units
    ddf.def$cutpoints = ecsas.ship.ctpoints
    ddf.def$distance.centers = ecsas.ship.distance.centers
  } else if (str_detect(nm, "Aerial")) {
    ddf.def$convert_units = ecsas.air.convert.units
    ddf.def$distance.centers = ecsas.air.distance.centers
  } else
    stop("set.def.df.spec.values: unrecognized ddf name.")
  ddf.def
}

# Parse a ddf name of the form "ECSAS_Ship_F" and return a list of the
# 3 elements.
parse.df.name <- function(nm){
  res <- str_split_1(nm, fixed("_"))
  list(dataset = res[1], platform = res[2], behav = res[3])
}


# Three functions to peel off specific aspects of a ddf spec name
get.dataset <- function(nm) {
  parse.df.name(nm)$dataset
}

get.platform <- function(nm) {
  parse.df.name(nm)$platform
}

get.behav <- function(nm) {
  parse.df.name(nm)$behav
}

# Process all ddf specs for a given species. Called from Generic_1_ddf_fitting.Rmd
do.det.fcn.specs <-
  function(df.specs,
           species,
           distdata,
           do.final.only,
           rerun,
           parallel,
           nCores,
           do.eda = TRUE,
           ...) {

  ret <- df.specs %>%
    imap(
      do.det.fcn.spec,
      species = species,
      distdata = distdata,
      do.final.only = do.final.only,
      rerun = rerun,
      parallel = parallel,
      nCores = nCores,
      do.eda = do.eda,
      ...
    )

  # resave ddf spec for this species which may now have final fitted model
  #  and fitted distdata objects
  save.ddf.specs(ret, species)

  ret
}

# Save ddf specs (multiple specs - one for each dataset, platform and behav)
# for a single species.
save.ddf.specs <- function(df.mod.specs, species) {
  filename <- file.path(RDataDir, paste0(species, dfModListSuffix))
  message(sprintf("Saving ddf specs for %s in '%s'", species, filename))
  save(df.mod.specs, file = filename)
}


# Do ddf fitting for a given df.spec (which will be for one
# dataset/platform/behav for 1 species), extracting the correct data, building
# either all candidate ddfs, or just the final one. If final.only, store the
# returned ddf in the df.mod.list structure.
#
# df.spec - the ddf specification (see my structure in analysis_settings.R)
#
# df.spec.nm - name of the df.spec - e.g., ECSAS_Ship_W
#
# species - char vector
#
# distdata - should only contain obs for given species, but we filter on species
#   anyway just in case.
#
# Typically called once for each ddf specification in each species in df.mod.list
do.det.fcn.spec <- function(df.spec,
                            df.spec.nm,
                            species,
                            distdata,
                            do.final.only,
                            rerun,
                            parallel,
                            nCores,
                            do.eda,
                            ...) {

  dataset <- get.dataset(df.spec.nm)
  platform <- get.platform(df.spec.nm)
  behav <- get.behav(df.spec.nm)

  message(sprintf(
    "Doing %s, %s, for (%s, %s, %s) ",
    species,
    ifelse(
      do.final.only,
      "final detection function",
      "all candidate det fcns"
    ),
    dataset,
    platform,
    behav
  ))

  # Folder to save ddf fitting results when doing all candidate ddfs.
  folder <-
    file.path(ResultsDir,
              species,
              paste("DF Summaries", dataset, platform, behav, sep = "_"))

  # Filter distdata based on whether this is a ddf for obs with no perp
  # distances, or a normal one. Record this fact in distdata for use by
  # create.dsm.data).
  if (behav == "S") {
    distdata <- filter(
      distdata,
      dataset == dataset,
      SurveyType == platform,
      Alpha %in% spec.grps[[species]],
      # only doing 300m strip transect, ignore others (only a few)
      !(DistanceCode %in% c("F", "G", "I", "K", "5")),
      # Weed out ones with NA distance but width of transect was more than 300m
      TransFarEdge <= 300 | is.na(TransFarEdge),
      (DistType != "Perp." | is.na(distance))
    )
  } else if (behav %in% c("W", "F")){
    # Normal ddf
    # Filter perp distances for this species, behav, platform, dataset
    distdata <- filter(
      distdata,
      dataset == dataset,
      SurveyType == platform,
      FlySwim == behav,
      Alpha %in% spec.grps[[species]],
      DistType == "Perp.",
      # Only use obs with perp. distances.
      !is.na(distance)
    )
  } else
    stop("do.det.fcn.spec: unknown behav: ", behav)

  # assign season
  distdata <- assign.season(distdata, seasons[[species]])

  # # XXXX should make dynamic instead of relying on a static list.
  if (isTRUE(df.spec$remove.fishing))
    distdata %<>% filter(!(CruiseID %in% fishingCruises))

  # EDA
  message("\nDistribution of distances:")
  print(table(distdata$distance, useNA = "always"))
  cat("\n")
  print(summary(distdata))
  # Look at covar plots
  if (isTRUE(do.eda) && nrow(distdata) > 0){
    ds.eda(
      distdata,
      species = species,
      do.final.only = do.final.only,
      df.spec = df.spec)
  }


  # Only need to call do.det.fcn() if we're doing final ddfs (no matter if strip
  # trnasect or normal), OR we're doing all possible ddfs and behav isn't
  # "S". In the latter case, we are doing all possible ddfs as a side effect
  # and ignoring the return value from do.ds() called in do.det.fcn. We only
  # bother to call it to fit all possible models if behav isn't "S" since if
  # it was S we already know what the final model will be and there is no
  # need to fit "all possible" models since there aren't any.
  #
  # XXX Need to think about moving much of the above code into this conditional
  # since most of it doesn't need to be run in do.final == TRUE.
  if(do.final.only == TRUE || behav != "S"){
    # fit det fcn(s)
    res <- do.det.fcn(distdata = distdata,
                      dsetname = dataset,
                      species = species,
                      df.model = df.spec,
                      df.model.nm = df.spec.nm,
                      folder = folder,
                      do.final.only = do.final.only,
                      rerun = rerun,
                      parallel = parallel,
                      nCores = nCores,
                      ...)
    # When do.final.only is FALSE, res is either a list of models or a dataframe
    # of model specifications (if runModels was FALSE in do.ds) so no point in
    # trying to get the df_final and distdata elements.
    if (do.final.only) {
      df.spec$fitted.model <- res$df_final
      df.spec$fitted.distdata <- res$distdata
      # ddftype combines info on platform (aerial, ship) and behav(W, F, S). Used
      # as factor in dsm models.
      df.spec$ddftype <- substr(platform, 1, 1) %>%
        paste0(behav) %>%
        factor(levels = ddftype.levels)
    }
  }

  # return the per-species df.model (with newly added fitted.model object if
  # do.final.only was true)
  df.spec
}

# Return a dummy_ddf() for various situations
create.strip.ddf <- function(distdata, dataset, platform) {
  # Set up appropriate one strip transect
  if (dataset == "ECSAS") {
    if (platform == "Aerial") {
      if (!is.numeric(distdata$object)) {
        warning(
          "create.strip.ddf: converting object id to numeric. This should be done in 00.01_Extract_data",
          immediate. = TRUE
        )
      }

      # XXX need to replace calculation of left with info from the aerial
      # transect for each observer??? Since it isn't really the smallest
      # distance if there are obs and it isn't really 0 if there are no obs.
      df_final <- dummy_ddf(
        distdata$object %>%
          str_replace("[a-zA-Z_]+", "") %>%
          as.numeric(),
        distdata$size,
        left = ifelse(nrow(distdata) == 0, 0, min(distdata$distbegin)),
        width = ecsas.air.right.trunc
      )
    } else if (platform == "Ship") {
      if (!is.numeric(distdata$object)) {
        warning(
          "create.strip.ddf: converting object id to numeric. This should be done in 00.01_Extract_data",
          immediate. = TRUE
        )
      }

      # For ECSAS ship we have defined cutpoints and no left trunc
      df_final <- dummy_ddf(
        distdata$object %>%
          str_replace("[a-zA-Z_]+", "") %>%
          as.numeric(),
        distdata$size,
        width = max(ecsas.ship.ctpoints)
      )
    } else
      # Platform == ...
      stop(sprintf(
        "create.strip.ddf: unrecognized survey platform type: %s",
        platform
      ))
  } else
    # dataset == ...
    # Add PQ stuff here later
    stop(sprintf(
      "create.strip.ddf: unrecognized survey dataset name: %s",
      dataset
    ))

  # replace the minimal data the dummy_ddf adds to df_final with our own with
  # more columns so that downstream processing can grab the distdata from df_final
  df_final$data <- distdata
  df_final
}

# Called to run final ddf for one specie and dataset combination. Run
# ds() analysis either for all possible models or for just the final model
# specified by df.model.
#
#
# distdata = the observation dataframe to be subsetted for the analysis
# species = char vector
# df.model = det fcn model specification used if do.final.only is true.
# dataset = chr string indicating which dataset (ECSAS, SOMEC, etc) we're analyzing.
# do.final.only = if true then the final chosen det fcn is fitted else all
#   candidate det fcns are run
# rerun = logical, if true any previously run candidate det fcns will be re run.
#   Previously run detection functions are those that have results in folder.
# folder = where to save the df summaries.
# parallel = logical, should parallel processing be used to run multiple
#   candidate det fcns at the same time.
# nCores = number of (virtual) cores to use for parallel processing
# ... = other args passed on to ds().
#
# Return list with df_final and distdata (augmented with detProb and adjSize)
do.det.fcn <- function(distdata,
                       species,
                       df.model = NULL,
                       df.model.nm = NULL,
                       dsetname = NULL,
                       do.final.only = TRUE,
                       cleanFolder = FALSE,
                       rerun = FALSE,
                       parallel = FALSE,
                       folder,
                       nCores = 1,
                       ...) {

  # set up species and spill specific pathnames
  # source(here::here("species settings.r"), echo = T, local = TRUE)

  dataset <- get.dataset(df.model.nm)
  platform <- get.platform(df.model.nm)

  # Set truncation distance if necessary for surveys without known set cutpoints
  if (dataset == "ECSAS") {
    if (platform == "Aerial") {
      truncation <-
        list(left = ifelse(nrow(distdata) == 0, 0, min(distdata$distbegin)),
             right = ecsas.air.right.trunc)
    } else if (platform == "Ship") {
      truncation <- max(df.model$cutpoints)
    } else
      stop(sprintf("do.det.fcn: unrecognized survey platform type: %s", platform))
  } else
    # Add PQ stuff here later
    stop(sprintf("do.det.fcn: unrecognized survey dataset name: %s", dataset))

  if (do.final.only) {
    # CleanFolder only applies when fitting all possible ddfs.
    if (cleanFolder)
      warning("do.det.fcn: do.final == TRUE, ignoring cleanFolder == TRUE.",
              immediate. = TRUE)

    # Doing final ddf
    if (df.model$strip) {
      df_final <- create.strip.ddf(distdata, dataset, platform)
    } else { # df.model$strip == FALSE

      if (rerun) {
        # (re-)fit the model - Doing a real ddf

        # remove unneeded cols, but be careful because some cols are needed
        # downstream by create.dsm.data) for example.
        # and remove redundant distance or distbegin/distend cols.
        distdata <-
          distdata %>%
          dplyr::select(object, size, distbegin, distend, distance, Season,
                        SurveyType, FlySwim, Sample.Label, WatchID, Alpha, Dataset,
                        LatStart, LongStart,
                        all_of(all.vars(df.model$final.formula))) %>%
          check.distdata.cols()

        ###### use do.ds machinery to re-fit final model
        df_final <- do.ds(
          data = distdata,
          runModels = TRUE,
          incl.adj = FALSE,
          parallel = parallel,
          folder = folder,
          rerun = TRUE,
          models = with(
            df.model,
            data.frame(
              label =  create.model.name(final.key, final.formula, final.adj),
              key = final.key,
              # Can't insert a formula into a dataframe and as.character
              # separates the ~ from the rest, so paste it back together
              # then turn it back into a formulat in do.ds
              form =  paste(as.character(final.formula), collapse = " "),
              adj = ifelse(is.null(final.adj), "", final.adj)
            )
          ),
          max_adjustments = df.model$max.adjustments,
          nadj = df.model$final.nadj,
          cutpoints = df.model$cutpoints,
          convert_units = df.model$convert.units,
          truncation = truncation
        ) %>%
          # Note df_final will be a list with one element and we need to peel it
          # off.
          pluck(1)

      } else { # rerun = FALSE
        # Look for file with results of running this model already

        # Sometimes I forget to add as.formula() to final model specs, so
        # check it here.
        stopifnot(class(df.model$final.formula) == "formula")

        nam <- with(df.model,
                    create.model.name(final.key, final.formula, final.adj))
        pat <- paste0("AIC.*_", nam, ".RData")
        file <- list.files(folder, pattern = pat, full.names = TRUE,)
        if (length(file) != 1) {
          message("do.det.fcn: do.final == TRUE and rerun == FALSE. Looking for existing ddf model results file. Found ", length(file), " file(s), should be 1. Quitting.")
          stop()
        }

        # get the model
        load(file)
        df_final <- model
      } # End rerun = FALSE
    } # End strip == FALSE

    # Check model and augment distdata with detprob and adjSize
    distdata <- check.det.fcn(df_final, species, df.model.nm)

    return(list(df_final = df_final, distdata = distdata))
  } else {  # do.final.only == FALSE, doing all possible models
    # Do all candidate ddfs
    # folder <- ifelse(is.null(folder), summaryDir, folder) # folder is required now.

    # Remove redundant distance or distbegin/distend cols.
    distdata <- check.distdata.cols(distdata)
    do.ds(
      data = distdata,
      covars = df.model$candidate.covars,
      cutpoints = df.model$cutpoints,
      folder = folder,
      convert_units = df.model$convert.units,
      max_adjustments = df.model$max.adjustments,
      truncation = truncation,
      parallel = parallel,
      nCores = nCores,
      rerun = rerun,
      cleanFolder = cleanFolder,
      ...)
  }
}

# Exploratory data analysis
# Produce plots distance vs various covariates
ds.eda <-
  function(distdata = NULL,
           species = NULL,
           do.final.only = T,
           df.spec = NULL) {

  if (nrow(distdata) == 0)
    return()

  if (do.final.only == TRUE) {
    vars <- all.vars(df.spec$final.formula)
  } else {
    vars <- df.spec$candidate.covars
  }

  # plot individual covars
  walk(vars, plot.covar, distdata = distdata, species = species, df.spec = df.spec)
}


# Make a plot of a covar vs distance. Note the use of get(covar) to handle non-standard
# evaluation and allows us to pass a column name in covar.
plot.covar <-
  function(covar = NULL,
           distdata = NULL,
           species = NULL,
           df.spec = NULL) {

  if (covar == "size") {
    p <- ggplot(data = as.data.frame(distdata), aes(cut(distance, df.spec$cutpoints, right = FALSE), size))
    p <- p + geom_boxplot(varwidth = TRUE)
    p <- p + labs(x = "Distance Category", y = "Size", title = species)
    print(p)
  } else if (is.factor(distdata[, covar][[1]])) {
    p <- ggplot(data = as.data.frame(distdata), aes(get(covar), distance))
    p <- p + geom_violin(draw_quantiles = c(.25, .5, .75), scale = "count")
    p <- p + scale_y_continuous(labels = as.character(df.spec$distance.centers),
                                breaks = df.spec$distance.centers)
    p <- p + labs(y = "Distance Category", x = covar, title = species)
    print(p)
  } else {
    p <- ggplot(data = as.data.frame(distdata), aes(get(covar), distance))
    p <- p + geom_point()
    p <- p + geom_smooth(method = "lm")
    p <-
      p + scale_y_continuous(
        labels = as.character(df.spec$distance.centers),
        breaks = df.spec$distance.centers
      )
    p <- p + labs(y = "Distance Category", x = covar, title = species)
    print(p)
  }
}

# The next 3 functions are used to convert a dataframe that contains
# coords of endpoints of a line (ie transect or watch) into an sf line object
#
# Make a linestring from from a matrix of coords for 2 points
make.line <- function(xy2){
  st_linestring(matrix(xy2, nrow=2, byrow=TRUE))
}

# convert dataframe of coords to lines
make.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  m = as.matrix(df[,names])
  lines = apply(m, 1, make.line, simplify=FALSE)
  st_sfc(lines, crs = crs)
}

# Convert a dataframe to lines. Used in Extract_data.Rmd
sf.pts.to.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  geom = make.lines(df, names, crs)
  df = st_sf(df, geometry=geom)
  df
}

# Assign season based on date ignoring year.  Works across new year boundary
#
# Inputs:
#   dat - object containing data to be assigned.
#   season.def - definition of the seasons.
#   datefield - Character string giving the name of the date field in dat.
#
# Outputs:
#   Input dataframe with a Season field added (if needed) set to  correct season.
assign.season <- function(dat, season.def, datefield = "Date"){

  if (is.null(season.def)) {
    stop("Could not find a season definition for this species!")
  }

  # Short circuit exit if no data
  if (nrow(dat) == 0)
    return(dat)

  # init
  dat$Season <- NA
  season.index <- rep(NA, nrow(dat))

  # convert dates to my format
  dates <- pull(dat, datefield)
  dat$monthday <- month(dates) * 100 + day(dates)

  # step through each season with cheesy for loop
  for (i in seq_along(season.def)) {

    # Get logical index of matching rows
    if (season.def[[i]]["from"] > season.def[[i]]["to"]) {
      # wraps around new year
      criteria <-
        dat$monthday >= season.def[[i]]["from"] |
        dat$monthday <= season.def[[i]]["to"]
    } else {
      criteria <-
        between(dat$monthday, season.def[[i]]["from"], season.def[[i]]["to"])
    }

    # assign the current season to the index of matching rows
    season.index[criteria] <- i
  }

  # Assign season
  dat$Season <- factor(
    names(season.def)[season.index],
    levels = season.names,
  )

  if(any(is.na(dat$Season)))
    warning(sprintf("assign.season: %d rows failed to have season assigned",
                    sum(is.na(dat$Season))), immediate. = TRUE)
  dat
}


# Convert a dataframe to sf and save as shapefile. Called from Extract_data.Rmd
df.to.shapefile <- function(df,
                            coords = c("LongStart", "LatStart"),
                            crs = st_crs(4326),
                            out.proj,
                            dsn = ShapeDir,
                            layer)
{
  st_as_sf(df,
           coords = coords,
           crs = crs,
           remove = FALSE) %>%
    st_transform(out.proj) %>%
    st_write(
      dsn = dsn,
      layer = layer,
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
}

# Get distances and SD between successive GPS points along a watch. Called from
# Extract_data.Rmd.
add.dist.sd <- function(watch) {
  if (nrow(watch$posns[[1]]) > 2) {
    watch$dists <-
      list(geodist::geodist(watch$posns[[1]], sequential = TRUE, measure = "geodesic"))
    watch$dists.sd <- round(sd(unlist(watch$dists)), 4)
  } else {
    watch$dists.sd <- -1
  }
  watch
}

# get the sum of lengths between GPS positions for a given watch. Return value
# is in km. Called from Extract_data.Rmd.
get.gps.length <- function(watch, posns) {
  # get positions in this watch
  posns <- filter(posns,
                  between(posns$datetime, watch$WatchStartTime, watch$WatchEndTime))

  # sum distances between the points. Note that the units of val will
  # depend on the projection of the coords in posns, but we assume here
  # they are in meters
  if (nrow(posns) >= 2) {
    val <-
      list(geodist::geodist(posns, sequential = TRUE, measure = "geodesic")) %>%
      unlist %>%
      sum
  } else {
    val <- 0
  }

  val/1000
}


# make a raster from an sf object element for a given season. Called from
# dsm.pred()
make.season.raster <- function(season, obj, variable) {
  ret <- filter(obj, Season == season) %>%
    select(all_of(variable), geometry) %>%
    stars::st_rasterize() %>%
    rast
  ret
}

# copy objects from one env to another. Useful after running a background job
# whose results are returned in an environment that you save in an .rda file.
# Reloading the .rda into a clean session will give one object (the environment)
# in the global env. So use:
#
# copy.env(name_of_loaded_env, .GlobalEnv)
copy.env <- function(src, dst) {
  for(n in ls(src, all.names=TRUE))
    assign(n, get(n, src), dst)
}

## Create initial segment data from watches
# inproj - proj4string or EPSG of input watches
# outproj - proj4string or EPSG of output segdata
# depth and depht.g rasters are created by predgrid creation and are passed in.
#
# Note internally all extraction of environmental values from rasters uses
# segProj and all saved rasters are in segProj. So, watch data coords will need
# to be in segProj to do the extractino
#
# Not using rxtractogon any more. Instead, netcdf files have been downloaded by hand
# and rasters extracted from them. Those rasters are utilized here.
create.segdata <- function(the.data,
                           study.area,
                           inproj,
                           outproj,
                           scale.factors,
                           verbose = FALSE) {

  # create initial segdata and reproject
  if (verbose) message("Creating initial segdata from watches")
  segdata <- the.data$watches %>%
    rename(Effort = WatchLenKm) %>%
    mutate(TransectID = as.character(TransectID),
           year = lubridate::year(Date),
           yday = lubridate::yday(Date),
           MonthYear = format(Date, "%Y-%m"),
           segment.area = Effort * TotalWidthKm) %>%
    st_as_sf(coords = c("LongStart", "LatStart"), crs = inproj, remove = FALSE) %>%
    st_transform(outproj) %>%
    cbind(st_coordinates(.)) %>%
    rename(x = X, y = Y)

  ###---------------------------------------------------------------------------
  #### Load rasters
  if (verbose) message("Loading environmental rasters")

  # get the dates of monthly rasters needed (ie. sst, sst.g, etc) so we can read
  # the needed files into a big SpatRaster
  dates.needed <- segdata$Date %>%
    as.character %>%
    str_sub(end = -4) %>%
    paste0("-16") %>%
    unique %>%
    str_sort

  ### Depth and other static rasters
  # Depth
  depth <- rast(file.path(predLayerStudyAreaDir, "depth.img"))

  # Depth gradient
  depth.g <- rast(file.path(predLayerStudyAreaDir, "depth.g.img"))

  # SST
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.", dates.needed, ".img"))
  sst <- rast(files)

  # SST gradient
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.g.", dates.needed, ".img"))
  sst.g <- rast(files)

  ###---------------------------------------------------------------------------
  ### Extract raster values at segdata locations

  ### Extract depth values
  if (verbose) message("Extracting depth at watch locations")
  segdata <-
    terra::extract(
      depth,
      vect(segdata),
      bind = TRUE
    ) %>%
    st_as_sf

  ### Extract depth gradient values
  if (verbose) message("Extracting depth gradient at watch locations")
  segdata <-
    terra::extract(
      depth.g,
      vect(segdata),
      bind = TRUE
    ) %>%
    st_as_sf

  ### Extract SST values
  if (verbose) message("Extracting sst at watch locations")
  needed.layers <- paste0("sst.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst,
      vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>%
    st_as_sf %>%
    rename(sst = value) %>%
    select(-layer) # Note that layer is off by one even though the sst values is correct

  ### Extract SST gradient values
  if (verbose) message("Extracting sst gradient at watch locations")
  needed.layers <- paste0("sst.g.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst.g,
      vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>%
    st_as_sf %>%
    rename(sst.g = value) %>%
    select(-layer) # Note that layer is off by one even though the sst values is correct

  ###---------------------------------------------------------------------------
  ## Add scaled versions of all preds.
  if (verbose) message("Scaling covars")

  segdata %<>%
    mutate(depth.sc = (depth - scale.factors$depth_mean)/scale.factors$depth_sd,
           depth.g.sc = (depth.g - scale.factors$depth.g_mean)/scale.factors$depth.g_sd,
           x.sc = (x - scale.factors$x_mean)/scale.factors$x_sd,
           y.sc = (y - scale.factors$y_mean)/scale.factors$y_sd,
           sst.sc = (sst - scale.factors$sst_mean)/scale.factors$sst_sd,
           sst.g.sc = (sst.g - scale.factors$sst.g_mean)/scale.factors$sst.g_sd,
    )

  ###---------------------------------------------------------------------------
  ## Save segdata
  if (verbose) message("Saving results")
  save(segdata, file = segdatloc)

  st_write(
    segdata,
    dsn = ShapeDir,
    layer = "segdata.shp",
    driver = "ESRI Shapefile",
    delete_layer = TRUE
  )

  if (verbose) print("Done.")
  segdata
}

# NOTE: no longer used due to issues with timeouts/lags/delays/etc.
# but kept for posterity in case I return to this approach.
#
#
# Extract remote sensing data given by dataset and parameter in a polygon given
# by xcoord, ycoord (decimal degrees) for date given by tcoord and convert to a
# raster with projection prj.
#
# note: if tcoord is a length 2 vector, rxtractogon will return a whole series
# of data from tcoord[1] to tcoord[2], however that is not the intended use
# here.
#
# NOTE: to get the right monthly raster, tcoord must be of form "year-month-15"
# (or so), b/c rxtractogon() will choose the monthly dataset with the
# closest data so specifying "year-month" or "year-month-01" may get the previous
# month, which will cause no end of tears when trying to debug....
#
#
# this is designed to be used as: map(tcoords, rxtractogon.rast, blah, blah,
# blah, ...) to extract one raster per tcoord
# #
# Params:
# tcoord - time coordinate for dataset to download. Not all datasets (e.g. ETOPO depth)
#     have a tcoord in which case it should be null.
#
# dataset - name of ERDDAP dataset to access (e.g. jplMURSST41mday)
#
# parameter - name of paramter to get from dataset (e.g. sst)
#
# xcoord, ycoord - outline of area required (it will extract the bounding box)
#
# plotit - plot downloaded data
#
# saveit - save downloaded data to RasterDir
#
# folder -
rxtractogon.rast <-
  function(tcoord,
           dataset,
           parameter,
           xcoord,
           ycoord,
           plotit = FALSE,
           saveit = FALSE,
           folder) {

  # Some dataset (ie depth) don't require a tcoord
  if (is.null(tcoord)) {
    print(sprintf("Extracting %s from %s", parameter, dataset))
  } else {
    if (length(tcoord) != 1)
      stop(sprintf("tcoord should have length 1 but has length %d", length(tcoord)))

    print(sprintf("Extracting %s from %s for %s", parameter, dataset, tcoord))
  }

  dat <-
    rxtractogon(
      rerddap::info(dataset),
      parameter =  parameter,
      xcoord = xcoord,
      ycoord = ycoord,
      tcoord = tcoord
    )
  layername <- names(dat)[1]
  dat <- pluck(dat, 1) # use 1 since it is always 1st element, but not always called same as value of parameter
  if (length(dim(dat)) > 2) # remove useless third dimension
    dat <- dat[,,1]

  rast <- dat %>%
    t %>%                   # rxtracto returns matrix in odd order with x and y transposed and south to north so: transpose
    .[nrow(.):1, ] %>%       # ... and reverse order of rows
    raster(
      xmn = min(xcoord),
      xmx = max(xcoord),
      ymn = min(ycoord),
      ymx = max(ycoord),
      crs = latlongproj
    )

  if (plotit)
    plot(rast, main = sprintf("%s from %s for %s", parameter, dataset, tcoord))

  if (saveit) {

    if (is.null(tcoord))
      filename <- file.path(folder, paste(layername, "img", sep = "."))
    else
      filename <- file.path(folder, paste(layername, tcoord, "img", sep = "."))

    message("Saving downloaded ERDDAP raster to ", filename)
    if (!dir.exists(folder))
      dir.create(folder, recursive = TRUE)
    writeRaster(rast, filename = filename, format = "HFA", overwrite = TRUE)
  }

  rast
}

# XXX No longer used??
# function stolen from http://r-sig-geo.2731867.n2.nabble.com/Run-focal-function-on-a-multi-layer-raster-td7589931.html
# to perform focal() on RasterStack
multi.focal <- function(x, w = matrix(1, nrow = 3, ncol = 3), ...) {

  if (is.character(x)) {
    x <- brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w = w, ...)
  }

  n <- seq(nlayers(x))
  list <- lapply(X = n, FUN = fun, x = x, w = w, ...)

  out <- stack(list)
  return(out)
}

# Reads all raster files in folder with names matching pattern and returns a raster
# stack.
#
# Useful for debugging and to avoid downloading from ERDDAP again.
#
# e.g. x <- recreate.sst.mnth.from.files(RasterDir, "sst.+img$")
recreate.sst.mnth.from.files <- function(folder, pattern){
  files <- list.files(folder, pattern = pattern, full.names = T)
  r <- map(files, raster) %>%
    stack
  names(r) <- basename(files) %>%
    str_replace(fixed(".img"), "")
  r
}

# Extract data from a 3d matrix of data retrieved from a netCDF file
# with ncvar_get() (arg. dat). Create a raster and name it according to the.date.
#
# Note: the data in the y-dimension need to be flipped since they are "upside down"
#    in netcdf data.
#
#
# dat - 3D matrix of values. Dimensions are (lon, lat, month) (may have lon and lat
#   reversed)
# the.date - the date (16 of the month for "sst")
# index -  index of which month to pick off
# datname - the name of the layer (e.g. "sst") that dat represents. only used
#    to label the plot.
# x - x coords of dat
# y - y coords of dat
# inproj - projection of dat
# outproj - projection of return raster dat
create.ncdf.rast <-
  function(the.date,
           index,
           dat,
           datname,
           x,
           y,
           inproj,
           outproj = NULL,
           to = NULL) {

  message(sprintf(
      "Creating %s raster %s",
      datname,
      ifelse(as.character(the.date) == "", "", as.character(the.date))
    ))

  the.data <- switch(length(dim(dat)), NULL, dat, dat[,, index])

  # Make sure data dimensionality is sensible
  if (is.null(the.data))
    stop(paste0("create.ncdf.rast: illegal data dimension: ", dims ))

  res <-
    raster::raster(
      t(the.data),
      xmn = min(x),
      xmx = max(x),
      ymn = min(y),
      ymx = max(y),
      crs = CRS(inproj)
    ) %>%
    raster::flip(direction = "y") %>%
    rast

  # is reprojection/resampling required
  if (!is.null(to)) {
    if (!is.null(outproj))
      warning("create.ncdf.rast: both 'outproj' and 'to' are provided, ignoring outproj",
              immediate. = TRUE)

    res <- terra::project(res, rast(to), threads = TRUE) %>%
      terra::mask(study.area)
  } else if (!is.null(outproj))
    res <- terra::project(res, outproj, threads = TRUE) %>%
      terra::mask(study.area)

  plot(res, main = paste(datname, the.date))
  res
}



# Extract data from one NetCDF file
# dataset - character string giving name of dataset to extract. (eg. "sst")
# XXX TODO: figure out how to replace this with stars::read_ncdf()
ncdf.to.raster <- function(filename,
                           dataset,
                           inproj,
                           outproj = NULL,
                           to = NULL) {

  message(sprintf("Extracting %s from netCDF file: %s", dataset, filename))

  x <- ncdf4::nc_open(filename)
  lon <- ncdf4::ncvar_get(x, "longitude")
  lat <- ncdf4::ncvar_get(x, "latitude")
  dat <- ncdf4::ncvar_get(x, dataset)

  # Get dimension descriptors
  dims <- tidync::tidync(filename) %>%
    tidync::activate(dataset) %>%
    tidync::hyper_dims()

  # netCDF may have multiple "layers" - one for each time step. Note there may
  # be only 1 time step.
  if ("time" %in% dims$name)
    dates <-
      ncdf4::ncvar_get(x, "time")  %>%  # returns seconds since 01/01/1970
      `/`(86400) %>%  # convert to days since 01/01/1970
      as_date(origin = lubridate::origin)
  else
    dates <- ""

  # get rasters for each date
  res <-
    map2(
      dates,
      seq_along(dates),
      create.ncdf.rast,
      dat = dat,
      datname = dataset,
      x = lon,
      y = lat,
      inproj = inproj,
      outproj = outproj,
      to = to
  )

  # If there was more than one date then stack 'em.
  # Otherwise, just peel off the single raster
  if(length(res) > 1){
    res <- rast(res)
    names(res) <- paste(dataset, dates, sep = ".")
  } else {
    res <- res[[1]]
    names(res) <- dataset
  }
  res
}

# Reads all raster files in folder with names matching pattern and returns a raster
# or rasterstack.Called from Extract_env_rasters.Rmd
#
# folder - where to find the files
# pattern - filename pattern to match e.g. "^jplMURSST41mday.+nc$")
# variable - the name of the variable in the NetCDF to extract, e.g. "sst"
# inproj - projection of the input netcdf
# outproj - projection to reproject to
# to - a raster to be used to reproject/resample/match extent to. If both
#     "outproj" and "to" are provided, "outproj" is ignored with a warning.
import.netCDF <-
  function(folder, pattern, variable, inproj, outproj = NULL, to = NULL) {

  files <- list.files(folder, pattern = pattern, full.names = T)

  if (length(files) == 0) {
    stop("import.netCDF: no matching filenames")
  } else {
    # reading one or more
    res <-
      map(
        files,
        ncdf.to.raster,
        dataset = variable,
        inproj = inproj,
        outproj = outproj,
        to = to
      )

    # If only 1 file, just return the first element of the list.
    # Note that if the single file had multiple layers (say years etc)
    # then the single element of res can still be a rasterbrick, but thats ok. I
    # just want to avoid returning a brick when a rasterlayer is expected.
    if (length(files) > 1) {
      res <- rast(res)
    } else
      res <- res[[1]]

  }
  res
}


# Function to take a SpatRaster of monthly rasters (all months), extract the monthly
# rasters of interest and take their mean. Usefule for creating, for example,
# mean jan, feb, etc monthly averages. Used by Extract_env_rasters.rmd.
rast.monthly.mean <- function(mnth, r){
  message("Getting monthly means for month ", mnth)
  r.mnths <- names(r) %>%
    str_split_fixed(fixed("."), n = Inf) %>%
    extract(, 2) %>%
    str_split_fixed(fixed("-"), n = Inf) %>%
    extract(, 2) %>%
    as.integer

  sel <- r[[which(r.mnths == mnth)]]
  mean(sel)
}

# Produce a dotchart of a dataframe column. Used by Generic_1.5_dsm_EDA.Rmd
do.dotchart <- function(varname, dat) {
  if (varname %in% names(dat))
    dotchart(dat[, varname],
             main = varname,
             xlab = "Values of variable",
             ylab = "Order of the data")
}

# Compute seasonal mean for dynamic variable "var" (and its scaled version) from
# monthly values stored in separate columns in dat.
# Called once for each season and dynamic variable by get.seas.mean()
# var - variable of interest (e.g., "sst")
# dat - normally predgrid containing columns of monthly values variable var
# start, end - start and end MONTHS to have values averaged across.
get.seas.mean.var <- function(var, dat, start, end) {
  var.names <- paste(var, start:end, sep = ".")
  var.names.sc <- paste0(var.names, "_sc")
  ret = list(dat %>%
               select(all_of(var.names)) %>%
               rowMeans,
             dat %>%
               select(all_of(var.names.sc)) %>%
               rowMeans)
  names(ret) <- c(var, paste0(var, "_sc"))
  ret
}


# Compute species-specific seasonal means for dynamic vars listed in dyn.vars
# mean of appropriate months for each dynamic variable (and their scaled
# verions) listed in dyn.vars. Called once for each season and species from
# create.seasonal.predgrid()
#
# seas - character string given season
# season.spec - species-specific season boundaries,
# dat - normally predgrid containing columns of monthly values for various dynamic
#     variable
# dyn.vars - char vector of names of dynamic variables (e.g. "sst", "sst.g")
get.seas.mean <- function(seas, season.spec, dat, dyn.vars) {
  start <- season.spec[[seas]]["from"] %/% 100
  end <- season.spec[[seas]]["to"] %/% 100
  dat <- filter(dat, as.character(Season) == seas)
  new.cols <- map(dyn.vars, get.seas.mean.var, dat = dat, start = start, end = end)
  cbind(dat, new.cols)
}


# Save current prediction htmls to folder. Useful for saving and comparing to
# subsequent improved iterations.
#
# XXX Perhaps not needed anymore now that i'm saving results htmls in git LFS?
save.prediction.htmls <- function(folder) {
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)

  names(spec.grps) %>%
    map(function(species) {
      filename <- file.path(ResultsDir, species, paste0(species, "_3_prediction.html"))
      message(sprintf("Copying '%' in '%s'", filename, folder))
      file.copy(filename, folder)
    })
}


# Reclassify values in a terra raster, r, according to class.arg, project to
# layer, to, and save in filename. Used by Extract_env_rasters.rmd
#
# Note: class.arg can be either a 3-col matrix (from, to, becomes), a 2-col
#     matrix (is, becomes) or a 1-col matrix (or vector) which specifies cut points.
#     normally used to convert Inf to NA with class.arg = cbind(Inf, NA)
reclassify.project.save <- function(r, class.arg, to, filename){
  outdir <- unique(dirname(filename))

  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE )

  r %>%
    classify(rcl = class.arg) %>%
    terra::project(y = to,
                   method = "bilinear",
                   threads = TRUE) %>%
    terra::mask(to) %>%
    terra::writeRaster(filename = filename,
                       overwrite = TRUE)
}

# extract names for dynamic variables (and their gradients if used) from
# env_covars table
dynamic.env.covar.names <- function(){
  dyn_vars <- filter(env_covars, var_type == "dynamic")
  c(dyn_vars$var_name, paste0(dyn_vars$var_name[x$do_gradient], ".g"))
}

# Create seasonal predgrid - used by Generic_3_prediction.rmd
#
# We will need 1 copy of the basic predgrid for each season. Seasonal dynamic
# predictors are formed as the average of monthly values within each season. Seasons
# are defined on a per-species basis.
#
# Once this seasonal predgrid is created, we will again need to make multiple copies
# of it (just like segdata - 1 for Fly and 1 for Swim).

create.seasonal.predgrid <- function(species, predgrid) {
  # Get species-specific season setting and create predgrid.all.seas with 4 seasons
  season.spec <- seasons[[species]]
  ret <-
    rbind(predgrid, predgrid, predgrid, predgrid) %>%
    mutate(Season = as.factor(rep(season.names, each = nrow(predgrid))))

  # Get seasonal means for dynamic variables
  match.re <- c("[0-9]$", "[0-9]_sc$")
  p.geom <- st_geometry(ret) # save geometry
  ret <- season.names %>%
    map_dfr(
      get.seas.mean,
      season.spec = season.spec,
      dat = st_drop_geometry(ret),
      dyn.vars = dynamic.env.covar.names()
    ) %>%
    # Remove monthly values from predgrid to make it smaller now that we have
    # seasonal means computed, add area, and make back into sf object.
    select(!matches(match.re)) %>%
    cbind(p.geom) %>% # add geometry back in
    st_sf

  # Make two copies - one for Fly and one for Swim (technically only needed by
  # factor and fs models but make copy anyways and run.dsm.pred will handle it
  # correctly for nofactor model)
  ret <- rbind(ret, ret) %>%
    mutate(FlySwim =
             factor(rep(c("F", "W"), each = nrow(
               ret
             )),
             c("F", "W")))

  # From multiddf paper code:
  # create an extra column to account for the variance propagation model
  # the variance propagation adds a random effect covariate named "XX"
  # which we can safely give the value 0 now the variance has been propagated.
  # pred$XX <- matrix(0, nrow=nrow(pred), ncol=3)
  ret
}

# Do prediction maps for all four seasons, for a single model, and species.
# Called from Generic_3_prediction.rmd
#
# dat - a dataframe with 3 sets of seasonal predicions: 1 set for flying,
#     1 set for swimming, and 1 set for combined.
#
# subs - which subset of predictions to plot: flying, swimming or combined
do.pred.maps <-
  function(dat,
           model,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {


  # Get data subset
  subs <- match.arg(subs)
  dat <- filter(dat, subset == subs)
  segdata <- model$data %>%
    st_as_sf

  message(sprintf("%s, %s: Doing %s abundance prediction map for",
                  species, modname, subs))


  ret <- season.names %>%
    map(do.pred.map, dat, segdata, modname, species, subs, ...)

  names(ret) <- season.names
  ret
}

# Produce a leaflet map of Density predictions from model modname for one
# species, season, and subset (Combined, F, or W). Called from do.pred.maps()
#
# dat - predgrid (sf polygons) with NHat and Dens for all seasons
# segdata - segdata for all seasons
# samp_n - number of polygons from dat to plot. If NA, plot all. Otherwise, draw a
#    sample of samp_n from the rows of dat (after filtering by season)
do.pred.map <-
  function(season,
           dat,
           segdata,
           modname,
           species,
           subset,
           samp_n = NA) {


  message(sprintf("\t%s",season))

  # Filter by season, and create log Density for potential mapping - not currently
  # used.
  dat <- dat %>%
    filter(Season == season) %>%
    st_transform(latlongproj) %>%
    mutate(lDens = case_when(Dens == 0 ~ 0,
                             TRUE ~ log(Dens))) %>%
    select(Dens, geometry) %>%
    ms_simplify()

  segdata <- segdata %>%
    filter(Season == season) %>%
    get.combined.segdata() %>%
    st_transform(latlongproj)

  # Plot only a sample of the polygons for efficiency? Typically used for
  # testing.
  if (!is.na(samp_n)) {
    index <- sample(1:nrow(dat), size = samp_n)
    dat <- dat[index,]
  }

  # Remove ridiculously large densities b/c they mess up the legend and swamp
  # everything else
  dat <- mutate(dat,
                Dens = case_when(Dens > MAX_DENS_VALUE ~ NA,
                                 TRUE ~ Dens))

  if ((n.na <- sum(is.na(dat$Dens))) > 0)
    message("Warning: ", n.na, " cells larger than ", MAX_DENS_VALUE, " were converted to NA")

  groups <- c("est abund", "Pred Dens")
  m <-
    leaflet(
      data = dat,
      options = leafletOptions(preferCanvas = TRUE)
    ) %>%
    # Options help to speed up rendering. NOTE - dont't use addProviderTiles
    # if you want to save the map and reload in a subsequent R session - it won't
    # work.
    addTiles(options = tileOptions(updateWhenZooming = FALSE,
                                 updateWhenIdle = FALSE)) %>%
    addMapPane("density", zIndex = 410) %>%
    addMapPane("abund", zIndex = 420) %>%
    # Predicted density
    addPolygons(
      fillColor = ~ pal_pred(Dens),
      color = ~ pal_pred(Dens),
      fillOpacity = 1.0,
      opacity = 1.0,
      weight = 1,
      group = "Pred Dens",
      options = pathOptions(pane = "density")
    ) %>%
    # 0 Abund
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      stroke = FALSE,
      radius = rep(1000, times = nrow(filter(segdata, estAbund == 0))),
      color = "white",
      fillColor = "white",
      fillOpacity = 0.1,
      opacity = 0.9,
      group = "est abund",
      data = filter(segdata, estAbund == 0),
      options = pathOptions(pane = "abund")
    ) %>%
    # Est Abund
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      radius = ~ estAbund * (max.circ.radius/max(estAbund)),#  scale so largest is 100km
      color = "black",
      weight = 1,
      popup = ~ htmlEscape(paste0("Abund ", round(estAbund, 3), ", raw ", round(rawCount, 3))),
      group = "est abund",
      data = filter(segdata, estAbund != 0),
      options = pathOptions(pane = "abund")
    ) %>%
    addLegend(
      pal = pal_pred,
      opacity = 1,
      values = ~ Dens,
      # values = class_intervals$brks, # only for discrete color scale
      title = paste(species, season)
    ) %>%
    addMouseCoordinates() %>%
    addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = FALSE)) %>%
    addLayersControl(overlayGroups = groups,
                     options = layersControlOptions(collapsed = FALSE)) %>%
    hideGroup(c("est abund"))

  m
}

# Find the species grp that a species Alpha code belongs to. Used by
# Extract_data.Rmd.
find.spec.grp <- Vectorize(function(species) {
  if (is.na(species))
    return(NA)

  res <- names(spec.grps[grepl(species, spec.grps)])

  # No match
  if (length(res) == 0)
    return(NA)

  # Too many matches
  if (length(res) > 1)
    stop(sprintf(
      "%s belongs to more than one species group: %s",
      species,
      paste(res, collapse = ", ")
    ))

  res
})



# Take segdata that has copies for each ddftype and combine them by summing
# estimated densities, abundances, counts, etc.
#
# Used by do.pred.map() and create.species.shapefiles()
#
# segdata - assumed to be in format returned by create.dsm.data
#
get.combined.segdata <- function(segdata){
  # Note that the species-specific segdata has multiple copies for each original
  # segment - currently up to 3 per aerial and 3 per ship-based segments: one
  # for W[ater], one for F[lying], as one for S[trip] with Sample.Label suffixes
  # as specified by ddftype.levels ins analysis_settings.R, so we first need to chop these suffixes to do the filter and
  # then sum the estDen[sities].
  segdata <-
    select(
      segdata,
      Sample.Label,
      estDens,
      estAbund,
      rawCount,
      SurveyType,
      x,
      y,
      LatStart,
      LongStart,
      geometry,
      Season,
      Date,
      ddftype
    ) %>%
    # Remove Sample_Label suffix (_xx)
    mutate(Sample.Label = str_sub(Sample.Label, 1, nchar(Sample.Label) -
                                    nchar(ddftype.levels[1]) - 1)) %>%
    arrange(Sample.Label)

  # There are now multiple consecutive rows for each segment.

  # Figure out how many copies of each segment there are. Complain if not all the
  # same. This relies on having a fixed symmetric setup where, for example,
  # the aerial segments are copied the same number of times as the ship ones
  # (ie once each for Water, Fly, and Strip)
  rl <- rle(segdata$Sample.Label)
  stopifnot(length(unique(rl$lengths)) == 1)
  ncopies <- rl$lengths[1]

  # Keep every nth row, summing values within each group of ncopies segdata rows,
  # and make it a column in  new data frame containing only every nth segdata row
  # in order to get one row per segment. Note that rollapply() uses the full
  # segdata (before selecting every nth row)

  keep <- rep(c(TRUE, rep(FALSE, times = ncopies - 1)),
              times = nrow(segdata) / ncopies)
  res <- segdata[keep,] %>%
    mutate(
      estDens = zoo::rollapply(segdata$estDens, ncopies, by = ncopies, sum),
      estAbund = zoo::rollapply(segdata$estAbund, ncopies, by = ncopies, sum),
      rawCount = zoo::rollapply(segdata$rawCount, ncopies, by = ncopies, sum)
    )

  res
}

# Sum fly and water detection-corrected segment densities for a given species in
# each segment (optionally limited to only those sample labels in sample.labs)
# and save as a set of seasonal shapefiles in folder. Called from Generic_2_dsm.Rmd
create.species.shapefiles <-
  function(spec,
           sample.labs = NULL,
           folder,
           segdata,
           distdata) {

    # Get segdata for spec with fly/water combined.
    segdata <- get.combined.segdata(segdata)

    # Restrict to certain sample labels (useful for cropping to only those samples
    # in a certain spatial area)
    if (!is.null(sample.labs))
      segdata <- filter(segdata, Sample.Label %in% sample.labs)

    # save as seasonal shapefiles if required
    for (seas in season.names) {
      layer <- paste(spec, seas, "segdata", sep = "_")

      # if shapefile doesn't exist or recreateSpecSegShapefiles is TRUE then
      # save the shapefile
      if (!file.exists(paste0(folder, "/", layer, ".shp")) ||
          recreateSpecSegShapefiles) {
        message(sprintf("Creating segdata shapefile for %s %s", spec, seas))

        filter(segdata, Season == seas) %>%
          st_write(
            dsn = folder,
            layer = layer,
            driver = "ESRI Shapefile",
            delete_layer = T
          )
      }
    }

    #### Now do same for distdata.
    # Distdata isn't seasonal b/c it is used in its entirety for the ddf (but
    # with season as covar if needed)
    layer <- paste(spec, "distdata", sep = "_")

    if (!file.exists(paste0(folder, "/", layer, ".shp")) ||
        recreateDistdataShapefile) {
      message(sprintf("Creating distdata shapefile for %s", spec))

      # Save distdata as a shapefile.
      distdata %>%
        # Needed since numeric ids can get too big for shapefile numbers
        mutate(object = as.character(object)) %>%
        st_as_sf(
          coords = c("LongStart", "LatStart"),
          crs = st_crs(4326),
          remove = FALSE
        ) %>%
        st_write(
          dsn = ShapeDir,
          layer = layer,
          driver = "ESRI Shapefile",
          delete_layer = TRUE
        )
    }
  }

# Return names of prediction raster file (and aux.xml file)
# for a single species/season pair. Used by 03b_Save_chosen_model_predictions.Rmd
get.final.prediction.name <- function(spec, season){

  # Get final model predictions raster filename
  modname <- final.dsm.models$dsm_final_name[final.dsm.models$species == spec]
  filenames <- vector(mode = "list", length = 2)
  filenames[[1]] <- sprintf("%s.%s.%s.%d_sqkm.img", spec, season, modname, predgridCellArea)
  filenames[[2]] <- sprintf("%s.%s.%s.%d_sqkm.img.aux.xml", spec, season, modname, predgridCellArea)
  filenames
}

# undebug all debugged functions. This is tricky and sometimes doesn't work.
undebug.all <- function(where=search()) {
  aa <- all_debugged(where)
  lapply(aa$env,undebug)
  ## now debug namespaces
  invisible(mapply(function(ns,fun) {
    undebug(getFromNamespace(fun,ns))
  },names(aa$ns),aa$ns))
}


# Create  2x2 table of maps - one panel for each season.
# XXX Not currently used
create.annual.map.grid <- function(maps) {
  if(length(maps) != 4)
    stop(sprintf("create_annual_map_grids: maps argument does contains %d maps - should be 4."),
         length(maps))

  res <-
    tagList(tags$table(
      style = "width:100%",
      tags$tr(tags$td(tagList(maps$Spring)),
              tags$td(tagList(maps$Summer))),
      tags$tr(tags$td(tagList(maps$Fall)),
              tags$td(tagList(maps$Winter)))
    ))

  res
}

# Save a list of maps  (typically 4 seasonal prediction maps) from model
# modname to a summary folder. Used by Generic_3_prediction.Rmd.
save.map <- function(maps, modname, species) {
  dirname <- here(ResultsDir, species, "Prediction summaries")
  if (!dir.exists(dirname))
    dir.create(dirname, recursive = TRUE)

  filename <- here(dirname, paste0(modname, ".html"))
  message(sprintf("%s, %s: Saving map in %s.", species, modname, filename))
  list(h2(paste0(species, "_", modname)),
       leafsync::sync(maps)) %>%
    tagList %>%
    save_html(file = filename)
}

# Make quick and dirty prediction maps with ggplot. Called from
# Generic_3_prediction.Rmd
do.pred.maps.ggplot <-
  function(dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {


    # Get data subset
    subs <- match.arg(subs)
    dat <- filter(dat, subset == subs)

    message(sprintf("%s, %s: Doing %s abundance prediction map for",
                    species, modname, subs))

    ret <- season.names %>%
      map(do.pred.map.ggplot, dat, modname, species, subs, ...)

    ret <- wrap_plots(ret) + plot_annotation(title = modname)
    ret
  }

# Plot a predicton ggplot for one model with 4 seasons. Called from
# do.pred.maps.ggplot()
do.pred.map.ggplot <-
  function(season,
           dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           samp_n = NA) {

    message(sprintf("\t%s",season))
    dat <- dat %>%
      filter(Season == season) %>%
      select(Dens, geometry) %>%
      ms_simplify()


    # Plot only a sample of the polygons for efficiency? Typically used for
    # testing.
    if (!is.na(samp_n)) {
      index <- sample(1:nrow(dat), size = samp_n)
      dat <- dat[index, ]
    }

    # Remove ridiculously large densities b/c they mess up the legend and swamp
    # everything else
    dat <- mutate(dat,
                  Dens = case_when(Dens > MAX_DENS_VALUE ~ NA,
                                   TRUE ~ Dens))

    ret <- ggplot(dat = dat) +
      geom_sf(aes(fill = Dens), color = NA)  +
      scale_fill_continuous(
        type = "viridis"
        # breaks = class_intervals$brks,
        # labels = round(class_intervals$brks, 4)
      ) +
      ggtitle(season)
    ret
  }


# Quick and dirty leaflet map for watches. Used when doing debugging of watch
# data in Extract_data.Rmd
watch.map <- function(dat) {
  leaflet(dat) %>%
    addTiles() %>%
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      popup = ~ htmlEscape(paste0("WatchID ", WatchID))
    )
}

# Recursively create a pathname if it doesn't exist. Return pathname invisibly
# so this function can be used in a pipeline.
create.dir.if.needed <- function(pathname){
  if (!dir.exists(pathname)){
    dir.create(pathname, recursive = TRUE)
    message(sprintf("Creating needed folder %s", pathname))
  }
  return(invisible(pathname))

}


# Create an initial generic ddf model list for each species group, optionaly
# writing it to dfModlistLoc (default false). Called from
# Create_final_ddf_model_specs.Rmd and Generic_1_ddf_fitting.Rmd
init.df.mod.list <- function(saveit = FALSE) {

  # Create one copy of the basic list for each species group
  df.mod.list <- rep(list(def.ddf.list), length(spec.grps))
  names(df.mod.list) <- names(spec.grps)

  # Set default convert_units value based on whether survey is aerial or ship.
  # This may get updated dynamically by later processing.
  df.mod.list %<>%
    map( ~ imap(., set.def.df.spec.values))

  # Sometimes we want to save it (ie if it didn't already exist) but other times
  # we're just called to return a initial structure that can be modified (eg in
  # 01.02_Create_final_ddf_model_specs.Rmd)
  if (saveit)
    save(df.mod.list, file = dfModlistLoc)

  df.mod.list
}

# Take a dataframe that contains at least DistMeth and FlySwim and assign
# the distance type (None, Perp., or Radial) by looking it up in lkpDistMeth
# in ECSAS database.
#
# Used by Extract_data.Rmd and create.survey.data()
assign.dist.type <- function(dat) {
  distmeth <- ECSAS.get.table(ecsas.path = ECSAS.Path, "lkpDistMeth")
  DistType <- left_join(dat, distmeth, by = c("DistMeth" = "DistMethCode")) %>%
    mutate(DistType = as.factor(
      case_when(
        FlySwim == "F" ~ DistMethFlyHow,
        FlySwim == "W" ~ DistMethWaterHow,
        TRUE ~ NA_character_
      )
    )) %>%
    pull(DistType)

  if (any(is.na(DistType)))
    warning("assign.dist.type: ",
            sum(is.na(DistType)),
            " rows could not be assigned a distance type",
            immediate. = TRUE)

  DistType
}

# Generate the cononical model name given a key, formula and adj term.
#
# Used by do.det.fcn().
#
create.model.name <- function(key, form, adj) {
  # Was there a null formula?
  if (form == as.formula(~ 1)) {
    # Grab adjustment, converting names to shorthand used in filenames
    form <-
      if (is.null(adj))
        NULL
    else
      c(cos = "cos",
        herm = "hp",
        poly = "sp")[adj]
  } else {
    # Convert formula to dotted notation
    form <-
      str_replace(as.character(form)[2], fixed(" + "), ".")
  }
  # Deal with adjustment that might be NULL with ifelse()
  nam <- paste0(key, ifelse(is.null(form), "", paste0(".", form)))
  nam
}


## A few functions to make backing up and cleaning up DSM summaries easier

# Make backup copy of current DSM summaries for species labeling them with their
# modification date/time.
backup.dsm.summary <- function(species){
  folder <- here(ResultsDir, "Backups", species, "DSM Summaries")
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)

  message("Backing up ", species, " DSM summaries to ", folder)

  # src files
  src <- list.files(path = here(ResultsDir, species, "DSM Summaries"),
                      pattern = "*.Rdata", full.names = T)

  if (length(src) == 0){
    message("\tNo files to backup - quitting.")
    return
  }

  ## create dst filenames

  # add modification dates to files
  mtimes <- str_replace_all(file.info(src)$mtime, " ", "_") %>%
    str_replace_all(":", "") %>%
    str_replace("\\..*$", "") # remove trailing milliseconds
  stopifnot(length(src) == length(mtimes))
  dst <- src %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    here(folder, .) %>%
    paste0("_", mtimes, ".", tools::file_ext(src))

  file.copy(src, dst, overwrite = TRUE, copy.date = TRUE) %>%
    invisible
}

# Do DSM summaries backup for all species
backup.dsm.summaries <- function(){
  names(spec.grps) %>%
    walk(backup.dsm.summary)
}


# Remove all current DSM summaries for all species
remove.dsm.summaries <- function(){
  names(spec.grps) %>%
    walk(\(species){
      message("Removing DSM summaries for ", species)
      files <- list.files(path = here(ResultsDir, species, "DSM Summaries"),
                 pattern = "*.Rdata", full.names = T)
      file.remove(files)
      })
}


# Get the dates of a dynamic env covariate given a variable name. Useful to see
# what files we need to download.
# Returns a dataframe with two columns: date, filename (typcially containing
# 12 rows - 1 for each month)
get.covar.netCDF.dates <- function(var.name){
  # list files in the folder here(predLayerDir, "NetCDF", var.name)
  # get dates associated with files, convert to char and return
  #
  files <- list.files(here(predLayerDir, "NetCDF", var.name),
                      ".*\\.nc$",
                      full.names = TRUE)

  suppressWarnings(
    dates <- map_dfr(
      files,
      \(filenm) {
        data.frame(
          date = stars::read_ncdf(filenm, var = "time", proxy = TRUE) %>%
            stars::st_get_dimension_values("time") %>%
            str_replace(" UTC", ""),
          filename = filenm
        )
  }))

  dups <- duplicated(dates$date)
  if(any(dups)) {
    warning(sprintf(c("get.covar.netCDF.dates: variable %s: the following dates",
                      " were found in multiple files: "), var.name), immediate. = TRUE)
    dates[dups,]
  }

  dates
}

# called once with each row of env_covars to get netcdf if needed, extract
# data from it, reproject and clip, etc.
get.env.covar <- function(env_covar_spec,
                          dates.needed,
                          verbose = TRUE
) {


  var_name <- env_covar_spec$var_name

  message("get.env.covar: getting covariate '", var_name, "'")

  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env_dat <- stars::read_ncdf(here(
      predLayerDir,
      "NetCDF",
      paste0(env_covar_spec$ERDDAP_dataset_name, ".nc")
    ), var = env_covar_spec$netcdf_vars)


  } else if (env_covar_spec$var_type == "dynamic") {
    # Dynamic involves making sure we have the dates we want.
    netcdf.dates.have  <- get.covar.netCDF.dates(var_name)
    dates.to.get <- setdiff(dates.needed, netcdf.dates.have$date)

    if (length(dates.to.get) > 0) {
      # Download needed files
      # XXXX TODO: need to use curl (or something ) to download the netcdf files
      # we need by constructing the correct URL.
    }

    # Now we should have all files we need.
    # Figure out which files to read, just in case we have files we don't need
    files <- filter(netcdf.dates.have, date %in% dates.needed) %>%
      pull(filename) %>%
      unique()

    # Read all needed netcdf files
    #
    # # This approach didn't work. Something to do with fact that some of the
    # .nc files have only one var in them ("sst") and some have both ("sst" and
    # "mask") due to different ways they were downloaded.
    #
    # res <- map(files, \(filenm) {
    #   stars::read_ncdf(filenm, var = var_name)
    # })
    #
    # # Collapse list of stars objects to single object with all dates combined.
    # env.dat <- Reduce(c, res)
    env_dat <- stars::read_stars(files, sub = var_name)
    if(is.na(st_crs(env_dat)))
      st_crs(env_dat) <- env_covar_spec$CRS
  } else
    stop("get.env.covar: variable ",
         var_name,
         ": illegal var_type '",
         env_covar_spec$var_type, "'")

  if (verbose)
    message("\tProjecting and clipping to study area...", appendLF = FALSE)

  # Re-project, clip, etc
  final <- st_transform(env_dat, segProj) %>%
    `[`(study.area) %>%
    setNames(var_name) # b/c previous processing steps loose the name

  if (verbose)
    message("done.")

  final
}



# called once with each row of env_covars to get netcdf if needed, extract
# data from it, reproject and clip, etc.
#
do.env.covar <- function(env_covar_spec,
                         segdata.dates.needed,
                         predgrid.dates.needed) {
  var_name <- env_covar_spec$var_name

  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env.dat <- stars::read_ncdf(here(
      predLayerDir,
      "NetCDF",
      paste0(env_covar_spec$ERDDAP_dataset_name, ".nc")
    ), var = env_covar_spec$netcdf_vars)


  } else if (env_covar_spec$var_type == "dynamic") {
    # Dynamic involves making sure we have the dates we want.
    netcdf.dates.have  <- get.covar.netCDF.dates(var_name)
    all.dates.needed <- union(segdata.dates.needed, predgrid.dates.needed)
    dates.to.get <- setdiff(all.dates.needed, netcdf.dates.have)

    if (length(dates.to.get) > 0) {
      # Download needed files
      # XXXX TODO: need to use curl (or something ) to download the netcdf files
      # we need by constructing the correct URL.
    }

    # Now we should have all files we need.
    # Figure out which files to read
    files <- filter(netcdf.dates.have, date %in% all.dates.needed) %>%
      pull(filename) %>%
      unique()

    # Read all needed netcdf files
    res <- map(files, \(filenm) {
      stars::read_ncdf(filenm, var = var_name)
    })

    # Collapse list of stars objects to single object with all dates combined.
    # xxx this doesn't work, returned object seems to only have "attr" as it's
    # values
    env.dat <- Reduce(c, res)
  } else
    stop("do.env.covar: variable ",
         var_name,
         ": illegal var_type '",
         env_covar_spec$var_type, "'")

  # Re-project, clip, etc and save - xxx this line doesn't work
  final <- st_transform(env.dat, segProj) %>%
    `[`(study.area)

  final
}

# Aggregate observations by watch and then reproject them to grid
# obs - observations
# obs.name - name for the observation data (e.g., "tbmu")
# watches - watches
# grid - grid to project to
# time.period - character string defining name of time.period the data is for.
#     Typically "year_round, "Spring", "Summer", "Fall", or "Winter". Used to
#     name the resulting shapefile
agg.by.grid <- function(obs,
                        obs.name,
                        watches,
                        grid,
                        time.period,
                        save.shapefile = FALSE,
                        save.RDS = FALSE) {

  message("Aggregating ", obs.name, " during ", time.period)

  # Join the watch and obs data
  dat <- watches %>%
    mutate(year = year(Date)) %>%
    left_join(obs, by = "WatchID")

  # sum number of birds and number of obs by watch.
  # Set tot_size and n_obs to 0 for watches with no obs.
  agg.data <- dat %>%
    group_by(WatchID) %>%
    summarise(year = unique(year), nbirds = sum(size), n_obs = n()) %>%
    mutate(nbirds = case_when(is.na(nbirds) ~ 0, .default = nbirds),
           n_obs = case_when(nbirds == 0 ~ 0, .default = n_obs)) %>%
    left_join(select(watches, WatchID, LatStart, LongStart), by = "WatchID") %>%
    st_as_sf(
      coords = c("LongStart", "LatStart"),
      crs = st_crs("EPSG:4326"),
      remove = FALSE
    ) %>%
    st_transform(proj) %>%
    vect()

  ### NOTE: choose final field names that will be shapefile friendly
  ### when combined with 4-letter species codes.

  # pack years into a single number (product of primes) for aggregating to
  # raster and then unpack as a point object.
  years <- rasterize(agg.data, grid, field = "year", fun = encode_years) %>%
    as.points() %>%
    st_as_sf() %>%
    mutate(!!paste0(obs.name, "_yr") := decode_years(year)) %>%
    select(-year)

  # rasterize nbird summing all points that fall in each cell and convert
  # to points
  nbirds <- rasterize(agg.data, grid, field = "nbirds", fun = sum) %>%
    as.points() %>%
    st_as_sf() %>%
    mutate(sum = as.integer(sum)) %>%
    rename(!!paste0(obs.name, "_tot") := sum)

  # rasterize summing all points that fall in each cell n_obs and convert to points
  nobs <- rasterize(agg.data, grid, field = "n_obs", fun = sum) %>%
    as.points() %>%
    st_as_sf() %>%
    mutate(sum = as.integer(sum)) %>%
    rename(!!paste0(obs.name, "_n") := sum)

  # combine and rename columns
  agg.data.point <- nbirds %>%
    cbind(nobs %>% st_drop_geometry(),
          years %>% st_drop_geometry()
    )

  # Save
  filename <- paste("ECSAS", obs.name, "abundance", time.period, sep = "_")

  if (save.shapefile)
    st_write(
      agg.data.point,
      # dsn = file.path(share_drive, "../products/shapefiles"),
      dsn = ShapeDir,
      layer = filename,
      driver = "ESRI Shapefile",
      delete_layer = T
    )

  if (save.RDS)
    saveRDS(agg.data.point,
            file = file.path(GenDataDir, paste0(filename, ".rds")))

  invisible(agg.data.point)
}

# Encode a vector of years into a single integer by collapsing the 4-digit years
# NOTE - this assumes all years >=  base and < base + 1000 (since decoding uses)
# the list of the first 1000 primes.
#
# Works by subtracting base from each unique year and adding 1 to get an index into
# the list of nth prime numbers. Then multiply these prime numbers together
# to get the answer. To decode simply get the prime factors of the number
# and use them to figure out what position (ie index) they are in the
# list of the first 1000 primes. Then un-standardize the index to recover the year.
# E.g. years = c(2020, 2005)
#   subtracting 2000 and adding 1 gives c(21, 6) therefore prms is a vector
#   containing the 21st and 6th primes (ie 73 and 13).
#   - multiplying these gives ans == 949
encode_years <- function(years, base = 2000){
  yrs <- unique(years)

  # standardize to 1 == base, 2 == base + 1, etc and use as index to choose
  # primes
  prms <- ((na.omit(yrs) - base) + 1) %>%
    primes::nth_prime() %>%
    as.numeric()

  # multiply primes together to get a single number. Easy to recover
  # the primes later since there is only 1 prime factorization of a number
  # resulting from mutliply all primes.
  ans <- Reduce(`*`, prms)

  # Just in case something funky has gone on
  if (is.na(ans) || is.null(ans))
    return(0)

  ans
}


# Create a vectorized version of gmp::factorize for use in decode_years
factorize.vec <- Vectorize(gmp::factorize)

# Convert number representing encoded years back to years
# Uses prime_factorization to get the prime factors of num and then finds
# their indices in the first 1000 primes then unstandardizes by subracting 1
# and adding base.
# E.g. num = 949
#   prime factors are 73 and 13
#   these are the 21st and 6th primes
#   unstandardizing recovers 2020 and 2005
#
# For datapoints with many years the value of num was too big for
# primes::prime_factors() to handle, so I had to switch to using
# a vectorized verseion of gmp::factorize()
#
# Return a character where each element is a comma separated string of years
decode_years <- function(num, base = 2000){

    zero.indices <- num == 0
    res <- vector("character", length(num))
    yrs <- num[!zero.indices] %>%
      factorize.vec() %>%
      map(\(x) match(x, primes::primes) - 1 + base) %>%
      map_chr(\(x) paste(x, collapse = ", "))

    # Just in case
    res[zero.indices] <- NA
    res[!zero.indices] <- yrs
    res
}

# Aggregate observations dat on watches for a particular taxa (grp) to
# the raster (grid) both seasonally and year-round. Returns a list of
# (typically) 5 sf point objects: 1 for year-round and one each for the
# (tyipcally) 4 seasons.
agg_species_group <- function(grp, dat, watches, grid){
  # get species of interest
  if (grp == "Sbrd")
    dat.filt <- filter(the.data$distdata, Seabird == -1)
  else if (grp == "Wbrd")
    dat.filt <- filter(the.data$distdata, Waterbird == -1)
  else
    dat.filt <- filter(dat, Alpha %in% spec.grps[[grp]])

  # Do year round
  yr <- list(agg.by.grid(dat.filt, grp, watches, grid, "year_round")) %>%
    setNames("year_round")

  # Now do seasonally
  seas <- watches %>%
    split(.$Season) %>%
    map(\(seas.watches) {
      agg.by.grid(dat.filt,
                  grp,
                  seas.watches,
                  grid,
                  unique(seas.watches$Season))
    })

  append(yr, seas)
}

# Do (possibly weighted) kernels
# dat - sf point object
# species - eg "tbmu"- just for progress message and file naming
# weights - vector of weights - one value per row in dat.
# bw - bandwidth in units of dat sf layer
# res - resolution of kernel grid in units of dat sf layer
# percent - percent volume contour to produce (e.g 50)
# time.period - time period covered by dat (e.g. "year_round", "spring", etc)
#    for output file naming
# shapedir - folder to save shapefiles in
do.kde <- function(dat,
                   species,
                   weights = NULL,
                   bw = NULL,
                   res = 1000,
                   percent,
                   time.period,
                   shapedir
) {
  wghted <- ifelse(is.null(weights), "unweighted", "weighted")
  message("Creating ",  wghted,
          " kernel for ", species, " during ", time.period,
          "\nbandwidth = ", bw,
          "\noutput grid resolution = ", res,
          "\nvolume percent(s) = ", paste(percent, collapse = ", ")
  )

  # Re-project
  dat <- st_transform(dat, proj)
  plot(st_geometry(dat), pch = ".", main = "raw data")

  # Do weighting if requested.
  if (!is.null(weights))
    dat <- dat %>%
    splitstackshape::expandRows(count = weights,
                                count.is.col = FALSE,
                                drop = FALSE)

  # Choose bandwidth automatically?
  if (is.null(bw))
    bw <- dat %>%
    st_geometry() %>%
    as.ppp() %>%
    bw.diggle()

  st_write(
    dat,
    # dsn = file.path(share_drive, "../products/shapefiles"),
    dsn = shapedir,
    layer = paste("ECSAS", species, time.period, "kernel points", sep = "_"),
    driver = "ESRI Shapefile",
    delete_layer = TRUE
  )

  # create adehabitatHR kernel
  dat.sp <- dat %>%
    mutate(ID=1) %>%
    as_Spatial() %>%
    subset(select = "ID" )

  k <- kernelUD(dat.sp, h = bw, grid = res)
  vud <- getvolumeUD(k)

  # Funciton to get percent polygon of vud for a single percentage
  create.vud.poly <- function(percent, vud) {
    message("Creating ", percent, " percent UD polygons...", appendLF = FALSE)
    poly <- getverticeshr(vud, percent)
    message("done.")
    poly
  }

  # Get percent polygons
  polys <- map(percent, \(perc) {
    create.vud.poly(perc, vud = vud) %>%
      st_as_sf()
  }) %>%
    setNames(as.character(percent))

  filename <- paste("ECSAS",
                    species,
                    wghted,
                    "kernUD",
                    time.period,
                    paste0("bw", bw),
                    paste0("res", res),
                    sep = "_")

  # Save each percentage polygon
  imap(polys, \(poly, pname) {
    st_write(
      poly,
      dsn = shapedir,
      layer = paste0(filename, "_", pname, "pct"),
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
  })

  saveRDS(polys,
          # file = file.path(share_drive, "ECSAS", paste0(filename, ".RDS")))
          file = file.path("raw-data", paste0(filename, ".RDS")))
  polys
}


# Function to combine like seasons across all species in dat.
# dat is a list of species where each species element is a list of seasons
# and each season is a dataframe for that species-season combo.
combine_species <- function(dat) {
  message("Combining species")
  # Check that all lists have the same names for the elements
  if (any(sapply(dat, function(x) !all(names(x) == names(dat[[1]]))))) {
    stop("All input lists must have the same structure (same names).")
  }

  # Combine the lists for each season (spring, summer, fall, winter)
  combined <- lapply(names(dat[[1]]), combine_season, dat)

  # Return a list with the combined values for each season
  names(combined) <- names(dat[[1]])
  return(combined)
}

# Combine a single season from all species.
# dat is a list of species where each species element is a list of seasons
# and each season is a dataframe for that species-season combo.
combine_season <- function(season, dat) {
  # Get the geometry for this season from 1st species in dat
  geom <- st_geometry(dat[[1]][[season]])

  # Extract given season from each species
  season_values <- sapply(dat, function(species){
    st_drop_geometry(species[[season]])
  }, simplify = FALSE)

  # Combine the columns for all species and reset the geometry that was removed
  # earlier
  ret <- Reduce(cbind, season_values) %>%
    st_set_geometry(geom)
  ret
}


# Take a named list of seasonal sf objects and save each one as a shapefile.
save_seasonal_shapefiles <- function(seas_list){
  imap(seas_list, \(dat, season) {
    message("Saving shapefile for ", season)
    st_write(dat,
             dsn = ShapeDir,
             layer = paste("ECSAS", season, "aggregated_obs", sep = "_"),
             driver = "ESRI Shapefile",
             delete_layer = TRUE
    )

  })

  invisible(seas_list)
}



# Function to create a nicely formatted HTML list with indentation
generate_html_species_list <- function(spec.grps) {
  # Create a list of HTML content
  html_content <- lapply(names(spec.grps), function(group_name) {
    group_items <- spec.grps[[group_name]]

    # Create an unordered list for each group
    group_html <- htmltools::tags$ul(
      # Group name with a smaller indentation
      htmltools::tags$li(style = "font-weight: bold; margin-left: 0;", group_name),  # Add group name as list header

      # List items with more indentation
      lapply(group_items, function(item) {
        htmltools::tags$li(style = "margin-left: 20px;", item)  # Add each item under the group with more indentation
      })
    )
    group_html
  })

  # Combine all lists into a single HTML page
  html_page <- htmltools::tags$html(
    htmltools::tags$head(),
    htmltools::tags$body(
      # htmltools::tags$h2("Multi-species groups"),
      do.call(htmltools::tags$div, html_content)  # Add all groups to the body
    )
  )

  # Return the HTML content as a browsable object
  htmltools::browsable(html_page)
}
