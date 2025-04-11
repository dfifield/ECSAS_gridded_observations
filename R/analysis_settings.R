### Global settings #############################################

# Generates a warning if partial object column name matching with $ occurs.
# Prevents things like accessing the return value of dht() with
# obj$individuals$N from matching obj$individuals$Nhat.by.sample if column N
# doesn't exist (which it won't if the area passed to dht() is 0)
options(warnPartialMatchDollar = TRUE)
options(tibble.width = Inf) # All tibble columns shown
options(todor_patterns = c("FIXME", "TODO", "xxx", "XXX"))

## Load the packages
library(todor)
library(terra) # do not change the order of this one
library(tidyverse)# do not change the order of this one
library(lubridate)
library(geosphere)
library(Distance)
library(dsm)
library(knitr)
library(sp)
library(viridis)
library(stringr)
library(ggplot2)
library(here)
library(furrr)
library(ECSASconnect)
library(rmarkdown)
library(htmltools)
library(sf)
library(checkmate)
library(dftools)

# stop summarise's new habit of telling us it has dropped the last grouping level (its default behaviour)
options(dplyr.summarise.inform = FALSE)

#### files and folders ####
(RDir <- here("R"))
(DataDir <- here("Data"))
(RawDataDir <- here("Data/Raw"))
(GenDataDir <- here("Data/Generated"))
(RDataDir <- GenDataDir)
(GISDir <- here("GIS"))
(RasterDir <- file.path(GISDir, "Rasters"))
(ShapeDir <- file.path(GISDir, "Shapefiles"))
(ResultsDir <- here("Results"))

# Data locations
# ECSAS database location
(ECSAS.Path <- file.path(RawDataDir, "Master_ECSAS_v_3.73.mdb"))
(the.data.file <- file.path(GenDataDir, "ECSAS_ship_the_data.rds"))
(study.area.file <- file.path(GenDataDir, "study_area.rds"))
(proj.file <- file.path(GenDataDir, "projection.rds"))
(grid.file <- file.path(GenDataDir, "grid.rds"))

# Data options
ECSAS.reextract <- FALSE # Requery the database or use existing saved data.
aerial.watch.len <- 30 # seconds

spec.grps <- list(
  ATPU = "ATPU",
  BLKI = "BLKI",
  COEI = "COEI",
  KIEI = "KIEI",
  DOVE = "DOVE",
  Gull = c(
    "BHGU",
    "BOGU",
    "BTGU",
    "CAGU",
    "FRGU",
    "GBBG",
    "GLGU",
    "GWGU",
    "HEEG",
    "HERG",
    "HYGU",
    "ICGU",
    "IVGU",
    "KUGU",
    "LAGU",
    "LBBG",
    "LIGU",
    "MEGU",
    "RBGU",
    "ROGU",
    "SAGU",
    "SBGU",
    "SMGU",
    "THGU",
    "UNGU",
    "UNLA",
    "UNWW",
    "VEGU",
    "WEGU",
    "YLGU"
  ),
  # LESP = "LESP",
  Murr = c("COMU", "TBMU", "UNMU"),
  NOFU = "NOFU",
  NOGA = "NOGA",
  Ptrl = c("LESP", "WISP", "UNSP"),
  RAZO = "RAZO",
  Swtr = c(
    "AUSH",
    "BASH",
    "COSH",
    "GRSH",
    "MASH",
    "SMSH",
    "SOSH",
    "TOSH",
    "UNCS",
    "UNSH",
    "YESH"
  )
)

spec.grps.names <- list(Gull = "Gulls", Murr = "Murres", Ptrl = "Storm-Petrels",
                        Swtr = "Shearwaters")

# messed up watches removed
# ECSAS.remove <- c(-1997427588, 1865804800)
ECSAS.ship.cruise.remove <- c(
  # 1695932151, # Only Nova Scotia industry cruise
  # -157502116  # Ewan Wakefield RRS discovery - problematic watchlenkm
)

ECSAS.ship.watch.remove <- c(
  # -1780495786, -883836309 # Two messed up watchlenkm
)


# date range of seasons: FORMAT is month * 100 + day. These are used to assign
# seasons to distdata and segdata.
season.names <- c("Fall", "Spring", "Summer", "Winter")
seasons <- list(Spring = c(from = 401, to = 531),
              Summer = c(from = 601, to = 831),
              Fall = c(from = 901, to = 1031),
              Winter = c(from = 1101, to = 331))

#### detection function info ####
# No longer used - set in Do_final_ddf_fitting.Rmd
# df.options <- list(
#   do.final.only = TRUE, # Should we fit all possible models or just the final chosen one. Useful for knitting final version - set to T.
#   rerun = FALSE, # rerun models that already have results in folder
#   parallel = TRUE, # turn off parallel processing if it gives trouble
#   nCores = 55
# )

# filename suffix for df model
(dfModSuffix <- "detfcn.results.rda") # XXX not used anymore??? was for storing ddf

# Distance bin cut points will differ depending on distance sampling bins
# All are in km.
#ecsas.air.ctpoints <- c(0, .05, .1, .2, .3, .5) # Unneeded since they vary?
#
# Note that these convert_units values are no longer actually used because we
# only use ds() to compute the detection function and not
# abundance. Also, it won't get used by dsm() (or more correctly dsm:::make_data())
# either b/c I only use the default value of convert.units ( ie. 1) in calls to
# rum.dsm.model (and thsu dsm()). See note in docs for create.dsm.data() to see
# how convert.units is handled in creation of segment.area.
ecsas.air.convert.units = 1
ecsas.air.distance.centers <- c(0.025, 0.075, 0.15, 0.25, 0.4)
ecsas.air.right.trunc <- 0.5 #

ecsas.ship.distance.centers <- c(0.025, 0.075, 0.15, 0.25)
ecsas.ship.ctpoints <- c(0, .05, .1, .2, .3)
ecsas.ship.convert.units = 0.5 # accounts for 1-sided transects

# Default Covariates for DF model if not specified
# dfCovars <- c("Season", "size", "Windspeed", "Visibility", "ObserverName")
dfCovars <- c("Season", "size")

# Detection function model def. This is used as a base structure for df.mod.list.
# There is one of this type structure for each ddf model (ie. "platform").
generic.ddf.spec <- list(
  transect_width = NA,
  remove.fishing = NA,
  strip = FALSE,
  candidate.keys = c("hr", "hn"),
  candidate.covars = dfCovars,
  max.adjustments = 2,
  final.key = NA,
  final.formula = as.formula("~1"),
  final.adj = NULL,
  final.nadj = NULL,
  convert_units = NA,
  cutpoints = NULL,
  distance.centers = NA,
  notes = NA
)

# special case for the strip transects for obs with no perp distance. only difference
# is that strip = TRUE
generic.strip.ddf.spec <- list(
  transect_width = NA,
  remove.fishing = NA,
  strip = TRUE,
  candidate.keys = c("hr", "hn"),
  candidate.covars = dfCovars,
  max.adjustments = 2,
  final.key = NA,
  final.formula = as.formula("~1"),
  final.adj = NULL,
  final.nadj = NULL,
  convert_units = NA,
  cutpoints = NULL,
  distance.centers = NA,
  notes = NA
)

# setup default ddf model specs These each represent different observation
# process and are the same as "platforms" in Miller et al. (2021). Initialize
# each one to be either a normal ddf spec or a strip transect one. Note: values
# in these specs are changed in Create_final_ddf_model_specs.Rmd
def.ddf.list <- list(
  ECSAS_Ship_W = generic.ddf.spec,
  ECSAS_Ship_F = generic.ddf.spec,
  ECSAS_Aerial_W = generic.ddf.spec,
  ECSAS_Aerial_F = generic.ddf.spec,
  ECSAS_Ship_S = generic.strip.ddf.spec,  # Strip transect - no perp. distances
  ECSAS_Aerial_S = generic.strip.ddf.spec # Strip transect - no perp. distances
)

# Different types of ddfs in use. Used to characterize segdata copies and by dsm
# as a factor in gams.
# SW = Ship Water
# SF = Ship FLy
# AW = Aerial Water
# AF = Aerial Fly
# SS = Ship strip
# AS = Aerial strip
#
# Note the order is important as this will be the order of factor level numbering
# for ddfobj in distdata and segdata.
ddftype.levels <- c("SW", "SF", "AW", "AF", "SS", "AS")

# beaufort conversion taken from ECSAS protocol Appendix III
beaufort_conversion <- tribble(
  ~beaufort, ~speed.kts,
  0, 0,
  1, 2,
  2, 5,
  3, 8.5,
  4, (11 + 16)/2,
  5, (17 + 21)/2,
  6, (22 + 27)/2,
  7, (28 + 33)/2,
  8, (24 + 40)/2,
  9, (41 + 47)/2,
  10, (48 + 55)/2,
  11, (56 + 63)/2
)

# Adjustment only models. Label column must not contain "_" since this is used
# to split output filenames on, to see which models have been run already in do.ds
adj.models <- tribble(
  ~label, ~key, ~adj, ~form,

  "hn.cos", "hn", "cos", "~1",
  "hn.sp", "hn", "poly", "~1",
  "hn.hp", "hn", "herm", "~1",

  "hr.cos", "hr", "cos", "~1",
  "hr.sp", "hr", "poly", "~1",
  "hr.hp", "hr", "herm", "~1",

  "unif.cos", "unif", "cos", "~1",
  "unif.sp", "unif", "poly", "~1",
  "unif.hp", "unif", "herm", "~1"
)

# some ddfs with very small number of obs (e.g. 9) with aerial truncation limits
# caused the fitting optimizer to never finish (with 1 cos adjustment). So if
# number of observations is less than SMALL_DDF_DATA_LIMIT it will
# be considered failed, and a "failed_model_xxx" file will be created in the DF
# Summaries folder. If all ddfs for a species fail, then just use a dummy_ddf
# instead.
SMALL_DDF_DATA_LIMIT <- 20


# location of DS utils
source(file.path(RDir, "ds_utils_0.6.r"))
source(file.path(RDir, "functions.R"))
