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
library(ECSASconnect)
library(rmarkdown)
library(htmltools)
library(sf)
library(checkmate)
library(dftools)

SubProject <- "Atl_OPP" # Atlantic OPP IMRP/response study area
# SubProject <- "All_ECSAS" # Entire ECSAS dataset including both NL and NS


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
(the.data.file <- file.path(GenDataDir, "ECSAS_the_data.rds"))
(study.area.file <- file.path(GenDataDir, "study_area.rds"))
(proj.file <- file.path(GenDataDir, "projection.rds"))
(grid.file <- file.path(GenDataDir, "grid.rds"))

# Data options
ECSAS.reextract <- FALSE # Requery the database or use existing saved data.
save.grid.polygons <- FALSE # Takes hours to save!!
aerial.watch.len <- 30 # seconds


if (SubProject == "All_ECSAS") {
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
} else if (SubProject == "Atl_OPP"){
  spec.grps <- list(
    ATPU = "ATPU",
    BLKI = "BLKI",
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
} else
  stop("analysis_settings.R: unknown SubProject ", SubProject)
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

# location of DS utils
source(file.path(RDir, "ds_utils_0.6.r"))
source(file.path(RDir, "functions.R"))
