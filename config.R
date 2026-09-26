# config.R — settings for the Paper 2 inversion pipeline ---------------------------------------
# Every path comes from an environment variable with a default that matches the Methane_AMMBEC
# folder layout on Priyanka's Mac (this repo sits next to Paper 1's methane_AMMBEC repo):
#
#   Methane_AMMBEC/
#     methane_AMMBEC/        Paper 1 repo (read-only from here: its R/ helpers are sourced)
#     methane_inversion/     this repo
#     instruments/Aircraft/  Twin Otter ICARTT files
#     inventories/GRA2PES/   GRA2PES trees (sectors/<ver>/<SECTOR>/<MONTH>/<daytype>, v2total, ch4only, allspec)
#     inventories/EPA_GHGI/  EPA gridded GHGI
#     met/hrrr/sub/          HRRR ARL subgrids (shared with Paper 1 script 64)
#     inversion_outputs/     everything this repo writes
#
# Paper 1 is never modified. Nothing here writes outside INV_OUT, MET_DIR and BASE/stilt (all on the Mac).
# -----------------------------------------------------------------------------------------------
.here <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), ".."), mustWork = FALSE) else normalizePath(".", mustWork = FALSE)
}
INV_ROOT <- normalizePath(Sys.getenv("METHANE_INV_ROOT", if (file.exists("config.R")) "." else .here()), mustWork = FALSE)
# Two roots (25 Sep 2026): input data are READ from the Elements drive; code and everything the
# pipeline WRITES (outputs, HRRR subgrids it downloads, STILT binaries, footprints) stay on the Mac.
BASE      <- normalizePath(Sys.getenv("METHANE_BASE", file.path(INV_ROOT, "..")), mustWork = FALSE)   # Mac: /Users/priyanka/Methane_AMMBEC
DATA_BASE <- path.expand(Sys.getenv("METHANE_DATA_BASE", "/Volumes/Elements/Methane_AMMBEC"))          # read-only inputs
env <- function(k, d) path.expand(Sys.getenv(k, d))

# R packages and CmdStan live in the Mac folder too (the system R library is not writable)
R_LIB       <- env("METHANE_R_LIB", file.path(BASE, "Rlib"))
CMDSTAN_DIR <- env("METHANE_CMDSTAN_DIR", file.path(BASE, "cmdstan"))
dir.create(R_LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(R_LIB, .libPaths()))
use_local_cmdstan <- function() {        # point cmdstanr at the newest CmdStan in CMDSTAN_DIR
  if (nzchar(Sys.getenv("CMDSTAN"))) return(cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN")))
  d <- sort(list.dirs(CMDSTAN_DIR, recursive = FALSE), decreasing = TRUE)
  d <- d[grepl("cmdstan-", basename(d))]
  if (length(d)) cmdstanr::set_cmdstan_path(d[1])
}

PAPER1_DIR <- env("METHANE_PAPER1_DIR", file.path(DATA_BASE, "methane_AMMBEC"))  # sourced, never modified
DATA_DIR   <- env("METHANE_DATA_DIR",   file.path(DATA_BASE, "instruments"))      # holds Aircraft/
GRA2PES_DIR <- env("METHANE_GRA2PES_DIR", file.path(DATA_BASE, "inventories", "GRA2PES"))
GHGI_FILE  <- env("METHANE_GHGI", file.path(DATA_BASE, "inventories", "EPA_GHGI", "Express_Extension_Gridded_GHGI_Methane_v2_2020.nc"))
MET_DIR    <- env("METHANE_HRRR_DIR",   file.path(BASE, "met", "hrrr"))           # Mac; sub/ links the Elements subgrids
MET_READ_DIR <- env("METHANE_HRRR_READ_DIR", file.path(DATA_BASE, "met", "hrrr", "sub"))
INV_OUT    <- env("METHANE_INV_OUT",    file.path(BASE, "inversion_outputs"))
# Transport configuration. Everything that depends on the footprints (footprints, Jacobians, OSSE,
# inversion results) goes in <INV_OUT>/runs/<TRANSPORT>/, so particle-number and mixed-layer
# sensitivity runs sit side by side. The first runs (24-26 Sep) are runs/np300_h10.
NUMPAR    <- as.integer(Sys.getenv("METHANE_INV_NUMPAR", "1000"))    # STILT v2 default (Fasoli et al. 2018)
N_HOURS   <- -abs(as.integer(Sys.getenv("METHANE_INV_NHOURS", "10")))
ZISCALE   <- as.numeric(Sys.getenv("METHANE_ZISCALE", "1"))          # mixed-layer height scaling (1 = HRRR as is)
TRANSPORT <- Sys.getenv("METHANE_TRANSPORT", sprintf("np%d_h%d%s", NUMPAR, abs(N_HOURS),
                        if (ZISCALE != 1) sprintf("_zi%.2f", ZISCALE) else ""))
TRANSPORT_DIR <- file.path(INV_OUT, "runs", TRANSPORT)                     # footprints + transport checks
FOOT_DIR  <- env("METHANE_INV_FOOT", file.path(TRANSPORT_DIR, "footprints"))
# Prior variant (26 Sep 2026; default "complex" = main analysis, "wwtp" = secondary analysis). "wwtp": the Metro Water Recovery plant is its own component at its own
# location (see FACILITIES below); "complex": the 24 Sep definition (plant + neighbouring v2 cell).
# The two variants keep separate priors and jacobian/posterior files, so neither overwrites the other:
#   complex -> priors in INV_OUT,                  results in runs/<TRANSPORT>/
#   wwtp    -> priors in INV_OUT/priors/metro_wwtp, results in runs/<TRANSPORT>/metro_wwtp/
METRO_COMPONENT <- match.arg(Sys.getenv("METHANE_METRO_COMPONENT", "complex"), c("wwtp", "complex"))
PRIOR_TAG <- if (METRO_COMPONENT == "complex") "" else "metro_wwtp"
PRIOR_DIR <- if (nzchar(PRIOR_TAG)) file.path(INV_OUT, "priors", PRIOR_TAG) else INV_OUT
RUN_DIR   <- if (nzchar(PRIOR_TAG)) file.path(TRANSPORT_DIR, PRIOR_TAG) else TRANSPORT_DIR
stopifnot(!startsWith(normalizePath(INV_OUT, mustWork = FALSE), "/Volumes/Elements"))

## ---- STILT / HYSPLIT engine --------------------------------------------------------------------
# hycs_std with STILT emulation. Priyanka's macOS HYSPLIT 5.3 build cannot read the land-use files
# or the HRRR header (see Paper 1 notes), so the default is STILT's own build (uataq/stilt bin/),
# which on macOS needs the gfortran runtime in /usr/local/gfortran. Set METHANE_STILT_EXE to any
# working hycs_std; METHANE_STILT_ENGINE=external only writes run directories (for run_stilt_linux.sh).
STILT_EXE    <- env("METHANE_STILT_EXE", file.path(BASE, "stilt", "bin", if (Sys.info()[["sysname"]] == "Darwin") "macos_x64" else "linux_x64", "hycs_std"))
XTRCT_EXE    <- env("METHANE_XTRCT_EXE", file.path(dirname(STILT_EXE), "xtrct_grid"))   # cuts HRRR to DOMAIN (02_met.R)
STILT_BDY    <- env("METHANE_HYSPLIT_BDY", "~/hysplit/bdyfiles")   # LANDUSE.ASC, ROUGLEN.ASC, TERRAIN.ASC, ASCDATA.CFG
STILT_ENGINE <- Sys.getenv("METHANE_STILT_ENGINE", "local")         # "local" or "external"
N_CORES      <- as.integer(Sys.getenv("METHANE_CORES", max(1, parallel::detectCores() - 1)))

## ---- domain and receptors ----------------------------------------------------------------------
DOMAIN   <- c(ymn = 37.5, xmn = -108.5, ymx = 42.5, xmx = -101.5)  # met subgrid and flux domain
OBS_BOX  <- c(ymn = 39.3, xmn = -105.5, ymx = 41.0, xmx = -103.7)  # observations used: Denver + DJB
URBAN_BOX <- list(lat_s = 39.50, lat_n = 39.95, lon_w = -105.20, lon_e = -104.55)   # Paper 1 box
BASIN_LAT <- 40.05
SEG_S     <- as.integer(Sys.getenv("METHANE_INV_SEG_S", "60"))   # observation segment length (s)
MIN_SEG_N <- 30L                                                 # valid 1-Hz samples per segment
MAX_AGL_M <- 3000                                                # keep segments below this height

## ---- source components -------------------------------------------------------------------------
# CH4 components the inversion scales. og_urban / og_basin split GRA2PES OG by URBAN_MASK
# (cells inside the Paper 1 box, grown by URBAN_BUFFER_KM). Facility components carve named waste
# cells out of the waste sector so that DADS, Tower Road and the Metro/Suncor complex each get
# their own scale factor (Paper 1 found their inventory errors have opposite signs).
URBAN_BUFFER_KM <- 8
# Point sources around the Commerce City industrial corridor, used by 08_metro_wwtp.R and (when
# METHANE_METRO_COMPONENT = "wwtp") by 04 to give the wastewater plant its own component.
# METRO_WWTP is APPROXIMATE (from the 6450 York St address: the plant lies between York St and the
# South Platte, north of the Cherokee power station). Paper 1's biogenic_sources.csv puts it at
# 39.792 N 104.938 W, south-east of Suncor, which looks wrong; check on a map and override with
# METHANE_METRO_LAT / METHANE_METRO_LON.
POINT_SOURCES <- data.frame(
  name = c("metro_wwtp", "suncor", "cherokee"),
  label = c("Metro Water Recovery (Hite)", "Suncor refinery", "Cherokee power station"),
  lat = c(as.numeric(Sys.getenv("METHANE_METRO_LAT", "39.812")), 39.802, 39.8074),
  lon = c(as.numeric(Sys.getenv("METHANE_METRO_LON", "-104.955")), -104.945, -104.9639),
  stringsAsFactors = FALSE)
FACILITIES <- data.frame(
  name = c("dads", "tower_road", if (METRO_COMPONENT == "wwtp") "metro_wwtp" else "metro_complex"),
  # centres = the GRA2PES v2.0beta waste hotspot cells found in Paper 1 (v2_pointsource_hotspots.csv);
  # 3 km picks the single 4-km cell, 5 km also takes the neighbouring cell holding the Metro plant
  lat  = c(39.6592, 39.8364, 39.8239),
  lon  = c(-104.692, -104.759, -104.946),
  radius_km = c(3, 3, 5), stringsAsFactors = FALSE)
# "wwtp": the waste methane the inventory puts within 5 km of the Metro hotspot (the plant, wherever
# the inventory placed it) is gathered and moved, hour by hour and mass-conserving, into the single
# cell holding the plant's own location (POINT_SOURCES). Its scale factor is then the plant's alone,
# and its footprint sensitivity is the plant's, not the neighbouring cell's.
RELOCATE <- if (METRO_COMPONENT == "wwtp")
  list(metro_wwtp = c(lat = POINT_SOURCES$lat[POINT_SOURCES$name == "metro_wwtp"],
                      lon = POINT_SOURCES$lon[POINT_SOURCES$name == "metro_wwtp"])) else list()

# Ethane:methane source ratios (mol/mol): prior mean and sd per component. Waste and livestock
# carry none. Values follow Paper 1 section 3.2 / S5 (Denver delivered gas 0.110, 0.095-0.141;
# DJB 0.061-0.102). "other" is mostly combustion and small; weakly constrained.
R_PRIOR <- list(
  og_basin      = c(mean = 0.080, sd = 0.015),
  og_urban      = c(mean = 0.110, sd = 0.012),
  postmeter     = c(mean = 0.110, sd = 0.012),
  waste         = c(mean = 0,     sd = 0),
  dads          = c(mean = 0,     sd = 0),
  tower_road    = c(mean = 0,     sd = 0),
  metro_complex = c(mean = 0,     sd = 0),
  metro_wwtp    = c(mean = 0,     sd = 0),
  ag            = c(mean = 0,     sd = 0),
  other         = c(mean = 0.020, sd = 0.020))

## ---- units --------------------------------------------------------------------------------------
MW <- c(CH4 = 16.04, CO = 28.01)
# GRA2PES flux: mol km-2 h-1 -> umol m-2 s-1 (x 1e6 umol/mol / 1e6 m2/km2 / 3600 s/h)
GRA2PES_TO_UMOL <- 1 / 3600
# EPA GHGI flux: molec cm-2 s-1 -> umol m-2 s-1
GHGI_TO_UMOL <- 1e4 / 6.02214076e23 * 1e6

dir.create(INV_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(INV_OUT, "figures"), showWarnings = FALSE)
dir.create(file.path(RUN_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(TRANSPORT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(PRIOR_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Paper 1 helpers (ICARTT reader, leg detection, enhancements) -------------------------------
p1 <- function(f) file.path(PAPER1_DIR, "R", f)
if (!file.exists(p1("read_icartt.R"))) stop("Paper 1 repo not found at ", PAPER1_DIR, " (set METHANE_PAPER1_DIR)")
for (f in c("read_icartt.R", "paths.R", "massbalance.R", "enhancements.R", "ratios.R")) source(p1(f))
