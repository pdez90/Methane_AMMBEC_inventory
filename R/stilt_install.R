# stilt_install.R — get STILT's hycs_std / xtrct_grid running on an Apple Silicon Mac, base R ------
# STILT ships Intel-only macOS binaries (uataq/stilt bin/macos_x64) linked against the gfortran
# runtime at /usr/local/gfortran/lib. On Apple Silicon they need:
#   1. Rosetta 2 (one-time, admin):  softwareupdate --install-rosetta --agree-to-license
#   2. the Intel gfortran runtime. Instead of installing it system-wide, the libraries are pulled
#      out of the gfortran-for-macOS Intel .dmg into <BASE>/stilt/bin/macos_x64/ and the binaries
#      are relinked to @executable_path (install_name_tool + ad-hoc codesign). Nothing outside
#      <BASE>/stilt is changed, and no admin rights are needed for this part.
# ----------------------------------------------------------------------------------------------
STILT_RAW  <- "https://github.com/uataq/stilt/raw/main/bin/macos_x64/"
GFORTRAN_DMG <- "https://github.com/fxcoudert/gfortran-for-macOS/releases/download/14.2-sequoia/gfortran-14.2-Intel-Sequoia.dmg"

rosetta_ok <- function() Sys.info()[["machine"]] != "arm64" ||
  suppressWarnings(system("arch -x86_64 /usr/bin/true > /dev/null 2>&1")) == 0

.deps <- function(f) {             # /usr/local/gfortran/... libraries a Mach-O file links against
  x <- suppressWarnings(system2("otool", c("-L", shQuote(f)), stdout = TRUE, stderr = FALSE))
  x <- trimws(sub("\\(compatibility.*", "", x[-1]))
  x[grepl("^/usr/local/gfortran/", x)]
}

.local_name <- function(b, bin_dir) {   # libgcc_s.1.dylib -> libgcc_s.1.1.dylib if that is what we have
  if (file.exists(file.path(bin_dir, b))) return(b)
  h <- list.files(bin_dir, pattern = "\\.dylib$"); h <- h[startsWith(h, sub("\\.dylib$", ".", b))]
  if (length(h)) h[1] else b
}

#' Download the STILT binaries into bin_dir, bring in the Intel gfortran runtime next to them and
#' relink everything to @executable_path / @loader_path. Returns TRUE when hycs_std starts.
install_stilt_mac <- function(bin_dir) {
  dir.create(bin_dir, FALSE, TRUE)
  exes <- c("hycs_std", "xtrct_grid")
  for (e in exes) {
    f <- file.path(bin_dir, e)
    if (!file.exists(f) || file.size(f) < 1e4) {
      cat("  downloading STILT", e, "\n")
      download.file(paste0(STILT_RAW, e), f, mode = "wb", quiet = TRUE)
    }
    Sys.chmod(f, "755")
  }
  need <- unique(unlist(lapply(file.path(bin_dir, exes), .deps)))
  have <- file.exists(file.path(bin_dir, vapply(basename(need), .local_name, "", bin_dir = bin_dir)))
  missing_libs <- need[!have & !file.exists(need)]
  if (length(missing_libs)) {
    cat("  fetching the Intel gfortran runtime (one time, ~", "300 MB download, only 3-4 small libraries kept)\n", sep = "")
    tmp <- tempfile("gf_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    dmg <- file.path(dirname(bin_dir), "gfortran-Intel.dmg")   # kept until the runtime is in place
    mnt <- file.path(tmp, "mnt"); exp <- file.path(tmp, "pkg")
    options(timeout = max(3600, getOption("timeout")))
    if (!file.exists(dmg) || file.size(dmg) < 1e8) download.file(GFORTRAN_DMG, dmg, mode = "wb", quiet = FALSE)
    stopifnot(system2("hdiutil", c("attach", "-nobrowse", "-readonly", "-mountpoint", shQuote(mnt), shQuote(dmg)),
                      stdout = FALSE, stderr = FALSE) == 0)
    pkg <- list.files(mnt, pattern = "\\.pkg$", full.names = TRUE, recursive = TRUE)[1]
    system2("pkgutil", c("--expand-full", shQuote(pkg), shQuote(exp)), stdout = FALSE)
    system2("hdiutil", c("detach", shQuote(mnt), "-quiet"))
    libs <- list.files(exp, pattern = "\\.dylib$", recursive = TRUE, full.names = TRUE)
    # closure: the libraries the binaries need, plus whatever those libraries need
    todo <- basename(need); done <- character()
    while (length(todo)) {
      b <- todo[1]; todo <- todo[-1]; if (b %in% done) next
      src <- libs[basename(libs) == b][1]
      # newer GCC ships libgcc_s.1.1.dylib (a compatible superset) instead of libgcc_s.1.dylib
      if (is.na(src)) src <- libs[startsWith(basename(libs), sub("\\.dylib$", ".", b))][1]
      if (is.na(src)) stop("could not find ", b, " in the gfortran package; it has: ",
                           paste(sort(unique(basename(libs))), collapse = " "))
      dst <- file.path(bin_dir, basename(src))   # keep the real name; links are pointed at it below
      file.copy(normalizePath(src), dst, overwrite = TRUE); Sys.chmod(dst, "755")
      done <- c(done, b, basename(src)); todo <- c(todo, setdiff(basename(.deps(dst)), done))
    }
    unlink(dmg)
  }
  # relink: executables -> @executable_path/<lib>, libraries -> @loader_path/<lib>
  for (f in list.files(bin_dir, full.names = TRUE)) {
    if (!(basename(f) %in% exes || grepl("\\.dylib$", f))) next
    d <- .deps(f); if (!length(d) && !grepl("\\.dylib$", f)) next
    pre <- if (grepl("\\.dylib$", f)) "@loader_path/" else "@executable_path/"
    tgt <- vapply(basename(d), .local_name, "", bin_dir = bin_dir)
    args <- if (length(d)) as.vector(rbind("-change", d, paste0(pre, tgt))) else character()
    if (grepl("\\.dylib$", f)) args <- c("-id", paste0("@loader_path/", basename(f)), args)
    if (!length(args)) next
    system2("install_name_tool", c(args, shQuote(f)), stderr = FALSE)
    system2("codesign", c("--force", "-s", "-", shQuote(f)), stdout = FALSE, stderr = FALSE)
  }
  system2("xattr", c("-dr", "com.apple.quarantine", shQuote(bin_dir)), stderr = FALSE)
  invisible(TRUE)
}

#' TRUE if the binary loads and executes: it then stops by itself for lack of input (hycs_std prints its
#' HYSPLIT banner; xtrct_grid hits end-of-file on stdin, which is a Fortran runtime message, not a
#' load failure). Load failures look like dyld / "Bad CPU type" / "Library not loaded".
stilt_starts <- function(exe) {
  wd <- tempfile("stilt_"); dir.create(wd); on.exit(unlink(wd, recursive = TRUE))
  suppressWarnings(system(paste("cd", shQuote(wd), "&&", shQuote(exe), "> log 2>&1 < /dev/null"), timeout = 20))
  lg <- paste(readLines(file.path(wd, "log"), warn = FALSE), collapse = " ")
  loaded <- !grepl("dyld|Bad CPU|Library not loaded|cannot execute|Killed|not permitted", lg, ignore.case = TRUE) &&
    grepl("HYSPLIT|CONTROL|STILT|Fortran runtime|xtrct", lg, ignore.case = TRUE)
  structure(loaded, log = substr(lg, 1, 300))
}
