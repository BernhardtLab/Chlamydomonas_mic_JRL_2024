# =============================================================================
# 01_data_upload.R
# Chlamydomonas x microbe stress-resilience experiment
# Jason R Laurich
#
# Assemble all raw plate-reader exports (RFU + OD600 time series) into a single
# tidy data frame, joined to the experimental design, for all four blocks.
#
# ---- What this script produces --------------------------------------------
#   Data-processed/01_raw_timeseries_data.csv   (133,320 observations)
#
# ---- Refactor note (2026) --------------------------------------------------
# This is a refactor of the original 01_data_upload.R. The DATA-HANDLING LOGIC
# is unchanged and reproduces the original output value-for-value; only the
# structure changed:
#   * the 36 copy-pasted temperature blocks -> one function + a parameter table
#   * rbind-in-loop -> collect in lists, bind once (bind_rows)
#   * plate/read/excel-range constants pulled into one place
#   * unused packages dropped
#   * NEW: a documented `carbon.added` column (the ONLY change to the output;
#     it is appended as the last column, so every original column is identical).
#
# Block 2 did not receive organic carbon in its media. It is retained here and
# flagged (carbon.added == FALSE) so the include/exclude decision can be made
# downstream; it will not feature in the main-text analyses.
# =============================================================================

library(tidyverse)   # dplyr / tibble / readr / %>%
library(readxl)      # read_excel()

# ---- Configuration ---------------------------------------------------------
# Paths are relative to the RStudio project root (open the .Rproj, or setwd()).
raw.dir  <- "Data-raw"          # plate .xlsx files AND the *-design.csv files
proc.dir <- "Data-processed"    # cleaned outputs land here (new home)
out.file <- "01_raw_timeseries_data.csv"

if (!dir.exists(proc.dir)) dir.create(proc.dir, recursive = TRUE)

# Experiment structure ------------------------------------------------------
blocks <- 1:4
temps  <- c(8, 14, 20, 25, 30, 33, 35, 39, 43)   # SOURCE ORDER — do not reorder
                                                 # (sets output row order)

design.files <- c("1" = "02-block-1-design.csv",
                  "2" = "04-block-2-design.csv",
                  "3" = "06-block-3-design.csv",
                  "4" = "08-block-4-design.csv")

# Plates per (temperature): 30 C carried the N + salt gradients (32 plates);
# every other temperature is a 2-plate thermal-gradient run.
plate.bound <- function(temp) if (temp == 30) 32L else 2L

# Highest read index attempted for a (block, temperature). Missing files are
# skipped, so these are upper bounds. They reproduce the original loop limits
# exactly — including block 3's warm runs, which stopped at read 12.
read.bound <- function(block, temp) {
  if (temp %in% c(8, 14, 20, 25)) return(20L)   # long cold runs
  if (temp == 30)                 return(15L)   # 30 C gradient plates
  if (block == 3)                 return(12L)   # block 3 warm runs (33-43 C)
  15L                                           # all other warm runs
}

# Fixed plate-reader export layout.
# The exports are full 96-well grids: well-column numbers in row 48, then
# plate-rows A-H each as an RFU row followed by an OD600 row (rows 49-64),
# numeric data in spreadsheet columns C-N. This experiment uses only the
# INTERIOR 6 x 10 block (plate-rows B-G, well-columns 2-11); the edge wells
# (row A / H, columns 1 / 12) are excluded by design.
#   D51:M62  -> that interior block: 6 well-rows x (RFU, OD600) x 10 columns
#   B7       -> run date;  B8 -> run time
# NOTE: the original script read D50:M62 with headers ON, which used row 50
# (row A's OD600 values) as a throwaway header. Reading D51:M62 with headers
# OFF is identical value-for-value but avoids readxl's "New names" renaming of
# those numeric header cells. Data begins at row 51 (verified against the
# reference output: well 1 = plate-row B, well-column 2 = cell D51).
data.range <- "D51:M62"   # 12 rows x 10 cols of pure numbers (no header row)
head.range <- "B7:B8"

# Accumulator for any file whose data block is not the expected 12 x 10.
problem.files <- list()

# ---- Load & annotate the design files -------------------------------------
design.list <- lapply(blocks, function(b) {
  d <- read.csv(file.path(raw.dir, design.files[[as.character(b)]]))
  d$Block     <- b
  d$unique.id <- paste0("b", b, ".t", d$Temperature.C,
                        ".p", d$Plate.at.T, ".w", d$Well.at.T)
  d
})
names(design.list) <- as.character(blocks)

# Flag any wells with hand notes (potential pipetting errors). These are NOT
# corrected here — that happens in 02_mu_estimation.R — but printing them keeps
# the original script's disclosure behaviour.
for (b in blocks) {
  notes <- design.list[[as.character(b)]] %>%
    filter(!is.na(Notes), trimws(Notes) != "")
  if (nrow(notes)) message("Block ", b, ": ", nrow(notes), " wells carry notes.")
}

# ---- Parse one plate-reader file ------------------------------------------
# Returns a 60-row data frame (one row per well) or NULL if the file / header
# cannot be read.
read_one_plate <- function(fp, plate, read) {

  # Read sheet 1 explicitly: a few exports carry a second sheet, and the
  # default (sheet 1) is the one used everywhere else.
  raw <- tryCatch(
    as.data.frame(read_excel(fp, sheet = 1, range = data.range, col_names = FALSE)),
    error = function(e) e)
  if (inherits(raw, "error")) {
    problem.files[[length(problem.files) + 1L]] <<- data.frame(
      file = basename(fp), reason = "read error", nrow = NA, ncol = NA)
    message("READ ERROR : ", basename(fp), " (", conditionMessage(raw), ") - skipping")
    return(NULL)
  }

  # Guard the layout: anything that is not the expected 12 x 10 numeric block
  # (e.g. a differently structured export) is logged and skipped rather than
  # crashing the whole run.
  if (nrow(raw) != 12L || ncol(raw) != 10L) {
    problem.files[[length(problem.files) + 1L]] <<- data.frame(
      file = basename(fp), reason = "unexpected dimensions",
      nrow = nrow(raw), ncol = ncol(raw))
    message("BAD LAYOUT : ", basename(fp), " -> ", nrow(raw), " x ", ncol(raw),
            " (expected 12 x 10) - skipping")
    return(NULL)
  }

  vars <- tryCatch(as.data.frame(read_excel(fp, sheet = 1, range = head.range,
                                            col_names = FALSE)),
                   error = function(e) NULL)
  if (is.null(vars)) {
    problem.files[[length(problem.files) + 1L]] <<- data.frame(
      file = basename(fp), reason = "header unreadable", nrow = 12L, ncol = 10L)
    message("NO HEADER  : ", basename(fp), " - skipping")
    return(NULL)
  }

  # 12 data rows: odd rows (1,3,..,11) = RFU per well-row; even rows = OD600.
  # 10 columns = well-columns. Flatten row-major (well-row outer, column inner)
  # to reproduce the original well ordering exactly.
  vals    <- raw[, 1:10]
  vals[]  <- lapply(vals, function(z) as.numeric(as.character(z)))
  m       <- as.matrix(vals)
  rfu     <- as.vector(t(m[seq(1, 11, by = 2), , drop = FALSE]))
  od600   <- as.vector(t(m[seq(2, 12, by = 2), , drop = FALSE]))

  # Date (B7) and datetime (B8) -> a single POSIXct plus the julian date the
  # original recorded (origin 2030-01-01; used only as a stored offset).
  d    <- as.Date(vars[1, 1], format = "%Y-%m-%d")
  tstr <- format(as.POSIXct(vars[2, 1], format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 "%H:%M:%S")

  data.frame(
    Row      = rep(1:6, each = 10L),
    Column   = rep(1:10, times = 6L),
    Plate    = plate,
    Read     = read,
    RFU      = rfu,
    OD600    = od600,
    date     = as.numeric(julian(d, origin = as.Date("2030-01-01"))),
    time     = tstr,
    datetime = as.POSIXct(paste(d, tstr), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Row       = as.integer(Row),
           Column    = as.integer(Column),
           Well.at.T = (Plate - 1L) * 60L + (Row - 1L) * 10L + Column)
}

# ---- Assemble one (block, temperature) run --------------------------------
read_block_temp <- function(block, temp, design.block) {

  design.bt <- design.block %>% filter(Temperature.C == temp)
  pmax <- plate.bound(temp)
  imax <- read.bound(block, temp)

  plate.frames <- list()
  for (p in 1:pmax) {
    design.tp <- design.bt %>% filter(Plate.at.T == p)

    read.frames <- list()
    for (i in 0:imax) {
      fp  <- file.path(raw.dir, sprintf("JRL_block%d_%d_%d_%d.xlsx", block, temp, p, i))
      if (!file.exists(fp)) next
      plt <- read_one_plate(fp, plate = p, read = i)
      if (is.null(plt)) next
      read.frames[[length(read.frames) + 1L]] <-
        left_join(plt, design.tp, by = "Well.at.T")
    }
    if (length(read.frames))
      plate.frames[[length(plate.frames) + 1L]] <- bind_rows(read.frames)
  }
  if (!length(plate.frames)) return(NULL)

  frame <- bind_rows(plate.frames)
  # Elapsed time since the first read OF THIS run (per block x temperature).
  frame$days <- as.numeric(difftime(frame$datetime,
                                    min(frame$datetime, na.rm = TRUE),
                                    units = "days"))
  frame
}

# ---- Run over every block x temperature (source order) --------------------
frames <- list()
for (b in blocks) {
  design.b <- design.list[[as.character(b)]]
  for (t in temps) {
    key <- paste0("b", b, ".t", t)
    message("Reading ", key, " ...")
    frames[[key]] <- read_block_temp(b, t, design.b)
  }
}

df.sum <- bind_rows(frames)

# ---- Document the carbon-free block ---------------------------------------
# Block 2 received no organic carbon; every other block did.
df.sum$carbon.added <- df.sum$Block != 2

# ---- Sanity checks ---------------------------------------------------------
cat("\n--- 01_data_upload summary ---\n")
cat("Total observations :", nrow(df.sum), "(original reference: 133320)\n")
cat("Plate-reads parsed  :", nrow(df.sum) / 60, "\n")
cat("NA design joins     :", sum(is.na(df.sum$Microbe)), "(expected 0)\n")
cat("Blocks present      :", paste(sort(unique(df.sum$Block)), collapse = ", "), "\n")
cat("carbon.added FALSE  :", sum(!df.sum$carbon.added), "rows (block 2)\n")

# Report any files that could not be parsed as a clean 12 x 10 block.
if (length(problem.files)) {
  probs <- bind_rows(problem.files)
  cat("\nNON-CONFORMING FILES (", nrow(probs), ") - skipped:\n", sep = "")
  print(probs)
  write.csv(probs, file.path(proc.dir, "01_problem_files.csv"), row.names = FALSE)
  cat("Logged to:", file.path(proc.dir, "01_problem_files.csv"), "\n")
} else {
  cat("\nNo non-conforming files: every plate parsed as 12 x 10.\n")
}

# ---- Write -----------------------------------------------------------------
write.csv(df.sum, file.path(proc.dir, out.file))
cat("\nWritten to:", file.path(proc.dir, out.file), "\n")
