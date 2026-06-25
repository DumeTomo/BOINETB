#############################################################################
# BOIN-ETB Code for simulation
#############################################################################
rm(list = ls())

# ---- Load packages ----
packages <- c("Iso", "UniIsoRegression", "mfp", "openxlsx",
              "tidyverse", "boinet", "grid", "gridExtra", "readxl")
lapply(packages, function(pkg) {
  message(paste("Loading package:", pkg))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "is not installed. Installing now..."))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

FOLDER <- "~/."

# Load function definitions
source(paste0(FOLDER, "/_FUNC/FUNC_CALC.SUM.R"))
source(paste0(FOLDER, "/_FUNC/FUNC_SIM.BOINETB.R"))

# ---- Simulation settings ----

K            <- 6                                        # Number of baskets
J            <- 4                                        # Number of dose levels
startdose    <- c(1, 5, 9, 13, 17, 21)                   # Starting dose index for each basket
nsim         <- 1000                                     # Number of simulations
ncohort.list <- 6                                        # List of maximum cohort numbers to evaluate
cohortsize   <- c(3, 3, 3, 3, 3, 3)                      # Cohort size per basket
mn           <- rep(99, K)                               # Minimum patients to stop each basket
mnodc        <- cohortsize * 2                           # Minimum patients required for ODC selection
nmethod      <- "BOIN-ETB"                               # Method name for output
tgt.t        <- c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3)          # Target toxicity probability per basket
tgt.e        <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)         # Target efficacy probability per basket
dev          <- 0.4                                      # Allowable deviation from the target
u            <- c(100, 40, 60, 0)                        # Utility scores (c1, c2, c3, c4)
method       <- paste0("BOIN-ETB_util", u[1], "-", u[2], "-", u[3], "-", u[4])
OUTFILE      <- paste0(FOLDER, "/Out/")

# Build output file name from simulation settings
fn_base  <- paste0("BOINETB_b", K, "d", J, "c", ncohort.list)
fn_tgtt  <- paste0("_tgtt", paste(tgt.t, collapse = "-"))
fn_tgte  <- paste0("_tgte", paste(tgt.e, collapse = "-"))
fn_dev   <- paste0("_dev", dev)
fn_util  <- paste0("_util", paste(u, collapse = "-"))
fn_csize <- paste0("_cohortsize", paste(cohortsize, collapse = "-"))
FILENAME <- paste0(fn_base, fn_tgtt, fn_tgte, fn_dev, fn_util, fn_csize)

# ---- Load scenario definitions ----
sc.ds   <- read_excel(paste0(FOLDER, "/SIM.BOINETB_scenario.xlsx"), sheet = "Example")
sc.name <- dplyr::distinct(sc.ds, sc)
sc.name <- sc.name$sc

# ---- Run simulations ----

# Loop over cohort numbers
for (inc in 1:length(ncohort.list)) {

  cat("\n")
  cat(paste0("### ncohort = ", ncohort.list[inc], "\n"))
  cat("\n")

  ncohort <- ncohort.list[inc]

  # Initialize result containers before scenario loop
  SIM.ALL1.isc <- data.frame()
  SIM.ALL2.isc <- data.frame()
  SC.TRUE.isc  <- data.frame()

  # Loop over scenarios
  for (isc in 1:length(sc.name)) {

    cat("\n")
    cat(paste0("### scenario: ", sc.name[isc], "\n"))
    cat("\n")

    # Extract true toxicity probabilities for this scenario
    pt.true <- dplyr::filter(sc.ds, sc %in% sc.name[isc]) %>%
      select(contains("pt"))
    pt.true <- as.vector(as.matrix(t(pt.true)))

    # Extract true efficacy probabilities for this scenario
    pe.true <- dplyr::filter(sc.ds, sc %in% sc.name[isc]) %>%
      select(contains("pe"))
    pe.true <- as.vector(as.matrix(t(pe.true)))

    # Run simulation for this scenario
    sim.sc <- SIM.BOINETB.ALG1(
      isc, pt.true, pe.true, startdose, ncohort, cohortsize,
      nsim, mn, mnodc, tgt.t, tgt.e, dev, K, J,
      paste0(FILENAME, "_SC", isc, "_nc", ncohort.list[inc])
    )

    # Build true probability and utility matrices for this scenario
    m.pt.true <- matrix(pt.true, ncol = J, byrow = TRUE)
    m.pe.true <- matrix(pe.true, ncol = J, byrow = TRUE)
    m.ut.true <- u[1] * (1 - m.pt.true) *      m.pe.true +
                 u[2] * (1 - m.pt.true) * (1 - m.pe.true) +
                 u[3] *      m.pt.true  *      m.pe.true +
                 u[4] *      m.pt.true  * (1 - m.pe.true)

    # Build scenario summary data frame
    m.sc  <- cbind(rep(sc.name[isc], K), 1:K, m.pt.true, m.pe.true, m.ut.true)
    df.sc <- data.frame(m.sc)
    colnames(df.sc) <- c("scenario", "basket",
                         paste0("ptox", 1:J),
                         paste0("peff", 1:J),
                         paste0("util", 1:J))

    # Accumulate results across scenarios
    SIM.ALL1.isc <- rbind(SIM.ALL1.isc, sim.sc$out.sum1)
    SIM.ALL2.isc <- rbind(SIM.ALL2.isc, sim.sc$out.sum2)
    SC.TRUE.isc  <- rbind(SC.TRUE.isc,  df.sc)
  }

  # Accumulate results across cohort numbers
  if (inc == 1) {
    SIM.ALL1 <- SIM.ALL1.isc
    SIM.ALL2 <- SIM.ALL2.isc
    SC.TRUE  <- SC.TRUE.isc
  } else {
    SIM.ALL1 <- rbind(SIM.ALL1, SIM.ALL1.isc)
    SIM.ALL2 <- rbind(SIM.ALL2, SIM.ALL2.isc)
  }

  # ---- Output to Excel file ----
  wb <- createWorkbook()
  addWorksheet(wb, "Scr")
  addWorksheet(wb, "Sum1")
  addWorksheet(wb, "Sum2")
  writeDataTable(wb, "Scr",  x = SC.TRUE)
  writeDataTable(wb, "Sum1", x = SIM.ALL1)
  writeDataTable(wb, "Sum2", x = SIM.ALL2)
  saveWorkbook(
    wb,
    paste0(OUTFILE, "summary_", FILENAME, "_ALL.xlsx"),
    overwrite = TRUE
  )
}

#############################################################################
# End of program
#############################################################################
