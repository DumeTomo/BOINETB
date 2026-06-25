####################################################################
# BOIN-ETB
# Code to define function CALC.SUM (simulation)
####################################################################
CALC.SUM <- function(
  OUTFILE,  # Path to the output folder
  FILENAME, # Output file name (without extension)
  mnodc,    # Minimum number of patients required for ODC selection per basket
  SC        # Scenario number to reference
) {

  # ---- Load simulation data ----
  data  <- read_csv(file.path(OUTFILE, paste0(FILENAME, "_data.csv")))
  sdata <- read_csv(file.path(OUTFILE, paste0(FILENAME, "_stop.csv")))
  DATA  <- as.data.frame(sapply(data[-2:-1, ], as.numeric))
  DOSE  <- data[2, ]
  BKT   <- data[1, ]

  # Number of dose levels (per basket, total)
  nd1 <- J; nd2 <- K; nd <- nd1 * nd2
  num <- 1:(J * K)

  # Number of studies without ODC
  non.ODC <- colSums(sdata)

  # ---- Define dose categories ----

  # Target dose level (max utility, not too toxic)
  ud  <- dplyr::filter(sc.ds, sc %in% sc.name[SC]) %>% select(udose)
  ud  <- as.vector(as.matrix(t(ud)))
  ud1 <- rep(ud, each = J)

  # Target dose level (closest to target efficacy, not too toxic)
  ed  <- dplyr::filter(sc.ds, sc %in% sc.name[SC]) %>% select(edose)
  ed  <- as.vector(as.matrix(t(ed)))
  ed1 <- rep(ed, each = J)

  # Too toxic dose level
  od  <- dplyr::filter(sc.ds, sc %in% sc.name[SC]) %>% select(odose)
  od  <- as.vector(as.matrix(t(od)))
  od1 <- rep(od, each = J)

  # Dose category based on utility (3 categories: under / target / over)
  c3udose <- rep(1, J * K)
  c3udose[ud] <- 2
  c3udose[which(num > ud1)] <- 3
  c3udose <- factor(c3udose, levels = 1:3)

  # Dose category based on utility (4 categories: under / target / over / tox)
  c4udose <- rep(1, J * K)
  c4udose[ud] <- 2
  c4udose[which(num > ud1 & num < od1)] <- 3
  c4udose[which(num > ud1 & num >= od1)] <- 4
  c4udose <- factor(c4udose, levels = 1:4)

  # Dose category based on efficacy (3 categories: under / target / over)
  c3edose <- rep(1, J * K)
  c3edose[ed] <- 2
  c3edose[which(num > ed1)] <- 3
  c3edose <- factor(c3edose, levels = 1:3)

  # Dose category based on efficacy (4 categories: under / target / over / tox)
  c4edose <- rep(1, J * K)
  c4edose[ed] <- 2
  c4edose[which(num > ed1 & num < od1)] <- 3
  c4edose[which(num > ed1 & num >= od1)] <- 4
  c4edose <- factor(c4edose, levels = 1:4)

  # ---- Main simulation loop ----

  for (i in 1:nsim) {

    dose  <- as.numeric(DOSE[1:nd])
    bkt   <- as.numeric(BKT[1:nd])
    npts  <- as.numeric(DATA[i, 1:nd])
    ntox  <- as.numeric(DATA[i, (nd + 1):(nd * 2)])
    neff  <- as.numeric(DATA[i, (nd * 2 + 1):(nd * 3)])
    ppts  <- as.numeric(DATA[i, (nd * 3 + 1):(nd * 4)])
    ptox  <- as.numeric(DATA[i, (nd * 4 + 1):(nd * 5)])
    peff  <- as.numeric(DATA[i, (nd * 5 + 1):(nd * 6)])
    nct1  <- as.numeric(DATA[i, (nd * 6 + 1):(nd * 7)])
    nct2  <- as.numeric(DATA[i, (nd * 7 + 1):(nd * 8)])
    nct3  <- as.numeric(DATA[i, (nd * 8 + 1):(nd * 9)])
    nct4  <- as.numeric(DATA[i, (nd * 9 + 1):(nd * 10)])
    pct1  <- as.numeric(DATA[i, (nd * 10 + 1):(nd * 11)])
    pct2  <- as.numeric(DATA[i, (nd * 11 + 1):(nd * 12)])
    pct3  <- as.numeric(DATA[i, (nd * 12 + 1):(nd * 13)])
    pct4  <- as.numeric(DATA[i, (nd * 13 + 1):(nd * 14)])

    # ---- Find ODC by basket ----

    for (ii in 1:nd2) {

      doseb <- dose[which(bkt == ii)]
      neffb <- neff[which(bkt == ii)]
      ntoxb <- ntox[which(bkt == ii)]
      nptsb <- npts[which(bkt == ii)]

      # Estimate toxicity rate with small prior (avoid zero division)
      ptoxb <- (ntoxb + 0.05) / (nptsb + 0.1)

      # Isotonic regression on toxicity estimates
      phatb <- pava(ptoxb)

      # Observed efficacy rate
      peffb <- neffb / nptsb

      # Find MTD: dose with observed toxicity closest to target
      # among doses with sufficient patients
      idx <- which(nptsb >= mnodc[ii])
      if (length(idx) > 0) {
        min_diff   <- min(abs(phatb[idx] - tgt.t[ii]))
        candidates <- idx[which(abs(phatb[idx] - tgt.t[ii]) == min_diff)]
        mtdb       <- max(candidates)  # Select the highest qualifying dose
      } else {
        mtdb <- NA  # No dose meets the minimum patient requirement
      }

      # Indicator of tolerable doses (at or below MTD)
      ind.tdb <- doseb <= mtdb

      # ---- Calculate utility ----

      # Utility based on observed outcome counts (considers MTD)
      nctb1  <- nct1[which(bkt == ii)]
      nctb2  <- nct2[which(bkt == ii)]
      nctb3  <- nct3[which(bkt == ii)]
      nctb4  <- nct4[which(bkt == ii)]
      x      <- (u[1] * nctb1 + u[2] * nctb2 + u[3] * nctb3 + u[4] * nctb4) / 100
      x[is.na(x)] <- 0
      utilb  <- x * ind.tdb   # Utility with MTD constraint
      utilbn <- x             # Utility without MTD constraint

      # Utility based on estimated probabilities
      xx <- u[1] * peffb * (1 - phatb) +
            u[2] * (1 - peffb) * (1 - phatb) +
            u[3] * peffb * phatb +
            u[4] * (1 - peffb) * phatb
      xx[is.na(xx)] <- 0
      utilc  <- xx * ind.tdb  # Utility with MTD constraint
      utilcn <- xx            # Utility without MTD constraint

      # ---- Find ODC among doses with sufficient patients ----
      # ODC definitions:
      # odc : Observed probability, with MTD constraint
      # odc1: Observed probability, without MTD constraint
      # odc2: Estimated probability, with MTD constraint    <- primary metric
      # odc3: Estimated probability, without MTD constraint

      utilb[nptsb < mnodc[ii]]  <- 0
      if (sum(utilb) == 0)  { odcb  <- 99
      } else {
        odcb <- max(which(utilb == max(utilb, na.rm = TRUE)))
      }

      utilbn[nptsb < mnodc[ii]] <- 0
      if (sum(utilbn) == 0) { odcb1 <- 99
      } else {
        odcb1 <- max(which(utilbn == max(utilbn, na.rm = TRUE)))
      }

      utilc[nptsb < mnodc[ii]]  <- 0
      if (sum(utilc) == 0)  { odcb2 <- 99
      } else {
        odcb2 <- max(which(utilc == max(utilc, na.rm = TRUE)))
      }

      utilcn[nptsb < mnodc[ii]] <- 0
      if (sum(utilcn) == 0) { odcb3 <- 99
      } else {
        odcb3 <- max(which(utilcn == max(utilcn, na.rm = TRUE)))
      }

      # Accumulate results across baskets
      if (ii == 1) {
        phat.tox <- phatb
        mtd  <- mtdb
        odc  <- odcb
        odc1 <- odcb1
        odc2 <- odcb2
        odc3 <- odcb3
        util <- utilb
      } else {
        phat.tox <- c(phat.tox, phatb)
        mtd  <- c(mtd,  odcb  + (ii - 1) * J)
        odc  <- c(odc,  odcb  + (ii - 1) * J)
        odc1 <- c(odc1, odcb1 + (ii - 1) * J)
        odc2 <- c(odc2, odcb2 + (ii - 1) * J)
        odc3 <- c(odc3, odcb3 + (ii - 1) * J)
        util <- cbind(util, utilb)
      }
    }

    # ---- Accumulate simulation results across iterations ----
    sim_result_vars <- list(
      NPTS = npts, NTOX = ntox, NEFF = neff,
      PPTS = ppts, PTOX = ptox, PEFF = peff,
      NCT1 = nct1, NCT2 = nct2, NCT3 = nct3, NCT4 = nct4,
      PCT1 = pct1, PCT2 = pct2, PCT3 = pct3, PCT4 = pct4,
      MTD  = mtd,
      ODC  = odc, ODC1 = odc1, ODC2 = odc2, ODC3 = odc3
    )
    if (i == 1) {
      SIM_RESULTS <- lapply(sim_result_vars, list)
    } else {
      SIM_RESULTS <- mapply(
        function(acc, val) append(acc, list(val)),
        SIM_RESULTS, sim_result_vars,
        SIMPLIFY = FALSE
      )
    }

    # Accumulate total sample size
    if (i == 1) { SUM.N <- sum(npts)
    } else      { SUM.N <- c(SUM.N, sum(npts)) }
  }

  # Unpack accumulated results for readability
  NPTS <- SIM_RESULTS$NPTS; NTOX <- SIM_RESULTS$NTOX; NEFF <- SIM_RESULTS$NEFF
  PPTS <- SIM_RESULTS$PPTS; PTOX <- SIM_RESULTS$PTOX; PEFF <- SIM_RESULTS$PEFF
  NCT1 <- SIM_RESULTS$NCT1; NCT2 <- SIM_RESULTS$NCT2
  NCT3 <- SIM_RESULTS$NCT3; NCT4 <- SIM_RESULTS$NCT4
  PCT1 <- SIM_RESULTS$PCT1; PCT2 <- SIM_RESULTS$PCT2
  PCT3 <- SIM_RESULTS$PCT3; PCT4 <- SIM_RESULTS$PCT4
  MTD  <- SIM_RESULTS$MTD
  ODC  <- SIM_RESULTS$ODC;  ODC1 <- SIM_RESULTS$ODC1
  ODC2 <- SIM_RESULTS$ODC2; ODC3 <- SIM_RESULTS$ODC3

  # ---- Summary statistics ----

  # Grid search to find the thresholds
  optim_list <- vector("list", length(tgt.t))
  for (i in seq_along(tgt.t)) {
    optim_list[[i]] <- gridoptim(
      pi     = rep(1/6, 6),
      phi    = tgt.t[i],
      phi1   = dev * tgt.t[i],
      phi2   = (1 + dev) * tgt.t[i],
      delta  = tgt.e[i],
      delta1 = (1 - dev) * tgt.e[i],
      n      = 100
    )
  }
  optim_df <- bind_rows(
    lapply(seq_along(optim_list), function(i) {
      df        <- as.data.frame(optim_list[[i]])
      df$basket <- paste0("basket", i)
      df$tgt.t  <- tgt.t[i]
      df$tgt.e  <- tgt.e[i]
      df
    })
  )

  # Mean counts and rates by dose level
  m.NPTS  <- round(Reduce('+', NPTS) / nsim, digits = 2)        # Mean patients per dose
  m.NTOX  <- round(Reduce('+', NTOX) / nsim, digits = 2)        # Mean toxicity count per dose
  m.NEFF  <- round(Reduce('+', NEFF) / nsim, digits = 2)        # Mean efficacy count per dose
  m.PPTS  <- round(Reduce('+', PPTS) / nsim * 100, digits = 1)  # % patients treated per dose
  m.pTOX  <- round(Reduce('+', NTOX) / Reduce('+', NPTS) * 100, digits = 1)  # Observed P(tox)
  m.pEFF  <- round(Reduce('+', NEFF) / Reduce('+', NPTS) * 100, digits = 1)  # Observed P(eff)

  # % selecting each dose as MTD or ODC
  pMTD  <- round(table(factor(Reduce(function(a, b) { c(a, b) }, MTD),  levels = 1:nd)) / nsim * 100, digits = 1)
  pODC  <- round(table(factor(Reduce(function(a, b) { c(a, b) }, ODC),  levels = 1:nd)) / nsim * 100, digits = 1)
  pODC1 <- round(table(factor(Reduce(function(a, b) { c(a, b) }, ODC1), levels = 1:nd)) / nsim * 100, digits = 1)
  pODC2 <- round(table(factor(Reduce(function(a, b) { c(a, b) }, ODC2), levels = 1:nd)) / nsim * 100, digits = 1)
  pODC3 <- round(table(factor(Reduce(function(a, b) { c(a, b) }, ODC3), levels = 1:nd)) / nsim * 100, digits = 1)

  # Summary by basket
  correct.pODC <- pODC2[ud1]   # % selecting correct ODC
  target.PPTS  <- m.PPTS[ud1]  # % patients assigned to target dose

  dat.bkt <- data.frame(
    bkt      = bkt,
    c3udose  = c3udose,
    c4udose  = c4udose,
    c3edose  = c3edose,
    c4edose  = c4edose,
    pMTD     = pMTD,
    pODC     = pODC,
    pODC1    = pODC1,
    pODC2    = pODC2,
    pODC3    = pODC3,
    m.NPTS   = m.NPTS,
    m.NTOX   = m.NTOX,
    m.NEFF   = m.NEFF,
    m.PPTS   = m.PPTS
  )

  # Summary by basket and dose category
  sum.bkt.c3udose <- dat.bkt %>%
    group_by(bkt, c3udose, .drop = FALSE) %>%
    summarize(c.pODC = sum(pODC2.Freq), c.NPTS = sum(m.NPTS), .groups = "drop") %>%
    group_by(bkt) %>%
    mutate(c.PPTS = round(c.NPTS / sum(c.NPTS) * 100, digits = 1)) %>%
    ungroup()

  sum.bkt.c4udose <- dat.bkt %>%
    group_by(bkt, c4udose, .drop = FALSE) %>%
    summarize(c.pODC = sum(pODC2.Freq), c.NPTS = sum(m.NPTS), .groups = "drop") %>%
    group_by(bkt) %>%
    mutate(c.PPTS = round(c.NPTS / sum(c.NPTS) * 100, digits = 1)) %>%
    ungroup()

  sum.bkt.c3edose <- dat.bkt %>%
    group_by(bkt, c3edose, .drop = FALSE) %>%
    summarize(c.pODC = sum(pODC2.Freq), c.NPTS = sum(m.NPTS), .groups = "drop") %>%
    group_by(bkt) %>%
    mutate(c.PPTS = round(c.NPTS / sum(c.NPTS) * 100, digits = 1)) %>%
    ungroup()

  sum.bkt.c4edose <- dat.bkt %>%
    group_by(bkt, c4edose, .drop = FALSE) %>%
    summarize(c.pODC = sum(pODC2.Freq), c.NPTS = sum(m.NPTS), .groups = "drop") %>%
    group_by(bkt) %>%
    mutate(c.PPTS = round(c.NPTS / sum(c.NPTS) * 100, digits = 1)) %>%
    ungroup()

  # True toxicity, efficacy, and utility matrices
  aaa <- matrix(pt.true, ncol = J, byrow = TRUE,
                dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))
  bbb <- matrix(pe.true, ncol = J, byrow = TRUE,
                dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))
  true.util <- u[1] * (1 - aaa) * bbb +
               u[2] * (1 - aaa) * (1 - bbb) +
               u[3] * aaa * bbb +
               u[4] * aaa * (1 - bbb)
  dimnames(true.util) <- list(paste0("basket", 1:K), paste0("dose", 1:J))

  # ---- Text summary output (start) ----
  sink(paste0(OUTFILE, FILENAME, "_sum.txt"), append = FALSE)

  cat("\n")
  cat(paste0("# Simulation results of ", nmethod, "\n"))
  cat("# Simulation Date:", format(Sys.time(), "%Y%m%d_%H%M"), "\n")
  cat("\n\n")
  cat("# Settings \n")
  cat("\n")
  cat(paste0("simulation times        = ", nsim, "\n"))
  cat(paste0("number of max cohort    = ", ncohort, "\n"))
  cat("thresholds of escalation/de-escalation \n")
  print(optim_df); cat("\n")
  cat(paste0("scenario = ", SC, "\n"))
  cat("cohortsize \n");   print(cohortsize)
  cat("target tox \n");   print(tgt.t)
  cat("target eff \n");   print(tgt.e)
  cat("true probability (tox) \n"); print(aaa); cat("\n")
  cat("true probability (eff) \n"); print(bbb); cat("\n")
  cat("parameter for utility calculation \n"); print(u)
  cat("utility by basket and dose under true tox/eff probability \n")
  print(true.util)
  cat("\n\n")
  cat("# Results \n")
  cat("\n")

  cat("# mean number of patients by dose level \n")
  print(matrix(m.NPTS, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage of number of patients treated by dose level \n")
  print(matrix(m.PPTS, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# mean number of tox by dose level \n")
  print(matrix(m.NTOX, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# mean number of eff by dose level \n")
  print(matrix(m.NEFF, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# observed P(tox) by dose level \n")
  print(matrix(m.pTOX, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# observed P(eff) by dose level \n")
  print(matrix(m.pEFF, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to select as MTD by dose level \n")
  print(matrix(pMTD, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to select as ODC by dose level (based on observed prob) \n")
  print(matrix(pODC, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to select as ODC by dose level (based on observed prob, not consider MTD) \n")
  print(matrix(pODC1, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to select as ODC by dose level (based on estimated prob) \n")
  print(matrix(pODC2, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to select as ODC by dose level (based on estimated prob, not consider MTD) \n")
  print(matrix(pODC3, ncol = J, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), paste0("dose", 1:J)))); cat("\n")

  cat("# percentage to stop early \n")
  print(round(non.ODC / nsim * 100, digits = 1)); cat("\n")

  cat("# percentage to select ODC by dose category (utility; 3 cat) \n")
  print(matrix(sum.bkt.c3udose$c.pODC, ncol = 3, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over")))); cat("\n")

  cat("# percentage of patients assigned by dose category (utility; 3 cat) \n")
  print(matrix(sum.bkt.c3udose$c.PPTS, ncol = 3, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over")))); cat("\n")

  cat("# percentage to select ODC by dose category (utility; 4 cat) \n")
  print(matrix(sum.bkt.c4udose$c.pODC, ncol = 4, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over", "tox")))); cat("\n")

  cat("# percentage of patients assigned by dose category (utility; 4 cat) \n")
  print(matrix(sum.bkt.c4udose$c.PPTS, ncol = 4, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over", "tox")))); cat("\n")

  cat("# percentage to select ODC by dose category (efficacy; 3 cat) \n")
  print(matrix(sum.bkt.c3edose$c.pODC, ncol = 3, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over")))); cat("\n")

  cat("# percentage of patients assigned by dose category (efficacy; 3 cat) \n")
  print(matrix(sum.bkt.c3edose$c.PPTS, ncol = 3, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over")))); cat("\n")

  cat("# percentage to select ODC by dose category (efficacy; 4 cat) \n")
  print(matrix(sum.bkt.c4edose$c.pODC, ncol = 4, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over", "tox")))); cat("\n")

  cat("# percentage of patients assigned by dose category (efficacy; 4 cat) \n")
  print(matrix(sum.bkt.c4edose$c.PPTS, ncol = 4, byrow = TRUE,
               dimnames = list(paste0("basket", 1:K), c("under", "target", "over", "tox")))); cat("\n")

  cat("# mean number of patients \n")
  print(summary(SUM.N))

  sink()
  # ---- Text summary output (end) ----

  # ---- Build output data frames ----

  # Basket-level base information
  out.bkt <- data.frame(
    method     = rep(method, K),
    ncohort    = rep(ncohort, K),
    cohortsize = cohortsize,
    scenario   = rep(sc.name[SC], K),
    basket     = 1:K
  )

  # Study-level base information
  out.sum <- data.frame(
    method     = method[1],
    ncohort    = ncohort[1],
    cohortsize = cohortsize[1],
    scenario   = sc.name[SC],
    basket     = "pooled"
  )

  # ODC selection rates by dose level (primary metric: estimated prob with MTD)
  out.pODC2       <- data.frame(matrix(pODC2,  ncol = J, byrow = TRUE))
  names(out.pODC2) <- paste0("pODC2_", 1:J)

  # ODC selection rates by dose level (observed prob with MTD)
  out.pODC        <- data.frame(matrix(pODC,   ncol = J, byrow = TRUE))
  names(out.pODC)  <- paste0("pODC_",  1:J)

  # ODC selection rates by dose level (observed prob without MTD)
  out.pODC1       <- data.frame(matrix(pODC1,  ncol = J, byrow = TRUE))
  names(out.pODC1) <- paste0("pODC1_", 1:J)

  # ODC selection rates by dose level (estimated prob without MTD)
  out.pODC3       <- data.frame(matrix(pODC3,  ncol = J, byrow = TRUE))
  names(out.pODC3) <- paste0("pODC3_", 1:J)

  # Mean patients and % patients by dose level
  out.mPPTS <- data.frame(matrix(m.PPTS, ncol = J, byrow = TRUE))
  names(out.mPPTS) <- paste0("m.PPTS", 1:J)

  out.mNPTS <- data.frame(matrix(m.NPTS, ncol = J, byrow = TRUE))
  names(out.mNPTS) <- paste0("m.NPTS", 1:J)

  # ODC and patient summaries by dose category (utility; 3 cat)
  out.c3u.pODC <- data.frame(matrix(sum.bkt.c3udose$c.pODC, ncol = 3, byrow = TRUE))
  names(out.c3u.pODC) <- paste0("c3u.pODC", 1:3)

  out.c3u.NPTS <- data.frame(matrix(sum.bkt.c3udose$c.NPTS, ncol = 3, byrow = TRUE))
  names(out.c3u.NPTS) <- paste0("c3u.NPTS", 1:3)

  out.c3u.PPTS <- data.frame(matrix(sum.bkt.c3udose$c.PPTS, ncol = 3, byrow = TRUE))
  names(out.c3u.PPTS) <- paste0("c3u.PPTS", 1:3)

  # Early stopping rate
  out.earlystop <- data.frame(round(non.ODC / nsim * 100, digits = 1))
  names(out.earlystop) <- "p.earlystop"

  # Total sample size statistics
  out.meanN   <- data.frame(round(mean(SUM.N),   digits = 1)); names(out.meanN)   <- "mean.npat"
  out.sdN     <- data.frame(round(sd(SUM.N),     digits = 2)); names(out.sdN)     <- "sd.npat"
  out.minN    <- data.frame(min(SUM.N));                        names(out.minN)    <- "min.npat"
  out.medN    <- data.frame(round(median(SUM.N), digits = 1)); names(out.medN)    <- "med.npat"
  out.maxN    <- data.frame(max(SUM.N));                        names(out.maxN)    <- "max.npat"

  # ODC and patient summaries by dose category (utility; 4 cat)
  out.c4u.pODC <- data.frame(matrix(sum.bkt.c4udose$c.pODC, ncol = 4, byrow = TRUE))
  names(out.c4u.pODC) <- paste0("c4u.pODC", 1:4)

  out.c4u.NPTS <- data.frame(matrix(sum.bkt.c4udose$c.NPTS, ncol = 4, byrow = TRUE))
  names(out.c4u.NPTS) <- paste0("c4u.NPTS", 1:4)

  out.c4u.PPTS <- data.frame(matrix(sum.bkt.c4udose$c.PPTS, ncol = 4, byrow = TRUE))
  names(out.c4u.PPTS) <- paste0("c4u.PPTS", 1:4)

  # ODC and patient summaries by dose category (efficacy; 3 cat)
  out.c3e.pODC <- data.frame(matrix(sum.bkt.c3edose$c.pODC, ncol = 3, byrow = TRUE))
  names(out.c3e.pODC) <- paste0("c3e.pODC", 1:3)

  out.c3e.NPTS <- data.frame(matrix(sum.bkt.c3edose$c.NPTS, ncol = 3, byrow = TRUE))
  names(out.c3e.NPTS) <- paste0("c3e.NPTS", 1:3)

  out.c3e.PPTS <- data.frame(matrix(sum.bkt.c3edose$c.PPTS, ncol = 3, byrow = TRUE))
  names(out.c3e.PPTS) <- paste0("c3e.PPTS", 1:3)

  # ODC and patient summaries by dose category (efficacy; 4 cat)
  out.c4e.pODC <- data.frame(matrix(sum.bkt.c4edose$c.pODC, ncol = 4, byrow = TRUE))
  names(out.c4e.pODC) <- paste0("c4e.pODC", 1:4)

  out.c4e.NPTS <- data.frame(matrix(sum.bkt.c4edose$c.NPTS, ncol = 4, byrow = TRUE))
  names(out.c4e.NPTS) <- paste0("c4e.NPTS", 1:4)

  out.c4e.PPTS <- data.frame(matrix(sum.bkt.c4edose$c.PPTS, ncol = 4, byrow = TRUE))
  names(out.c4e.PPTS) <- paste0("c4e.PPTS", 1:4)

  # Combine into summary data frames
  out.sum1 <- cbind(out.bkt,
                    out.pODC2, out.mNPTS,
                    out.c3u.pODC, out.c3u.NPTS, out.c3u.PPTS,
                    out.earlystop,
                    out.c4u.pODC, out.c4u.NPTS, out.c4u.PPTS,
                    out.c3e.pODC, out.c3e.NPTS, out.c3e.PPTS,
                    out.c4e.pODC, out.c4e.NPTS, out.c4e.PPTS)

  out.sum2 <- cbind(out.sum,
                    out.meanN, out.sdN, out.minN, out.medN, out.maxN)

  # Return results
  list(
    pt.true      = aaa,
    pe.true      = bbb,
    ut.true      = true.util,
    out.sum1     = out.sum1,
    out.sum2     = out.sum2,
    correct.pODC = correct.pODC
  )
}

##################################################################
# End of program
##################################################################
