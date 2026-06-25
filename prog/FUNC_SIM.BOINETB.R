####################################################################
# BOIN-ETB
# Code to define function SIM.BOINETB (simulation)
####################################################################
SIM.BOINETB <- function(
  SC,          # Scenario number
  pt.true,     # True toxicity probability vector (all doses x baskets)
  pe.true,     # True efficacy probability vector (all doses x baskets)
  startdose,   # Starting dose index for each basket
  ncohort,     # Maximum number of cohorts
  cohortsize,  # Cohort size per basket
  nsim,        # Number of simulations
  mn,          # Minimum number of patients to stop each basket
  mnodc,       # Minimum number of patients required for ODC selection
  tgt.t,       # Target toxicity probability per basket
  tgt.e,       # Target efficacy probability per basket
  dev,         # Allowable deviation from the target
  K,           # Number of baskets
  J,           # Number of dose levels
  FILENAME     # Output file name (without extension)
) {

  ### Derive parameters
  numdose   <- 1:(K * J)           # Sequential index of all dose levels
  doselevel <- rep(1:J, K)         # Dose level corresponding to each index
  basket    <- rep(1:K, each = J)  # Basket number corresponding to each index

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

  # Extract thresholds
  lambda.t1 <- optim_df$lambda1  # Toxicity: escalation threshold
  lambda.t2 <- optim_df$lambda2  # Toxicity: de-escalation threshold
  lambda.e  <- optim_df$eta1     # Efficacy: threshold

  ### Prepare OUT.STUDY
  OUT.STUDY <- rbind(
    rep(rep(1:K, each = J), 14),
    rep(rep(1:J, K), 14)
  )
  colnames(OUT.STUDY) <- paste(
    rep(c("npts", "ntox", "neff",
          "ppts", "ptox", "peff",
          "nct1", "nct2", "nct3", "nct4",
          "pct1", "pct2", "pct3", "pct4"),
        each = K * J),
    rep(rep(1:K, each = J), 14),
    rep(rep(1:J, K), 14),
    sep = ""
  )

  ### Run simulations
  set.seed(20230502)

  for (nnn in 1:nsim) {

    ndose <- length(numdose)

    # Initialize variables
    # y.t   : Number of toxicity events per dose
    # y.e   : Number of efficacy events per dose
    # y.c1  : No toxicity & efficacy (category 1)
    # y.c2  : No toxicity & no efficacy (category 2)
    # y.c3  : Toxicity & efficacy (category 3)
    # y.c4  : Toxicity & no efficacy (category 4)
    # n     : Number of patients treated per dose
    # elm.bd: Early termination flag per dose level (1 = terminated)
    # elm.b : Early termination flag per basket (1 = terminated)
    y.t <- y.e <- y.c1 <- y.c2 <- y.c3 <- y.c4 <- n <- elm.bd <- rep(0, J * K)
    elm.b <- numeric(K)

    d  <- startdose  # Current dose index for each basket
    dd <- d[1]       # Dose index currently being processed

    for (icohort in 1:ncohort) {

      ### Simulate toxicity and efficacy data by basket
      for (ibkt in 1:K) {

        dd <- d[ibkt]          # Current dose for basket ibkt
        nn <- cohortsize[ibkt] # Cohort size for basket ibkt

        if (elm.b[ibkt] == 0) {  # Skip if basket is already terminated

          if (elm.bd[dd] == 0) {  # Skip if dose is already terminated

            # Simulate toxicity and efficacy events
            r.y.t <- runif(nn) < pt.true[dd]
            r.y.e <- runif(nn) < pe.true[dd]

            y.t[dd] <- y.t[dd] + sum(r.y.t)  # Accumulate toxicity count
            y.e[dd] <- y.e[dd] + sum(r.y.e)  # Accumulate efficacy count
            n[dd]   <- n[dd] + nn             # Accumulate patient count

            # Accumulate outcome category counts
            y.c1[dd] <- y.c1[dd] + sum(r.y.t == 0 & r.y.e == 1)  # No tox, eff
            y.c2[dd] <- y.c2[dd] + sum(r.y.t == 0 & r.y.e == 0)  # No tox, no eff
            y.c3[dd] <- y.c3[dd] + sum(r.y.t == 1 & r.y.e == 1)  # Tox, eff
            y.c4[dd] <- y.c4[dd] + sum(r.y.t == 1 & r.y.e == 0)  # Tox, no eff

            # Early termination criteria for current dose level in basket ibkt
            # Condition: toxicity too high, or efficacy negligible
            if (pbeta(tgt.t[ibkt], y.t[dd] + 1, n[dd] - y.t[dd] + 1) < 0.05 ||
                pbeta(lambda.e[ibkt], y.e[dd] + 1, n[dd] - y.e[dd] + 1) > 0.95) {
              elm.bd[dd] <- 1
            }
          }

          # Summarize toxicity counts by dose level across all baskets
          t.dl <- data.frame(
            doselevel = doselevel,
            y.t       = y.t,
            n         = n
          )
          sum.t.dl <- t.dl %>%
            group_by(doselevel) %>%
            summarize(n = sum(n), y.t = sum(y.t))

          # Observed toxicity and efficacy rates
          p.t.dl   <- sum.t.dl$y.t / sum.t.dl$n  # Toxicity rate by dose level
          p.e.list <- y.e / n                      # Efficacy rate by dose and basket

          # Early termination criteria for each dose level across all baskets
          for (jj in 1:J) {
            if (pbeta(tgt.t[ibkt],
                      sum.t.dl$y.t[jj] + 1,
                      sum.t.dl$n[jj] - sum.t.dl$y.t[jj] + 1) < 0.05) {
              elm.bd[jj + J * (0:(K - 1))] <- 1
            }
          }
        }
      }

      ### Select next dose by basket
      for (ibkt in 1:K) {

        dd <- d[ibkt]        # Current dose index
        dl <- doselevel[dd]  # Current dose level

        if (elm.b[ibkt] == 0) {  # Skip if basket is already terminated

          ### Select admissible doses
          if (p.t.dl[dl] <= lambda.t1[ibkt] &
              p.e.list[dd] <= lambda.e[ibkt]) {
            # Low toxicity, low efficacy -> escalate by one level
            if (doselevel[dd] == J) { ad.ibkt <- dd
            } else                  { ad.ibkt <- dd + 1 }

          } else if (lambda.t1[ibkt] < p.t.dl[dl] &
                     p.t.dl[dl] < lambda.t2[ibkt] &
                     p.e.list[dd] <= lambda.e[ibkt]) {
            # Moderate toxicity, low efficacy -> consider neighboring doses
            if      (doselevel[dd] == 1) { ad.ibkt <- c(dd, dd + 1)
            } else if (doselevel[dd] == J) { ad.ibkt <- c(dd - 1, dd)
            } else                          { ad.ibkt <- c(dd - 1, dd, dd + 1) }

          } else if (p.t.dl[dl] < lambda.t2[ibkt] &
                     p.e.list[dd] > lambda.e[ibkt]) {
            # Toxicity within acceptable range, high efficacy -> stay at current dose
            ad.ibkt <- dd

          } else if (p.t.dl[dl] >= lambda.t2[ibkt]) {
            # High toxicity -> de-escalate, or stop basket if at lowest dose
            if (doselevel[dd] == 1) {
              ad.ibkt <- dd
              elm.bd[dl + J * (0:(K - 1))] <- 1
              elm.b[ibkt] <- 1
            } else {
              ad.ibkt <- dd - 1
            }
          }

          ### Select next dose among admissible doses
          if (length(ad.ibkt) == 1) {
            d[ibkt] <- ad.ibkt

          } else if ((dd + 1) %in% ad.ibkt & n[dd + 1] == 0) {
            # Prioritize the next untreated dose if available
            d[ibkt] <- dd + 1

          } else {
            # Select the dose with the highest observed efficacy rate
            # (random selection if tied)
            ddd <- numdose[
              p.e.list == max(p.e.list[ad.ibkt]) &
              numdose %in% ad.ibkt &
              basket == ibkt
            ]
            if (length(ddd) == 1) {
              d[ibkt] <- ddd
            } else {
              d[ibkt] <- sample(ddd, size = 1)
            }
          }

          ### Stop basket if number of patients at next dose reaches minimum
          if (n[d[ibkt]] >= mn[ibkt]) { elm.b[ibkt] <- 1 }
        }
      }
    }

    # Store simulation results grouped by category
    out.study <- c(
      n,                                     # Number of patients treated
      y.t, y.e,                              # Toxicity and efficacy counts
      round(n / sum(n),        digits = 3),  # Proportion of patients
      round(y.t / n,           digits = 3),  # Observed toxicity rate
      round(y.e / n,           digits = 3),  # Observed efficacy rate
      y.c1, y.c2, y.c3, y.c4,               # Outcome category counts
      round(y.c1 / n,          digits = 3),  # Outcome category proportions
      round(y.c2 / n,          digits = 3),
      round(y.c3 / n,          digits = 3),
      round(y.c4 / n,          digits = 3)
    )
    OUT.STUDY <- rbind(OUT.STUDY, out.study)

    # Record early termination flags by basket
    if (nnn == 1) {
      STOP.STUDY <- elm.b
      names(STOP.STUDY) <- paste0("basket", 1:K)
    } else {
      STOP.STUDY <- rbind(STOP.STUDY, elm.b)
    }
  }

  # Output results to CSV files
  write.csv(
    x         = OUT.STUDY,
    file      = paste(OUTFILE, FILENAME, "_data.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    x         = STOP.STUDY,
    file      = paste(OUTFILE, FILENAME, "_stop.csv", sep = ""),
    row.names = FALSE
  )

  # Calculate summary statistics
  aaa <- CALC.SUM(OUTFILE, FILENAME, mnodc, SC)

  pt.true      <- aaa$pt.true
  pe.true      <- aaa$pe.true
  ut.true      <- aaa$ut.true
  out.sum1.all <- aaa$out.sum1
  out.sum2.all <- aaa$out.sum2

  ### Output to Excel file ###
  wb <- createWorkbook()
  addWorksheet(wb, "Sum1")
  addWorksheet(wb, "Sum2")
  writeDataTable(wb, "Sum1", x = out.sum1.all)
  writeDataTable(wb, "Sum2", x = out.sum2.all)
  saveWorkbook(
    wb,
    paste(OUTFILE, "summary_", FILENAME, ".xlsx", sep = ""),
    overwrite = TRUE
  )

  # Return results
  list(
    pt.true    = pt.true,
    pe.true    = pe.true,
    ut.true    = ut.true,
    out.sum1   = out.sum1.all,
    out.sum2   = out.sum2.all,
    OUT.STUDY  = OUT.STUDY,
    STOP.STUDY = STOP.STUDY
  )
}

##################################################################
# End of program
##################################################################
