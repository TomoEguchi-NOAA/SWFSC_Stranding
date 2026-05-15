# SWFSC_Stranding_fcns.R
# 
# Functions common to many scripts in the SWFSC_Stranding project
library(tidyverse)
library(readr)

connection.string <- function(database){
  return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=",
                database, ";Trusted_Connection=yes;TrustServerCertificate=yes;"))
  
  # return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-ut.nmfs.local;Database=",
  #               database, ";Trusted_Connection=yes;Port=1433;TrustServerCertificate=yes;"))

  }

compute.LOOIC <- function(loglik.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  n.samples <- dim(loglik.array)[1]
  
  # Create an empty list to hold the log likelihood values

  loglik.list <- lapply(1:n.samples, function(s) {
    # Extract the slice for the current sample 's'
    slice <- loglik.array[s, , , ]
    # Return the non-NA values from that slice
    slice[!is.na(slice)]
  })
  
  loglik.mat <- do.call(rbind, loglik.list)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = MCMC.params$n.chains)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = MCMC.params$n.chains, 
                           k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}

# Compute the "rank-normalized R-hat" by Vehtari et al. (2021) from jagsUI
# output.
# Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667–718.
# https://doi.org/10.1214/20-BA1221
# 
# The first input is MCMC samples. If the jagsUI output is jm, this is jm$samples.
# The second input is a string of regular expression. This is a bit
# complicated. For example, to select "BF.Fixed" and all K parameters, which 
# are indexed, Use "^BF\\.Fixed|^K\\[" A '^' specifies that the following letter
# is the beginning of a string. '\\.' specifies a literal period, which needs
# to be "escaped" by two backslashes (\\). A square bracket needs to be escaped
# with two backslashes as well. The pipe (|) indicates 'or'.  
# 
rank.normalized.R.hat <- function(samples, params, MCMC.params){
  library(posterior)
  library(coda)
  
  col.names <- grep(params, varnames(samples), value = TRUE, perl = TRUE)
  subset.mcmc.samples <- samples[, col.names]
  
  subset.mcmc.array <- as_draws_array(subset.mcmc.samples, .nchains = MCMC.params$n.chains)
  
  rhat.values <- apply(subset.mcmc.array, 
                       MARGIN = 3, 
                       FUN = posterior::rhat)
  
  return(rhat.values)
}

summary.logistic <- function(jags.out,
                             jags.data,
                             xvar = "Age",
                             sp.names){
  # 1. Extract the posterior means for the B matrix
  # Assuming 'samples' is your mcmc.list object
  post_summary <- tidybayes::summarise_draws(jags.out$samples) %>%
    filter(grepl("B\\[", variable)) %>%
    select(variable, mean)
  
  # 1. Corrected Coefficients Data Frame (as before)
  plot_coeffs <- data.frame(
    species_idx = rep(1:3, times = 2),
    param = rep(c("intercept", "slope"), each = 3),
    value = post_summary$mean
  ) %>%
    tidyr::pivot_wider(names_from = param, values_from = value)
  
  # 2. Generate smooth curves for the logistic fit
  # We'll create a grid of ages for each species
  X_range <- seq(min(jags.data$X), max(jags.data$X), length.out = 100)
  pred_data <- plot_coeffs %>%
    group_by(species_idx) %>%
    do(data.frame(
      X = X_range,
      prob = 1 / (1 + exp(-(.$intercept + .$slope * X_range)))
    )) %>%
    mutate(species = factor(species_idx, levels = 1:3, labels = sp.names))
  
  X_grid <- seq(min(jags.data$X), max(jags.data$X), length.out = 50)
  
  # 1. Extract all posterior draws for the B matrix
  # This creates a long-format data frame with columns: .chain, .iteration, .draw, B[row,col]
  draws_df <- jags.out$samples %>%
    tidybayes::spread_draws(B[spec_idx, param_idx]) %>% 
    ungroup() %>%
    # Pivot to get intercept and slope in separate columns
    mutate(param_name = ifelse(param_idx == 1, "intercept", "slope")) %>%
    select(-param_idx) %>%
    pivot_wider(names_from = param_name, values_from = B)
  
  # 2. Calculate the curves for every draw
  # Note: This can be computationally heavy if you have many draws; 
  # we'll use a subset of the grid for efficiency
  fitted_draws <- draws_df %>%
    rename(species_idx = spec_idx) %>%
    # Crossing joins every draw with every age in our grid
    crossing(X = X_grid) %>%
    mutate(
      # The inverse logit calculation for every single MCMC draw
      prob = 1 / (1 + exp(-(intercept + slope * X)))
    ) %>%
    # 3. Summarize the distribution at each age point per species
    group_by(species_idx, X) %>%
    summarise(
      low = quantile(prob, 0.025),
      high = quantile(prob, 0.975),
      median_prob = median(prob),
      .groups = "drop"
    ) %>%
    # Apply species names
    mutate(species = factor(species_idx, levels = 1:3, labels = sp.names))
  
  # 3. Prepare observed data and apply labels
  obs_data <- data.frame(
    X = jags.data$X,
    mature = jags.data$mature,
    species_idx = jags.data$species_idx
  ) %>%
    mutate(species = factor(species_idx, levels = 1:3, labels = sp.names))
  
  a50_summary <- draws_df %>%
    mutate(a50_draw = -intercept / slope) %>%
    group_by(spec_idx) %>%
    summarise(
      median_a50 = median(a50_draw),
      low_a50 = quantile(a50_draw, 0.025),
      high_a50 = quantile(a50_draw, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      species = factor(spec_idx, levels = 1:3, labels = sp.names),
      # Coordinates for the segments
      x_start = median_a50, x_end = median_a50, y_start = 0, y_end = 0.5, # Vertical
      x_h_start = 0, x_h_end = median_a50, y_h_start = 0.5, y_h_end = 0.5  # Horizontal
    )
  
  # 4. Final Plot with Named Facets
  # 5. Plot with Ribbon
  p.logistic <- ggplot(fitted_draws, aes(x = X)) +
    # 95% Credible Interval Ribbon
    geom_ribbon(aes(ymin = low, ymax = high, fill = species), alpha = 0.2) +
    # Median Line
    geom_line(aes(y = median_prob, color = species), size = 1) +
    # Horizontal line at 50% probability
    # 3. Add Segments (Vertical then Horizontal)
    geom_segment(data = a50_summary, 
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linetype = "dashed", color = "gray40") +
    #geom_segment(data = a50_summary, 
    #             aes(x = x_h_start, xend = x_h_end, y = y_h_start, yend = y_h_end),
    #             linetype = "dashed", color = "gray40") +
    # Text Label for A50 value
    geom_text(data = a50_summary, 
              aes(x = median_a50, y = 0.55, label = paste0("A50: ", round(median_a50, 2))),
              hjust = -0.1, color = "firebrick", fontface = "bold") +
    # Observed Data Points
    geom_jitter(data = obs_data, aes(x = X, y = mature), 
                height = 0.03, width = 0, alpha = 0.4, color = "gray30") +
    facet_wrap(~species) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = "Maturity Curves with 95% Credible Intervals",
         subtitle = "Shaded areas represent posterior uncertainty",
         x = xvar,
         y = "Probability of Maturity") +
    theme_bw() +
    theme(legend.position = "none")
  
  return(list(plot = p.logistic,
              A50.summary = a50_summary,
              data = obs_data,
              draws = fitted_draws,
              pred.data = pred_data))
}


species.col.types <- cols(ID = col_integer(),
                          TaxanomicOrder = col_character(),
                          SubOrder = col_character(),
                          InfraOrder = col_character(),
                          Family = col_character(),
                          FamilyCommonName = col_character(),
                          Genus = col_character(),
                          Species = col_character(),
                          Subspecies = col_character(),
                          CommonName = col_character(),
                          SpType = col_character(),
                          SpName = col_character(),
                          SpCode = col_character(),
                          AerialFisheryCode = col_character(),
                          BirderCode = col_character(),
                          AlternateCode1 = col_character(),
                          AlternateCode2 = col_character(),
                          EditDate = col_datetime(),
                          EditUser = col_character(),
                          RecordCreationDate = col_datetime())

morph.col.types <- cols(Specimen = col_character(),
                        TotalLength_LAB = col_double(),
                        IsStandardTL_LAB = col_factor(levels = c("Y", "N", "y", "n")),
                        TotalLength_FIELD = col_double(),
                        IsStandardTL_FIELD = col_factor(levels = c("Y", "N", "y", "n")),
                        IsAltMeasurementDeviceUsed = col_factor(levels = c("Y", "N", "y", "n")),
                        Spotter_Color = col_integer(),
                        Spinner_Morph = col_integer(),
                        Spinner_Cape = col_integer(),
                        Spinner_Fin = col_integer(),
                        Spinner_Belly = col_integer(),
                        STOANUS = col_double(),
                        STOGENSLIT = col_double(),
                        STOUMBIL = col_double(),
                        STOTHRGROO = col_double(),
                        STODOFINTIP = col_double(),
                        STOANTDOR = col_double(),
                        STOFLIPPER = col_double(),
                        STOEAR = col_double(),
                        STOEYE = col_double(),
                        STOGAPE = col_double(),
                        STOBLOHOLE = col_double(),
                        STOMELAPEX = col_double(),
                        ETOEAR = col_double(),
                        ETOGAPE = col_double(),
                        ETOBLOHL_L = col_double(),
                        ETOBLOHL_R = col_double(),
                        BLOHL_LGTH = col_double(),
                        BLOHLWDTH = col_double(),
                        DIAM_EAR = col_character(),
                        HEAD_DIAM = col_double(),
                        LGTH_EOP = col_double(),
                        ROSTWIDTH = col_double(),
                        PROJECTUP = col_character(),
                        PROJECTLOW = col_character(),
                        THROATGROOVE_N = col_integer(),
                        LGTH_GROO = col_double(),
                        FLIPLGTH_A = col_double(),
                        FLIPLGTH_P = col_double(),
                        FLIPWIDTH = col_double(),
                        LGTHMAMS_R = col_double(),
                        LGTHMAMS_L = col_double(),
                        MAMMSLIT_N = col_integer(),
                        LGTHGENSLI = col_double(),
                        LGTHANASLI = col_double(),
                        PERILGTH = col_double(),
                        FLUKWDTH = col_double(),
                        FLUKDPTH_L = col_double(),
                        FLUKDPTH_N = col_double(),
                        FLUKNTDPTH = col_double(),
                        DORFNHGT = col_double(),
                        DORFNBLGTH = col_double(),
                        GATEYE = col_double(),
                        GAAXILLA = col_double(),
                        GIRTHMAX = col_double(),
                        GATANUS = col_double(),
                        GMIDANTONT = col_double(),
                        HGTSMPL = col_double(),
                        THICKSMPL = col_double(),
                        BLUBTHIK_D = col_double(),
                        BLUBTHIK_L = col_double(),
                        BLUBTHIK_V = col_double(),
                        BLUBTHIK_C = col_double(),
                        BlubberCompData = col_factor(levels = c("Y", "N", "y", "n")),
                        Pinn_CurvilinearLength = col_double(),
                        Pinn_FlipperLn_ForeAnt = col_double(),
                        Pinn_FlipperLn_HindAnt = col_double(),
                        EditDate = col_datetime(),
                        EditUser = col_character(),
                        RecordCreationDate = col_datetime())

animal.col.types <- cols(Specimen = col_character(),
                         OriginID = col_integer(),
                         IsSWFSC = col_factor(levels = c("Y", "N", "y", "n")),
                         IsDataSheet = col_factor(levels = c("Y", "N", "y", "n")),
                         Cruise = col_character(),
                         Cruise_Set = col_character(),
                         Year = col_integer(),
                         Month = col_integer(),
                         Day = col_integer(),
                         Latitude = col_double(),
                         Latitude_Precision = col_double(),
                         Latitude_Precision_Unit = col_character(),
                         Longitude = col_double(),
                         Longitude_Precision = col_double(),
                         Longitude_Precision_Unit = col_character(),
                         CityID = col_integer(),
                         CountyID = col_integer(),
                         StateID = col_integer(),
                         CountryID = col_integer(),
                         SpeciesID = col_character(),
                         SpDeterminationID = col_character(),
                         Sex = col_factor(levels = c("F", "M", "U", "f", "m", "u")),
                         Adrenals_Coll = col_character(),
                         Genetics_Biopsy_Coll = col_character(),
                         Blood_Coll = col_character(),
                         Blubber_Coll = col_character(),
                         Brain_Coll = col_character(),
                         Carcass_Coll = col_character(),
                         Feces_Coll = col_character(),
                         Fetus_Coll = col_character(),
                         FetusGenetics_Coll = col_character(),
                         Gonad_Coll = col_character(),
                         Head_Coll = col_character(),
                         Histo_Coll = col_character(),
                         Kidney_Coll = col_character(),
                         Liver_Coll = col_character(),
                         Lung_Coll = col_character(),
                         Morphometry_Coll = col_character(),
                         Muscle_Coll = col_character(),
                         Photos_Coll = col_character(),
                         Radiology_Coll = col_character(),
                         Skeleton_Coll = col_character(),
                         Spleen_Coll = col_character(),
                         Stomach_Coll = col_character(),
                         Teeth_Coll = col_character(),
                         Urine_Coll = col_character(),
                         Other_Coll = col_character(),
                         EditDate = col_datetime(),
                         EditUser = col_character(),
                         RecordCreationDate = col_datetime())

age.col.types <- cols(ID = col_integer(),
                      Specimen = col_character(),
                      Age = col_double(),
                      IsAnalysisQuality = col_factor(levels = c("Y", "N")),
                      EstimationMethod = col_factor(),
                      AgeReader1 = col_double(), 
                      ReaderID1 = col_integer(),
                      AgeReader2 = col_double(),
                      ReaderID2 = col_integer(),
                      AgeReader3 = col_double(),
                      ReaderID3 = col_integer(),
                      AgeReader4 = col_double(),
                      ReaderID4 = col_integer(),
                      EditDate = col_datetime(),
                      EditUser = col_character(),
                      RecordCreationDate= col_datetime())

repro.col.types <- cols(Specimen = col_character(),
                        IsSideKnown = col_factor(levels = c("Y", "N", "y", "n")),
                        IsMature = col_factor(levels = c("Y", "N", "U", "y", "n", "u")),
                        MaturityID = col_integer(),
                        IsLactating = col_factor(levels = c("Y", "N", "y", "n")),
                        IsPregnant = col_factor(levels = c("Y", "N", "y", "n")),
                        Follicle_Diam = col_double(),
                        OvaryWeight_R = col_double(),
                        OvaryWeight_L = col_double(),
                        OvaryLength_R = col_double(),
                        OvaryWidth_R = col_double(),
                        OvaryDepth_R = col_double(),
                        OvaryLength_L = col_double(),
                        OvaryWidth_L = col_double(),
                        OvaryDepth_L = col_double(),
                        CL_Diam1 = col_integer(),
                        CL_Diam2 = col_integer(),
                        CL_Diam3 = col_integer(),
                        CL_InternalDiam1 = col_integer(),
                        CL_InternalDiam2 = col_integer(),
                        CA1_R = col_integer(),
                        CA2_R = col_integer(),
                        CA3_R = col_integer(),
                        CA4_R = col_integer(),
                        CA5_R = col_integer(),
                        CA6_R = col_integer(),
                        CA_RIGHT = col_integer(),
                        CA1_L = col_integer(),
                        CA2_L = col_integer(),
                        CA3_L = col_integer(),
                        CA4_L = col_integer(),
                        CA5_L = col_integer(),
                        CA6_L = col_integer(),
                        CA_LEFT = col_integer(),
                        TotalCorpora = col_integer(),
                        CL_LocationID = col_integer(),
                        FetusLength_Standard = col_double(),
                        FetusLength_Curvilinear = col_double(),
                        FetusSex = col_factor(levels = c("M", "F", "U", "m", "f", "u")),
                        FetusWeight = col_double(),
                        WeightWEpi_L = col_double(),
                        WeightWEpi_R = col_double(),
                        WeightWOEPI_R = col_double(),
                        WeightWOEPI_L = col_double(),
                        TestisLength_R = col_double(),
                        TestisWidth_R = col_double(),
                        TestisDepth_R = col_double(),
                        TestisLength_L = col_double(),
                        TestisWidth_L = col_double(),
                        TestisDepth_L = col_double(),
                        EditDate = col_datetime(),
                        EditUser = col_character(),
                        RecordCreationDate = col_datetime())

weight.col.types <- cols(Specimen = col_character(),
                         Adrenal_1 = col_double(),
                         AdrenalSide_1 = col_character(),
                         Adrenal_2 = col_double(),
                         AdrenalSide_2 = col_character(),
                         Blubber = col_double(),
                         Brain = col_double(),
                         Carcass_Intact = col_double(),
                         Epaxial = col_double(),
                         Heart = col_double(),
                         Hypaxial = col_double(),
                         Intestines = col_double(),
                         Intestine_Length = col_double(),
                         Kidney_R = col_double(),
                         Kidney_L = col_double(),
                         Liver = col_double(),
                         Lung_R = col_double(),
                         Lung_L = col_double(),
                         Misc = col_double(),
                         Muscle_Total = col_double(),
                         Pancreas = col_double(),
                         Spleen = col_double(),
                         Stomach_Full = col_double(),
                         Stomach_Empty = col_double(),
                         Thymus = col_double(),
                         Viscera = col_double(),
                         EditDate = col_datetime(),
                         EditUser = col_character(),
                         RecordCreationDate = col_datetime())