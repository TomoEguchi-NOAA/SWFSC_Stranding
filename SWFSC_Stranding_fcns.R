# SWFSC_Stranding_fcns.R
# 
# Functions common to many scripts in the SWFSC_Stranding project
library(tidyverse)
library(readr)

jags.logistic <- function(jags.data,
                          MCMC.params,
                          out.filename,
                          jags.model,
                          jags.params,
                          inits = NULL,
                          xvar,
                          yvar,
                          sp.names){
  if (!file.exists(out.filename)){
    
    tic <- Sys.time()
    jm. <- jagsUI::jags(jags.data,
                        inits = inits,
                        parameters.to.save= jags.params,
                        model.file = jags.model,
                        n.chains = MCMC.params$n.chains,
                        n.burnin = MCMC.params$n.burnin,
                        n.thin = MCMC.params$n.thin,
                        n.iter = MCMC.params$n.samples,
                        DIC = T,
                        parallel=T)
    
    toc <- Sys.time()
    Rhat. <- rank.normalized.R.hat(jm.$samples,
                                   params = "^mu_|^L|^slope",
                                   MCMC.params = MCMC.params)
    
    post <- posterior::as_draws(jm.$samples)
    summary.posterior <- posterior::summarise_draws(post)
    
    jm.summary <- summary.logistic.no.B(jags.out = jm.,
                                        jags.data = jags.data, 
                                        xvar = xvar,
                                        yvar = yvar,
                                        sp.names = sp.names)
    
    jags.out <- list(jags.model = jags.model,
                     jags.data = jags.data,
                     MCMC.params = MCMC.params,
                     jags.out = jm.,
                     Rhat = Rhat.,
                     posterior.summary = summary.posterior,
                     logistic.summary = jm.summary,
                     Run.Date = Sys.Date(),
                     Run.Time = toc - tic,
                     System = Sys.getenv())
    
    saveRDS(jags.out, 
            file = out.filename)
  } else {
    jags.out <- read_rds(out.filename)
  } 
  return(jags.out)
}
  
jags.Laird.growth <- function(LAB.data,
                              length.age.data,
                              MCMC.params,
                              out.filename,
                              jags.model,
                              jags.params,
                              inits.fcn){
  
  if (!file.exists(out.filename)){
    
    # Define an identity matrix for the Wishart scale matrix
    #R_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
    
    jags.data <- list(
      # Logistic Regression Data -
      N_post = nrow(LAB.data),
      is_postnatal = LAB.data$Status,       # Must be 0 (prenatal) or 1 (postnatal)
      length_post = LAB.data$Length,
      sp_post = LAB.data$Species.int,
      
      # Growth Curve Data
      N_growth = nrow(length.age.data),
      age_growth = length.age.data$Age,
      length_growth = length.age.data$TotalLength_FIELD,
      sp_growth = length.age.data$species_idx,
      sex_growth = length.age.data$Sex_idx,
      
      # Shared Constants
      N_species = max(length.age.data$species_idx),
      N_sexes = 2
    )
    
    tic <- Sys.time()
    jm. <- jagsUI::jags(jags.data,
                        inits = lapply(1:MCMC.params$n.chains, 
                                       function(i) inits.fcn()),
                        parameters.to.save= jags.params,
                        model.file = jags.model,
                        n.chains = MCMC.params$n.chains,
                        n.burnin = MCMC.params$n.burnin,
                        n.thin = MCMC.params$n.thin,
                        n.iter = MCMC.params$n.samples,
                        DIC = T,
                        parallel=T)
    toc <- Sys.time()
    
    Rhat. <- rank.normalized.R.hat(jm.$samples,
                                   params = "^B1|^L0|^a|^tc|^s",
                                   MCMC.params = MCMC.params)
    
    post <- posterior::as_draws(jm.$samples)
    summary.posterior <- posterior::summarise_draws(post)
    
    jags.out. <- list(jags.model = jags.model,
                      jags.data = jags.data,
                      MCMC.params = MCMC.params,
                      jags.out = jm.,
                      Rhat = Rhat.,
                      posterior.summary = summary.posterior,
                      Run.Date = Sys.Date(),
                      Run.Time = toc - tic,
                      System = Sys.getenv())
    
    saveRDS(jags.out.age.length, 
            file = jags.out.filename.age.length.sex)
  } else {
    jags.out. <- read_rds(jags.out.filename.age.length)
  }   
  
  return(jags.out)
}



connection.string <- function(database){
  # return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=",
  #               database, ";Trusted_Connection=yes;TrustServerCertificate=yes;"))
  
  # return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-ut.nmfs.local;Database=",
  #               database, ";Trusted_Connection=yes;Port=1433;TrustServerCertificate=yes;"))
  
  return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-g.nmfs.local;Database=",
                database, ";Trusted_Connection=yes;Port=1433;TrustServerCertificate=yes;"))

}

# read tables that were extracted by sQL_data_extraction.R
read.tables <- function(data.extraction.date){
  # All column types are stored below
  # There were some extensive editing (mostly deleting entries) for the Species
  # table. So, I'm not going to change the data file from 2026-05-12
  table.Species <- read_csv(file = paste0("Data/tblSpecies.csv"),
                            col_types = species.col.types) %>%
    rename(SpeciesID = SpCode) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  table.Morph <- read_csv(file = paste0("Data/tbl_Morphology_", 
                                        data.extraction.date, ".csv"),
                          col_types = morph.col.types) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  table.Animal <- read_csv(file = paste0("Data/tbl_Animal_", 
                                         data.extraction.date, ".csv"),
                           col_types = animal.col.types) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  table.Age <- read_csv(file = paste0("Data/tbl_Age_", 
                                      data.extraction.date, ".csv"),
                        col_types = age.col.types) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  table.Age %>%
    filter(IsAnalysisQuality == "Y") %>%
    left_join(table.Animal, by = "Specimen") -> table.Age.1
  
  table.Repro <- read_csv(file = paste0("Data/tbl_Reproduction_", 
                                        data.extraction.date, ".csv"),
                          col_types = repro.col.types) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  table.Weight <- read_csv(file = paste0("Data/tbl_Weight_", 
                                         data.extraction.date, ".csv"),
                           col_types = weight.col.types) %>%
    select(-c(EditDate, EditUser, RecordCreationDate))
  
  # Remove a fin whale that had very small Lab measurement
  # The filter removes if ratio = NA. Need to explicitly keep NAs. 
  table.Morph  %>%
    select(Specimen, TotalLength_LAB, TotalLength_FIELD, IsStandardTL_FIELD,
           IsStandardTL_LAB, GIRTHMAX) %>%
    mutate(ratio = TotalLength_FIELD/TotalLength_LAB) %>%
    filter(ratio < 2 | is.na(ratio)) -> table.Morph.1
  
  Field.Lab.lm <- lm(TotalLength_LAB ~ TotalLength_FIELD, data = table.Morph.1)
  
  table.Morph.1 %>%
    left_join(table.Animal %>% 
                select(Specimen, SpeciesID, Year, Latitude, Longitude),
              by = "Specimen") %>% 
    left_join(table.Species %>% 
                select(SpeciesID, Genus, Species, CommonName, SpName), 
              by = "SpeciesID") %>% 
    filter(Genus != "Unid") %>%
    mutate(Genus.f = as.factor(Genus)) %>% 
    left_join(table.Weight %>%
                select(Specimen, Carcass_Intact, Heart, Kidney_R, Kidney_L, Liver), 
              by = "Specimen") %>%
    left_join(table.Age.1 %>%
                select(Specimen, Age, IsAnalysisQuality, EstimationMethod), 
              by = "Specimen") %>%
    left_join(table.Repro %>% 
                select(Specimen, IsMature, MaturityID, IsLactating, IsPregnant,
                       CA_LEFT, CA_RIGHT, TotalCorpora, FetusLength_Standard,
                       FetusWeight), 
              by = "Specimen") -> table.Morph.2  

  return(list(Morph.2 = table.Morph.2,
              Morph.1 = table.Morph.1,
              Morph = table.Morph,
              Weight = table.Weight,
              Repro = table.Repro,
              Age.1 = table.Age.1,
              Age = table.Age,
              Animal = table.Animal,
              Species = table.Species))  
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
                             yvar = "Probability of Maturity",
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
    y = jags.data$y,
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
    geom_line(aes(y = median_prob, color = species), linewidth = 1) +
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
    geom_jitter(data = obs_data, aes(x = X, y = y), 
                height = 0.03, width = 0, alpha = 0.4, color = "gray30") +
    facet_wrap(~species) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = "Maturity Curves with 95% Credible Intervals",
         subtitle = "Shaded areas represent posterior uncertainty",
         x = xvar,
         y = yvar) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(list(plot = p.logistic,
              A50.summary = a50_summary,
              data = obs_data,
              draws = fitted_draws,
              pred.data = pred_data))
}

# no B matrix logistic model
summary.logistic.no.B <- function(jags.out,
                                  jags.data,
                                  xvar = "Age",
                                  yvar = "Probability of Maturity",
                                  sp.names){
  
  # 1. Extract the MCMC samples into a matrix
  # Assuming your coda.samples output is named 'mcmc_logistic'
  post_matrix <- as.matrix(jags.out$samples)
  
  # 2. Create a smooth sequence of lengths for the x-axis
  # Adjust the min and max lengths based on your actual data
  min_len <- min(jags.data$X, na.rm = TRUE)
  max_len <- max(jags.data$X, na.rm = TRUE)
  length_seq <- seq(min_len, max_len, length.out = 100)
  
  # 3. Create an empty dataframe to store the calculated curve data
  logistic_plot_data <- data.frame()
  L0_summary_data <- data.frame()
  
  # 4. Loop through the 3 species to build the prediction intervals
  for (s in 1:3) {
    
    # Extract the posteriors for the current species
    L0_samp <- post_matrix[, paste0("L0[", s, "]")]
    slope_samp <- post_matrix[, paste0("slope[", s, "]")]
    
    # Create a matrix to hold the probability predictions
    # (Rows = MCMC iterations, Columns = length points)
    p_matrix <- matrix(NA, nrow = nrow(post_matrix), ncol = length(length_seq))
    
    # Calculate the probability for every length and every MCMC sample
    for (i in 1:length(length_seq)) {
      X <- length_seq[i]
      
      # Apply the inverse-logit math formula directly
      p_matrix[, i] <- 1 / (1 + exp(-slope_samp * (X - L0_samp)))
    }
    
    # 5. Summarize the thousands of predictions into median and 95% intervals
    species_df <- data.frame(
      species = s,
      X = length_seq,
      fit = apply(p_matrix, 2, quantile, probs = 0.500),
      lwr = apply(p_matrix, 2, quantile, probs = 0.025),
      upr = apply(p_matrix, 2, quantile, probs = 0.975)
    )
    
    logistic_plot_data <- rbind(logistic_plot_data, species_df)
    
    L0_df <- data.frame(
      species = s,
      L0_med = median(L0_samp)
    )
    L0_summary_data <- rbind(L0_summary_data, L0_df)
    
  }
  
  # 6. Apply factor labels for clean ggplot facet titles
  #species_names <- c("S. attenuata", "S. coeruleoalba", "S. longirostris")
  
  logistic_plot_data$species_label <- factor(logistic_plot_data$species, 
                                             levels = 1:3, 
                                             labels = sp.names)
  
  L0_summary_data$species_label <- factor(L0_summary_data$species, 
                                          levels = 1:3, 
                                          labels = sp.names)
  
  # 7. Prepare observed data and apply labels
  obs_data <- data.frame(
    X = jags.data$X,
    y = jags.data$y,
    species_idx = jags.data$species_idx
  ) %>%
    mutate(species = factor(species_idx, levels = 1:3, labels = sp.names))

   # 8. Generate the plot
  p.logistic <- ggplot() +
    # Add the 95% credible interval ribbon
    geom_ribbon(data = logistic_plot_data, 
                aes(x = X, ymin = lwr, ymax = upr), 
                fill = "seagreen", alpha = 0.3) +
    
    # Add the median fitted logistic curve
    geom_line(data = logistic_plot_data, 
              aes(x = X, y = fit), 
              color = "darkgreen", linewidth = 1) +
    
    # Add the raw observation points (0 = Prenatal/Immature, 1 = Postnatal/Mature)
    # Using geom_jitter slightly separates overlapping points on the 0 and 1 lines
    geom_jitter(data = obs_data, 
                aes(x = X, y = y), 
                width = 0, height = 0.02, alpha = 0.3, size = 1.5, color = "black") +
    
    # Add a dotted line across the 50% mark to visualize L0
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    
    geom_vline(data = L0_summary_data, 
               aes(xintercept = L0_med), 
               linetype = "dotted", color = "firebrick", linewidth = 0.8) +
    geom_text(data = L0_summary_data, 
              aes(x = L0_med, y = 0.25, label = paste0("p50 = ", round(L0_med, 1))), 
              angle = 90, vjust = -0.7, color = "firebrick", size = 4) +
    # Create a separate panel for each species
    facet_wrap(~species_label, scales = "free_x") +
    
    
    # Formatting and labels
    theme_minimal() +
    labs(
      title = "Logistic Regression of Maturity/Postnatal Status",
      x = xvar,
      y = yvar,
      subtitle = "Solid line represents median fit; shaded region represents 95% credible intervals"
    ) +
    theme(
      strip.text = element_text(face = "italic", size = 12),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(list(plot = p.logistic,
              data = obs_data,
              p_matrix = p_matrix,
              plot.data = logistic_plot_data))
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