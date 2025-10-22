## -------------------------------------------------------------------------- ##
# Script associated with Steger et al. (2025):
# "Impact-based early warning of mass movements – A dynamic spatial modelling approach for the Alpine region"
# doi: https://doi.org/10.5194/egusphere-2025-4940
## -------------------------------------------------------------------------- ##

# The R code uses data from the repository to fit and analyze the models and to
# visualize data and results.
# It also contains code snippets on how certain anlayses
# (e.g. basin-based cross validation, variable importance) were performed.
#
# Structure of the script
# - Step 1: Load data and basic visualizations
# - Step 2: Fit the models
# - Step 3: Calculate fitting performance and plot ROC curves
#           (Fig. 4a)
# - Step 4: Cross validation: Snippet code for basin-based cross validation and
#           visualizing final performance (Fig. 4a-d)
# - Step 5: Variable importance plots (incl. snippet for calculating feature importance)
#           (Fig. 5)
# - Step 6: Visualizing partial effects from fitted models
#           (Figs. 6, Fig. 7, Fig. 8)

## -------------------------------------------------------------------------- ##
## Setup ----
## -------------------------------------------------------------------------- ##

logger::log_info("S0 » Running setup and config")

# clean workspace and run garbage collection to free memory
rm(list = ls())
gc()

# Load packages
suppressPackageStartupMessages({
  library("mgcv") # fitting GAMs (bam/gam)
  library("dplyr") # data manipulation
  library("purrr") # functional programming
  library("sf") # spatial vector data (simple features)
  library("ggplot2") # plotting
  library("patchwork") # compose ggplots
  library("tmap") # thematic mapping for sf objects
  library("sperrorest") # spatial error estimation
  library("pROC") # ROC / AUC calculations
  library("vip") # variable importance computation
  library("gratia") # diagnostic & plotting helpers for GAMs
  library("logger") # logging
})

## -------------------------------------------------------------------------- ##
## Config ----
## -------------------------------------------------------------------------- ##

# Colors for plots
cols <- "#56b899" # color slide-type models (SL)
cold <- "#9052b6" # color flows-type  models (DF)
colr <- "#e06c3e" # color fall-type models (RF)

ncores <- 8L

# Custom ggplot theme ====
common_theme <- theme_minimal(base_size = 12) + # base theme for ggplot
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # subtle border for clarity
  )

## -------------------------------------------------------------------------- ##
## Step 1: Load data ----
## -------------------------------------------------------------------------- ##

log_info("S1 » Loading data")
# Load preprocessed input data

basins <- readRDS("data/basins_compressed.rds") # Basins with static data (sf) -- polygons with attributes
AS <- readRDS("data/Alpine_Space_Boundary_4326.rds") # Alpine Space outline (sf)
df_slides <- readRDS("data/df_slides.rds") # Training data for slide-type (SL) model
df_flows <- readRDS("data/df_flows.rds") # Training data for flow-type (DF) model
df_falls <- readRDS("data/df_falls.rds") # Training data for fall-type (RF) model
SLcv <- readRDS("data/mySLcv_5f_10r.rds") # Cross validation results for slide type (SL) model
DFcv <- readRDS("data/myDFcv_5f_10r.rds") # Cross validation results for flow type (DF) model
RFcv <- readRDS("data/myRFcv_5f_10r.rds") # Cross validation results for fall type (RF) model
varimpo_raw <- readRDS("data/varimpo_raw.rds") # Results (raw) from variable importance assessment (raw permute outputs)
varimpoSL <- readRDS("data/varimpoSL.rds") # Plot-ready variable importance data for slide-type (SL)
varimpoDF <- readRDS("data/varimpoDF.rds") # Plot-ready variable importance data for flow-type (DF)
varimpoRF <- readRDS("data/varimpoRF.rds") # Plot-ready variable importance data for fall-type (RF)

# Example: Alpine Space outline map
log_info("-- » Creating Alpine space outline map")
tm_shape(AS) +
  tm_fill(col = "lightgrey") +
  tm_borders(col = "black", lwd = 2) +
  tm_title("Alpine Space Boundary")

# Example plot for basins:
#     - Mean slope angles in the respective PPAs
#     - "Missing" (NA) means no PPA pixel in the entire basin;
#     - Facets show slope maps for slides (SL), flows (DF) and falls (RF), respectively
#     - Basin polygons are filled without outlines; use tm_polygons() for visualizing borders
log_info("-- » Creating example plot for basins")
tm_shape(basins) +
  tm_fill(c("SL_Slope", "DF_Slope", "RF_Slope"),
    fill.scale = tm_scale(values = c("#d7f4be", "#f0d076", "#e73030")),
    fill.legend = tm_legend(title = "Mean slope in PPA"), fill.free = FALSE
  ) +
  tm_facets(ncol = 1)


## -------------------------------------------------------------------------- ##
## Step 2: Fit the models ----
## -------------------------------------------------------------------------- ##

## Config ====

# k-factor to restrict maximum flexibility of smooth terms (keeps smooths conservative/generalized)
maxk <- 4

## -------------------------------------------------------------------------- ##
## Helper functions ====
## -------------------------------------------------------------------------- ##

#' Fit a Generalized Additive Mixed Model (GAMM)
#'
#' This function fits a Generalized Additive Mixed Model (GAMM) using the `mgcv::bam` function.
#' It provides options to print a summary of the model and generate exploratory partial effect plots.
#'
#' @param formula A formula specifying the model structure.
#' @param data A data frame containing the variables used in the model.
#' @param summary Logical; if `TRUE` (default), prints a summary of the fitted model.
#' @param plot Logical; if `TRUE` (default), generates quick exploratory partial effect plots.
#'   This is for exploratory purposes only. For publication-quality plots, consider using `gratia::draw()`.
#'
#' @return The fitted GAMM model object returned by `mgcv::bam`.
#'   If `summary = TRUE`, the model summary is printed to the console.
#'   If `plot = TRUE`, exploratory partial effect plots are displayed.
#'
#' @details
#' The function uses the `mgcv::bam` function to fit a GAMM with the following settings:
#' - Family: Binomial
#' - Discrete: 100 (discretize covariates for storage and efficiency reasons)
#' - Method: "fREML" (fast REML)
#' - Select: `TRUE` (perform model selection by adding selection penalties to smooth effects)
fit_gamm <- function(formula, data, summary = TRUE, plot = TRUE, ...) {
  bam_model <- mgcv::bam(formula, data = data, family = binomial, discrete = 100, method = "fREML", select = TRUE, ...)
  if (summary) summary(bam_model)
  if (plot) plot(bam_model, pages = 1)
  return(bam_model)
}

## -------------------------------------------------------------------------- ##

# Slide-type (SL) model ====
log_info("S2 » Fitting slide type model")

# Define formula for SL ####
fo_slides <-
  SL01 ~ # Binary response variable for slide-type (SL) -- dependent variable in training data
  s(tp_2, k = maxk) + # Short-term precipitation (smoothed)
  s(tp_30, k = maxk) + # Antecedent precipitation: 30-day (smoothed)
  s(tp_su_mean_annual_mean, k = maxk) + # Mean annual precipitation (smoothed)
  s(SL_Slope, k = maxk) + # Slope angle in slide-type PPA (smoothed)
  SL_dominant_landcover + # Dominant land cover class in slide-type PPA (parametric categorical)
  SL_dominant_lithology + # Dominant lithology class in slide-type PPA (parametric categorical)
  s(log_B_SL_n, k = maxk) + # Log Number of buildings in slide-type PPA (smoothed)
  s(SL_roads_12, k = maxk) + # Portion of “road” cells in slide PPA (smoothed)
  s(doy, bs = "cc", k = maxk) + # Day-of-Year as circular variable (cyclic spline)
  s(temp_mean, k = maxk) + # Mean daily temperature (smoothed)
  temp_sub0 + # Binary: Day below 0 °C (yes/no) (parametric)
  s(cat, bs = "re") + # Random effect: Sampling location (Basin-ID) using bs="re"
  s(year, bs = "re") # Random effect: Sampling Year using bs="re"

# Fit the slide type model ####
fit_slides <- fit_gamm(formula = fo_slides, data = df_slides, nthreads = ncores)

## -------------------------------------------------------------------------- ##

# Flow-type (DF) model ====
log_info("-- Fitting flow type model")

# Define formula for DF ####
fo_flows <-
  DF01 ~ # Response for flow-type (DF) (binary)
  s(tp_2, k = maxk) + # Short term precipitation (smoothed)
  s(tp_21, k = maxk) + # Antecedent precipitation (21-day) (smoothed)
  s(tp_su_mean_annual_mean, k = maxk) + # Mean annual precipitation (smoothed)
  s(DF_Slope, k = maxk) + # Slope angle in flow-type PPA (smoothed)
  s(DF_CI_mean, k = maxk) + # Convergence index in flow-type PPA (smoothed)
  DF_dominant_landcover + # Dominant land cover class in flow-type PPA (parametric categorical)
  DF_dominant_lithology + # Dominant lithology class in flow-type PPA (parametric categorical)
  s(log_B_DF_n, k = maxk) + # Log Number of buildings in flow-type PPA (smoothed)
  s(doy, bs = "cc", k = maxk) + # Day-of-Year as circular variable (cyclic spline)
  s(temp_mean, k = maxk) + # Mean daily temperature (smoothed)
  temp_sub0 + # Binary: Day below 0 °C (yes/no) (parametric)
  s(cat, bs = "re") + # Random effect: Sampling location (Basin-ID) using bs="re"
  s(year, bs = "re") # Random effect: Sampling Year using bs="re"

# Fit the flow type model ####
fit_flows <- fit_gamm(formula = fo_flows, data = df_flows, nthreads = ncores)

## -------------------------------------------------------------------------- ##

# Fall-type (RF) model ====

# Define formula for (RF) ####
log_info("-- » Fitting fall-type model")
fo_falls <-
  RF01 ~ # Response for fall-type (RF) (binary)
  s(tp_2, k = maxk) + # Short term precipitation (smoothed)
  s(tp_14, k = maxk) + # Antecedent precipitation (14-day) (smoothed)
  s(tp_su_mean_annual_mean, k = maxk) + # Mean annual precipitation (smoothed)
  s(RF_Slope, k = 3) + # Slope angle in PPA; k=3 chosen for less wiggliness
  s(RF_CI_mean, k = maxk) + # Convergence index in fall-type PPA (smoothed)
  s(RF_bare.other, k = maxk) + # Portion of bare surface in fall-type PPA (smoothed)
  RF_dominant_landcover + # Dominant land cover class in fall-type PPA (parametric categorical)
  RF_dominant_lithology + # Dominant lithology class in fall-type PPA (parametric categorical)
  s(log_B_RF_n, k = maxk) + # Log Number of buildings in fall-type PPA (smoothed)
  s(RF_roads_12, k = maxk) + # Portion of “road” cells in fall PPA (smoothed)
  temp_mixed + # Binary: Day crossing 0 °C (yes/no) (parametric)
  s(cat, bs = "re") + # Random effect: Sampling location (Basin-ID) using bs="re"
  s(year, bs = "re") # Random effect: Sampling Year using bs="re"

# Fit the fall type model ####
fit_falls <- fit_gamm(formula = fo_falls, data = df_falls, nthreads = ncores)


## -------------------------------------------------------------------------- ##
## Step 3: Calculate fitting performance and plot ROC curves ----
##         (Fig. 4a in publication)
## -------------------------------------------------------------------------- ##

log_info("S3 » Assessing model performance")

## -------------------------------------------------------------------------- ##
## Helper functions ====
## -------------------------------------------------------------------------- ##

#' Performance evaluation for GAMM ####
#'
#' This function evaluates the performance of a predictive model by computing
#' the predicted probabilities, generating a Receiver Operating Characteristic
#' (ROC) curve, calculating the Area Under the Curve (AUC), and identifying the
#' best threshold for classification.
#'
#' @param model A fitted model object. Only tested for `mgcv::bam()`.
#' @param data A data frame containing the data used for prediction.
#'   This should include the response variable and any predictors required by the model.
#' @param response_col A string specifying the name of the column in `data` that
#'   contains the binary response variable (e.g., 0/1 or TRUE/FALSE).
#' @param prob_col A string specifying the name of the column where the
#'   predicted probabilities will be stored in `data`.
#' @param exclude_terms A character vector specifying terms to average out during prediction.
#'   Uses the random random effects (`bs = "re"`) `s(cat)` and `s(year)` as default.
#'
#' @return A list containing the following elements:
#' - roc: A `pROC::roc` object representing the ROC curve for the model's predictions.
#' - auc: A numeric value representing the Area Under the Curve (AUC) of the ROC curve.
#' - best_threshold: A data frame containing the best threshold for classification, along with its sensitivity and specificity.
evaluate_performance <- function(model, data, response_col, prob_col, exclude_terms = c("s(cat)", "s(year)")) {
  # Predict probabilities for training data w/ random effects excluded
  data[[prob_col]] <- predict(model, type = "response", newdata = data, exclude = exclude_terms)

  # Compute ROC and AUC
  roc_obj <- roc(response = data[[response_col]], predictor = data[[prob_col]], auc = TRUE, ci = TRUE)
  auc_value <- roc_obj$auc

  # Identify "best" threshold (c method = closest to (0,1))
  best_threshold <- pROC::coords(roc_obj, x = "best", input = "threshold", best.method = "c")

  # Return results as a list
  return(list(roc = roc_obj, auc = auc_value, best_threshold = best_threshold))
}

## -------------------------------------------------------------------------- ##

#' Prepare ROC Data for Plotting ####
#'
#' This function converts a `pROC::roc` object into a data frame suitable for plotting,
#' including the False Positive Rate (FPR), True Positive Rate (TPR), and a process type label.
#'
#' @param roc_obj A `pROC::roc` object representing the Receiver Operating Characteristic (ROC) curve.
#'   This object is typically generated using the `pROC::roc` function.
#' @param process_name A string specifying the name of the process type associated with the ROC curve.
#'
#' @return A data frame with the following columns:
#' - FPR: False Positive Rate (fall-out, 1 - specificity.
#' - TPR: True Positive Rate (sensitivity).
#' - Process: A string indicating the process type.
prepare_roc_data <- function(roc_obj, process_name) {
  data.frame(
    FPR = 1 - roc_obj$specificities, # False positive rate
    TPR = roc_obj$sensitivities, # True positive rate
    Process = process_name
  )
}

## -------------------------------------------------------------------------- ##

#' Extract Cutpoints for ROC Analysis ####
#'
#' This function creates a data frame containing the False Positive Rate (FPR),
#' True Positive Rate (TPR), and process type name based on the best threshold
#' identified in a ROC analysis.
#'
#' @param best_threshold A data frame or list containing the best threshold information,
#'   typically obtained using the `pROC::coords` function. It should include `specificity` and `sensitivity` values.
#' @param process_name A string specifying the name of the process associated with the ROC curve.
#'
#' @return A data frame with the following columns:
#' - FPR: False Positive Rate (fall-out, 1 - specificity.
#' - TPR: True Positive Rate (sensitivity).
#' - Process: A string indicating the process type.
extract_cutpoints <- function(best_threshold, process_name) {
  data.frame(
    FPR = 1 - best_threshold$specificity, # FPR at best thresholds
    TPR = best_threshold$sensitivity, # TPR at best thresholds
    Process = process_name
  )
}

## -------------------------------------------------------------------------- ##

# Define a list of processes with their respective parameters ====
processes <- tibble(
  model = list(fit_slides, fit_flows, fit_falls),
  data = list(df_slides, df_flows, df_falls),
  response_col = c("SL01", "DF01", "RF01"),
  prob_col = c("probSL", "probDF", "probRF"),
  process_name = c("Slide-type", "Flow-type", "Fall-type")
)

# Map fit_performance across all three processes ====
results <- processes |>
  mutate(
    fit_result = pmap(
      list(model, data, response_col, prob_col),
      ~ evaluate_performance(..1, ..2, ..3, ..4)
    ),
    roc_data = map2(fit_result, process_name, ~ prepare_roc_data(.x$roc, .y)),
    cutpoint_data = map2(fit_result, process_name, ~ extract_cutpoints(.x$best_threshold, .y))
  )

# Extract AUC values and labels ====
auc_labels <- results |>
  mutate(
    Process = process_name,
    AUC = map_dbl(fit_result, ~ round(.x$auc, 3)),
    x = 0.4,
    y = c(0.15, 0.11, 0.07),
    .keep = "none"
  )

log_info("-- » RESULT: AUC per process type:")
auc_vals <- setNames(
  paste0(auc_labels$Process, " (AUC = ", round(auc_labels$AUC, 2), ")"),
  auc_labels$Process
) |>
  print()

# Prepare ROC data and cutpoints for plotting ====
roc_all <- bind_rows(results$roc_data) |>
  mutate(Process = factor(Process, levels = auc_labels$Process))

log_info("-- » RESULT: Optimal cutpoints per process type:")
cutpoints <- bind_rows(results$cutpoint_data) |>
  mutate(Process = factor(Process, levels = auc_labels$Process)) |>
  print()

# Build ROC plot (Fig 4a) ====
p_roc <- ggplot(roc_all, aes(x = FPR, y = TPR, color = Process)) +
  geom_line(linewidth = 1) + # ROC curves
  geom_abline(linetype = "dashed", color = "gray60") + # no-skill line
  geom_point(data = cutpoints, aes(x = FPR, y = TPR, color = Process), size = 2) + # best-threshold points
  scale_color_manual(
    values = c("Slide-type" = cols, "Flow-type" = cold, "Fall-type" = colr),
    labels = auc_vals,
    breaks = c("Slide-type", "Flow-type", "Fall-type")
  ) +
  labs(
    title = "Fitting performance",
    x = "False Positive Rate (FPR)",
    y = "True Positive Rate (TPR)",
    color = NULL
  ) + # removes legend title (we use labels instead)
  coord_equal() + # equal aspect ratio to avoid distortion
  common_theme +
  theme(
    legend.position = c(0.95, 0.05), # place legend inside plot at bottom-right
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA), # translucent bg for legend
    legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12)
  )
log_info("-- » Creating ROC plot")
p_roc


## -------------------------------------------------------------------------- ##
## Step 4: Cross validation ----
##         Snippet code for basin-based CV and visualizing final performance
##         (Fig. 4a-d in publication)
## -------------------------------------------------------------------------- ##

log_info("S4 » Performing cross validation")

# This section
#    - describes how basin-based cross-validation was performed and
#    - provides data and plots.
#
# Description:
# - The function "sperrorest::partition_factor_cv()" was used to create
#   train/test splits based on basin IDs ("cat").
# - Five folds and ten repetitions were applied.
# - The process was repeated for each model (slide, flow, fall) and for multiple
#   antecedent precipitation time windows (tp_null, tp_7, tp_14, tp_21, tp_30).
#
# Example (commented) usage:

DONTRUN <- TRUE

# styler: off
if (!DONTRUN) {

  nreps <- 10                                   # Number of repetitions (10x repeated CV)
  nfolds <- 5                                   # Number of folds (5-fold CV)
  parti <- partition_factor_cv(                 # Function to generate basin-based cross-validation splits
    df_slides,                                  # Input data frame (here slide-type data)
    nfold = nfolds,                             # Number of folds (5-fold CV)
    repetition = nreps,                         # Number of repetitions (10x repeated CV)
    seed1 = 666,                                # Random seed for reproducibility
    fac = "cat"                                 # Factor variable ("cat") defining basin IDs
  )                                             # Creates a list of partition objects for repeated 5-fold CV

  # Example iteration (commented snippet):
  results_df_slides <- data.frame(
    # Create an empty df to store results
    Formula = character(),
    Repetition = integer(),
    Fold = integer(),
    AUROC = numeric(),
    stringsAsFactors = FALSE
  )
  for (formula_name in names(formula_list)) {   # Iterate over each formula (formula_list contains gam-formulas (cf. Step 2) with varying antec. precip. time windows)
    fo <- formula_list[[formula_name]]          # Extract current formula by name
    print(fo)                                   # Print formula to console for tracking progress
    # Loop through repetitions and folds
    for (j in 1:nreps) {                        # Outer loop over number of repetitions
      partiloop <- parti[[j]]                   # Access the j-th partition (one repetition)
      for (i in 1:nfolds) {                     # Inner loop over folds within repetition
        first <- partiloop[[i]][[2]]            # Extract test indices for the i-th fold
        test <- df_slides[first, ]
        ntesti <- nrow(test)                    # Create test set and count its rows
        train <- df_slides[-first, ]            # Define training data by excluding test indices
        # Fit model and calculate AUROC
        myfit <- mgcv::bam(                     # Fit GAM model using mgcv::bam (fast REML method)
          fo,                                   # Model formula
          data = train,                         # Training data subset
          family = binomial,                    # Binary family for presence/absence
          method = "fREML",                     # Fitting method: fast REML
          discrete = 100                        # Use discrete approximation (speed optimization)
        )
        test$prob <- predict.gam(               # Predict probabilities on the test set
          myfit,                                # Fitted GAM model
          type = "response",                    # Output predicted probabilities (0–1)
          newdata = test,                       # Data for prediction
          exclude = c("s(cat)", "s(year)")      # Exclude random effects (as defined above)
        )
        p_roc <- roc(                           # Compute ROC using pROC::roc
          response = test$SL01,                 # Binary response
          predictor = test$prob,                # Predicted probabilities
          auc = TRUE                            # Return AUC value
        )
        auroc <- round(p_roc$auc, 5)            # Round AUC
        print(auroc)                            # Print current AUC to console
        results_df_slides <- rbind(             # Append results to a data frame
          results_df_slides,                    # Existing results table
          data.frame(                           # Add new row with results for current fold
            Formula = formula_name,             # Formula identifier
            Repetition = j,                     # Current repetition number
            Fold = i,                           # Current fold number
            AUROC = auroc                       # Calculated AUC value
          )
        )
      }
    }
  }
  print(results_df_slides)                      # Display full table of AUROC results
  mySLcv_5f_10r <- results_df_slides            # Save cross-validation results object for slides (cf. data loaded in Step 1)
}
# styler: on

# End of cross-validation example.
# Repeat the same procedure for flows (DF) and falls (RF)
# Below: summarizing final results for all process types.

## -------------------------------------------------------------------------- ##

# extract stats from cross-validation outputs (slide / flow / fall) ====
log_info("-- » Assessing CV results")

#' Summarize Cross-Validation Results
#'
#' This function calculates summary statistics (mean, median, min, max, and IQR)
#' for the Area Under the Receiver Operating Characteristic Curve (AUROC)
#' from cross-validation results, grouped by formula.
#'
#' @param data A data frame containing cross-validation results.
#'   It must include a `Formula` column (grouping variable) and
#'   an `AUROC` column (numeric values to summarize).
#'
#' @return A summarized data frame with the following columns:
#' - Formula: The formula used in the model (grouping variable).
#' - Mean_AUROC: The mean AUROC value for each formula.
#' - Median_AUROC: The median AUROC value for each formula.
#' - Min_AUROC: The minimum AUROC value for each formula.
#' - Max_AUROC: The maximum AUROC value for each formula.
#' - IQR_AUROC: The interquartile range (IQR) of AUROC values for each formula.
convert_formula <- function(data, levels, labels) {
  # Convert Formula to factor with custom labels for consistent x-axis labels
  data |>
    mutate(Formula = factor(Formula, levels = levels, labels = labels))
}
summarize_cv <- function(data, custom_order, custom_labels) {
  data |>
    convert_formula(levels = custom_order, labels = custom_labels) |>
    group_by(Formula) |>
    summarise(
      Mean_AUROC = mean(AUROC),
      Median_AUROC = median(AUROC),
      Min_AUROC = min(AUROC),
      Max_AUROC = max(AUROC),
      IQR_AUROC = IQR(AUROC),
      .groups = "drop" # Avoid grouped output
    )
}

# custom ordering / labels for antecedent precipitation windows ====
custom_order <- c("fo_null", "fo_7", "fo_14", "fo_21", "fo_30")
custom_labels <- c("Null", "7 days", "14 days", "21 days", "30 days")

# Summarize cross-validation results for slides, flows, and falls ====
log_info("-- » RESULT: CV AUROC for slides, flows and falls:")
summarize_cv(SLcv, custom_order, custom_labels)
summarize_cv(DFcv, custom_order, custom_labels)
summarize_cv(RFcv, custom_order, custom_labels)

# Compute median AUROC per group ====

SLcv <- convert_formula(SLcv, levels = custom_order, labels = custom_labels)
SLcv_medians <- SLcv |>
  group_by(Formula) |>
  summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianSL <- max(SLcv_medians$median_AUROC) # used for title annotation

DFcv <- convert_formula(DFcv, levels = custom_order, labels = custom_labels)
DFcv_medians <- DFcv |>
  group_by(Formula) |>
  summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianDF <- max(DFcv_medians$median_AUROC)

RFcv <- convert_formula(RFcv, levels = custom_order, labels = custom_labels)
RFcv_medians <- RFcv |>
  group_by(Formula) |>
  summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianRF <- max(RFcv_medians$median_AUROC)

# IQR metrics used in text or figure annotations ====
log_info("-- » RESULT: IQR for SL30, DF21 and RF14:")
(iqr_SL_30 <- IQR(SLcv$AUROC[SLcv$Formula == "30 days"], na.rm = TRUE)) # Slide-type IQR for 30d formula
(iqr_DF_21 <- IQR(DFcv$AUROC[DFcv$Formula == "21 days"], na.rm = TRUE)) # Flow-type IQR for 21d
(iqr_RF_14 <- IQR(RFcv$AUROC[RFcv$Formula == "14 days"], na.rm = TRUE)) # Fall-type IQR for 14d

## -------------------------------------------------------------------------- ##
# Create boxplots for each process type with unified theme (Fig. 4) ====

# SL
p1 <- ggplot(SLcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = cols) +
  labs(
    title = paste0("Slide-type: Predictive performance\nBest: 30 days (median: ", round(max_medianSL, 2), ")"),
    x = "Antecedent precipitation time window",
    y = "AUROC"
  ) +
  coord_cartesian(ylim = c(0.725, 0.925)) + # y-limits tuned for visibility across plots
  common_theme

# DF
p2 <- ggplot(DFcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = cold) +
  labs(
    title = paste0("Flow-type: Predictive performance\nBest: 21 days (median: ", round(max_medianDF, 2), ")"),
    x = "Antecedent precipitation time window",
    y = "AUROC"
  ) +
  coord_cartesian(ylim = c(0.725, 0.925)) +
  common_theme

# RF
p3 <- ggplot(RFcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = colr) +
  labs(
    title = paste0("Fall-type: Predictive performance\nBest: 14 days (median: ", round(max_medianRF, 2), ")"),
    x = "Antecedent precipitation time window",
    y = "AUROC"
  ) +
  coord_cartesian(ylim = c(0.725, 0.925)) +
  common_theme

# Plot composition for fitting and predictive performance ====

# Define margin (top, right, bottom, left) for consistent spacing
mymargin <- theme(plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm"))

# Apply margin to ROC plot and CV plots
p_roc <- p_roc + mymargin
p1 <- p1 + mymargin
p2 <- p2 + mymargin
p3 <- p3 + mymargin

# Add tag letters for figure panels
p_roc_tagged <- p_roc + labs(tag = "a)")
p1_tagged <- p1 + labs(tag = "b)")
p2_tagged <- p2 + labs(tag = "c)")
p3_tagged <- p3 + labs(tag = "d)")

# Combine plots and keep tags using patchwork
# 2x2 layout: ROC + p1 (top row), p2 + p3 (bottom row)
combined_plot <- (p_roc_tagged + p1_tagged) / (p2_tagged + p3_tagged) +
  # keep tags on all panels
  plot_layout(tag_level = "keep") &
  # style tags
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1))

log_info("-- » Creating composite performance plot (Fig. 4)")
print(combined_plot)


## -------------------------------------------------------------------------- ##
## Step 5: Variable importance plots ----
##         (incl. snippet code for calculating feature importance)
##         (Fig. 5)
## -------------------------------------------------------------------------- ##

log_info("S5 » Assessing variable importance")

# Snippet code for slide-type models (SL) -- permutation importance via vip::vi_permute (for fast demo nsim small)
mynsim <- 5 # set the number of permutations per variable (set to 100 in original publication -> slow)

set.seed(1) # reproducible permutations
# Select the relevant columns from fitted model
cn_s <- colnames(fit_slides$model)
train_sel_s <- df_slides |>
  dplyr::select(any_of(cn_s)) |> # select columns used in model
  dplyr::select(-SL01) # drop response from training set used for vi_permute
# Define the target variable
target_s <- as.factor(df_slides$SL01) # binary response

# Calculate variable importance using permutation (vi_permute returns Importance & StDev)
result_s <- vi_permute(
  nsim = mynsim, # number of permutations (small here for speed)
  object = fit_slides, # fitted mgcv::bam object
  train = train_sel_s, # training predictors
  target = target_s, # target vector
  event_level = "second", # which factor level is event (depending on factor coding)
  metric = "roc_auc", # metric to evaluate (AUROC)
  pred_wrapper = predict.bam # use predict.bam as wrapper for consistency
) |>
  dplyr::arrange(-Importance) |> # sort by decreasing importance
  dplyr::mutate(process = "Slide-type") # annotate process type

# Optionally drop random effect variables from importance table
result_s <- result_s |>
  dplyr::filter(!(Variable %in% c("cat", "year"))) # remove random effect vars since permuting them is not meaningful
# repeat for flow-types and fall-types and create "varimpo" data (increase mynsim beforehand)

# Print raw and processed importance outputs ====
log_info("-- » RESULT: raw VI dataframe")
print(varimpo_raw, n = nrow(varimpo_raw))

# sorted and renamed importance results for all three process models
log_info("-- » RESULT: VI of slide-type model")
print(varimpoSL)

log_info("-- » RESULT: VI of flow-type model")
print(varimpoDF)

log_info("-- » RESULT: VI of fall-type model")
print(varimpoRF)

# create plots for variable importance w/ ggplot2 (Fig. 5) ====
p_vi_1 <- ggplot(varimpoSL, aes(x = Importance, y = Variable)) +
  geom_point(color = cols, size = 3) + # point for importance
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev),
    width = 0.3, color = cols, size = 1
  ) + # horizontal error bars represent permute SD
  theme_minimal(base_size = 14) + # readable baseline font size
  labs(title = "Slide-type", x = "Importance", y = "Variable") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

p_vi_2 <- ggplot(varimpoDF, aes(x = Importance, y = Variable)) +
  geom_point(color = cold, size = 3) +
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev),
    width = 0.3, color = cold, size = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Flow-type", x = "Importance", y = "Variable") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

p_vi_3 <- ggplot(varimpoRF, aes(x = Importance, y = Variable)) +
  geom_point(color = colr, size = 3) +
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev),
    width = 0.3, color = colr, size = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Fall-type", x = "Importance", y = "Variable") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

# Add tags for subfigures to variable importance plots ====
p_vi_1_tagged <- p_vi_1 + labs(tag = "a)")
p_vi_2_tagged <- p_vi_2 + labs(tag = "b)")
p_vi_3_tagged <- p_vi_3 + labs(tag = "c)")

# Arrange plots vertically, keeping tags (patchwork) ====
combined_plot <- (p_vi_1_tagged) / (p_vi_2_tagged) / (p_vi_3_tagged) +
  plot_layout(tag_level = "keep") &
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1))

log_info("-- » Creating variable importance plot (Fig. 5)")
print(combined_plot)


## -------------------------------------------------------------------------- ##
## Step 6: Visualizing partial effects from fitted models ----
##         (Fig. 6, Fig. 7, Fig. 8)
## -------------------------------------------------------------------------- ##

# Partial effects created separately for SL, DF, RF using gratia::draw()
log_info("S6 » Assessing partial effects")

# Partial effect plot for slide-type from fitted model `fit_slides` (Fig. 6) ####
log_info("-- » Estimating partial effects for slide-type model")

# get smooths
mys <- smooths(fit_slides)
# indices of smooths to extract (custom selection, first eight elements)
myselect <- mys[1:8]

# Define custom labels for each plot
custom_titles <- c(
  "Short-term precipitation (mm)",
  "Antecedent precipitation (mm) — 30 days",
  "Mean annual precipitation (mm)",
  "Slope angle (°)",
  "Number of buildings (log scale)",
  "Transport infrastructure (coverage %)",
  "Day-of-Year",
  "Mean daily temperature (°C)",
  "Land cover",
  "Lithology",
  "Day below 0 °C"
)

# Extract individual plots instead of combining into patchwork directly (draw() from gratia)
plots <- draw(fit_slides,
  data = df_slides, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = cols, smooth_col = "#070000",
  return = "list"
) # return list of ggplot objects for fine-grained customization

# Post-process plots: apply consistent styling and handle discrete variables ====
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) + # set title & x-label
    theme_bw(base_size = 10) + # classical background for publication
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(), # titles handled separately in arrangement
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  # Manually relabel the ticks for the 11th plot (discrete variable) if present
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  } # nicer labels for binary factor
  return(p)
})

tags <- letters[1:length(plots)] # generates sequence "a", "b", "c", ...
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) + # tag inserted for panel lettering
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1), # style tag (letter)
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  }
  return(p)
})

# Combine plots with tags preserved (wrap_plots from patchwork)
p_SL <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = "keep") # arrange into grid, 4 columns

log_info("-- » Creating partial effects plot for slide-type model (Fig. 6)")
print(p_SL)

## -------------------------------------------------------------------------- ##

# Partial effect plot for flow-type from fitted model `fit_flows` (Fig. 7)
log_info("-- » Estimating partial effects for slide-type model")

# get smooths
mys <- smooths(fit_flows)
# model diagnostics reminder
myselect <- mys[1:8]

# Define custom labels for each plot for flows
custom_titles <- c(
  "Short-term precipitation (mm)", "Antecedent precipitation (mm) — 21 days", "Mean annual precipitation (mm)", "Slope angle (°)", "Convergence index", "Number of buildings (log scale)",
  "Day-of-Year", "Mean daily temperature (°C)", "Land cover", "Lithology", "Day below 0 °C"
)

# Extract individual plots for flow model using gratia::draw()
plots <- draw(fit_flows,
  data = df_flows, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = cold, smooth_col = "#070000",
  return = "list"
)

# Post-process flow plots: styling and discrete label handling
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  # discrete variable relabeling if necessary
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  }
  return(p)
})

# Tag letters for flow-type plots and styling
tags <- letters[1:length(plots)]
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  }
  return(p)
})

# Combine and keep tags for flow plots
p_DF <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = "keep") # 4 columns layout

log_info("-- » Creating partial effects plot for flow-type model (Fig. 7)")
print(p_DF)

## -------------------------------------------------------------------------- ##

# Partial effect plot for fall-type from fitted model `fit_falls` (Fig. 8)
log_info("-- » Estimating partial effects for fall-type model")

# get smooths
mys <- smooths(fit_falls)
# selection indices for fall model
myselect <- mys[1:8]

# Define custom labels for fall plots
custom_titles <- c(
  "Short-term precipitation (mm)",
  "Antecedent precipitation (mm) — 14 days",
  "Mean annual precipitation (mm)",
  "Slope angle (°)",
  "Convergence index",
  "Bare surface (coverage %)",
  "Number of buildings (log scale)",
  "Transport infrastructure (coverage %)",
  "Land cover",
  "Lithology",
  "Day crossing 0 °C"
)

# Extract partial plots for fall-type model
plots <- draw(fit_falls,
  data = df_falls, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = colr, smooth_col = "#070000",
  return = "list"
)

# Style and discrete handling for fall plots
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  # Manually relabel the ticks for the 11th plot (discrete variable)
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  }
  return(p)
})

# Tag letters for Fall-type plots
tags <- letters[1:length(plots)]
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100)
    )
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes"))
  }
  return(p)
})

# Combine and keep tags for fall-type partials (4 columns grid)
p_RF <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = "keep")

log_info("-- » Creating partial effects plot for slide-type model (Fig. 8)")
print(p_RF)

## End of script ------------------------------------------------------------ ##
