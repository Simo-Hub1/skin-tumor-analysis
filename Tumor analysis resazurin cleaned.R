suppressPackageStartupMessages({
  library(ggplot2)   # 3.4.4
  library(dplyr)     # 1.1.3
  library(tidyr)     # 1.3.0
  library(pwr)       # 1.3-0
})

# ============================================================
# 1. INPUT / OUTPUT
# Paths are relative to the project root — open the .Rproj file
# in RStudio before running, or set your working directory manually
# with setwd("path/to/project").
# ============================================================
CSV_PATH   <- "data/all_plates_combined.csv"
output_dir <- "output/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 2. SETTINGS
# ============================================================

# Cell lines — TT (4532M3) excluded from the paper figure
CELL_ORDER       <- c("4532F1", "3847F1")
CELL_LABELS      <- c("4532F1 (T/WT)", "3847F1 (WT)")

# Raw condition names as they appear in the CSV
COND_ORDER  <- c("UV light media", "TT media", "TWT media", "4OHT", "Fresh media")

# Display names used in plots
COND_LABELS <- c(
  "Apoptotic media",
  "Moderate CIN media",
  "Very low CIN media",
  "Control media",
  "Fresh media"
)

COND_COLORS <- c(
  "Apoptotic media"    = "#4600CC",
  "Moderate CIN media" = "goldenrod3",
  "Very low CIN media" = "grey46",
  "Control media"      = "grey20",
  "Fresh media"        = "springgreen4"
)

# Reference (own) media per cell line — comparisons are always vs Control media
OWN_MEDIA <- c(
  "4532F1 (T/WT)" = "Very low CIN media",
  "3847F1 (WT)"   = "Control media"
)

# Outliers to exclude (Cell_type uses original CSV names)
OUTLIERS <- data.frame(
  Cell_type = c("4532F1",        "3847F1"),
  Condition = c("UV light media", "Fresh media"),
  Replicate = c(1L,               4L),
  stringsAsFactors = FALSE
)

ERROR_TYPE <- "sem"   # "sem" or "sd"

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Load CSV, recode cell type names, filter to relevant days and cell lines
load_and_prep <- function(day_filter = NULL) {
  df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
  df <- df[df$Cell_type != "SKIP", ]
  df$Cell_type <- dplyr::recode(df$Cell_type,
                                "T/WT1" = "4532F1", "T/WT2" = "4532F1",
                                "WT"    = "3847F1", "TT"    = "4532M3")
  df <- df[df$Cell_type %in% CELL_ORDER, ]
  if (!is.null(day_filter)) df <- df[df$Day %in% day_filter, ]
  df$Condition <- factor(df$Condition, levels = COND_ORDER,   labels = COND_LABELS)
  df$Cell_type <- factor(df$Cell_type, levels = CELL_ORDER,   labels = CELL_LABELS)
  df$Day       <- factor(df$Day,       levels = day_filter %||% c("Day 1", "Day 3", "Day 5"))
  df
}

# Remove pre-specified outlier replicates
remove_outliers <- function() {
  df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
  df <- df[df$Cell_type != "SKIP", ]
  df$Cell_type <- dplyr::recode(df$Cell_type,
                                "T/WT1" = "4532F1", "T/WT2" = "4532F1",
                                "WT"    = "3847F1", "TT"    = "4532M3")
  df <- df[df$Cell_type %in% CELL_ORDER, ]
  df$Condition <- factor(df$Condition, levels = COND_ORDER, labels = COND_LABELS)
  df$Cell_type <- factor(df$Cell_type, levels = CELL_ORDER, labels = CELL_LABELS)
  for (i in seq_len(nrow(OUTLIERS))) {
    excl      <- OUTLIERS[i, ]
    cl_label  <- CELL_LABELS[match(excl$Cell_type, CELL_ORDER)]
    cond_label <- COND_LABELS[match(excl$Condition, COND_ORDER)]
    mask <- df$Cell_type == cl_label &
      df$Condition == cond_label &
      df$Replicate == excl$Replicate
    df <- df[!mask, ]
  }
  df
}

# Compute mean, SD, SEM per group
summarise_df <- function(df) {
  df %>%
    dplyr::group_by(Cell_type, Condition, Day) %>%
    dplyr::summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD   = sd(Value,   na.rm = TRUE),
      N    = sum(!is.na(Value)),
      SEM  = SD / sqrt(N),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Err = if (ERROR_TYPE == "sem") SEM else SD)
}

# Compute per-replicate slopes between two days
compute_slopes <- function(day_filter, day_start, day_end) {
  df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
  df <- df[df$Cell_type != "SKIP", ]
  df$Cell_type <- dplyr::recode(df$Cell_type,
                                "T/WT1" = "4532F1", "T/WT2" = "4532F1",
                                "WT"    = "3847F1", "TT"    = "4532M3")
  df <- df[df$Cell_type %in% CELL_ORDER, ]
  
  df_bio <- df %>%
    dplyr::group_by(Cell_type, Condition, Day, Replicate) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  d1 <- df_bio %>% dplyr::filter(Day == day_start) %>%
    dplyr::rename(V1 = Value) %>% dplyr::select(-Day)
  d2 <- df_bio %>% dplyr::filter(Day == day_end) %>%
    dplyr::rename(V2 = Value) %>% dplyr::select(-Day)
  
  n_days <- as.numeric(gsub("Day ", "", day_end)) -
    as.numeric(gsub("Day ", "", day_start))
  
  d1 %>%
    dplyr::inner_join(d2, by = c("Cell_type", "Condition", "Replicate")) %>%
    dplyr::mutate(
      Slope     = (V2 - V1) / n_days,
      Condition = factor(Condition, levels = COND_ORDER,  labels = COND_LABELS),
      Cell_type = factor(Cell_type, levels = CELL_ORDER,  labels = CELL_LABELS)
    )
}

# Cohen's d (pooled SD)
cohens_d_fn <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  ps <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (ps == 0) return(NA_real_)
  abs((mean(x) - mean(y)) / ps)
}

# Cohen's d magnitude label
d_mag <- function(d) {
  ad <- abs(d)
  if (is.na(ad))   return("nd")
  if (ad >= 0.8)   return("large")
  if (ad >= 0.5)   return("medium")
  if (ad >= 0.2)   return("small")
  return("negligible")
}

# Significance label
sig_fn <- function(p) {
  if (is.na(p))    return("nd")
  if (p < 0.001)   return("***")
  if (p < 0.01)    return("**")
  if (p < 0.05)    return("*")
  if (p < 0.10)    return(sprintf("p=%.3f", p))
  return("ns")
}

# ============================================================
# 4. MAIN ANALYSIS FUNCTION
#    Produces the annotated slope plot used in the paper
# ============================================================
run_stats_and_plots <- function(slopes_rep, suffix) {
  
  # Summary stats per condition per cell line
  sum_sl <- slopes_rep %>%
    dplyr::group_by(Cell_type, Condition) %>%
    dplyr::summarise(
      Mean_slope = mean(Slope,  na.rm = TRUE),
      SD_slope   = sd(Slope,    na.rm = TRUE),
      N          = dplyr::n(),
      SEM_slope  = SD_slope / sqrt(N),
      .groups    = "drop"
    )
  
  # One-way ANOVA + Tukey HSD + Cohen's d vs own media per cell line
  anova_res  <- list()
  effect_res <- list()
  
  for (ct_label in CELL_LABELS) {
    sub       <- slopes_rep %>% dplyr::filter(Cell_type == ct_label)
    own_cond  <- OWN_MEDIA[ct_label]
    own_vals  <- sub %>% dplyr::filter(Condition == own_cond) %>% dplyr::pull(Slope)
    aov_fit   <- aov(Slope ~ Condition, data = sub)
    aov_s     <- summary(aov_fit)[[1]]
    f_val     <- aov_s[["F value"]][1]
    p_val     <- aov_s[["Pr(>F)"]][1]
    tukey     <- TukeyHSD(aov_fit)$Condition
    
    anova_res[[ct_label]] <- data.frame(
      Cell_type = ct_label, F_value = round(f_val, 3),
      p_ANOVA   = round(p_val, 4), Sig = sig_fn(p_val),
      stringsAsFactors = FALSE
    )
    
    message(sprintf("  %s ANOVA: F=%.3f p=%.4f %s", ct_label, f_val, p_val, sig_fn(p_val)))
    
    for (cond in levels(sub$Condition)[levels(sub$Condition) != own_cond]) {
      ov    <- sub %>% dplyr::filter(Condition == cond) %>% dplyr::pull(Slope)
      d     <- cohens_d_fn(own_vals, ov)
      pair1 <- paste0(own_cond, "-", cond)
      pair2 <- paste0(cond, "-", own_cond)
      tp    <- NA
      if (pair1 %in% rownames(tukey)) tp <- tukey[pair1, "p adj"]
      if (pair2 %in% rownames(tukey)) tp <- tukey[pair2, "p adj"]
      effect_res[[length(effect_res) + 1]] <- data.frame(
        Cell_type  = ct_label, Own_media = own_cond,
        Comparison = as.character(cond),
        Cohens_d   = round(d, 3), Magnitude = d_mag(d),
        Tukey_p    = round(tp, 4), Label = sig_fn(tp),
        stringsAsFactors = FALSE
      )
    }
  }
  
  anova_df  <- do.call(rbind, anova_res)
  effect_df <- do.call(rbind, effect_res)
  
  write.csv(anova_df,  file.path(output_dir, paste0("stats_anova_",    suffix, ".csv")), row.names = FALSE)
  write.csv(effect_df, file.path(output_dir, paste0("stats_cohens_d_", suffix, ".csv")), row.names = FALSE)
  
  # Annotation data for the plot
  annot_df <- effect_df %>%
    dplyr::left_join(
      sum_sl %>%
        dplyr::rename(Comparison = Condition) %>%
        dplyr::mutate(Comparison = as.character(Comparison)),
      by = c("Cell_type", "Comparison")
    ) %>%
    dplyr::mutate(
      Condition = factor(Comparison, levels = COND_LABELS),
      Cell_type = factor(Cell_type,  levels = CELL_LABELS)
    )
  
  y_range <- sum_sl %>%
    dplyr::group_by(Cell_type) %>%
    dplyr::summarise(y_max = max(Mean_slope + SEM_slope, na.rm = TRUE), .groups = "drop")
  
  annot_df <- annot_df %>%
    dplyr::left_join(y_range, by = "Cell_type") %>%
    dplyr::mutate(
      y_star = y_max * 1.05,
      y_d    = y_max * 1.18,
      annot  = paste0(Label, "\nd=", sprintf("%.2f", Cohens_d), " (", Magnitude, ")")
    )
  
  ref_df <- data.frame(
    Cell_type = factor(CELL_LABELS, levels = CELL_LABELS),
    Condition = factor(unname(OWN_MEDIA), levels = COND_LABELS)
  ) %>%
    dplyr::left_join(y_range, by = "Cell_type") %>%
    dplyr::mutate(y_ref = y_max * 1.05)
  
  # Figure: annotated slope plot (used in paper)
  slope_annot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_errorbar(
      data = sum_sl,
      aes(x = Condition, ymin = Mean_slope - SEM_slope,
          ymax = Mean_slope + SEM_slope, colour = Condition),
      width = 0.25, linewidth = 0.8
    ) +
    geom_crossbar(
      data = sum_sl,
      aes(x = Condition, y = Mean_slope, ymin = Mean_slope, ymax = Mean_slope),
      width = 0.28, linewidth = 0.7, colour = "black", fill = NA
    ) +
    geom_jitter(
      data = slopes_rep,
      aes(x = Condition, y = Slope, colour = Condition),
      width = 0.12, size = 3.5, alpha = 0.8, shape = 16
    ) +
    geom_text(
      data = annot_df,
      aes(x = Condition, y = y_star, label = annot),
      size = 2.8, colour = "black", vjust = 0, lineheight = 1.3
    ) +
    geom_text(
      data = ref_df,
      aes(x = Condition, y = y_ref, label = "reference"),
      size = 2.8, colour = "grey50", vjust = 0, fontface = "italic"
    ) +
    facet_wrap(~ Cell_type, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = COND_COLORS) +
    scale_y_continuous(
      labels = function(x) paste0(round(x / 1000), "k"),
      expand = expansion(mult = c(0.05, 0.45))
    ) +
    scale_x_discrete(labels = COND_LABELS) +
    labs(
      title    = paste0("Growth rate + stats \u2014 ", suffix),
      subtitle = "Tukey p vs own media  |  horizontal line = mean  |  error bars = mean \u00b1 SEM",
      x        = NULL,
      y        = "Slope (FI / day)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle     = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      strip.text        = element_text(face = "bold", size = 10),
      strip.background  = element_rect(fill = "grey95", colour = "grey70"),
      legend.position   = "none",
      axis.text.x       = element_text(size = 9, angle = 35, hjust = 1),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
      panel.spacing     = unit(1, "lines")
    )
  
  ggsave(
    file.path(output_dir, paste0("resazurin_slopes_stats_", suffix, ".png")),
    slope_annot, width = 9, height = 6.5, dpi = 150
  )
  ggsave(
    file.path(output_dir, paste0("resazurin_slopes_stats_", suffix, ".pdf")),
    slope_annot, width = 9, height = 6.5
  )
  
  message(sprintf("  Saved: resazurin_slopes_stats_%s.png", suffix))
}

# ============================================================
# 5. RUN ANALYSIS — DAY 1 & DAY 3
# ============================================================
DAY_FILTER <- c("Day 1", "Day 3")
SUFFIX     <- "day1_day3"

if (!file.exists(CSV_PATH)) {
  stop("Input CSV not found: ", CSV_PATH)
}

message("Computing slopes Day 1 -> Day 3 ...")
slopes_d3 <- compute_slopes(DAY_FILTER, "Day 1", "Day 3")
write.csv(slopes_d3,
          file.path(output_dir, "slopes_per_replicate_day1_day3.csv"),
          row.names = FALSE)

message("Running stats and generating plot ...")
run_stats_and_plots(slopes_d3, SUFFIX)

message("Done. Output saved to: ", output_dir)