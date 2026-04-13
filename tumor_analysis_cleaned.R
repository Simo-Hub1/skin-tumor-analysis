suppressPackageStartupMessages({
  library(stringr)      # 1.5.0
  library(ggplot2)      # 3.4.4
  library(grid)         # 4.3.1
  library(ggbeeswarm)   # 0.7.2
  library(openxlsx)     # 4.2.5.2
  library(readr)        # 2.1.4
  library(tidyr)        # 1.3.0
  library(dplyr)        # 1.1.3
  library(emmeans)      # 1.8.9
  library(ggh4x)        # 0.2.6
  library(ggpubr)       # 0.6.0
  library(rstatix)      # 0.7.2
  library(ggsignif)     # 0.6.4
  library(MASS)         # 7.3-60
  library(dunn.test)    # 1.3.5
  library(nlme)         # 3.1-163
})

# Prevent masking
select    <- dplyr::select
filter    <- dplyr::filter
summarise <- dplyr::summarise
mutate    <- dplyr::mutate
rename    <- dplyr::rename

# ============================================================
# 1. INPUT / OUTPUT
# Paths are relative to the project root — open the .Rproj file
# in RStudio before running, or set your working directory manually
# with setwd("path/to/project").
# ============================================================
input_file <- "data/tumor_sizes_raw.csv.xlsx"
output_dir <- "output/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 2. HELPERS / SETTINGS
# ============================================================
genotype_levels <- c("WT/WT", "T/WT", "D/WT", "T/T", "T/D", "D/D")
size_levels     <- c("<10", "10-100", "100-500", ">500")

genotype_palette <- c(
  "WT/WT" = "grey4",
  "T/WT"  = "grey46",
  "D/WT"  = "grey78",
  "T/T"   = "goldenrod3",
  "T/D"   = "darkorange3",
  "D/D"   = "red2"
)

size_palette <- c(
  "<10"     = "steelblue2",
  "10-100"  = "dodgerblue4",
  "100-500" = "deeppink3",
  ">500"    = "deeppink4"
)

# Kruskal-Wallis + Dunn post-hoc vs WT/WT reference group
perform_kw_dunn <- function(data, value_col, group_col, p_thresh = 0.05, ref_group = "WT/WT") {
  kw_test <- kruskal.test(as.formula(paste(value_col, "~", group_col)), data = data)
  if (kw_test$p.value < p_thresh) {
    dunn_res <- rstatix::dunn_test(
      data = data,
      formula = as.formula(paste(value_col, "~", group_col)),
      p.adjust.method = "BH"
    ) %>%
      dplyr::filter(group1 == ref_group | group2 == ref_group) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(g1 = ref_group, g2 = ifelse(group1 == ref_group, group2, group1)) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(
        group1 = g1, group2 = g2,
        p = p, p.adj = p.adj, p.adj.signif = p.adj.signif,
        p_label = dplyr::case_when(
          p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
          p.adj < 0.05  ~ "*",   TRUE ~ "ns"
        )
      ) %>%
      dplyr::filter(p.adj < p_thresh)
    if (nrow(dunn_res) == 0) dunn_res <- NULL
    return(list(kw = kw_test, dunn = dunn_res))
  } else {
    return(list(kw = kw_test, dunn = NULL))
  }
}

# Add significance brackets to an existing ggplot
add_significance_brackets <- function(p, comparisons, y_max, y_min = NULL,
                                      step_fraction = 0.08, textsize = 4) {
  if (is.null(comparisons) || nrow(comparisons) == 0) return(p)
  if (is.null(y_min)) y_min <- 0
  y_range <- max(y_max - y_min, 1e-6)
  comparisons <- comparisons %>%
    dplyr::mutate(y_pos = y_max + seq_len(nrow(.)) * (y_range * step_fraction))
  p + ggsignif::geom_signif(
    data = comparisons,
    aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p_label),
    manual = TRUE, textsize = textsize, vjust = 0.3,
    tip_length = 0.01, inherit.aes = FALSE
  )
}

# Build Poisson or negative binomial rate model with emmeans post-hoc
build_rate_model <- function(data, count_col, exposure_col, output_prefix, output_dir) {
  full_form <- as.formula(
    paste0(count_col, " ~ Genotype + offset(log(", exposure_col, "))")
  )
  null_form <- as.formula(
    paste0(count_col, " ~ 1 + offset(log(", exposure_col, "))")
  )
  
  pois_model    <- glm(full_form, data = data, family = poisson(link = "log"))
  overdisp_ratio <- pois_model$deviance / pois_model$df.residual
  
  nb_model <- tryCatch(
    MASS::glm.nb(full_form, data = data),
    error = function(e) NULL
  )
  
  final_model <- if (!is.null(nb_model) && overdisp_ratio > 1.5) nb_model else pois_model
  model_name  <- if (!is.null(nb_model) && overdisp_ratio > 1.5) "negative binomial" else "poisson"
  
  omnibus_p <- tryCatch({
    if (model_name == "negative binomial") {
      null_model <- MASS::glm.nb(null_form, data = data)
      an_tab     <- anova(null_model, final_model)
      as.numeric(an_tab$`Pr(Chi)`[2])
    } else {
      null_model <- glm(null_form, data = data, family = poisson(link = "log"))
      an_tab     <- anova(null_model, final_model, test = "Chisq")
      as.numeric(an_tab$`Pr(>Chi)`[2])
    }
  }, error = function(e) NA_real_)
  
  emm <- emmeans::emmeans(final_model, ~ Genotype)
  
  pairwise_table <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    dplyr::mutate(
      IRR     = exp(estimate),
      CI_low  = exp(estimate - 1.96 * SE),
      CI_high = exp(estimate + 1.96 * SE),
      p_label = dplyr::case_when(
        p.value < 0.0001 ~ "****", p.value < 0.001 ~ "***",
        p.value < 0.01   ~ "**",   p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  readr::write_csv(
    pairwise_table,
    file.path(output_dir, paste0(output_prefix, "_pairwise_all_genotypes.csv"))
  )
  
  list(
    final_model    = final_model,
    model_name     = model_name,
    overdisp_ratio = overdisp_ratio,
    omnibus_p      = omnibus_p,
    emmeans        = emm,
    pairwise_table = pairwise_table
  )
}

# Annotate stacked bar plots with a triangle where mouse count drops below 3
build_nmouse_annotations <- function(avg_df, x_col, y_pad_frac = 0.05, y_pad_min = 0.05) {
  totals <- avg_df %>%
    dplyr::group_by(Genotype, x = .data[[x_col]]) %>%
    dplyr::summarise(
      total_height = sum(Mean_Tumor_Count, na.rm = TRUE),
      N_Mice       = dplyr::first(N_Mice),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Genotype, x)
  
  global_max <- max(totals$total_height, na.rm = TRUE)
  y_pad      <- max(global_max * y_pad_frac, y_pad_min)
  
  totals %>%
    dplyr::group_by(Genotype) %>%
    dplyr::mutate(
      run_start = dplyr::row_number() == 1 | N_Mice != dplyr::lag(N_Mice),
      run_id    = cumsum(run_start)
    ) %>%
    dplyr::group_by(Genotype, run_id) %>%
    dplyr::summarise(
      x_start      = min(x, na.rm = TRUE),
      N_Mice       = dplyr::first(N_Mice),
      first_height = total_height[which.min(x)],
      .groups = "drop"
    ) %>%
    dplyr::group_by(Genotype) %>%
    dplyr::arrange(x_start, .by_group = TRUE) %>%
    dplyr::mutate(
      prev_N = dplyr::lag(N_Mice),
      y      = first_height + y_pad,
      label  = "\u25B2"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(N_Mice < 3, is.na(prev_N) | prev_N >= 3)
}

# Kruskal-Wallis + Dunn post-hoc vs WT/WT, with labelled p-values for plotting
run_kruskal_dunn_wt <- function(data, value_col, group_col = "Genotype",
                                label = "", output_dir = NULL) {
  vals   <- data[[value_col]]
  groups <- data[[group_col]]
  kw     <- kruskal.test(vals ~ groups)
  
  dunn_tbl <- tryCatch({
    rstatix::dunn_test(
      data    = data,
      formula = as.formula(paste(value_col, "~", group_col)),
      p.adjust.method = "BH"
    ) %>%
      dplyr::filter(group1 == "WT/WT" | group2 == "WT/WT") %>%
      dplyr::mutate(
        g1_new = ifelse(group2 == "WT/WT", group2, group1),
        g2_new = ifelse(group2 == "WT/WT", group1, group2)
      ) %>%
      dplyr::mutate(group1 = g1_new, group2 = g2_new) %>%
      dplyr::select(-g1_new, -g2_new) %>%
      dplyr::mutate(
        stars   = dplyr::case_when(
          p.adj < 0.0001 ~ "****", p.adj < 0.001 ~ "***",
          p.adj < 0.01   ~ "**",   p.adj < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_label = paste0(stars, "\np=", signif(p.adj, 2))
      ) %>%
      dplyr::filter(stars != "ns")
  }, error = function(e) { message("Dunn failed: ", e$message); tibble::tibble() })
  
  if (!is.null(output_dir) && nchar(label) > 0) {
    capture.output(kw, file = file.path(output_dir, paste0(label, "_kruskal.txt")))
    if (nrow(dunn_tbl) > 0)
      readr::write_csv(dunn_tbl,
                       file.path(output_dir, paste0(label, "_dunn_wt_comparisons.csv")))
  }
  list(kruskal = kw, dunn = dunn_tbl)
}

# Assign y-positions for significance brackets above the data range
make_y_positions <- function(pv_tbl, y_max, step = 0.10) {
  if (nrow(pv_tbl) == 0) return(pv_tbl)
  pv_tbl %>%
    dplyr::arrange(group2) %>%
    dplyr::mutate(
      y.position = y_max * seq(1.08, 1.08 + step * (dplyr::n() - 1), by = step)
    )
}

make_y_positions_all <- function(pv_tbl, y_max, y_min = NULL, step_fraction = 0.08) {
  if (nrow(pv_tbl) == 0) return(pv_tbl)
  if (is.null(y_min)) y_min <- 0
  y_range <- max(y_max - y_min, 1e-6)
  pv_tbl %>%
    dplyr::mutate(
      y.position = y_max + seq_len(dplyr::n()) * (y_range * step_fraction)
    )
}

# ============================================================
# 3. READ FILE
# ============================================================
# Input: one Excel sheet where rows are timepoints and columns are individual tumors.
# Column 1 = week label, column 2 = timepoint (day), remaining columns = tumor sizes.
df_raw <- openxlsx::read.xlsx(input_file, sheet = 1)
colnames(df_raw) <- colnames(df_raw) %>% stringr::str_trim()
colnames(df_raw)[1] <- "Week"
colnames(df_raw)[2] <- "Timepoint"
tumor_cols         <- colnames(df_raw)[-(1:2)]
colnames(df_raw)   <- c("Week", "Timepoint", make.unique(tumor_cols, sep = "_tumor"))

df_raw <- df_raw %>%
  dplyr::mutate(
    Week      = as.character(Week),
    Timepoint = readr::parse_number(as.character(Timepoint),
                                    locale = readr::locale(decimal_mark = "."))
  ) %>%
  dplyr::mutate(dplyr::across(
    -(Week:Timepoint),
    ~ {
      x <- as.character(.x)
      x <- stringr::str_replace_all(x, "(?<=\\d),(?=\\d)", ".")
      readr::parse_number(x, locale = readr::locale(decimal_mark = "."))
    }
  ))

# ============================================================
# 4. LONG FORMAT + MAP TUMORS TO MOUSE / GENOTYPE
# ============================================================
df_long <- df_raw %>%
  tidyr::pivot_longer(cols = -(Week:Timepoint),
                      names_to = "Tumor_Column", values_to = "Tumor_Size")

tumor_map <- tibble::tibble(Tumor_Column = unique(df_long$Tumor_Column)) %>%
  dplyr::mutate(
    Tumor_Column_Base = Tumor_Column %>%
      stringr::str_remove("\\.\\.\\.[0-9]+$") %>%
      stringr::str_remove("_tumor\\.[0-9]+$|_tumor[0-9]*$") %>%
      stringr::str_squish(),
    Mouse_ID = stringr::str_extract(Tumor_Column_Base, "^[0-9]+"),
    Genotype = stringr::str_extract(Tumor_Column_Base, "WT/WT|T/WT|D/WT|T/T|T/D|D/D")
  ) %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::mutate(Tumor_Number = dplyr::row_number()) %>%
  dplyr::ungroup()

df_long <- df_long %>%
  dplyr::mutate(
    Tumor_Column_Base = Tumor_Column %>%
      stringr::str_remove("\\.\\.\\.[0-9]+$") %>%
      stringr::str_remove("_tumor\\.[0-9]+$|_tumor[0-9]*$") %>%
      stringr::str_squish()
  ) %>%
  dplyr::left_join(tumor_map, by = c("Tumor_Column", "Tumor_Column_Base")) %>%
  dplyr::select(Week, Timepoint, Mouse_ID, Genotype, Tumor_Number, Tumor_Size) %>%
  dplyr::mutate(
    Genotype      = factor(Genotype, levels = genotype_levels),
    Size_Category = dplyr::case_when(
      is.na(Tumor_Size) | Tumor_Size <= 0 ~ NA_character_,
      Tumor_Size < 10   ~ "<10",
      Tumor_Size < 100  ~ "10-100",
      Tumor_Size < 500  ~ "100-500",
      Tumor_Size >= 500 ~ ">500"
    ),
    Size_Category = factor(Size_Category, levels = size_levels),
    Week_Num      = readr::parse_number(Week),
    Tumor_ID      = paste(Mouse_ID, Genotype, Tumor_Number, sep = "_")
  )

readr::write_csv(df_long, file.path(output_dir, "tumor_sizes_long.csv"))

# ============================================================
# 4b. FOLLOW-UP WINDOW PER MOUSE
#     True_Last_Timepoint = last day with any measurement
#     Mouse_First_Tumor_Timepoint = day of first observed tumor
# ============================================================
study_start_timepoint <- min(df_long$Timepoint, na.rm = TRUE)

mouse_true_last <- df_long %>%
  dplyr::filter(!is.na(Mouse_ID), !is.na(Genotype)) %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(
    True_Last_Timepoint         = max(Timepoint[!is.na(Tumor_Size)], na.rm = TRUE),
    Mouse_First_Tumor_Timepoint = min(Timepoint[!is.na(Tumor_Size)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Study_Start_Timepoint = study_start_timepoint,
    Global_Last           = max(True_Last_Timepoint, na.rm = TRUE),
    Status                = ifelse(True_Last_Timepoint < Global_Last, 1L, 0L),
    Followup_From_Start   = as.numeric(True_Last_Timepoint - Study_Start_Timepoint),
    Followup_From_Onset   = as.numeric(True_Last_Timepoint - Mouse_First_Tumor_Timepoint),
    Mouse_ID              = factor(Mouse_ID),
    Genotype              = factor(Genotype, levels = genotype_levels)
  ) %>%
  dplyr::select(-Global_Last)

# ============================================================
# 5. ONSET VARIABLES
# ============================================================
# Onset_Day = first timepoint at which each individual tumor was observed
tumor_onset <- df_long %>%
  dplyr::filter(Tumor_Size > 0) %>%
  dplyr::group_by(Tumor_ID) %>%
  dplyr::summarise(
    Onset_Day  = min(Timepoint, na.rm = TRUE),
    Onset_Week = min(Week_Num,  na.rm = TRUE),
    .groups    = "drop"
  )

df_long <- df_long %>%
  dplyr::left_join(tumor_onset, by = "Tumor_ID") %>%
  dplyr::mutate(
    Days_From_Onset  = Timepoint - Onset_Day,
    Weeks_From_Onset = Week_Num  - Onset_Week
  )

readr::write_csv(df_long, file.path(output_dir, "tumor_sizes_long_with_onset.csv"))
df_onset <- df_long %>% dplyr::filter(!is.na(Days_From_Onset), Days_From_Onset >= 0)

# ============================================================
# 6. AVERAGE TUMOR COUNTS OVER TIME (absolute days)
#    Zero-fill only at actually observed mouse-timepoints
# ============================================================
mouse_observed_days <- df_long %>%
  dplyr::filter(!is.na(Genotype), !is.na(Mouse_ID), !is.na(Timepoint)) %>%
  dplyr::group_by(Timepoint, Genotype, Mouse_ID) %>%
  dplyr::summarise(Observed = any(!is.na(Tumor_Size)), .groups = "drop") %>%
  dplyr::filter(Observed) %>%
  dplyr::select(Timepoint, Genotype, Mouse_ID)

all_mouse_day_size <- mouse_observed_days %>%
  tidyr::crossing(Size_Category = factor(size_levels, levels = size_levels))

tumor_counts_mouse_size <- df_long %>%
  dplyr::filter(!is.na(Genotype), !is.na(Mouse_ID), !is.na(Size_Category), Tumor_Size > 0) %>%
  dplyr::group_by(Timepoint, Genotype, Mouse_ID, Size_Category) %>%
  dplyr::summarise(Tumor_Count = dplyr::n(), .groups = "drop")

tumor_counts_mouse_size <- all_mouse_day_size %>%
  dplyr::left_join(tumor_counts_mouse_size,
                   by = c("Timepoint", "Genotype", "Mouse_ID", "Size_Category")) %>%
  dplyr::mutate(Tumor_Count = ifelse(is.na(Tumor_Count), 0, Tumor_Count))

avg_tumor_counts_genotype_size <- tumor_counts_mouse_size %>%
  dplyr::group_by(Timepoint, Genotype, Size_Category) %>%
  dplyr::summarise(
    Mean_Tumor_Count = mean(Tumor_Count, na.rm = TRUE),
    SD_Tumor_Count   = sd(Tumor_Count,   na.rm = TRUE),
    N_Mice           = dplyr::n(),
    SE_Tumor_Count   = SD_Tumor_Count / sqrt(N_Mice),
    .groups          = "drop"
  ) %>%
  dplyr::arrange(Timepoint, Genotype, Size_Category)

day_breaks         <- sort(unique(avg_tumor_counts_genotype_size$Timepoint))
day_breaks_to_show <- day_breaks[seq(1, length(day_breaks), by = 8)]

nmouse_ann_abs <- build_nmouse_annotations(
  avg_df = avg_tumor_counts_genotype_size,
  x_col  = "Timepoint"
)

readr::write_csv(
  avg_tumor_counts_genotype_size,
  file.path(output_dir, "average_tumor_counts_per_day_genotype_sizeclass_including_zeros.csv")
)

# Figure: average tumor counts over absolute time (lower panel, used in paper)
p_avg_counts_bar_lower <- ggplot(
  avg_tumor_counts_genotype_size,
  aes(x = Timepoint, y = Mean_Tumor_Count, fill = Size_Category)
) +
  geom_col(width = 6.5, position = position_stack(reverse = TRUE)) +
  geom_text(
    data        = nmouse_ann_abs,
    aes(x = x_start, y = y, label = label),
    inherit.aes = FALSE, hjust = 0.5, vjust = 0.5, size = 2.5, color = "black"
  ) +
  ggh4x::facet_wrap2(~ Genotype, ncol = 3, scales = "fixed",
                     axes = "all", remove_labels = "y") +
  scale_fill_manual(values = size_palette, breaks = size_levels, drop = FALSE) +
  scale_x_continuous(breaks = day_breaks_to_show) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x        = element_text(size = 10),
    axis.ticks.x       = element_line(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    panel.spacing.x    = unit(0.2, "lines"),
    panel.spacing.y    = unit(0.2, "lines"),
    plot.margin        = margin(10, 10, 18, 10)
  ) +
  labs(
    title = "Average tumor counts over time by genotype and size class",
    x     = "Day",
    y     = "Average number of tumors",
    fill  = "Tumor size"
  )

ggsave(
  file.path(output_dir, "average_tumor_counts_per_day_genotype_sizeclass_barplot_lower.png"),
  p_avg_counts_bar_lower, width = 10, height = 5, dpi = 300
)

# ============================================================
# 7. AVERAGE TUMOR COUNTS FROM MOUSE TUMOR ONSET (weekly)
#    Zero-fill only at actually observed mouse-weeks
# ============================================================
mouse_onset_week <- df_long %>%
  dplyr::filter(Tumor_Size > 0, !is.na(Mouse_ID), !is.na(Genotype)) %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(Mouse_Onset_Week = min(Week_Num, na.rm = TRUE), .groups = "drop")

df_mouse_onset_week <- df_long %>%
  dplyr::left_join(mouse_onset_week, by = c("Mouse_ID", "Genotype")) %>%
  dplyr::mutate(Weeks_From_Mouse_Onset = Week_Num - Mouse_Onset_Week) %>%
  dplyr::filter(!is.na(Weeks_From_Mouse_Onset), Weeks_From_Mouse_Onset >= 0)

mouse_observed_weeks_from_onset <- df_mouse_onset_week %>%
  dplyr::filter(!is.na(Genotype), !is.na(Mouse_ID)) %>%
  dplyr::group_by(Weeks_From_Mouse_Onset, Genotype, Mouse_ID) %>%
  dplyr::summarise(Observed = any(!is.na(Tumor_Size)), .groups = "drop") %>%
  dplyr::filter(Observed) %>%
  dplyr::select(Weeks_From_Mouse_Onset, Genotype, Mouse_ID)

all_mouse_week_size <- mouse_observed_weeks_from_onset %>%
  tidyr::crossing(Size_Category = factor(size_levels, levels = size_levels))

tumor_counts_mouse_size_mouse_onset <- df_mouse_onset_week %>%
  dplyr::filter(!is.na(Genotype), !is.na(Mouse_ID), !is.na(Size_Category), Tumor_Size > 0) %>%
  dplyr::group_by(Weeks_From_Mouse_Onset, Genotype, Mouse_ID, Size_Category) %>%
  dplyr::summarise(Tumor_Count = dplyr::n(), .groups = "drop")

tumor_counts_mouse_size_mouse_onset <- all_mouse_week_size %>%
  dplyr::left_join(
    tumor_counts_mouse_size_mouse_onset,
    by = c("Weeks_From_Mouse_Onset", "Genotype", "Mouse_ID", "Size_Category")
  ) %>%
  dplyr::mutate(Tumor_Count = ifelse(is.na(Tumor_Count), 0, Tumor_Count))

avg_tumor_counts_mouse_onset <- tumor_counts_mouse_size_mouse_onset %>%
  dplyr::group_by(Weeks_From_Mouse_Onset, Genotype, Size_Category) %>%
  dplyr::summarise(
    Mean_Tumor_Count = mean(Tumor_Count, na.rm = TRUE),
    SD_Tumor_Count   = sd(Tumor_Count, na.rm = TRUE),
    N_Mice           = dplyr::n(),
    SE_Tumor_Count   = SD_Tumor_Count / sqrt(N_Mice),
    .groups          = "drop"
  ) %>%
  dplyr::arrange(Weeks_From_Mouse_Onset, Genotype, Size_Category)

week_breaks_to_show_lower <- seq(0, 32, by = 6)
week_breaks_mouse_onset   <- seq(0, 32, by = 4)

nmouse_ann_onset <- build_nmouse_annotations(
  avg_df = avg_tumor_counts_mouse_onset,
  x_col  = "Weeks_From_Mouse_Onset"
)

readr::write_csv(
  avg_tumor_counts_mouse_onset,
  file.path(output_dir, "average_tumor_counts_per_week_from_mouse_onset_genotype_sizeclass.csv")
)

# Figure: average tumor counts from mouse onset (lower panel, used in paper)
p_avg_counts_bar_mouse_onset_lower <- ggplot(
  avg_tumor_counts_mouse_onset,
  aes(x = Weeks_From_Mouse_Onset, y = Mean_Tumor_Count, fill = Size_Category)
) +
  geom_col(width = 0.9, position = position_stack(reverse = TRUE)) +
  geom_text(
    data        = nmouse_ann_onset,
    aes(x = x_start, y = y, label = label),
    inherit.aes = FALSE, hjust = 0.5, vjust = 0.5, size = 2.5, color = "black"
  ) +
  ggh4x::facet_wrap2(~ Genotype, ncol = 3, scales = "fixed",
                     axes = "all", remove_labels = "y") +
  scale_fill_manual(values = size_palette, breaks = size_levels, drop = FALSE) +
  scale_x_continuous(breaks = week_breaks_to_show_lower) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(xlim = c(0, 32), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 9),
    axis.ticks.x     = element_line(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(0.15, "lines"),
    panel.spacing.y  = unit(0.15, "lines"),
    plot.margin      = margin(10, 10, 18, 10)
  ) +
  labs(
    title = "Average tumor counts from mouse tumor onset by genotype and size class",
    x     = "Weeks from tumor onset per mouse",
    y     = "Mean tumor count per mouse",
    fill  = "Tumor size"
  )

ggsave(
  file.path(output_dir, "average_tumor_counts_by_genotype_sizeclass_weeks_from_mouse_onset_lower.png"),
  p_avg_counts_bar_mouse_onset_lower, width = 10, height = 5, dpi = 300
)

# ============================================================
# 8. TUMOR GROWTH CURVES
# ============================================================
tumor_growth_onset_df <- df_long %>%
  filter(!is.na(Genotype), !is.na(Tumor_ID), !is.na(Tumor_Size),
         Tumor_Size > 0, !is.na(Days_From_Onset), Days_From_Onset >= 0)

onset_day_breaks_growth         <- sort(unique(tumor_growth_onset_df$Days_From_Onset))
onset_day_breaks_growth_to_show <- onset_day_breaks_growth[
  seq(1, length(onset_day_breaks_growth), by = 10)]

# Figure: per-tumor growth curves by genotype (upper panel, no x-axis labels; used in paper)
p_tumor_growth_onset_upper <- ggplot(
  tumor_growth_onset_df,
  aes(x = Days_From_Onset, y = Tumor_Size, group = Tumor_ID, color = Genotype)
) +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_point(alpha = 0.5, size = 1.2) +
  ggh4x::facet_wrap2(~ Genotype, ncol = 3, scales = "fixed",
                     axes = "all", remove_labels = "all") +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  scale_x_continuous(breaks = onset_day_breaks_growth_to_show) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.ticks.x     = element_line(),
    axis.text.x      = element_text(size = 9, angle = 45, hjust = 1),
    axis.title.x     = element_blank()
  ) +
  labs(title = "Tumor growth from tumor onset by genotype", y = "Tumor size")

ggsave(
  file.path(output_dir, "tumor_growth_curves_per_tumor_by_genotype_days_from_onset_upper.png"),
  p_tumor_growth_onset_upper, width = 12, height = 7, dpi = 300
)

# Figure: growth curves by genotype and maximum size class (lower panel, used in paper)
tumor_max_size <- df_long %>%
  dplyr::filter(!is.na(Tumor_ID), !is.na(Tumor_Size), Tumor_Size > 0) %>%
  dplyr::group_by(Tumor_ID, Genotype) %>%
  dplyr::summarise(Max_Size = max(Tumor_Size, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    Max_Size_Category = dplyr::case_when(
      Max_Size < 10   ~ "<10",
      Max_Size < 100  ~ "10-100",
      Max_Size < 500  ~ "100-500",
      Max_Size >= 500 ~ ">500"
    ),
    Max_Size_Category = factor(Max_Size_Category, levels = size_levels)
  )

tumor_growth_onset_fixedcat_df <- df_long %>%
  dplyr::filter(!is.na(Genotype), !is.na(Tumor_ID), !is.na(Tumor_Size),
                Tumor_Size > 0, !is.na(Days_From_Onset), Days_From_Onset >= 0) %>%
  dplyr::left_join(
    tumor_max_size %>% dplyr::select(Tumor_ID, Max_Size_Category),
    by = "Tumor_ID"
  )

p_tumor_growth_onset_fixedcat_lower <- ggplot(
  tumor_growth_onset_fixedcat_df,
  aes(x = Days_From_Onset, y = Tumor_Size, group = Tumor_ID, color = Max_Size_Category)
) +
  geom_line(alpha = 0.45, linewidth = 0.6) +
  geom_point(alpha = 0.5, size = 1.2) +
  ggh4x::facet_grid2(Max_Size_Category ~ Genotype, scales = "free_y",
                     drop = FALSE, axes = "all", remove_labels = "all") +
  scale_color_manual(values = size_palette, breaks = size_levels, drop = FALSE) +
  scale_x_continuous(breaks = onset_day_breaks_growth_to_show) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.ticks.x     = element_line(),
    axis.text.x      = element_text(size = 9, angle = 45, hjust = 1)
  ) +
  labs(
    title = "Tumor growth from tumor onset by genotype and maximum size class",
    x     = "Days from tumor onset",
    y     = "Tumor size"
  )

ggsave(
  file.path(output_dir, "tumor_growth_curves_by_genotype_and_max_sizeclass_days_from_onset_lower.png"),
  p_tumor_growth_onset_fixedcat_lower, width = 14, height = 10, dpi = 300
)

# ============================================================
# 9. TUMOR ONSET PER GENOTYPE
# ============================================================
tumor_lag_df <- df_long %>%
  dplyr::filter(!is.na(Tumor_ID), !is.na(Genotype), !is.na(Onset_Day)) %>%
  dplyr::distinct(Tumor_ID, Genotype, Onset_Day)

readr::write_csv(tumor_lag_df, file.path(output_dir, "tumor_lag_times_per_tumor.csv"))

lag_stats <- perform_kw_dunn(tumor_lag_df, "Onset_Day", "Genotype")
capture.output(lag_stats$kw, file = file.path(output_dir, "lag_times_kruskal_test.txt"))
if (!is.null(lag_stats$dunn))
  readr::write_csv(lag_stats$dunn, file.path(output_dir, "lag_times_dunn_test.csv"))

p_tumor_onset <- ggplot(
  tumor_lag_df, aes(x = Genotype, y = Onset_Day)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = NA, color = "black") +
  geom_jitter(aes(color = Genotype), width = 0.18, alpha = 0.7, size = 2) +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(title = "Tumor onset per genotype", x = "Genotype", y = "Tumor onset (day)")

if (!is.null(lag_stats$dunn) && nrow(lag_stats$dunn) > 0) {
  y_max         <- max(tumor_lag_df$Onset_Day, na.rm = TRUE)
  y_min         <- min(tumor_lag_df$Onset_Day, na.rm = TRUE)
  p_tumor_onset <- add_significance_brackets(p_tumor_onset, lag_stats$dunn, y_max, y_min)
}

ggsave(
  file.path(output_dir, "tumor_onset_per_genotype.png"),
  p_tumor_onset, width = 5, height = 4, dpi = 300
)

# ============================================================
# 10. CROSS-SECTIONAL TUMOR SIZE AT DAY 145
# ============================================================
target_day   <- 145
available_tp <- sort(unique(df_long$Timepoint))
closest_tp   <- available_tp[which.min(abs(available_tp - target_day))]
message("Target day ", target_day, " -> using: day ", closest_tp)

tp_data <- df_long %>%
  dplyr::filter(Timepoint == closest_tp, !is.na(Mouse_ID), !is.na(Genotype)) %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(
    Mean_Tumor_Size = ifelse(any(Tumor_Size > 0, na.rm = TRUE),
                             mean(Tumor_Size[Tumor_Size > 0], na.rm = TRUE), NA_real_),
    Max_Tumor_Size  = ifelse(any(Tumor_Size > 0, na.rm = TRUE),
                             max(Tumor_Size[Tumor_Size > 0], na.rm = TRUE), NA_real_),
    N_tumors        = sum(Tumor_Size > 0, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  dplyr::mutate(Genotype = factor(Genotype, levels = genotype_levels))

tp_data_tumors <- tp_data %>% dplyr::filter(!is.na(Mean_Tumor_Size))
readr::write_csv(tp_data, file.path(output_dir, paste0("tumor_size_crosssection_day", closest_tp, ".csv")))

tp_stats <- if (length(unique(tp_data_tumors$Genotype)) >= 2) {
  res <- perform_kw_dunn(tp_data_tumors, "Mean_Tumor_Size", "Genotype")
  capture.output(res$kw,
                 file = file.path(output_dir, paste0("tumor_size_crosssection_day", closest_tp, "_kruskal.txt")))
  if (!is.null(res$dunn))
    readr::write_csv(res$dunn,
                     file.path(output_dir, paste0("tumor_size_crosssection_day", closest_tp, "_dunn.csv")))
  res
} else {
  list(kw = list(p.value = NA), dunn = NULL)
}

p_tp <- ggplot(tp_data_tumors, aes(x = Genotype, y = Mean_Tumor_Size, color = Genotype)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 3) +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title    = paste0("Mean tumor size per mouse at day ", closest_tp),
    subtitle = paste0("Only mice with \u2265 1 tumor | KW p=", signif(tp_stats$kw$p.value, 3)),
    x        = "Genotype",
    y        = "Mean tumor size"
  )

if (!is.null(tp_stats$dunn) && nrow(tp_stats$dunn) > 0) {
  y_max  <- max(tp_data_tumors$Mean_Tumor_Size, na.rm = TRUE)
  y_min  <- min(tp_data_tumors$Mean_Tumor_Size, na.rm = TRUE)
  p_tp   <- add_significance_brackets(p_tp, tp_stats$dunn, y_max, y_min)
}

ggsave(
  file.path(output_dir, paste0("tumor_size_crosssection_day", closest_tp, ".png")),
  p_tp, width = 9, height = 7, dpi = 300, bg = "white"
)

# ============================================================
# 11. TUMOR EXPONENTIAL GROWTH RATES
#     Per-tumor linear model on log(size) ~ days from onset
#     Requires >= 3 observations per tumor (n = 2 gives unreliable slope estimates)
#     Mouse-level random effect via linear mixed model
#     Note: includes tumors that plateau or shrink; growth rate
#     estimates for those are not reliable exponential rates.
# ============================================================
tumor_growth_fit_df <- df_long %>%
  dplyr::filter(!is.na(Tumor_ID), !is.na(Genotype), !is.na(Days_From_Onset),
                Days_From_Onset >= 0, !is.na(Tumor_Size), Tumor_Size > 0) %>%
  dplyr::group_by(Tumor_ID, Genotype) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::group_modify(~ {
    fit <- lm(log(Tumor_Size) ~ Days_From_Onset, data = .x)
    tibble::tibble(
      Growth_Rate = coef(fit)[["Days_From_Onset"]],
      Intercept   = coef(fit)[["(Intercept)"]],
      N_points    = nrow(.x)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Mouse_ID = stringr::str_extract(Tumor_ID, "^[^_]+"),
    Genotype = factor(Genotype, levels = genotype_levels),
    Mouse_ID = factor(Mouse_ID)
  )

readr::write_csv(
  tumor_growth_fit_df,
  file.path(output_dir, "tumor_exponential_growth_rates.csv")
)

growth_rate_lmm <- nlme::lme(
  fixed     = Growth_Rate ~ Genotype,
  random    = ~ 1 | Mouse_ID,
  data      = tumor_growth_fit_df,
  method    = "REML",
  na.action = na.omit,
  control   = nlme::lmeControl(returnObject = TRUE)
)

capture.output(
  summary(growth_rate_lmm),
  file = file.path(output_dir, "growth_rates_mixed_model_summary.txt")
)

growth_rate_anova <- anova(growth_rate_lmm)
capture.output(
  growth_rate_anova,
  file = file.path(output_dir, "growth_rates_mixed_model_anova.txt")
)

growth_rate_emm <- emmeans::emmeans(growth_rate_lmm, ~ Genotype)

growth_rate_pairs <- emmeans::contrast(
  growth_rate_emm, method = "pairwise", adjust = "tukey"
) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  dplyr::mutate(
    p_label = dplyr::case_when(
      p.value < 0.0001 ~ "****", p.value < 0.001 ~ "***",
      p.value < 0.01   ~ "**",   p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

readr::write_csv(
  growth_rate_pairs,
  file.path(output_dir, "growth_rates_mixed_model_posthoc_all_pairs.csv")
)

# Mixed model results are saved to file; for plot annotations use KW/Dunn
# so that group names match the factor levels on the x-axis
kw_gr <- run_kruskal_dunn_wt(
  tumor_growth_fit_df, "Growth_Rate",
  label = "growth_rates", output_dir = output_dir
)

pv_gr <- kw_gr$dunn

# Figure: growth rate boxplot (used in paper)
y_max_gr <- max(tumor_growth_fit_df$Growth_Rate, na.rm = TRUE)

p_growth_rates_box <- ggplot(
  tumor_growth_fit_df,
  aes(x = Genotype, y = Growth_Rate, color = Genotype)
) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 1.8) +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title    = "Tumor exponential growth rates per genotype",
    subtitle = paste0(
      "LMM omnibus p = ", signif(growth_rate_anova$`p-value`[2], 3),
      "; KW p = ", signif(kw_gr$kruskal$p.value, 3)
    ),
    x = "Genotype",
    y = "Growth rate (r)"
  )

if (!is.null(pv_gr) && nrow(pv_gr) > 0) {
  # Map genotype names to numeric x positions
  geno_pos <- setNames(seq_along(genotype_levels), genotype_levels)
  pv_gr <- pv_gr %>%
    dplyr::arrange(group2) %>%
    dplyr::mutate(
      x_start  = geno_pos[group1],
      x_end    = geno_pos[group2],
      y_pos    = y_max_gr + seq(0.06, by = 0.07, length.out = dplyr::n())
    )
  
  for (i in seq_len(nrow(pv_gr))) {
    p_growth_rates_box <- p_growth_rates_box +
      ggsignif::geom_signif(
        xmin        = pv_gr$x_start[i],
        xmax        = pv_gr$x_end[i],
        y_position  = pv_gr$y_pos[i],
        annotation  = pv_gr$stars[i],
        tip_length  = 0.005,
        textsize    = 3.5,
        vjust       = 0.3,
        color       = "black"
      )
  }
}

ggsave(
  file.path(output_dir, "tumor_growth_rates_boxplot.png"),
  p_growth_rates_box, width = 7, height = 5, dpi = 300
)

# ============================================================
# 12. INCIDENCE RATE MODELS
# ============================================================

# --- 12a. Overall incidence rate (from experiment start) ---
mouse_obs_time <- mouse_true_last %>%
  dplyr::mutate(Obs_Weeks = Followup_From_Start / 7) %>%
  dplyr::filter(Obs_Weeks > 0)

mouse_tumor_counts_total <- df_long %>%
  dplyr::filter(!is.na(Mouse_ID), !is.na(Genotype), Tumor_Size > 0, !is.na(Tumor_ID)) %>%
  dplyr::distinct(Mouse_ID, Genotype, Tumor_ID) %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(N_Tumors = dplyr::n(), .groups = "drop")

mouse_ir_input <- mouse_obs_time %>%
  dplyr::left_join(mouse_tumor_counts_total, by = c("Mouse_ID", "Genotype")) %>%
  dplyr::mutate(
    N_Tumors  = ifelse(is.na(N_Tumors), 0L, N_Tumors),
    Obs_Weeks = as.numeric(Obs_Weeks),
    Genotype  = factor(Genotype, levels = genotype_levels)
  ) %>%
  dplyr::filter(Obs_Weeks > 0)

ir_start_stats <- build_rate_model(
  mouse_ir_input, "N_Tumors", "Obs_Weeks",
  "incidence_rate_from_start", output_dir
)

# --- 12b. Incidence rate from first tumor onset ---
mouse_obs_time_onset <- mouse_true_last %>%
  dplyr::mutate(Obs_Weeks_From_Onset = Followup_From_Onset / 7) %>%
  dplyr::filter(Obs_Weeks_From_Onset > 0)

new_tumors_after_onset <- df_long %>%
  dplyr::filter(!is.na(Tumor_ID), !is.na(Mouse_ID), !is.na(Genotype), Tumor_Size > 0) %>%
  dplyr::group_by(Tumor_ID, Mouse_ID, Genotype) %>%
  dplyr::summarise(Onset_Day = min(Timepoint, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    df_long %>%
      dplyr::filter(Tumor_Size > 0, !is.na(Mouse_ID), !is.na(Genotype)) %>%
      dplyr::group_by(Mouse_ID, Genotype) %>%
      dplyr::summarise(First_Onset_Day = min(Timepoint, na.rm = TRUE), .groups = "drop"),
    by = c("Mouse_ID", "Genotype")
  ) %>%
  dplyr::filter(Onset_Day > First_Onset_Day)

mouse_new_tumors_after_onset <- new_tumors_after_onset %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(N_New_Tumors = dplyr::n(), .groups = "drop")

mouse_ir_onset_input <- mouse_obs_time_onset %>%
  dplyr::left_join(mouse_new_tumors_after_onset, by = c("Mouse_ID", "Genotype")) %>%
  dplyr::mutate(
    N_New_Tumors         = ifelse(is.na(N_New_Tumors), 0L, N_New_Tumors),
    Obs_Weeks_From_Onset = as.numeric(Obs_Weeks_From_Onset),
    Genotype             = factor(Genotype, levels = genotype_levels)
  ) %>%
  dplyr::filter(Obs_Weeks_From_Onset > 0)

ir_onset_stats <- build_rate_model(
  mouse_ir_onset_input, "N_New_Tumors", "Obs_Weeks_From_Onset",
  "incidence_rate_from_onset", output_dir
)

# ============================================================
# 13. INCIDENCE RATE PLOTS WITH STATISTICS
# ============================================================

# --- Build per-mouse incidence rate (from experiment start) ---
mouse_obs_time_plot <- mouse_true_last %>%
  dplyr::mutate(Obs_Weeks = Followup_From_Start / 7) %>%
  dplyr::filter(Obs_Weeks > 0)

new_tumors_plot <- df_long %>%
  dplyr::filter(
    !is.na(Tumor_ID), !is.na(Mouse_ID), !is.na(Genotype),
    !is.na(Timepoint), is.finite(Timepoint), Tumor_Size > 0
  ) %>%
  dplyr::group_by(Tumor_ID, Mouse_ID, Genotype) %>%
  dplyr::summarise(First_Appearance = min(Timepoint, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    mouse_true_last %>%
      dplyr::select(Mouse_ID, Mouse_First_Tumor_Timepoint, True_Last_Timepoint),
    by = "Mouse_ID"
  ) %>%
  dplyr::filter(First_Appearance <= True_Last_Timepoint) %>%
  dplyr::mutate(
    Days_From_Onset = First_Appearance - Mouse_First_Tumor_Timepoint,
    Genotype        = factor(Genotype, levels = genotype_levels)
  )

mouse_ir_plot <- new_tumors_plot %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(N_New_Tumors = dplyr::n(), .groups = "drop") %>%
  dplyr::right_join(
    mouse_obs_time_plot %>% dplyr::select(Mouse_ID, Genotype, Obs_Weeks),
    by = c("Mouse_ID", "Genotype")
  ) %>%
  dplyr::mutate(
    N_New_Tumors = ifelse(is.na(N_New_Tumors), 0, N_New_Tumors),
    IR_mouse     = N_New_Tumors / Obs_Weeks,
    Genotype     = factor(Genotype, levels = genotype_levels)
  )

# --- Build per-mouse incidence rate (from first tumor onset) ---
mouse_obs_time_onset_plot <- mouse_true_last %>%
  dplyr::mutate(Obs_Weeks_From_Onset = Followup_From_Onset / 7) %>%
  dplyr::filter(Obs_Weeks_From_Onset > 0)

new_tumors_after_onset_plot <- new_tumors_plot %>% dplyr::filter(Days_From_Onset > 0)

mouse_ir_onset_plot <- new_tumors_after_onset_plot %>%
  dplyr::group_by(Mouse_ID, Genotype) %>%
  dplyr::summarise(N_New_Tumors_Onset = dplyr::n(), .groups = "drop") %>%
  dplyr::right_join(
    mouse_obs_time_onset_plot %>% dplyr::select(Mouse_ID, Genotype, Obs_Weeks_From_Onset),
    by = c("Mouse_ID", "Genotype")
  ) %>%
  dplyr::mutate(
    N_New_Tumors_Onset = ifelse(is.na(N_New_Tumors_Onset), 0, N_New_Tumors_Onset),
    IR_mouse_onset     = N_New_Tumors_Onset / Obs_Weeks_From_Onset,
    Genotype           = factor(Genotype, levels = genotype_levels)
  ) %>%
  dplyr::filter(is.finite(IR_mouse_onset))

# --- Kruskal-Wallis + Dunn stats vs WT/WT ---
# Using KW/Dunn (mouse as experimental unit) so group names match the factor levels
kw_ir_overall <- run_kruskal_dunn_wt(
  mouse_ir_plot, "IR_mouse",
  label = "incidence_rate_overall", output_dir = output_dir
)

kw_ir_onset <- run_kruskal_dunn_wt(
  mouse_ir_onset_plot, "IR_mouse_onset",
  label = "incidence_rate_from_onset", output_dir = output_dir
)

pv_ir_overall <- make_y_positions(
  kw_ir_overall$dunn,
  y_max = max(mouse_ir_plot$IR_mouse, na.rm = TRUE)
)

pv_ir_onset <- make_y_positions(
  kw_ir_onset$dunn,
  y_max = max(mouse_ir_onset_plot$IR_mouse_onset, na.rm = TRUE)
)

# Figure: overall incidence rate with stats (used in paper)
p_ir_overall_stats <- ggplot(
  mouse_ir_plot, aes(x = Genotype, y = IR_mouse, color = Genotype)
) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title    = "Tumor incidence rate per mouse per week",
    subtitle = paste0(
      "Kruskal-Wallis p=", signif(kw_ir_overall$kruskal$p.value, 3),
      " \u2014 significant vs WT/WT shown"
    ),
    x = "Genotype",
    y = "Incidence rate (tumors / mouse / week)"
  )

if (nrow(pv_ir_overall) > 0) {
  p_ir_overall_stats <- p_ir_overall_stats +
    ggpubr::stat_pvalue_manual(
      pv_ir_overall, label = "p_label", tip.length = 0.01, size = 3.2
    )
}

ggsave(
  file.path(output_dir, "incidence_rate_overall_stats.png"),
  p_ir_overall_stats, width = 8, height = 6, dpi = 300, bg = "white"
)

# Figure: incidence rate from onset with stats (used in paper)
p_ir_onset_stats <- ggplot(
  mouse_ir_onset_plot, aes(x = Genotype, y = IR_mouse_onset, color = Genotype)
) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_color_manual(values = genotype_palette, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title    = "Tumor incidence rate per mouse per week (from onset)",
    subtitle = paste0(
      "Kruskal-Wallis p=", signif(kw_ir_onset$kruskal$p.value, 3),
      " \u2014 significant vs WT/WT shown"
    ),
    x = "Genotype",
    y = "Incidence rate (tumors / mouse / week)"
  )

if (nrow(pv_ir_onset) > 0) {
  p_ir_onset_stats <- p_ir_onset_stats +
    ggpubr::stat_pvalue_manual(
      pv_ir_onset, label = "p_label", tip.length = 0.01, size = 3.2
    )
}

ggsave(
  file.path(output_dir, "incidence_rate_from_onset_stats.png"),
  p_ir_onset_stats, width = 8, height = 6, dpi = 300, bg = "white"
)

# ============================================================
# 14. SUMMARY STATISTICS TABLE
# ============================================================
summary_stats <- dplyr::bind_rows(
  tibble::tibble(Analysis = "Incidence rate (overall)",    Model_p_value = kw_ir_overall$kruskal$p.value),
  tibble::tibble(Analysis = "Incidence rate (from onset)", Model_p_value = kw_ir_onset$kruskal$p.value),
  tibble::tibble(Analysis = "Growth rate (mixed model omnibus)", Model_p_value = growth_rate_anova$`p-value`[2]),
  tibble::tibble(Analysis = "Growth rate (KW vs WT/WT)",         Model_p_value = kw_gr$kruskal$p.value),
  tibble::tibble(Analysis = "Tumor onset day (KW)",        Model_p_value = lag_stats$kw$p.value),
  tibble::tibble(Analysis = "Tumor size at day 145 (KW)",  Model_p_value = tp_stats$kw$p.value)
) %>%
  dplyr::mutate(
    Significant = dplyr::case_when(
      is.na(Model_p_value) ~ "insufficient data",
      Model_p_value < 0.05 ~ "Yes",
      TRUE ~ "No"
    )
  )

readr::write_csv(summary_stats, file.path(output_dir, "statistics_summary_all_analyses.csv"))
message("=== Summary of all main tests ===")
print(summary_stats)