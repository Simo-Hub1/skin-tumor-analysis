suppressPackageStartupMessages({
  library(readr)      # 2.1.4
  library(dplyr)      # 1.1.3
  library(tidyr)      # 1.3.0
  library(stringr)    # 1.5.0
  library(ggplot2)    # 3.4.4
  library(readxl)     # 1.4.3
  library(rstatix)    # 0.7.2
  library(ggpubr)     # 0.6.0
  library(ggh4x)      # 0.2.6
  library(ggsignif)   # 0.6.4
})

# Prevent masking
select    <- dplyr::select
filter    <- dplyr::filter
summarise <- dplyr::summarise
mutate    <- dplyr::mutate

# ============================================================
# 1. INPUT / OUTPUT
# Paths are relative to the project root — open the .Rproj file
# in RStudio before running, or set your working directory manually
# with setwd("path/to/project").
# ============================================================
input_file <- "C:/Users/s.lafirenze/Desktop/Skin paper/Figure carprofen/carprofen data.xlsx"
output_dir <- "C:/Users/s.lafirenze/Desktop/Skin paper/Figure carprofen/output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 2. SETTINGS
# ============================================================
genotype_levels  <- c("WT/WT", "TA/TA", "TA/KD")
treatment_levels <- c("Control", "Carprofen")
condition_levels <- c(
  "WT/WT Control", "WT/WT Carprofen",
  "TA/TA Control", "TA/TA Carprofen",
  "TA/KD Control", "TA/KD Carprofen"
)
size_levels <- c("<10", "10-100", "100-500", ">500")

condition_palette <- c(
  "WT/WT Control"   = "grey20",
  "WT/WT Carprofen" = "grey60",
  "TA/TA Control"   = "goldenrod4",
  "TA/TA Carprofen" = "goldenrod1",
  "TA/KD Control"   = "darkorange4",
  "TA/KD Carprofen" = "darkorange1"
)

size_palette <- c(
  "<10"     = "steelblue2",
  "10-100"  = "dodgerblue4",
  "100-500" = "deeppink3",
  ">500"    = "deeppink4"
)

# ============================================================
# 3. READ FILE
# ============================================================
# Input: Excel file where row 1 is a header to skip, row 2 onward is data.
# Column 1 = label (dropped), column 2 = timepoint (day),
# remaining columns = tumor sizes in mm³.
# Column headers encode mouse ID, genotype, and treatment.
df_raw <- readxl::read_excel(
  input_file, col_names = TRUE, skip = 0, .name_repair = "minimal"
)

df_raw <- df_raw[-1, ]
colnames(df_raw) <- make.unique(colnames(df_raw) %>% str_trim(), sep = "_tumor")
colnames(df_raw)[1] <- "Label"
colnames(df_raw)[2] <- "Timepoint"
df_raw <- df_raw %>% dplyr::select(-Label)

df_raw <- df_raw %>%
  mutate(Timepoint = as.numeric(Timepoint)) %>%
  filter(!is.na(Timepoint)) %>%
  mutate(across(-Timepoint, ~ as.numeric(.x)))

message("Rows: ", nrow(df_raw), " | Columns: ", ncol(df_raw))

# ============================================================
# 4. LONG FORMAT + MAP TUMORS TO MOUSE / GENOTYPE / TREATMENT
# ============================================================
df_long <- df_raw %>%
  pivot_longer(cols = -Timepoint, names_to = "Tumor_Column", values_to = "Tumor_Size")

tumor_map <- tibble(Tumor_Column = unique(df_long$Tumor_Column)) %>%
  mutate(
    Tumor_Column_Base = Tumor_Column %>%
      str_remove("_tumor[0-9]*$") %>%
      str_squish(),
    Treatment = ifelse(str_detect(Tumor_Column_Base, regex("carprofen", ignore_case = TRUE)),
                       "Carprofen", "Control"),
    Genotype  = str_extract(Tumor_Column_Base, "WT/WT|TA/TA|TA/KD"),
    Mouse_ID  = str_extract(Tumor_Column_Base, "^[^_]+")
  ) %>%
  group_by(Mouse_ID, Genotype, Treatment) %>%
  mutate(Tumor_Number = row_number()) %>%
  ungroup()

df_long <- df_long %>%
  mutate(
    Tumor_Column_Base = Tumor_Column %>%
      str_remove("_tumor[0-9]*$") %>%
      str_squish()
  ) %>%
  left_join(tumor_map, by = c("Tumor_Column", "Tumor_Column_Base")) %>%
  dplyr::select(Timepoint, Mouse_ID, Genotype, Treatment, Tumor_Number, Tumor_Size) %>%
  mutate(
    Genotype      = factor(Genotype,  levels = genotype_levels),
    Treatment     = factor(Treatment, levels = treatment_levels),
    Condition     = factor(paste(Genotype, Treatment), levels = condition_levels),
    Size_Category = case_when(
      is.na(Tumor_Size) | Tumor_Size <= 0 ~ NA_character_,
      Tumor_Size < 10   ~ "<10",
      Tumor_Size < 100  ~ "10-100",
      Tumor_Size < 500  ~ "100-500",
      TRUE              ~ ">500"
    ),
    Size_Category = factor(Size_Category, levels = size_levels),
    Tumor_ID      = paste(Mouse_ID, Genotype, Treatment, Tumor_Number, sep = "_")
  )

readr::write_csv(df_long, file.path(output_dir, "carprofen_tumor_sizes_long.csv"))

# ============================================================
# 4b. FOLLOW-UP WINDOW PER MOUSE
# ============================================================
mouse_true_last <- df_long %>%
  filter(!is.na(Mouse_ID), !is.na(Genotype), !is.na(Timepoint)) %>%
  group_by(Mouse_ID, Genotype, Treatment) %>%
  summarise(
    True_Last_Timepoint  = max(Timepoint[!is.na(Tumor_Size)], na.rm = TRUE),
    True_First_Timepoint = min(Timepoint[!is.na(Tumor_Size)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Global_Last   = max(True_Last_Timepoint),
    Status        = ifelse(True_Last_Timepoint < Global_Last, 1L, 0L),
    Survival_Time = as.numeric(True_Last_Timepoint - True_First_Timepoint),
    Mouse_ID      = factor(Mouse_ID),
    Condition     = factor(paste(Genotype, Treatment), levels = condition_levels)
  ) %>%
  dplyr::select(-Global_Last)

readr::write_csv(mouse_true_last, file.path(output_dir, "mouse_true_last_observation.csv"))

# ============================================================
# 5. ONSET VARIABLES
# ============================================================
# Onset_Day = first timepoint at which each individual tumor was observed
tumor_onset <- df_long %>%
  filter(Tumor_Size > 0) %>%
  group_by(Tumor_ID) %>%
  summarise(Onset_Day = min(Timepoint, na.rm = TRUE), .groups = "drop")

df_long <- df_long %>%
  left_join(tumor_onset, by = "Tumor_ID") %>%
  mutate(Days_From_Onset = Timepoint - Onset_Day)

readr::write_csv(df_long, file.path(output_dir, "carprofen_tumor_sizes_long_with_onset.csv"))

# ============================================================
# 6. AVERAGE TUMOR COUNTS PER MOUSE OVER TIME
#    Zero-fill: all mouse-timepoint combinations
# ============================================================
all_mouse_days <- df_long %>%
  filter(!is.na(Condition), !is.na(Mouse_ID)) %>%
  distinct(Timepoint, Condition, Genotype, Treatment, Mouse_ID)

all_mouse_day_size <- tidyr::crossing(
  all_mouse_days,
  tibble(Size_Category = factor(size_levels, levels = size_levels))
)

tumor_counts_mouse <- df_long %>%
  filter(!is.na(Condition), !is.na(Mouse_ID), !is.na(Size_Category), Tumor_Size > 0) %>%
  group_by(Timepoint, Condition, Mouse_ID, Size_Category) %>%
  summarise(Tumor_Count = n(), .groups = "drop")

avg_counts <- all_mouse_day_size %>%
  left_join(tumor_counts_mouse,
            by = c("Timepoint", "Condition", "Mouse_ID", "Size_Category")) %>%
  mutate(Tumor_Count = ifelse(is.na(Tumor_Count), 0, Tumor_Count)) %>%
  group_by(Timepoint, Condition, Genotype, Treatment, Size_Category) %>%
  summarise(Mean_Count = mean(Tumor_Count, na.rm = TRUE), .groups = "drop")

readr::write_csv(avg_counts,
                 file.path(output_dir, "average_tumor_counts_per_condition_over_time.csv"))

day_breaks         <- sort(unique(avg_counts$Timepoint))
day_breaks_to_show <- day_breaks[seq(1, length(day_breaks),
                                     by = max(1, floor(length(day_breaks) / 10)))]
bar_width <- min(diff(sort(unique(avg_counts$Timepoint)))) * 1

# Start x-axis just before first tumor appears
x_start <- df_long %>%
  filter(Tumor_Size > 0, !is.na(Timepoint)) %>%
  summarise(min_day = min(Timepoint, na.rm = TRUE)) %>%
  pull(min_day) - 10

# Fixed y-axis limits: WT/WT panels share one scale, CIN panels share another
wt_y_max <- avg_counts %>%
  filter(str_detect(as.character(Condition), "WT/WT")) %>%
  group_by(Timepoint, Condition) %>%
  summarise(Total = sum(Mean_Count, na.rm = TRUE), .groups = "drop") %>%
  pull(Total) %>% max(na.rm = TRUE) * 1.05

# Figure: average tumor counts over time by condition (used in paper)
p_avg_counts <- ggplot(
  avg_counts, aes(x = Timepoint, y = Mean_Count, fill = Size_Category)
) +
  geom_col(width = bar_width, position = position_stack(reverse = TRUE), na.rm = TRUE) +
  ggh4x::facet_wrap2(~ Condition, ncol = 2, scales = "free_y", axes = "all") +
  scale_fill_manual(values = size_palette, breaks = size_levels, drop = FALSE) +
  scale_x_continuous(breaks = day_breaks_to_show, limits = c(x_start, NA)) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, wt_y_max)),
      scale_y_continuous(limits = c(0, wt_y_max)),
      scale_y_continuous(limits = c(0, 1.3)),
      scale_y_continuous(limits = c(0, 1.3)),
      scale_y_continuous(limits = c(0, 1.3)),
      scale_y_continuous(limits = c(0, 1.3))
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 8, angle = 45, hjust = 1),
    axis.ticks.x     = element_line(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(0.3, "lines"),
    panel.spacing.y  = unit(0.3, "lines")
  ) +
  labs(
    title = "Average tumor counts per mouse over time by condition",
    x     = "Day",
    y     = "Mean tumor count per mouse",
    fill  = "Tumor size"
  )

ggsave(file.path(output_dir, "average_tumor_counts_per_condition_over_time.png"),
       p_avg_counts, width = 10, height = 8, dpi = 300)

# ============================================================
# 7. INCIDENCE RATE PER MOUSE WITH TREATMENT STATS
# ============================================================

# Build per-mouse observation time
mouse_obs_time <- mouse_true_last %>%
  mutate(Obs_Weeks = as.numeric(True_Last_Timepoint - True_First_Timepoint) / 7) %>%
  filter(Obs_Weeks > 0)

# Count all tumors per mouse
all_tumors_condition <- df_long %>%
  filter(!is.na(Tumor_ID), !is.na(Genotype), Tumor_Size > 0) %>%
  distinct(Tumor_ID, Mouse_ID, Genotype, Treatment, Condition)

mouse_tumor_counts <- all_tumors_condition %>%
  group_by(Mouse_ID, Genotype, Treatment) %>%
  summarise(N_Tumors = n(), .groups = "drop")

# Per-mouse incidence rate
mouse_ir <- mouse_obs_time %>%
  left_join(mouse_tumor_counts, by = c("Mouse_ID", "Genotype", "Treatment")) %>%
  mutate(
    N_Tumors  = ifelse(is.na(N_Tumors), 0L, N_Tumors),
    IR_mouse  = N_Tumors / Obs_Weeks,
    Genotype  = factor(Genotype,  levels = genotype_levels),
    Treatment = factor(Treatment, levels = treatment_levels),
    Condition = factor(paste(Genotype, Treatment), levels = condition_levels)
  )

readr::write_csv(mouse_ir, file.path(output_dir, "incidence_rate_per_mouse.csv"))

# Wilcoxon test comparing Control vs Carprofen within each genotype
wilcox_treatment <- mouse_ir %>%
  group_by(Genotype) %>%
  group_modify(~ {
    test <- tryCatch(
      wilcox.test(IR_mouse ~ Treatment, data = .x),
      error = function(e) NULL
    )
    if (is.null(test)) return(tibble::tibble(p_value = NA_real_))
    tibble::tibble(p_value = test$p.value)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE ~ "ns"
    )
  )

readr::write_csv(wilcox_treatment,
                 file.path(output_dir, "incidence_rate_treatment_wilcox.csv"))
message("=== Wilcoxon treatment effect on incidence rate per genotype ===")
print(wilcox_treatment)

# Build significance brackets for the plot
# Map condition names to numeric x positions within each facet
sig_brackets <- wilcox_treatment %>%
  filter(stars != "ns") %>%
  mutate(
    group1 = paste(Genotype, "Control"),
    group2 = paste(Genotype, "Carprofen")
  )

# Figure: incidence rate per mouse with treatment stats (used in paper)
p_ir_stats <- ggplot(
  mouse_ir,
  aes(x = Condition, y = IR_mouse, color = Condition)
) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_color_manual(values = condition_palette, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~ Genotype, ncol = 3, scales = "free_x") +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.x     = element_blank()
  ) +
  labs(
    title    = "Tumor incidence rate per mouse per week",
    subtitle = "Wilcoxon treatment effect shown where significant",
    y        = "Incidence rate (tumors / mouse / week)",
    x        = "Condition"
  )

# Add significance brackets where treatment effect is significant
if (nrow(sig_brackets) > 0) {
  p_ir_stats <- p_ir_stats +
    ggsignif::geom_signif(
      data        = sig_brackets,
      aes(xmin = group1, xmax = group2,
          y_position = max(mouse_ir$IR_mouse, na.rm = TRUE) * 1.05,
          annotations = stars),
      manual      = TRUE,
      textsize    = 4,
      vjust       = 0.3,
      tip_length  = 0.01,
      inherit.aes = FALSE
    )
}

ggsave(file.path(output_dir, "incidence_rate_per_mouse_stats.png"),
       p_ir_stats, width = 10, height = 6, dpi = 300, bg = "white")