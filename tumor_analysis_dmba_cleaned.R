suppressPackageStartupMessages({
  library(readr)     # 2.1.4
  library(dplyr)     # 1.1.3
  library(tidyr)     # 1.3.0
  library(stringr)   # 1.5.0
  library(ggplot2)   # 3.4.4
  library(readxl)    # 1.4.3
  library(ggh4x)     # 0.2.6
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
input_file <- "C:/Users/s.lafirenze/Desktop/Skin paper/Figure dmba_CIN/DMBA_CIN experiment.xlsx"
output_dir <- "C:/Users/s.lafirenze/Desktop/Skin paper/Figure dmba_CIN/DMBA_output"


if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 2. SETTINGS
# ============================================================
genotype_levels <- c("WT/WT", "WT/TA", "KD/WT", "TA/TA", "TA/KD", "KD/KD")
size_levels     <- c("<10", "10-100", "100-500", ">500")

genotype_palette <- c(
  "WT/WT" = "grey4",
  "WT/TA" = "grey46",
  "KD/WT" = "grey78",
  "TA/TA" = "goldenrod3",
  "TA/KD" = "darkorange3",
  "KD/KD" = "red2"
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
# Input: one Excel sheet where rows are timepoints and columns are individual tumors.
# Column 1 = timepoint (day), remaining columns = tumor sizes in mm³.
# Column headers encode mouse ID and genotype (e.g. "7080001_WT/WT").
df_raw <- readxl::read_excel(input_file, col_names = TRUE, .name_repair = "minimal")
colnames(df_raw) <- colnames(df_raw) %>% str_trim()
colnames(df_raw)[1] <- "Timepoint"
tumor_cols         <- colnames(df_raw)[-1]
colnames(df_raw)   <- c("Timepoint", make.unique(tumor_cols, sep = "_tumor"))

df_raw <- df_raw %>%
  mutate(
    Timepoint = as.numeric(Timepoint),
    Week      = paste0("week ", ceiling(Timepoint / 7))
  ) %>%
  dplyr::select(Week, Timepoint, everything()) %>%
  mutate(across(-c(1, 2), ~ as.numeric(.x)))

message("Columns: ", ncol(df_raw), " | Rows: ", nrow(df_raw))

# ============================================================
# 4. LONG FORMAT + MAP TUMORS TO MOUSE / GENOTYPE
# ============================================================
df_long <- df_raw %>%
  pivot_longer(cols = -c(1, 2), names_to = "Tumor_Column", values_to = "Tumor_Size")

tumor_map <- tibble(Tumor_Column = unique(df_long$Tumor_Column)) %>%
  mutate(
    Tumor_Column_Base = Tumor_Column %>%
      str_remove("_tumor[0-9]*$") %>%
      str_squish(),
    Mouse_ID = str_extract(Tumor_Column_Base, "^[^_ ]+"),
    Genotype = str_extract(Tumor_Column_Base, "WT/WT|WT/TA|W/TA|KD/WT|TA/TA|TA/KD|KD/KD") %>%
      str_replace("^W/TA$", "WT/TA")   # fix occasional typo in column headers
  ) %>%
  group_by(Mouse_ID, Genotype) %>%
  mutate(Tumor_Number = row_number()) %>%
  ungroup()

df_long <- df_long %>%
  mutate(
    Tumor_Column_Base = Tumor_Column %>%
      str_remove("_tumor[0-9]*$") %>%
      str_squish()
  ) %>%
  left_join(tumor_map, by = c("Tumor_Column", "Tumor_Column_Base")) %>%
  dplyr::select(Week, Timepoint, Mouse_ID, Genotype, Tumor_Number, Tumor_Size) %>%
  mutate(
    Genotype      = factor(Genotype, levels = genotype_levels),
    Size_Category = case_when(
      is.na(Tumor_Size) | Tumor_Size <= 0 ~ NA_character_,
      Tumor_Size < 10   ~ "<10",
      Tumor_Size < 100  ~ "10-100",
      Tumor_Size < 500  ~ "100-500",
      TRUE              ~ ">500"
    ),
    Size_Category = factor(Size_Category, levels = size_levels),
    Week_Num      = readr::parse_number(Week),
    Tumor_ID      = paste(Mouse_ID, Genotype, Tumor_Number, sep = "_")
  )

readr::write_csv(df_long, file.path(output_dir, "tumor_sizes_long.csv"))

# ============================================================
# 4b. FOLLOW-UP WINDOW PER MOUSE
#     True_Last_Timepoint  = last day with any measurement
#     True_First_Timepoint = first day with any measurement
# ============================================================
mouse_true_last <- df_long %>%
  filter(!is.na(Mouse_ID), !is.na(Genotype), !is.na(Timepoint)) %>%
  group_by(Mouse_ID, Genotype) %>%
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
    Genotype      = factor(Genotype, levels = genotype_levels)
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
  summarise(
    Onset_Day  = min(Timepoint, na.rm = TRUE),
    Onset_Week = min(Week_Num,  na.rm = TRUE),
    .groups    = "drop"
  )

df_long <- df_long %>%
  left_join(tumor_onset, by = "Tumor_ID") %>%
  mutate(
    Days_From_Onset  = Timepoint - Onset_Day,
    Weeks_From_Onset = Week_Num  - Onset_Week
  )

readr::write_csv(df_long, file.path(output_dir, "tumor_sizes_long_with_onset.csv"))

# ============================================================
# 6. AVERAGE TUMOR COUNTS OVER TIME (absolute days)
#    Zero-fill: all mouse-timepoint combinations (includes mice
#    with 0 tumors at any given timepoint)
# ============================================================
all_mouse_days <- df_long %>%
  filter(!is.na(Genotype), !is.na(Mouse_ID)) %>%
  distinct(Timepoint, Genotype, Mouse_ID)

all_mouse_day_size <- tidyr::crossing(
  all_mouse_days,
  tibble(Size_Category = factor(size_levels, levels = size_levels))
)

tumor_counts_mouse_size <- df_long %>%
  filter(!is.na(Genotype), !is.na(Mouse_ID), !is.na(Size_Category), Tumor_Size > 0) %>%
  group_by(Timepoint, Genotype, Mouse_ID, Size_Category) %>%
  summarise(Tumor_Count = n(), .groups = "drop")

tumor_counts_mouse_size <- all_mouse_day_size %>%
  left_join(tumor_counts_mouse_size,
            by = c("Timepoint", "Genotype", "Mouse_ID", "Size_Category")) %>%
  mutate(Tumor_Count = ifelse(is.na(Tumor_Count), 0, Tumor_Count))

avg_tumor_counts_genotype_size <- tumor_counts_mouse_size %>%
  group_by(Timepoint, Genotype, Size_Category) %>%
  summarise(
    Mean_Tumor_Count = mean(Tumor_Count, na.rm = TRUE),
    SD_Tumor_Count   = sd(Tumor_Count,   na.rm = TRUE),
    N_Mice           = n(),
    SE_Tumor_Count   = SD_Tumor_Count / sqrt(N_Mice),
    .groups          = "drop"
  ) %>%
  arrange(Timepoint, Genotype, Size_Category)

readr::write_csv(
  avg_tumor_counts_genotype_size,
  file.path(output_dir, "average_tumor_counts_per_day_genotype_sizeclass_including_zeros.csv")
)

day_breaks         <- sort(unique(avg_tumor_counts_genotype_size$Timepoint))
day_breaks_to_show <- day_breaks[
  seq(1, length(day_breaks), by = max(1, floor(length(day_breaks) / 10)))
]
bar_width_days <- min(diff(sort(unique(avg_tumor_counts_genotype_size$Timepoint)))) * 3

# Figure: average tumor counts over absolute time (lower panel, used in paper)
p_avg_counts_bar_lower <- ggplot(
  avg_tumor_counts_genotype_size,
  aes(x = Timepoint, y = Mean_Tumor_Count, fill = Size_Category)
) +
  geom_col(width = bar_width_days, position = position_stack(reverse = TRUE), na.rm = TRUE) +
  ggh4x::facet_wrap2(~ Genotype, ncol = 3, scales = "fixed",
                     axes = "all", remove_labels = "y") +
  scale_fill_manual(values = size_palette, breaks = size_levels, drop = FALSE) +
  scale_x_continuous(breaks = day_breaks_to_show, labels = day_breaks_to_show) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(size = 8, angle = 45, hjust = 1),
    axis.ticks.x     = element_line(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(0.2, "lines"),
    panel.spacing.y  = unit(0.2, "lines")
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