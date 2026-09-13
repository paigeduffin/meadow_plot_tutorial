#!/usr/bin/env Rscript

# Command-line version of the meadow plot and candidate peak tutorial.
#
# Usage:
#   Rscript scripts/meadow_plot_cli.R PROJECT_DIR \
#     [FST_FILE] [PI_FILE] [COMPARISON_MAP_FILE] [CHROMOSOME_MAP_FILE] \
#     [MAX_INTERVENING_NONOUTLIERS]
#
# File arguments may be absolute paths or paths relative to PROJECT_DIR.

usage <- paste(
  "Usage:",
  "  Rscript scripts/meadow_plot_cli.R PROJECT_DIR",
  "    [FST_FILE] [PI_FILE] [COMPARISON_MAP_FILE] [CHROMOSOME_MAP_FILE]",
  "    [MAX_INTERVENING_NONOUTLIERS]",
  "",
  "Defaults:",
  "  FST_FILE: input/example_pixy_fst_data.csv",
  "  PI_FILE: input/example_pixy_pi_data.csv",
  "  COMPARISON_MAP_FILE: input/example_comparison_map.txt",
  "  CHROMOSOME_MAP_FILE: input/example_chromosome_map.txt",
  "  MAX_INTERVENING_NONOUTLIERS: 3",
  sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(status = if (length(args) == 0) 1 else 0)
}
if (length(args) > 6) stop(usage)

required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "tibble", "gridExtra",
  "ggtext", "scales", "gtable"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

project_dir <- normalizePath(args[1], mustWork = TRUE)

resolve_project_path <- function(path) {
  is_absolute <- grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)
  if (is_absolute) path else file.path(project_dir, path)
}

fst_file <- resolve_project_path(if (length(args) >= 2) args[2] else "input/example_pixy_fst_data.csv")
pi_file <- resolve_project_path(if (length(args) >= 3) args[3] else "input/example_pixy_pi_data.csv")
comparison_map_file <- resolve_project_path(if (length(args) >= 4) args[4] else "input/example_comparison_map.txt")
chromo_map_file <- resolve_project_path(if (length(args) >= 5) args[5] else "input/example_chromosome_map.txt")
max_intervening_nonoutliers <- if (length(args) >= 6) suppressWarnings(as.numeric(args[6])) else 3

if (length(max_intervening_nonoutliers) != 1 ||
    is.na(max_intervening_nonoutliers) ||
    max_intervening_nonoutliers < 0 ||
    max_intervening_nonoutliers %% 1 != 0) {
  stop("MAX_INTERVENING_NONOUTLIERS must be one non-negative whole number.")
}

wrangled_dir <- file.path(project_dir, "wrangled_data")
plot_dir <- file.path(project_dir, "plots")
dir.create(wrangled_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

input_files <- c(
  "Fst input" = fst_file,
  "Pi input" = pi_file,
  "Pairwise comparison map" = comparison_map_file,
  "Chromosome map" = chromo_map_file
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop(
    "The following input files were not found:\n",
    paste(names(missing_inputs), missing_inputs, sep = ": ", collapse = "\n")
  )
}

read_delimited <- function(path) {
  if (tolower(tools::file_ext(path)) == "csv") {
    read_csv(path, show_col_types = FALSE)
  } else {
    read_tsv(path, show_col_types = FALSE)
  }
}

message("Reading input files...")
fst_df <- read_delimited(fst_file)
pi_df <- read_delimited(pi_file)
comparison_map <- read_delimited(comparison_map_file)
chromosome_map <- read_delimited(chromo_map_file)

validate_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
}

validate_columns(
  fst_df,
  c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2", "avg_wc_fst"),
  "The Fst file"
)
validate_columns(
  pi_df,
  c("pop", "chromosome", "window_pos_1", "window_pos_2", "avg_pi"),
  "The pi file"
)
validate_columns(
  comparison_map,
  c("comp", "population_1", "population_2", "plot_order"),
  "The comparison map"
)
validate_columns(
  chromosome_map,
  c("chromosome", "chromo.num"),
  "The chromosome map"
)

if (anyDuplicated(comparison_map$comp)) {
  stop("Each value in comparison_map$comp must be unique.")
}
if (anyDuplicated(comparison_map$plot_order)) {
  stop("Each value in comparison_map$plot_order must be unique.")
}
if (anyDuplicated(chromosome_map$chromosome) || anyDuplicated(chromosome_map$chromo.num)) {
  stop("Each chromosome and chromo.num value in the chromosome map must be unique.")
}

chromosome_accessions <- chromosome_map %>%
  arrange(chromo.num) %>%
  pull(chromosome)
comps <- comparison_map$comp
plot_order <- comparison_map %>% arrange(plot_order) %>% pull(comp)

pair_key <- function(pop_a, pop_b) {
  paste(pmin(pop_a, pop_b), pmax(pop_a, pop_b), sep = "__")
}

comparison_map <- comparison_map %>%
  mutate(pair_key = pair_key(population_1, population_2))

fst_df_prepared <- fst_df %>%
  mutate(
    chromo.num = match(chromosome, chromosome_accessions),
    pair_key = pair_key(pop1, pop2)
  ) %>%
  filter(!is.na(chromo.num)) %>%
  inner_join(comparison_map %>% select(comp, pair_key), by = "pair_key") %>%
  select(chromosome, chromo.num, window_pos_1, window_pos_2, comp, avg_wc_fst)

missing_fst_comparisons <- setdiff(comps, unique(fst_df_prepared$comp))
if (length(missing_fst_comparisons) > 0) {
  stop("No Fst rows were found for: ", paste(missing_fst_comparisons, collapse = ", "))
}

pi_df_wide <- pi_df %>%
  select(pop, chromosome, window_pos_1, window_pos_2, avg_pi) %>%
  pivot_wider(names_from = pop, values_from = avg_pi)

expected_populations <- unique(c(comparison_map$population_1, comparison_map$population_2))
missing_pi_populations <- setdiff(expected_populations, names(pi_df_wide))
if (length(missing_pi_populations) > 0) {
  stop("The pi file contains no data for: ", paste(missing_pi_populations, collapse = ", "))
}

pi_df_ratios <- lapply(seq_len(nrow(comparison_map)), function(i) {
  comparison_name <- comparison_map$comp[[i]]
  population_1 <- comparison_map$population_1[[i]]
  population_2 <- comparison_map$population_2[[i]]

  pi_df_wide %>%
    transmute(
      chromosome, window_pos_1, window_pos_2,
      comp = comparison_name,
      population_1 = population_1,
      population_2 = population_2,
      pi_ratio = .data[[population_1]] / .data[[population_2]]
    )
}) %>%
  bind_rows() %>%
  left_join(chromosome_map, by = "chromosome") %>%
  mutate(log2_pi.ratio = log2(pi_ratio)) %>%
  filter(!is.na(chromo.num))

pi_fst_df_merge_clean <- pi_df_ratios %>%
  inner_join(
    fst_df_prepared,
    by = c("chromosome", "chromo.num", "window_pos_1", "window_pos_2", "comp")
  ) %>%
  filter(!is.na(avg_wc_fst), !is.na(log2_pi.ratio)) %>%
  mutate(
    chr_pos1_pos2_comp = paste(
      chromo.num, window_pos_1, window_pos_2, comp, sep = "_"
    )
  ) %>%
  select(
    chr_pos1_pos2_comp, chromo.num, window_pos_1, window_pos_2,
    comp, avg_wc_fst, pi_ratio, log2_pi.ratio
  )

finite_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  unname(quantile(x, probs = probability, na.rm = TRUE))
}

pi_fst_merge_no.neg.fst <- pi_fst_df_merge_clean %>%
  filter(avg_wc_fst >= 0)

global_pi_lower <- finite_quantile(pi_fst_merge_no.neg.fst$log2_pi.ratio, 0.01)
global_pi_upper <- finite_quantile(pi_fst_merge_no.neg.fst$log2_pi.ratio, 0.99)
global_fst_upper <- unname(quantile(
  pi_fst_merge_no.neg.fst$avg_wc_fst,
  probs = 0.99,
  na.rm = TRUE
))

pi_fst_merge_no.neg.fst <- pi_fst_merge_no.neg.fst %>%
  mutate(
    outlier.labels_GLOBAL = if_else(
      avg_wc_fst >= global_fst_upper &
        (log2_pi.ratio <= global_pi_lower | log2_pi.ratio >= global_pi_upper),
      "outlier",
      "not.outlier"
    )
  )

global_outlier_summary <- pi_fst_merge_no.neg.fst %>%
  group_by(comp) %>%
  summarise(
    n_outliers = sum(outlier.labels_GLOBAL == "outlier"),
    n_total = n(),
    proport_outliers = n_outliers / n_total,
    .groups = "drop"
  ) %>%
  arrange(match(comp, plot_order))

pi_fst_merge_no.neg.fst_merged_individ.comp <- pi_fst_merge_no.neg.fst %>%
  group_by(comp) %>%
  mutate(
    pi_lower.thresh = finite_quantile(log2_pi.ratio, 0.01),
    pi_high.thresh = finite_quantile(log2_pi.ratio, 0.99),
    fst_thresh = unname(quantile(avg_wc_fst, probs = 0.99, na.rm = TRUE)),
    outlier.labels_individ.comp = if_else(
      avg_wc_fst >= fst_thresh &
        (log2_pi.ratio <= pi_lower.thresh | log2_pi.ratio >= pi_high.thresh),
      "outlier",
      "not.outlier"
    )
  ) %>%
  ungroup()

individual_outlier_summary <- pi_fst_merge_no.neg.fst_merged_individ.comp %>%
  group_by(comp) %>%
  summarise(
    n_outliers = sum(outlier.labels_individ.comp == "outlier"),
    n_total = n(),
    proport_outliers = n_outliers / n_total,
    .groups = "drop"
  ) %>%
  arrange(match(comp, plot_order))

main_df_clean <- pi_fst_merge_no.neg.fst_merged_individ.comp %>%
  mutate(
    log2_pi.ratio_inf2NA = if_else(is.finite(log2_pi.ratio), log2_pi.ratio, NA_real_),
    outlier.labels_inf.labels = case_when(
      outlier.labels_individ.comp == "outlier" & log2_pi.ratio == -Inf ~ "outlier_pop1.pi.0",
      outlier.labels_individ.comp == "outlier" & log2_pi.ratio == Inf ~ "outlier_pop2.pi.0",
      outlier.labels_individ.comp == "outlier" ~ "outlier_not.inf",
      log2_pi.ratio == -Inf ~ "not.out_pop1.pi.0",
      log2_pi.ratio == Inf ~ "not.out_pop2.pi.0",
      TRUE ~ "not.out_not.inf"
    ),
    alpha.for.asterisk = case_when(
      outlier.labels_inf.labels %in% c("outlier_pop1.pi.0", "outlier_pop2.pi.0") ~ 1,
      outlier.labels_inf.labels %in% c("not.out_pop1.pi.0", "not.out_pop2.pi.0") ~ 0.6,
      TRUE ~ 0
    ),
    alpha.for.main.pts = case_when(
      outlier.labels_inf.labels == "outlier_not.inf" ~ 1,
      outlier.labels_inf.labels == "not.out_not.inf" ~ 0.6,
      TRUE ~ 0
    ),
    outlier.labels_individ.comp = factor(
      outlier.labels_individ.comp,
      levels = c("not.outlier", "outlier")
    )
  )

outlier_counts <- main_df_clean %>%
  filter(outlier.labels_individ.comp == "outlier") %>%
  mutate(
    outlier_population = if_else(log2_pi.ratio < 0, "population_1", "population_2")
  ) %>%
  count(comp, outlier_population, name = "n_outlier_windows") %>%
  pivot_wider(
    names_from = outlier_population,
    values_from = n_outlier_windows,
    values_fill = 0,
    names_prefix = "n_outliers_"
  ) %>%
  rename(comparison = comp)

message("Saving prepared data and outlier summaries...")
write_csv(
  pi_fst_merge_no.neg.fst_merged_individ.comp,
  file.path(wrangled_dir, "pi_fst_outlier_windows_with_infinite_values.csv")
)
write_csv(main_df_clean, file.path(wrangled_dir, "pi_fst_outlier_windows_final.csv"))
write_csv(
  main_df_clean %>%
    filter(!is.finite(log2_pi.ratio)) %>%
    select(chr_pos1_pos2_comp),
  file.path(wrangled_dir, "infinite_pi_ratio_windows.csv")
)
write_csv(outlier_counts, file.path(wrangled_dir, "outlier_counts_by_comparison.csv"))
write_csv(global_outlier_summary, file.path(wrangled_dir, "global_outlier_summary.csv"))
write_csv(
  individual_outlier_summary,
  file.path(wrangled_dir, "comparison_specific_outlier_summary.csv")
)

transparent_theme <- theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.position = "none"
)

auto_pi_scale <- function(x, title = "Legend", n = 50) {
  x <- x[is.finite(x)]
  limits <- if (length(x) == 0) c(-1, 1) else range(c(x, 0), na.rm = TRUE)
  if (diff(limits) == 0) limits <- limits + c(-1e-9, 1e-9)

  negative_colors <- c(
    "#3E207B", "#4C3194", "#5A42AD", "#5E63CC", "#5970D2",
    "#5C89DC", "#94BBEC", "#A3CCD0", "#9CCD9E", "#92D050"
  )
  positive_colors <- c(
    "#92D050", "#CAD961", "#FFE171", "#FACA42", "#F3A532",
    "#EB8E2D", "#E2662A", "#BE3D28", "#B12F27", "#752A24"
  )

  zero_position <- scales::rescale(0, from = limits)
  colors <- c(
    grDevices::colorRampPalette(negative_colors)(ceiling(n / 2)),
    grDevices::colorRampPalette(positive_colors)(ceiling(n / 2))[-1]
  )
  values <- c(
    seq(0, zero_position, length.out = ceiling(n / 2)),
    seq(zero_position, 1, length.out = ceiling(n / 2))[-1]
  )

  scale_color_gradientn(
    name = title,
    colors = colors,
    values = values,
    limits = limits,
    oob = scales::squish
  )
}

make_meadow_components <- function(comparison_name) {
  df <- main_df_clean %>% filter(comp == comparison_name)
  y_max <- max(df$avg_wc_fst, na.rm = TRUE)
  y_limit <- if (y_max > 0) y_max * 1.05 else 1

  pop_names <- strsplit(comparison_name, ".v.", fixed = TRUE)[[1]]
  if (length(pop_names) != 2) {
    stop("Comparison names must follow population_1.v.population_2: ", comparison_name)
  }
  legend_title <- paste0(
    "\u03c0 ratio<br>",
    "<span style='color:#E44E20; font-weight:bold'>", pop_names[2], "</span><br>",
    "vs<br>",
    "<span style='color:#545FBD; font-weight:bold'>", pop_names[1], "</span>"
  )

  common_facet <- facet_grid(
    ~chromo.num,
    scales = "free_x",
    space = "free_x",
    switch = "x"
  )
  common_theme <- theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.spacing.x = grid::unit(0.03, "lines"),
    plot.margin = margin(t = 2.75, r = 5.5, b = 2.75, l = 5.5, unit = "pt")
  )

  main_points <- ggplot(df, aes(window_pos_1, avg_wc_fst)) +
    geom_point(
      data = filter(df, is.finite(log2_pi.ratio)),
      aes(
        color = log2_pi.ratio,
        size = outlier.labels_individ.comp,
        alpha = alpha.for.main.pts
      )
    ) +
    geom_point(
      data = filter(df, log2_pi.ratio == -Inf),
      aes(size = outlier.labels_individ.comp, alpha = alpha.for.asterisk),
      shape = 8,
      color = "#6579D7"
    ) +
    geom_point(
      data = filter(df, log2_pi.ratio == Inf),
      aes(size = outlier.labels_individ.comp, alpha = alpha.for.asterisk),
      shape = 8,
      color = "#D55C2C"
    ) +
    auto_pi_scale(df$log2_pi.ratio, title = legend_title) +
    scale_alpha_identity() +
    scale_size_manual(
      values = c("not.outlier" = 0.3, "outlier" = 1.6),
      guide = "none"
    ) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0, 0.01))
    ) +
    common_facet +
    transparent_theme +
    common_theme +
    theme(
      legend.position = "right",
      legend.title = ggtext::element_markdown(hjust = 0, size = 7, lineheight = 1),
      legend.text = element_text(
        size = 7,
        color = "black",
        margin = margin(l = 3, unit = "pt")
      ),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.height = grid::unit(0.8, "cm"),
      legend.key.width = grid::unit(0.3, "cm")
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(0.4, "cm"),
        frame.colour = "black",
        frame.linewidth = 0.3,
        ticks.colour = "black",
        ticks.linewidth = 0.15
      )
    )

  list(main_points = main_points)
}

message("Creating meadow plots...")
meadow_plots <- setNames(lapply(comps, make_meadow_components), comps)

add_combined_panel_border <- function(plot) {
  grob <- ggplotGrob(plot)
  panels <- grob$layout[grepl("^panel", grob$layout$name), ]

  gtable::gtable_add_grob(
    grob,
    grid::rectGrob(gp = grid::gpar(color = "black", fill = NA, lwd = 1)),
    t = min(panels$t),
    b = max(panels$b),
    l = min(panels$l),
    r = max(panels$r),
    z = Inf,
    clip = "off"
  )
}

main_point_panels <- lapply(seq_along(plot_order), function(i) {
  plot <- meadow_plots[[plot_order[i]]]$main_points + labs(x = "Chromosome")

  if (i < length(plot_order)) {
    plot <- plot + theme(
      strip.background = element_blank(),
      strip.text.x = element_text(color = "transparent", size = 1),
      axis.title.x = element_text(color = "transparent", size = 1)
    )
  } else {
    plot <- plot + theme(
      strip.background = element_blank(),
      strip.text.x = element_text(),
      axis.title.x = element_text(size = 9, color = "black")
    )
  }

  add_combined_panel_border(plot)
})

window_size_kb <- median(fst_df$window_pos_2 - fst_df$window_pos_1 + 1, na.rm = TRUE) / 1000
window_steps <- fst_df %>%
  distinct(chromosome, window_pos_1) %>%
  arrange(chromosome, window_pos_1) %>%
  group_by(chromosome) %>%
  mutate(step = window_pos_1 - lag(window_pos_1)) %>%
  ungroup() %>%
  filter(!is.na(step), step > 0) %>%
  pull(step)
window_step_kb <- median(window_steps, na.rm = TRUE) / 1000

format_kb <- function(x) format(x, trim = TRUE, scientific = FALSE)
shared_y_title <- grid::textGrob(
  bquote(paste(
    "Average ", italic(F)[italic(ST)], " (",
    .(format_kb(window_size_kb)), "-kb windows, ",
    .(format_kb(window_step_kb)), "-kb step)"
  )),
  rot = 90,
  gp = grid::gpar(fontsize = 9)
)

meadow_plots_combined <- gridExtra::arrangeGrob(
  grobs = main_point_panels,
  ncol = 1,
  left = shared_y_title
)

ggsave(
  file.path(plot_dir, "meadow_plots.png"),
  meadow_plots_combined,
  width = 20,
  height = min(8 * nrow(comparison_map), 60),
  units = "cm",
  dpi = 600,
  bg = "transparent"
)

message("Grouping candidate outlier windows into peaks...")
indexed_windows <- main_df_clean %>%
  arrange(comp, chromo.num, window_pos_1, window_pos_2) %>%
  group_by(comp, chromo.num) %>%
  mutate(retained_window_index = row_number()) %>%
  ungroup()

candidate_peaks <- indexed_windows %>%
  filter(outlier.labels_individ.comp == "outlier") %>%
  group_by(comp, chromo.num) %>%
  arrange(retained_window_index, .by_group = TRUE) %>%
  mutate(
    new_peak = is.na(lag(retained_window_index)) |
      retained_window_index - lag(retained_window_index) >
        max_intervening_nonoutliers + 1,
    peak_number = cumsum(new_peak)
  ) %>%
  group_by(comp, chromo.num, peak_number) %>%
  summarise(
    first_window_start = min(window_pos_1),
    last_window_end = max(window_pos_2),
    n_outlier_windows = n(),
    n_intervening_retained_nonoutliers =
      max(retained_window_index) - min(retained_window_index) + 1 - n(),
    total_genomic_span_kb =
      (max(window_pos_2) - min(window_pos_1) + 1) / 1000,
    .groups = "drop"
  ) %>%
  arrange(match(comp, plot_order), chromo.num, first_window_start)

if (nrow(candidate_peaks) == 0) {
  stop("No candidate outlier windows were available for peak analysis.")
}

comparison_peak_summary <- candidate_peaks %>%
  group_by(comp) %>%
  summarise(
    n_candidate_peaks = n(),
    n_outlier_windows = sum(n_outlier_windows),
    mean_outlier_windows_per_peak = mean(n_outlier_windows),
    max_outlier_windows_per_peak = max(n_outlier_windows),
    .groups = "drop"
  ) %>%
  arrange(match(comp, plot_order))

peak_size_summary <- candidate_peaks %>%
  count(comp, n_outlier_windows, name = "n_candidate_peaks") %>%
  arrange(match(comp, plot_order), n_outlier_windows)

write_csv(candidate_peaks, file.path(wrangled_dir, "candidate_peaks.csv"))
write_csv(
  comparison_peak_summary,
  file.path(wrangled_dir, "candidate_peak_summary_by_comparison.csv")
)
write_csv(
  peak_size_summary,
  file.path(wrangled_dir, "candidate_peak_size_counts.csv")
)

candidate_peaks_for_plot <- candidate_peaks %>%
  mutate(comp = factor(comp, levels = plot_order))

candidate_peak_size_plot <- ggplot(
  candidate_peaks_for_plot,
  aes(comp, n_outlier_windows, group = comp)
) +
  geom_violin(
    width = 1,
    color = "#797979",
    fill = "#797979",
    alpha = 0.25
  ) +
  geom_jitter(
    aes(size = n_outlier_windows),
    color = "black",
    fill = "#797979",
    alpha = 0.9,
    width = 0.25,
    shape = 21
  ) +
  scale_size(range = c(1.5, 3), guide = "none") +
  scale_x_discrete(
    drop = FALSE,
    labels = function(x) gsub(".v.", "\nvs\n", x, fixed = TRUE)
  ) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  labs(
    x = "Population comparison",
    y = "Outlier windows per candidate peak"
  ) +
  theme_classic() +
  transparent_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  file.path(plot_dir, "candidate_peak_size_distributions.png"),
  candidate_peak_size_plot,
  width = min(5 * nrow(comparison_map), 60),
  height = 15,
  units = "cm",
  dpi = 600,
  bg = "transparent"
)

writeLines(capture.output(sessionInfo()), file.path(wrangled_dir, "session_info.txt"))

message("Analysis complete.")
message("Prepared data: ", wrangled_dir)
message("Plots: ", plot_dir)
