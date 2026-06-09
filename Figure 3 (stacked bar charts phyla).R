
# ---- Packages --------------------------------------------------------------
# install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ---- Paths -----------------------------------------------------------------
DATA_DIR <- "/Users/mfoley2/Desktop/msystems resubmit/taxonomic groups"
OUT_DIR  <- DATA_DIR

# Bacterial sources
bact_tax_file <- file.path(DATA_DIR, "16S_qsip_results_18O_outliers_removed_with_taxonomy.csv")
bact_strategy <- list(
  "Rhizosphere"             = "rhizosphere.microbes.csv",
  "Detritusphere"           = "detritusphere.microbes.csv",
  "Rhizosphere + detritus"  = "rhizosphere.detritusphere.microbes.types.csv"
)

# Fungal sources
fung_tax_file <- file.path(DATA_DIR, "ITS_qsip_results_18O_outliers_removed_with_taxonomy.csv")
fung_strategy <- list(
  "Rhizosphere"             = "rhizosphere.taxa.ITS.csv",
  "Detritusphere"           = "detritusphere.taxa.ITS.csv",
  "Rhizosphere + detritus"  = "rhizosphere.detritusphere.taxa.ITS.csv"
)

# Phyla with fewer than MIN_ASVS unique ASVs are dropped from each panel.
MIN_ASVS <- 5

# ---- Load + clean taxonomies ----------------------------------------------
bact_tax <- read_csv(bact_tax_file, show_col_types = FALSE)
bact_tax <- bact_tax[, c("feature_id", "phylum")]
bact_tax <- bact_tax[!duplicated(bact_tax$feature_id), ]
names(bact_tax)[2] <- "phylum"

fung_tax <- read_csv(fung_tax_file, show_col_types = FALSE)
fung_tax <- fung_tax[, c("feature_id", "Phylum")]
fung_tax <- fung_tax[!duplicated(fung_tax$feature_id), ]
fung_tax$Phylum <- sub("^p__", "", fung_tax$Phylum)
fung_tax$Phylum[is.na(fung_tax$Phylum) | fung_tax$Phylum == "" |
                fung_tax$Phylum == "unassigned"] <- "unassigned"
names(fung_tax)[2] <- "phylum"

# ---- Per-panel summariser --------------------------------------------------
summarise_panel <- function(strategy_path, tax_df, min_asvs = MIN_ASVS) {
  strat <- read_csv(strategy_path, show_col_types = FALSE)
  strat <- strat[, c("feature_id", "type")]

  df <- merge(strat, tax_df, by = "feature_id")

  phylum_counts <- as.data.frame(table(df$phylum), stringsAsFactors = FALSE)
  names(phylum_counts) <- c("phylum", "n_asv")
  keep_phyla <- phylum_counts$phylum[phylum_counts$n_asv >= min_asvs]
  df <- df[df$phylum %in% keep_phyla, , drop = FALSE]

  all_phyla <- unique(df$phylum)
  all_types <- c("specialist", "generalist")
  grid <- expand.grid(phylum = all_phyla, type = all_types,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

  obs <- as.data.frame(table(df$phylum, df$type), stringsAsFactors = FALSE)
  names(obs) <- c("phylum", "type", "n_asv")
  out <- merge(grid, obs, by = c("phylum", "type"), all.x = TRUE)
  out$n_asv[is.na(out$n_asv)] <- 0

  totals <- aggregate(n_asv ~ phylum, data = out, FUN = sum)
  names(totals)[2] <- "total_asv"
  out <- merge(out, totals, by = "phylum", all.x = TRUE)
  out$fraction <- out$n_asv / out$total_asv

  totals_ord <- totals[order(-totals$total_asv), ]
  out$phylum <- factor(out$phylum, levels = totals_ord$phylum)
  out$type   <- factor(out$type,   levels = c("specialist", "generalist"))

  list(summary = out, totals = totals[order(-totals$total_asv), ])
}

# ---- Per-panel plot builder -----------------------------------------------
build_panel <- function(summary_df, totals_df, panel_title,
                        show_y_left = TRUE, show_y_right = TRUE,
                        show_x_text = TRUE) {

  max_count <- max(totals_df$total_asv)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1
  count_df  <- totals_df
  count_df$y_scaled <- count_df$total_asv / max_count
  count_df$phylum   <- factor(count_df$phylum, levels = levels(summary_df$phylum))

  p <- ggplot(summary_df, aes(x = phylum, y = fraction, fill = type)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.3) +
    geom_line(
      data        = count_df,
      mapping     = aes(x = phylum, y = y_scaled, group = 1),
      inherit.aes = FALSE,
      linetype    = "solid",
      colour      = "black",
      linewidth   = 0.8
    ) +
    geom_point(
      data        = count_df,
      mapping     = aes(x = phylum, y = y_scaled),
      inherit.aes = FALSE,
      shape       = 21,
      size        = 3.6,
      fill        = "black",
      colour      = "white",
      stroke      = 0.9
    ) +
    scale_y_continuous(
      name     = if (show_y_left)  "Fraction of ASVs" else NULL,
      labels   = percent_format(accuracy = 1),
      expand   = expansion(mult = c(0, 0.02)),
      sec.axis = sec_axis(
        trans = ~ . * max_count,
        name  = if (show_y_right) "Total ASVs" else NULL
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(
      values = c(specialist = "#D55E00", generalist = "#0072B2"),
      name   = "Strategy"
    ) +
    labs(title = panel_title, x = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x         = if (show_x_text)
                              element_text(angle = 45, hjust = 1, size = 11,
                                           colour = "black")
                            else element_blank(),
      axis.text.y         = element_text(size = 11, colour = "black"),
      axis.title.y        = element_text(size = 12),
      axis.title.y.right  = element_text(size = 12),
      axis.ticks.x        = if (show_x_text) element_line() else element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.y  = element_line(colour = "grey90", linewidth = 0.3),
      legend.position     = "none",
      plot.title          = element_text(size = 13, face = "bold"),
      plot.title.position = "plot"
    )

  p
}

# ---- Build all 6 panels ---------------------------------------------------
habitats <- c("Rhizosphere", "Detritusphere", "Rhizosphere + detritus")

panels <- list()
for (h in habitats) {
  for (king in c("Bacteria", "Fungi")) {
    if (king == "Bacteria") {
      s <- summarise_panel(file.path(DATA_DIR, bact_strategy[[h]]), bact_tax)
    } else {
      s <- summarise_panel(file.path(DATA_DIR, fung_strategy[[h]]), fung_tax)
    }
    panels[[paste(h, king, sep = " | ")]] <- build_panel(
      summary_df   = s$summary,
      totals_df    = s$totals,
      panel_title  = paste0(king, " — ", h),
      show_y_left  = (king == "Bacteria"),
      show_y_right = (king == "Fungi"),
      show_x_text  = TRUE  # keep on every row; phyla differ between rows
    )
  }
}

# ---- Assemble with patchwork ----------------------------------------------
# Row 1: rhizosphere (bact, fungi)
# Row 2: detritusphere (bact, fungi)
# Row 3: rhizo + detritus (bact, fungi)
# Bacteria panels have far more phyla than fungi (~14 vs ~5), so make the
# left column wider.

combined <- (panels[["Rhizosphere | Bacteria"]]            | panels[["Rhizosphere | Fungi"]]            ) /
            (panels[["Detritusphere | Bacteria"]]          | panels[["Detritusphere | Fungi"]]          ) /
            (panels[["Rhizosphere + detritus | Bacteria"]] | panels[["Rhizosphere + detritus | Fungi"]] )

# A single shared legend across the whole figure
combined <- combined +
  plot_layout(widths = c(2.4, 1), guides = "collect") +
  plot_annotation(
    title    = "Strategy composition by phylum",
    subtitle = sprintf("Phyla with < %d ASVs excluded. Stacked bars: generalist vs. specialist fraction. Dots/line: total ASVs (right axis).", MIN_ASVS),
    theme    = theme(plot.title    = element_text(face = "bold", size = 16),
                     plot.subtitle = element_text(size = 11, colour = "grey25"))
  ) &
  theme(
    legend.position = "top",
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

# Smaller canvas + bigger fonts so the figure is legible when Word scales it
# down to ~6.5 inches wide.
out_path <- file.path(OUT_DIR, "strategy_composition_combined.png")
ggsave(out_path, combined,
       width = 9.5, height = 10, dpi = 400, bg = "white")

message("Saved: ", out_path)

# ============================================================
# Also export each row separately (for Affinity Designer / manual layout)
# Same fonts/points/line as the combined figure.
# ============================================================

row_files <- list(
  "Rhizosphere"            = "strategy_row_rhizosphere.png",
  "Detritusphere"          = "strategy_row_detritusphere.png",
  "Rhizosphere + detritus" = "strategy_row_rhizo_plus_detritus.png"
)

for (h in habitats) {
  row_plot <- (panels[[paste(h, "Bacteria", sep = " | ")]] |
               panels[[paste(h, "Fungi",    sep = " | ")]]) +
    plot_layout(widths = c(2.4, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 12, face = "bold"),
      legend.key.size = unit(0.6, "cm")
    )

  out_row <- file.path(OUT_DIR, row_files[[h]])
  ggsave(out_row, row_plot,
         width = 9.5, height = 4.2, dpi = 400, bg = "white")
  message("Saved row: ", out_row)
}
