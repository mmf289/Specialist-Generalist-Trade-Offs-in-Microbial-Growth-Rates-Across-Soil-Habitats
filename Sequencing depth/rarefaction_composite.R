# =============================================================================
# Rarefaction curves on composite (tube-level) samples
#
# Builds composite samples by SUMMING ASV counts across all fractions within
# each qSIP tube (sipID), then computes rarefaction curves and plots
# mean +/- SD per habitat group for Normal-moisture samples only.
#
# Habitats:  Rhizo  /  Detritus  /  Rhizo_Detritus
# Amplicons: 16S    /  ITS  (processed independently)
# =============================================================================

# -------------------- USER: edit these paths if needed -----------------------
data_dir <- "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/calculate sequencing depth"
out_dir  <- "/Users/mfoley2/Documents/Claude/Projects/Sequencing depth/rarefaction_output"

asv_paths <- list(
  "16S" = file.path(data_dir, "DRIP16S_table_2x2.txt"),
  "ITS" = file.path(data_dir, "DRIPITS_table_2x2.txt")
)
metadata_path <- file.path(data_dir, "18O.lab.data.density.corrected.csv")

# Rarefaction step size (reads). Smaller = smoother curve, slower.
rar_step  <- 2000
# Require at least this many tubes contributing at a given depth to plot the
# group mean/SD there (curves trim back as fewer tubes have that depth).
min_n_for_group <- 3

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Packages ----------------------------------------------
required <- c("vegan", "ggplot2", "dplyr", "tidyr", "scales")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# -------------------- Load metadata ------------------------------------------
md <- read.csv(metadata_path, stringsAsFactors = FALSE)
md <- md[md$water.trt == "Normal", ]
md$tube <- md$sipID

cat("Normal-moisture fractions in metadata:", nrow(md), "\n")
cat("Unique tubes per habitat (metadata):\n")
print(table(unique(md[, c("tube","c.trt")])$c.trt))

# -------------------- Helper: build composite (tube-level) table ------------
build_composites <- function(asv_path, md) {
  message("Reading ", basename(asv_path), " ...")
  asv <- read.table(asv_path, header = TRUE, sep = "\t",
                    check.names = FALSE, row.names = 1,
                    comment.char = "", quote = "")

  shared <- intersect(colnames(asv), md$sample.name)
  asv    <- asv[, shared, drop = FALSE]
  md_sub <- md[md$sample.name %in% shared, ]

  sample2tube <- setNames(md_sub$tube, md_sub$sample.name)
  tubes <- sort(unique(sample2tube))

  # Sum fraction columns within each tube
  comp <- sapply(tubes, function(t) {
    cols <- names(sample2tube)[sample2tube == t]
    cols <- intersect(cols, colnames(asv))
    if (length(cols) == 1) asv[, cols] else rowSums(asv[, cols, drop = FALSE])
  })
  comp <- t(comp)                              # rows = tubes, cols = ASVs
  rownames(comp) <- as.character(tubes)
  storage.mode(comp) <- "integer"

  # Habitat lookup per tube
  tube2hab <- unique(md_sub[, c("tube", "c.trt")])
  habitat  <- setNames(tube2hab$c.trt, as.character(tube2hab$tube))[rownames(comp)]

  list(mat = comp, habitat = habitat)
}

# -------------------- Helper: compute rarefaction curves --------------------
rarefaction_curves <- function(comp_mat, step = 2000) {
  totals <- rowSums(comp_mat)
  global_max <- max(totals)
  sizes_all  <- unique(c(1, seq(step, global_max, by = step), global_max))

  do.call(rbind, lapply(seq_len(nrow(comp_mat)), function(i) {
    s <- sizes_all[sizes_all <= totals[i]]
    if (length(s) == 0) return(NULL)
    rich <- vegan::rarefy(comp_mat[i, ], sample = s)
    data.frame(
      tube     = rownames(comp_mat)[i],
      depth    = s,
      richness = as.numeric(rich),
      stringsAsFactors = FALSE
    )
  }))
}

# -------------------- Helper: group summary (mean +/- SD) -------------------
summarize_curves <- function(curves, habitat_map, min_n = 3) {
  curves$habitat <- habitat_map[as.character(curves$tube)]
  curves %>%
    group_by(habitat, depth) %>%
    summarise(
      n             = dplyr::n(),
      mean_richness = mean(richness),
      sd_richness   = stats::sd(richness),
      .groups       = "drop"
    ) %>%
    dplyr::filter(n >= min_n)
}

# -------------------- Helper: plot ------------------------------------------
plot_rarefaction <- function(summary_df, per_tube_df, habitat_map, title) {
  per_tube_df <- per_tube_df
  per_tube_df$habitat <- habitat_map[as.character(per_tube_df$tube)]

  ggplot() +
    # Individual tube curves: thin background lines
    geom_line(data = per_tube_df,
              aes(x = depth, y = richness, group = tube, color = habitat),
              alpha = 0.25, linewidth = 0.3) +
    # Group mean +/- SD ribbon
    geom_ribbon(data = summary_df,
                aes(x = depth,
                    ymin = pmax(0, mean_richness - sd_richness),
                    ymax = mean_richness + sd_richness,
                    fill = habitat),
                alpha = 0.25, color = NA) +
    # Group mean line
    geom_line(data = summary_df,
              aes(x = depth, y = mean_richness, color = habitat),
              linewidth = 1.1) +
    scale_x_continuous(labels = scales::comma) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      x = "Sequencing depth (reads per composite tube)",
      y = "ASV richness (rarefied)",
      color = "Habitat",
      fill  = "Habitat",
      title = title,
      subtitle = paste("Composite samples (fractions summed per tube),",
                       "Normal moisture only.  Ribbon = mean +/- 1 SD.")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title.position = "plot")
}

# -------------------- Run for each amplicon ---------------------------------
results <- list()
for (tag in names(asv_paths)) {
  cat("\n========== Processing", tag, "==========\n")
  comp_obj <- build_composites(asv_paths[[tag]], md)
  cat("Composite tubes per habitat:\n"); print(table(comp_obj$habitat))
  cat("Composite read totals (per tube):\n")
  print(summary(rowSums(comp_obj$mat)))

  curves  <- rarefaction_curves(comp_obj$mat, step = rar_step)
  summary <- summarize_curves(curves, comp_obj$habitat, min_n = min_n_for_group)

  p <- plot_rarefaction(summary, curves, comp_obj$habitat,
                        paste0(tag, " rarefaction (composite tubes, Normal moisture)"))

  ggsave(file.path(out_dir, paste0("rarefaction_", tag, ".pdf")),
         p, width = 8, height = 5.5)
  ggsave(file.path(out_dir, paste0("rarefaction_", tag, ".png")),
         p, width = 8, height = 5.5, dpi = 300)

  write.csv(curves,  file.path(out_dir, paste0("rarefaction_", tag, "_per_tube.csv")),
            row.names = FALSE)
  write.csv(summary, file.path(out_dir, paste0("rarefaction_", tag, "_summary.csv")),
            row.names = FALSE)

  results[[tag]] <- list(plot = p, curves = curves, summary = summary,
                         comp = comp_obj)
}

cat("\nDone. Outputs in:", out_dir, "\n")
cat("Files:\n")
print(list.files(out_dir))
