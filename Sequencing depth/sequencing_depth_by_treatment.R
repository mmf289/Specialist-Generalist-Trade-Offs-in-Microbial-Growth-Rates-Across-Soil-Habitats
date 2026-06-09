# sequencing_depth_by_treatment.R
# ------------------------------------------------------------------
# Compute min / median / max TOTAL_READS (read pairs) per sample for
# three treatments under water.trt = "Normal":
#   - Rhizo
#   - Detritus
#   - Rhizo_Detritus
#
# Inputs:
#   - 18O.lab.data.density.corrected.csv
#   - filtered/16S/fastq_info_results_run*_filtered.txt
#   - filtered/ITS/fastq_info_results_run*_filtered.txt
#     (produced by subset_fastq_info.R)
#
# One value per sample per run is used: the Forward row (PAIR == "F"),
# since F and R rows report identical TOTAL_READS for paired-end data.
# If a sample appears in more than one run within an assay, read
# counts are SUMMED across runs. Stats are reported per assay
# (16S and ITS are different libraries and shouldn't be pooled).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

## ---- 1. Paths ----------------------------------------------------
base_dir  <- "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/calculate sequencing depth"
meta_file <- file.path(base_dir, "18O.lab.data.density.corrected.csv")
in_dir    <- file.path(base_dir, "filtered")    # output of subset_fastq_info.R

target_ctrts <- c("Rhizo", "Detritus", "Rhizo_Detritus")
target_water <- "Normal"

## ---- 2. Read metadata (treatment per sample) ---------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta_small <- meta %>%
  transmute(sample = trimws(as.character(sample.name)),
            c.trt, water.trt) %>%
  distinct(sample, .keep_all = TRUE)

## ---- 3. Helper: read a filtered fastq_info file ------------------
read_fastq_info <- function(fp) {
  # quote="" because RUN_INFO contains single quotes; comment.char="#"
  # skips the preserved "# ..." header lines at the top of the file.
  read.delim(
    fp, header = TRUE, sep = "\t", quote = "",
    comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE
  )
}

## ---- 4. Gather F-row read counts per sample per run --------------
all_rows <- list()
for (assay in c("16S", "ITS")) {
  fps <- list.files(
    file.path(in_dir, assay),
    pattern = "^fastq_info_results_run.*_filtered\\.txt$",
    full.names = TRUE
  )
  for (fp in fps) {
    run <- sub(".*results_(run\\d+)_filtered\\.txt$", "\\1", basename(fp))
    df  <- read_fastq_info(fp)
    df_F <- df %>%
      filter(PAIR == "F") %>%
      transmute(assay = assay,
                run   = run,
                sample = SAMPLE,
                total_reads = as.integer(TOTAL_READS))
    all_rows[[length(all_rows) + 1]] <- df_F
  }
}
reads <- bind_rows(all_rows)

## ---- 5. Sum reads per sample across runs within an assay ---------
per_sample <- reads %>%
  group_by(assay, sample) %>%
  summarise(total_reads = sum(total_reads), .groups = "drop") %>%
  left_join(meta_small, by = "sample")

## ---- 6. Depth summary for the three target treatments ------------
depth_summary <- per_sample %>%
  filter(water.trt == target_water,
         c.trt %in% target_ctrts) %>%
  group_by(assay, c.trt, water.trt) %>%
  summarise(
    n_samples    = n(),
    min_reads    = min(total_reads),
    median_reads = median(total_reads),
    max_reads    = max(total_reads),
    .groups      = "drop"
  ) %>%
  arrange(assay, factor(c.trt, levels = target_ctrts))

cat("\nSequencing depth (read pairs) by treatment:\n")
print(depth_summary, n = Inf)

## ---- 7. Save outputs --------------------------------------------
write.csv(depth_summary,
          file.path(in_dir, "sequencing_depth_by_treatment.csv"),
          row.names = FALSE)

per_sample_out <- per_sample %>%
  filter(water.trt == target_water, c.trt %in% target_ctrts) %>%
  arrange(assay, c.trt, sample)

write.csv(per_sample_out,
          file.path(in_dir, "sequencing_depth_per_sample.csv"),
          row.names = FALSE)

cat("\nSaved:\n  ", file.path(in_dir, "sequencing_depth_by_treatment.csv"),
    "\n  ", file.path(in_dir, "sequencing_depth_per_sample.csv"), "\n")

## ---- 8. Statistical test: does depth differ among habitats? -----
# Habitat here = c.trt (Rhizo, Detritus, Rhizo_Detritus), restricted
# to water.trt = "Normal". Read counts are typically right-skewed, so
# the primary test is non-parametric (Kruskal-Wallis), with pairwise
# Wilcoxon follow-ups (BH-adjusted). ANOVA on log-transformed reads
# is included as a parametric cross-check, with a Shapiro-Wilk
# normality check per group and a Levene-style variance check.

run_depth_tests <- function(df, assay_name) {
  cat(sprintf("\n================ %s ================\n", assay_name))
  df$c.trt <- factor(df$c.trt, levels = target_ctrts)

  # --- Group normality (Shapiro-Wilk; informational) ---
  cat("\nShapiro-Wilk normality per group (raw reads):\n")
  for (g in levels(df$c.trt)) {
    x <- df$total_reads[df$c.trt == g]
    if (length(x) >= 3 && length(x) <= 5000) {
      st <- shapiro.test(x)
      cat(sprintf("  %-18s n=%3d  W=%.3f  p=%.3g\n",
                  g, length(x), st$statistic, st$p.value))
    } else {
      cat(sprintf("  %-18s n=%3d  (skipped: n outside 3-5000)\n",
                  g, length(x)))
    }
  }

  # --- Equality-of-variance check (Bartlett, also sensitive to normality) ---
  bt <- bartlett.test(total_reads ~ c.trt, data = df)
  cat(sprintf("\nBartlett K^2=%.3f  df=%d  p=%.3g\n",
              bt$statistic, bt$parameter, bt$p.value))

  # --- PRIMARY: Kruskal-Wallis ---
  kw <- kruskal.test(total_reads ~ c.trt, data = df)
  cat(sprintf("\nKruskal-Wallis (primary): chi^2=%.3f  df=%d  p=%.4g\n",
              kw$statistic, kw$parameter, kw$p.value))

  # --- Post-hoc: pairwise Wilcoxon, BH-adjusted ---
  pw <- pairwise.wilcox.test(df$total_reads, df$c.trt,
                             p.adjust.method = "BH", exact = FALSE)
  cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
  print(pw$p.value)

  # --- Parametric cross-check: ANOVA on log10(reads) + Tukey HSD ---
  df$log_reads <- log10(df$total_reads)
  aov_fit <- aov(log_reads ~ c.trt, data = df)
  cat("\nANOVA on log10(reads):\n")
  print(summary(aov_fit))
  cat("\nTukey HSD on log10(reads):\n")
  print(TukeyHSD(aov_fit))

  # Return tidy summaries for saving
  list(
    assay        = assay_name,
    kruskal      = data.frame(assay = assay_name,
                              test = "Kruskal-Wallis",
                              statistic = unname(kw$statistic),
                              df = unname(kw$parameter),
                              p_value = kw$p.value),
    pairwise_wx  = {
      m <- pw$p.value
      out <- expand.grid(group1 = rownames(m), group2 = colnames(m),
                         stringsAsFactors = FALSE)
      out$p_adj <- mapply(function(r, c) m[r, c], out$group1, out$group2)
      out <- out[!is.na(out$p_adj), , drop = FALSE]
      out$assay <- assay_name
      out$method <- "Pairwise Wilcoxon, BH"
      out[, c("assay", "method", "group1", "group2", "p_adj")]
    },
    anova_log10  = {
      s <- summary(aov_fit)[[1]]
      data.frame(assay = assay_name,
                 test  = "ANOVA (log10 reads)",
                 F_statistic = s$`F value`[1],
                 df1 = s$Df[1], df2 = s$Df[2],
                 p_value = s$`Pr(>F)`[1])
    },
    tukey_log10  = {
      tk <- as.data.frame(TukeyHSD(aov_fit)$c.trt)
      tk$contrast <- rownames(tk)
      tk$assay <- assay_name
      tk[, c("assay", "contrast", "diff", "lwr", "upr", "p adj")]
    }
  )
}

# Data to feed into the tests: per-sample read counts (summed across
# runs within an assay), restricted to the three target habitats/Normal
test_data <- per_sample %>%
  filter(water.trt == target_water, c.trt %in% target_ctrts)

results_16S <- run_depth_tests(filter(test_data, assay == "16S"), "16S")
results_ITS <- run_depth_tests(filter(test_data, assay == "ITS"), "ITS")

## ---- 9. Save test results ---------------------------------------
write.csv(
  rbind(results_16S$kruskal, results_ITS$kruskal,
        results_16S$anova_log10, results_ITS$anova_log10),
  file.path(in_dir, "depth_tests_omnibus.csv"),
  row.names = FALSE
)
write.csv(
  rbind(results_16S$pairwise_wx, results_ITS$pairwise_wx),
  file.path(in_dir, "depth_tests_pairwise_wilcoxon.csv"),
  row.names = FALSE
)
write.csv(
  rbind(results_16S$tukey_log10, results_ITS$tukey_log10),
  file.path(in_dir, "depth_tests_tukey_log10.csv"),
  row.names = FALSE
)

cat("\nSaved test results:\n  ",
    file.path(in_dir, "depth_tests_omnibus.csv"), "\n  ",
    file.path(in_dir, "depth_tests_pairwise_wilcoxon.csv"), "\n  ",
    file.path(in_dir, "depth_tests_tukey_log10.csv"), "\n")
