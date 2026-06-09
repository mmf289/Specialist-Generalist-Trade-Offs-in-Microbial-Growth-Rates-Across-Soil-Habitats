
# --- packages ---------------------------------------------------------------
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2"))
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- settings ---------------------------------------------------------------
infile  <- "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/github/NRI:NTI/out.xlsx"             # path to the Excel file
outfile <- "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/github/NRI:NTI/nri_nti_figure.png"   # output figure

# Error bars: choose what the vertical lines represent.
#   "sd"    = mean +/- standard deviation (default; matches small-n spread)
#   "se"    = mean +/- standard error
#   "range" = min to max
errbar_type <- "sd"

# Sheet -> (metric, marker) mapping
sheet_map <- tibble::tribble(
  ~sheet,      ~metric, ~marker,
  "MPD 16S",   "NRI",   "16S",
  "MPD ITS",   "NRI",   "ITS",
  "MNTD 16S",  "NTI",   "16S",
  "MNTD ITS",  "NTI",   "ITS"
)

# --- read & combine all four sheets ----------------------------------------
raw <- sheet_map %>%
  rowwise() %>%
  mutate(data = list(read_excel(infile, sheet = sheet))) %>%
  ungroup() %>%
  select(metric, marker, data) %>%
  unnest(data)

# --- tidy / recode ----------------------------------------------------------
hab_levels  <- c("detritus", "rhizo", "rhizo det")
type_levels <- c("all taxa - amplicon sequencing", "specialist taxa", "generalist taxa")

dat <- raw %>%
  mutate(
    value   = -`mpd.obs.z`,                       # SES -> NRI / NTI
    habitat = recode(habitat,
                     "detritusphere"               = "detritus",
                     "rhizosphere"                 = "rhizo",
                     "rhizosphere + detritusphere" = "rhizo det"),
    taxa    = recode(type,
                     "whole community" = "all taxa - amplicon sequencing",
                     "specialist"      = "specialist taxa",
                     "generalist"      = "generalist taxa"),
    habitat = factor(habitat, levels = hab_levels),
    taxa    = factor(taxa,    levels = type_levels),
    metric  = factor(metric,  levels = c("NRI", "NTI")),
    marker  = factor(marker,  levels = c("16S", "ITS"))
  ) %>%
  filter(!is.na(value))

# --- summarise per group ----------------------------------------------------
summ <- dat %>%
  group_by(metric, marker, habitat, taxa) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    se   = sd(value) / sqrt(n()),
    lo   = min(value),
    hi   = max(value),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = dplyr::case_when(
      errbar_type == "sd"    ~ mean - sd,
      errbar_type == "se"    ~ mean - se,
      errbar_type == "range" ~ lo
    ),
    ymax = dplyr::case_when(
      errbar_type == "sd"    ~ mean + sd,
      errbar_type == "se"    ~ mean + se,
      errbar_type == "range" ~ hi
    )
  )

# --- colours ----------------------------------------------------------------
pal <- c(
  "all taxa - amplicon sequencing" = "#E1604E",  # coral / red
  "specialist taxa"                = "#7FC5CB",  # light teal
  "generalist taxa"                = "#2E8B8B"   # dark teal
)

dodge <- position_dodge(width = 0.6)

# --- plot -------------------------------------------------------------------
p <- ggplot(summ, aes(x = habitat, y = mean, colour = taxa)) +
  # non-significant zone (|SES| < 2) and significance threshold
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2,
           fill = "grey88", colour = NA) +
  geom_hline(yintercept = 2, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  geom_linerange(aes(ymin = ymin, ymax = ymax),
                 position = dodge, linewidth = 0.6) +
  geom_point(position = dodge, size = 2.2) +
  facet_grid(metric ~ marker, switch = "y") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_continuous(breaks = seq(-2, 6, 2)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold"),
    strip.text.x     = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    legend.key       = element_blank()
  )

# --- save -------------------------------------------------------------------
ggsave(outfile, p, width = 7.5, height = 5.2, dpi = 300)
message("Wrote ", outfile)
print(p)
