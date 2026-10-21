# ============================================================
# Regional Environmental Change (Springer Nature)
# Chapter 2 - publication-ready figures and tables
# Period: 1985-2022
#
# FINAL REEC VERSION v5:
# - Uses the local folder C:/Cap 2
# - American English throughout the artwork
# - Figure 5 finalized as a 2 x 2 PCA multipanel:
#     (a) Cachoeira River basin
#     (b) Dense Ombrophilous Forest (DOF)
#     (c) Semideciduous Seasonal Forest (SSF)
#     (d) Deciduous Seasonal Forest (DSF)
# - Final editorial refinement after inspection of the v4 PDF outputs
# - Plotmath notation avoids font/glyph problems in beta, R-squared, and yr^-1
# - DSF 1985-1986 annual-change audit printed in the console
# - High-resolution TIFF + vector PDF/EPS
# - Figure captions are NOT embedded in the artwork
#
# Springer Nature-oriented artwork settings:
# - final width: 174 mm (double-column)
# - lettering: 8-12 pt at final size
# - vector PDF/EPS preferred for plots
# - TIFF: 600 dpi for color/combination artwork
# - lowercase panel labels
# - color-blind-friendly viridis time scale
# - visual information remains interpretable without relying only on color
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------
setwd("C:/Cap 2")

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------
pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr", "janitor",
  "lmtest", "sandwich", "Kendall", "trend",
  "ggplot2", "patchwork", "ggrepel", "openxlsx",
  "scales", "viridisLite", "grid", "stringr"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------
# Accepts both the original names and timestamped copies, e.g.
# Bacia_percentual(20261006-135709).xlsx
find_one <- function(pattern) {
  x <- list.files(
    path = ".",
    pattern = pattern,
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(x) == 0) {
    stop("File not found for pattern: ", pattern,
         "\nWorking directory: ", normalizePath(getwd()))
  }
  if (length(x) > 1) {
    info <- file.info(x)
    x <- x[order(info$mtime, decreasing = TRUE)]
    message("Multiple files found. Using the most recently modified: ",
            basename(x[1]))
  }
  x[1]
}

arq_bacia <- find_one("^Bacia_percentual.*\\.xlsx$")
arq_fod   <- find_one("^fod_percentual.*\\.xlsx$")
arq_fes   <- find_one("^fes_percentual.*\\.xlsx$")
arq_fed   <- find_one("^fed_percentual.*\\.xlsx$")

dir_out <- file.path(getwd(), "REEC_outputs")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. Reading helpers
# ------------------------------------------------------------
detect_year_col <- function(nms) {
  cand <- nms[grepl("^ano$|ano_|_ano$|year|^t$|^time$", nms)]
  if (length(cand) == 0) NA_character_ else cand[1]
}

read_landuse_file <- function(path, label) {
  df <- readxl::read_xlsx(path, sheet = 1) |>
    janitor::clean_names()

  ycol <- detect_year_col(names(df))
  if (is.na(ycol)) stop("Year column not found in ", path)

  df <- df |>
    mutate(
      year = suppressWarnings(
        as.numeric(gsub(",", ".", as.character(.data[[ycol]])))
      )
    ) |>
    filter(is.finite(year)) |>
    arrange(year)

  fcol <- names(df)[grepl("^floresta$|florest|forest", names(df))][1]
  if (is.na(fcol)) stop("Forest column not found in ", path)

  # Convert all candidate land-use columns to numeric when possible.
  for (nm in setdiff(names(df), c(ycol, "year"))) {
    if (!is.numeric(df[[nm]])) {
      z <- suppressWarnings(as.numeric(gsub(",", ".", as.character(df[[nm]]))))
      if (sum(is.finite(z)) > 0) df[[nm]] <- z
    }
  }

  df <- df |>
    mutate(
      forest = suppressWarnings(
        as.numeric(gsub(",", ".", as.character(.data[[fcol]])))
      )
    )

  if (nrow(df) != 38 || min(df$year) != 1985 || max(df$year) != 2022) {
    warning(label, ": expected 38 annual observations from 1985 to 2022.")
  }

  if (anyDuplicated(df$year)) {
    stop(label, ": duplicated years were detected.")
  }

  attr(df, "label") <- label
  df
}

dat <- list(
  Basin = read_landuse_file(arq_bacia, "Cachoeira River basin"),
  DOF   = read_landuse_file(arq_fod,   "Dense Ombrophilous Forest"),
  SSF   = read_landuse_file(arq_fes,   "Semideciduous Seasonal Forest"),
  DSF   = read_landuse_file(arq_fed,   "Deciduous Seasonal Forest")
)

# ------------------------------------------------------------
# 4. Trend analysis
# ------------------------------------------------------------
analyze_series <- function(df, label) {
  mod <- lm(forest ~ year, data = df)
  n <- nrow(df)

  beta <- unname(coef(mod)[["year"]])
  r2 <- summary(mod)$r.squared

  lag_nw <- max(1, floor(4 * (n / 100)^(2 / 9)))
  ct_hac <- lmtest::coeftest(
    mod,
    vcov. = sandwich::NeweyWest(
      mod, lag = lag_nw, prewhite = FALSE, adjust = TRUE
    )
  )
  p_hac <- unname(ct_hac["year", "Pr(>|t|)"])

  dw <- tryCatch(
    unname(lmtest::dwtest(mod)$statistic),
    error = function(e) NA_real_
  )

  mk <- Kendall::MannKendall(df$forest)
  tau_mk <- unname(mk$tau)
  p_mk <- unname(mk$sl)

  sen <- trend::sens.slope(as.numeric(df$forest))
  sen_slope <- unname(sen$estimates)

  df <- df |>
    mutate(
      annual_change = 100 * (forest - lag(forest)) / lag(forest)
    )

  list(
    data = df,
    model = mod,
    summary = tibble(
      `Spatial unit` = label,
      Period = paste0(min(df$year), "-", max(df$year)),
      n = n,
      `beta (percentage points yr-1)` = beta,
      `HAC p-value` = p_hac,
      `Durbin-Watson` = dw,
      `Mann-Kendall tau` = tau_mk,
      `Mann-Kendall p-value` = p_mk,
      `Sen slope (percentage points yr-1)` = sen_slope,
      `R-squared` = r2
    )
  )
}

ana <- list(
  Basin = analyze_series(dat$Basin, "Basin"),
  DOF   = analyze_series(dat$DOF,   "DOF"),
  SSF   = analyze_series(dat$SSF,   "SSF"),
  DSF   = analyze_series(dat$DSF,   "DSF")
)

table1 <- purrr::map_dfr(ana, "summary") |>
  mutate(
    `Spatial unit` = factor(
      `Spatial unit`,
      levels = c("Basin", "DOF", "SSF", "DSF")
    )
  ) |>
  arrange(`Spatial unit`)

# ------------------------------------------------------------
# 5. REEC figure theme
# ------------------------------------------------------------
# Arial is normally available on Windows.
# If Arial is unavailable on another computer, change base_family to "sans".
theme_reec <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.30, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(5, 7, 5, 5),
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    panel.grid = element_blank()
  )

fmt_p <- function(p) {
  ifelse(
    p < 0.001,
    format(p, scientific = TRUE, digits = 3),
    formatC(p, format = "f", digits = 3)
  )
}

# ------------------------------------------------------------
# 6. Figure 3 - forest-cover trajectories
# ------------------------------------------------------------
make_trend_panel <- function(obj, panel_label) {
  df <- obj$data
  s <- obj$summary

  beta_txt <- sprintf("%.3f", s$`beta (percentage points yr-1)`)
  r2_txt   <- sprintf("%.3f", s$`R-squared`)
  p_txt <- ifelse(
    s$`HAC p-value` < 0.001,
    "HAC~italic(p)<0.001",
    paste0("HAC~italic(p)==", fmt_p(s$`HAC p-value`))
  )
  ann <- paste0(
    "beta==", beta_txt, "~'pp'~yr^{-1}",
    "*','~~R^2==", r2_txt,
    "*';'~~", p_txt
  )

  ggplot(df, aes(year, forest)) +
    geom_line(linewidth = 0.45, color = "grey25") +
    geom_point(size = 1.25, shape = 21, fill = "white",
               color = "black", stroke = 0.30) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      linetype = "dashed", linewidth = 0.50, color = "black"
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = panel_label,
      hjust = -0.04, vjust = 1.10, family = "Arial",
      fontface = "bold", size = 2.85
    ) +
    annotate(
      "text", x = Inf, y = Inf, label = ann,
      hjust = 1.04, vjust = 2.25, family = "Arial",
      size = 2.35, parse = TRUE
    ) +
    labs(x = "Year", y = "Forest cover (%)") +
    scale_x_continuous(
      breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)
    ) +
    theme_reec
}

p3 <- (
  make_trend_panel(ana$Basin, "(a) Cachoeira River basin") +
  make_trend_panel(ana$DOF,   "(b) Dense Ombrophilous Forest") +
  make_trend_panel(ana$SSF,   "(c) Semideciduous Seasonal Forest") +
  make_trend_panel(ana$DSF,   "(d) Deciduous Seasonal Forest")
) + patchwork::plot_layout(ncol = 2)

# ------------------------------------------------------------
# 7. Figure 4 - annual forest-cover change
# ------------------------------------------------------------
make_change_panel <- function(obj, panel_label) {
  df <- obj$data |> filter(!is.na(annual_change))

  ggplot(df, aes(year, annual_change)) +
    geom_hline(
      yintercept = 0, linetype = "dotted",
      linewidth = 0.40, color = "grey35"
    ) +
    geom_line(linewidth = 0.45, color = "grey25") +
    geom_point(
      size = 1.25, shape = 21, fill = "white",
      color = "black", stroke = 0.30
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = panel_label,
      hjust = -0.04, vjust = 1.10, family = "Arial",
      fontface = "bold", size = 2.85
    ) +
    labs(x = "Year", y = "Annual change (%)") +
    scale_x_continuous(
      breaks = c(1990, 1995, 2000, 2005, 2010, 2015, 2020)
    ) +
    theme_reec
}

p4 <- (
  make_change_panel(ana$Basin, "(a) Cachoeira River basin") +
  make_change_panel(ana$DOF,   "(b) Dense Ombrophilous Forest") +
  make_change_panel(ana$SSF,   "(c) Semideciduous Seasonal Forest") +
  make_change_panel(ana$DSF,   "(d) Deciduous Seasonal Forest")
) + patchwork::plot_layout(ncol = 2)

# ------------------------------------------------------------
# 8. PCA
# ------------------------------------------------------------
# American English labels used in Figure 5 and supplementary loadings.
translate_class <- c(
  agropecuaria = "Agricultural mosaic",
  apicum = "Salt flat",
  aquicultura = "Aquaculture",
  area_urbanizada = "Urban area",
  cafe = "Coffee",
  campo_alagado_e_area_pantanosa = "Wetland",
  floresta = "Forest",
  mangue = "Mangrove",
  mineracao = "Mining",
  outras_areas_nao_vegetadas = "Other non-vegetated areas",
  outras_formacoes_nao_florestais = "Other non-forest formations",
  outras_lavouras_perenes = "Other perennial crops",
  pastagem = "Pasture",
  praia_duna_e_areal = "Beach/dune/sand",
  restinga_herbacea = "Herbaceous restinga",
  rio_lago_e_oceano = "Water bodies",
  savana = "Savannic formation",
  silvicultura = "Silviculture"
)

run_pca <- function(df, unit_name) {
  # Keep the original land-use percentages.
  # Remove only the year aliases and the auxiliary 'forest' copy.
  x <- df |>
    select(-any_of(c("ano", "year", "forest"))) |>
    select(where(is.numeric))

  # Remove variables with all missing values or zero variance.
  keep <- vapply(
    x,
    function(z) {
      z <- z[is.finite(z)]
      length(z) >= 2 && sd(z) > 0
    },
    logical(1)
  )
  x <- x[, keep, drop = FALSE]

  if (ncol(x) < 3) {
    stop(unit_name, ": too few non-constant land-use variables for PCA.")
  }

  # PCA requires complete values. With these source spreadsheets the annual
  # percentages are expected to be complete; stop rather than silently impute.
  if (anyNA(x)) {
    bad <- names(x)[vapply(x, anyNA, logical(1))]
    stop(
      unit_name,
      ": missing values detected in PCA variables: ",
      paste(bad, collapse = ", "),
      ". Review the source spreadsheet before running PCA."
    )
  }

  pc <- prcomp(x, center = TRUE, scale. = TRUE)

  scores <- as.data.frame(pc$x[, 1:2, drop = FALSE])
  scores$year <- df$year

  load <- as.data.frame(pc$rotation[, 1:2, drop = FALSE])
  load$class <- rownames(load)

  # PCA signs are arbitrary. Orient PC1 so the time trajectory runs,
  # as far as possible, from earlier years on the left to later years
  # on the right. This improves cross-panel visual interpretation
  # without changing distances, variance, or relationships.
  early_pc1 <- mean(scores$PC1[scores$year <= 1990], na.rm = TRUE)
  late_pc1  <- mean(scores$PC1[scores$year >= 2017], na.rm = TRUE)

  if (early_pc1 > late_pc1) {
    scores$PC1 <- -scores$PC1
    load$PC1 <- -load$PC1
  }

  var_exp <- 100 * pc$sdev^2 / sum(pc$sdev^2)

  list(
    scores = scores,
    loadings = load,
    variance = var_exp,
    model = pc
  )
}

pca_all <- list(
  Basin = run_pca(dat$Basin, "Basin"),
  DOF   = run_pca(dat$DOF,   "DOF"),
  SSF   = run_pca(dat$SSF,   "SSF"),
  DSF   = run_pca(dat$DSF,   "DSF")
)

table2 <- tibble(
  `Spatial unit` = c("Basin", "DOF", "SSF", "DSF"),
  `PC1 variance (%)` = sapply(
    pca_all[c("Basin", "DOF", "SSF", "DSF")],
    function(x) x$variance[1]
  ),
  `PC2 variance (%)` = sapply(
    pca_all[c("Basin", "DOF", "SSF", "DSF")],
    function(x) x$variance[2]
  )
) |>
  mutate(
    `Cumulative variance (%)` =
      `PC1 variance (%)` + `PC2 variance (%)`
  )

# ------------------------------------------------------------
# 9. FIGURE 5 - MULTIPANEL PCA FOR REEC
# ------------------------------------------------------------
# Design rationale:
# 1) Four panels directly compare the basin with the three forest formations.
# 2) All 38 years are shown as a connected temporal trajectory.
# 3) A color-blind-friendly continuous scale represents time.
# 4) Selected years are labeled to avoid excessive text.
# 5) 1985 and 2022 receive distinct endpoint symbols, so the temporal
#    direction remains interpretable even in grayscale.
# 6) Only the strongest loading vectors are labeled in each panel.
# 7) No figure title or caption is embedded in the artwork.

# Years that provide temporal landmarks without overcrowding the panels.
years_to_label <- c(1985, 1990, 2000, 2010, 2020, 2022)

# Number of loading vectors displayed in each panel.
# Four provides a good compromise between interpretation and legibility.
n_loadings <- 4

prepare_pca_panel <- function(pca_obj, n_vectors = 4) {
  scores <- pca_obj$scores

  load <- pca_obj$loadings |>
    mutate(
      loading_norm = sqrt(PC1^2 + PC2^2),
      label = dplyr::recode(
        class,
        !!!translate_class,
        .default = stringr::str_to_sentence(gsub("_", " ", class))
      )
    ) |>
    slice_max(order_by = loading_norm, n = n_vectors, with_ties = FALSE)

  # Scale loading arrows relative to each panel's score cloud.
  # The factor is intentionally conservative to leave room for labels.
  pc1_span <- diff(range(scores$PC1, na.rm = TRUE))
  pc2_span <- diff(range(scores$PC2, na.rm = TRUE))

  denom1 <- max(abs(load$PC1), na.rm = TRUE)
  denom2 <- max(abs(load$PC2), na.rm = TRUE)

  arrow_mult <- 0.30 * min(
    pc1_span / denom1,
    pc2_span / denom2
  )

  load <- load |>
    mutate(
      xend = PC1 * arrow_mult,
      yend = PC2 * arrow_mult
    )

  list(scores = scores, load = load)
}

make_pca_panel <- function(pca_obj, panel_label, panel_name) {
  z <- prepare_pca_panel(pca_obj, n_vectors = n_loadings)
  scores <- z$scores
  load <- z$load

  # Endpoint data are plotted separately for print/grayscale redundancy.
  endpoints <- scores |>
    filter(year %in% c(1985, 2022)) |>
    mutate(
      endpoint = factor(
        year,
        levels = c(1985, 2022),
        labels = c("1985", "2022")
      )
    )

  ggplot(scores, aes(PC1, PC2)) +
    # Reference axes
    geom_hline(
      yintercept = 0, linetype = "dotted",
      linewidth = 0.32, color = "grey70"
    ) +
    geom_vline(
      xintercept = 0, linetype = "dotted",
      linewidth = 0.32, color = "grey70"
    ) +

    # Temporal trajectory: thin neutral line beneath the points
    geom_path(
      aes(group = 1),
      linewidth = 0.40,
      color = "grey55",
      alpha = 0.85
    ) +

    # Annual scores
    geom_point(
      aes(color = year),
      size = 1.65,
      alpha = 0.92
    ) +

    # Distinct first/last year markers preserve interpretation in grayscale
    geom_point(
      data = endpoints,
      mapping = aes(x = PC1, y = PC2, shape = endpoint),
      inherit.aes = FALSE,
      size = 2.55,
      stroke = 0.65,
      fill = "white",
      color = "black"
    ) +

    # Selected year labels
    ggrepel::geom_text_repel(
      data = subset(scores, year %in% years_to_label),
      aes(label = year),
      family = "Arial",
      size = 2.45,
      color = "grey15",
      box.padding = 0.24,
      point.padding = 0.12,
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = "grey55",
      max.overlaps = Inf,
      seed = 2026,
      show.legend = FALSE
    ) +

    # Loading vectors
    geom_segment(
      data = load,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = grid::arrow(
        length = grid::unit(1.45, "mm"),
        type = "closed"
      ),
      linewidth = 0.48,
      color = "grey15"
    ) +

    # Loading labels
    ggrepel::geom_text_repel(
      data = load,
      aes(x = xend, y = yend, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "italic",
      size = 2.10,
      color = "black",
      box.padding = 0.55,
      point.padding = 0.28,
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = "grey45",
      max.overlaps = Inf,
      force = 1.8,
      force_pull = 0.65,
      max.time = 4,
      seed = 2026
    ) +

    # Accessible sequential palette for time
    scale_color_viridis_c(
      option = "C",
      direction = 1,
      begin = 0.08,
      end = 0.92,
      breaks = c(1985, 2000, 2010, 2022),
      name = "Year",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(38, "mm"),
        barheight = grid::unit(3.0, "mm"),
        ticks.colour = "black",
        frame.colour = "grey35"
      )
    ) +

    scale_shape_manual(
      values = c("1985" = 21, "2022" = 24),
      name = "Endpoint"
    ) +

    labs(
      x = sprintf("PC1 (%.1f%%)", pca_obj$variance[1]),
      y = sprintf("PC2 (%.1f%%)", pca_obj$variance[2]),
      title = paste0(panel_label, " ", panel_name)
    ) +

    coord_cartesian(clip = "off") +
    theme_reec +
    theme(
      plot.title = element_text(
        family = "Arial",
        size = 8.5,
        face = "bold",
        margin = margin(b = 4)
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.key.width = grid::unit(8, "mm"),
      legend.spacing.x = grid::unit(2.5, "mm"),
      legend.spacing.y = grid::unit(1.5, "mm"),
      plot.margin = margin(8, 13, 8, 10)
    )
}

p5a <- make_pca_panel(
  pca_all$Basin,
  "(a)",
  "Cachoeira River basin"
)

p5b <- make_pca_panel(
  pca_all$DOF,
  "(b)",
  "Dense Ombrophilous Forest"
)

p5c <- make_pca_panel(
  pca_all$SSF,
  "(c)",
  "Semideciduous Seasonal Forest"
)

p5d <- make_pca_panel(
  pca_all$DSF,
  "(d)",
  "Deciduous Seasonal Forest"
)

# Collect the common legends once, beneath the four panels.
p5 <- (
  (p5a + p5b) /
  (p5c + p5d)
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# ------------------------------------------------------------
# 9.1. Audit of the first DSF annual-change observation
# ------------------------------------------------------------
# The first DSF annual-change value is visually large because the
# denominator is the previous year's forest-cover percentage. It is
# retained in Figure 4. This block prints the 1985 and 1986 source
# values and independently recalculates the percentage change so the
# analyst can verify the observation before submission.
dsf_1985_1986 <- ana$DSF$data |>
  filter(year %in% c(1985, 1986)) |>
  select(year, forest, annual_change)

if (all(c(1985, 1986) %in% dsf_1985_1986$year)) {
  dsf_1985 <- dsf_1985_1986$forest[dsf_1985_1986$year == 1985]
  dsf_1986 <- dsf_1985_1986$forest[dsf_1985_1986$year == 1986]
  dsf_change_check <- 100 * (dsf_1986 - dsf_1985) / dsf_1985

  cat("\nDSF first annual-change audit:\n")
  print(dsf_1985_1986)
  cat(
    sprintf(
      "Independent check for 1985-1986: %.6f%%\n",
      dsf_change_check
    )
  )
}

# ------------------------------------------------------------
# FINAL v5 publication note
# ------------------------------------------------------------
# Figure 4 data and scaling are intentionally unchanged from v4.
# The large first DSF annual-change observation is retained and audited
# from the source values. Figure 3 changes are typographic only.
# Figure 5 changes affect label placement only; PCA scores, loadings,
# explained variance, selected loading vectors, and temporal trajectory
# are unchanged.
#
# ------------------------------------------------------------
# 10. Export figures
# ------------------------------------------------------------
# Figure 3 and Figure 4: line-dominant artwork.
ggsave(
  file.path(dir_out, "Fig3_ForestCover_Trends_REEC_v5_FINAL.tif"),
  p3, width = 174, height = 145, units = "mm",
  dpi = 1200, compression = "lzw"
)
ggsave(
  file.path(dir_out, "Fig3_ForestCover_Trends_REEC_v5_FINAL.pdf"),
  p3, width = 174, height = 145, units = "mm",
  device = cairo_pdf
)
try(
  ggsave(
    file.path(dir_out, "Fig3_ForestCover_Trends_REEC_v5_FINAL.eps"),
    p3, width = 174, height = 145, units = "mm",
    device = cairo_ps
  ),
  silent = TRUE
)

ggsave(
  file.path(dir_out, "Fig4_Annual_Forest_Change_REEC_v5_FINAL.tif"),
  p4, width = 174, height = 145, units = "mm",
  dpi = 1200, compression = "lzw"
)
ggsave(
  file.path(dir_out, "Fig4_Annual_Forest_Change_REEC_v5_FINAL.pdf"),
  p4, width = 174, height = 145, units = "mm",
  device = cairo_pdf
)
try(
  ggsave(
    file.path(dir_out, "Fig4_Annual_Forest_Change_REEC_v5_FINAL.eps"),
    p4, width = 174, height = 145, units = "mm",
    device = cairo_ps
  ),
  silent = TRUE
)

# Figure 5: color + line artwork.
# 174 x 165 mm gives each PCA panel enough space for labels while
# remaining suitable for a two-column journal page.
ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL.tif"),
  p5,
  width = 174, height = 172, units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Vector PDF is the preferred working/submission version when accepted.
ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL.pdf"),
  p5,
  width = 174, height = 172, units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# EPS vector copy. On some Windows installations cairo_ps may be
# unavailable; the script continues and keeps TIFF/PDF if EPS fails.
try(
  ggsave(
    file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL.eps"),
    p5,
    width = 174, height = 172, units = "mm",
    device = cairo_ps,
    bg = "white"
  ),
  silent = TRUE
)

# High-resolution PNG only for quick visual checking, not the primary
# submission format.
ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL_preview.png"),
  p5,
  width = 174, height = 172, units = "mm",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 11. Export journal tables and PCA loadings
# ------------------------------------------------------------
table1_out <- table1 |>
  mutate(
    `beta (percentage points yr-1)` =
      round(`beta (percentage points yr-1)`, 3),
    `HAC p-value` = signif(`HAC p-value`, 3),
    `Durbin-Watson` = round(`Durbin-Watson`, 3),
    `Mann-Kendall tau` = round(`Mann-Kendall tau`, 3),
    `Mann-Kendall p-value` =
      signif(`Mann-Kendall p-value`, 3),
    `Sen slope (percentage points yr-1)` =
      round(`Sen slope (percentage points yr-1)`, 3),
    `R-squared` = round(`R-squared`, 3)
  )

table2_out <- table2 |>
  mutate(
    across(
      c(`PC1 variance (%)`,
        `PC2 variance (%)`,
        `Cumulative variance (%)`),
      ~ round(.x, 1)
    )
  )

write.csv(
  table1_out,
  file.path(dir_out, "Table1_Trend_REEC.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  table2_out,
  file.path(dir_out, "Table2_PCA_REEC.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Table 1 Trends")
openxlsx::writeData(wb, "Table 1 Trends", table1_out)

openxlsx::addWorksheet(wb, "Table 2 PCA")
openxlsx::writeData(wb, "Table 2 PCA", table2_out)

# Supplementary PCA loadings, including loading magnitude and rank.
for (nm in c("Basin", "DOF", "SSF", "DSF")) {
  sh <- paste0("PCA loadings ", nm)
  openxlsx::addWorksheet(wb, sh)

  tmp <- pca_all[[nm]]$loadings |>
    mutate(
      Class = dplyr::recode(
        class,
        !!!translate_class,
        .default = stringr::str_to_sentence(gsub("_", " ", class))
      ),
      `Loading magnitude` = sqrt(PC1^2 + PC2^2)
    ) |>
    arrange(desc(`Loading magnitude`)) |>
    mutate(Rank = row_number()) |>
    select(Rank, Class, PC1, PC2, `Loading magnitude`)

  openxlsx::writeData(wb, sh, tmp)
}

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sh in names(wb)) {
  dat_sh <- openxlsx::readWorkbook(wb, sh)

  openxlsx::addStyle(
    wb, sh, header_style,
    rows = 1,
    cols = seq_len(ncol(dat_sh)),
    gridExpand = TRUE
  )

  openxlsx::setColWidths(
    wb, sh,
    cols = seq_len(ncol(dat_sh)),
    widths = "auto"
  )

  openxlsx::freezePane(wb, sh, firstRow = TRUE)
}

openxlsx::saveWorkbook(
  wb,
  file.path(dir_out, "REEC_Tables_PCA_Loadings.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------------------
# 12. Console validation report
# ------------------------------------------------------------
cat("\n============================================================\n")
cat("REEC Chapter 2 processing completed.\n")
cat("Output folder:\n", normalizePath(dir_out), "\n\n", sep = "")

cat("Input files used:\n")
cat("Basin: ", basename(arq_bacia), "\n", sep = "")
cat("DOF:   ", basename(arq_fod), "\n", sep = "")
cat("SSF:   ", basename(arq_fes), "\n", sep = "")
cat("DSF:   ", basename(arq_fed), "\n\n", sep = "")

cat("Figure 5 PCA variance explained:\n")
print(table2_out)

cat("\nStrongest loading vectors displayed in Figure 5:\n")
for (nm in c("Basin", "DOF", "SSF", "DSF")) {
  z <- prepare_pca_panel(pca_all[[nm]], n_vectors = n_loadings)
  cat("\n", nm, ":\n", sep = "")
  print(
    z$load |>
      select(label, PC1, PC2, loading_norm) |>
      arrange(desc(loading_norm))
  )
}

cat("\nFigure 5 files:\n")
cat("  Fig5_PCA_Multipanel_REEC_v5_FINAL.tif  (600 dpi)\n")
cat("  Fig5_PCA_Multipanel_REEC_v5_FINAL.pdf  (vector)\n")
cat("  Fig5_PCA_Multipanel_REEC_v5_FINAL.eps  (vector, if Cairo EPS is available)\n")
cat("  Fig5_PCA_Multipanel_REEC_v5_FINAL_preview.png (300 dpi preview)\n")

cat("\nAbbreviations used in the manuscript:\n")
cat("DOF = Dense Ombrophilous Forest\n")
cat("SSF = Semideciduous Seasonal Forest\n")
cat("DSF = Deciduous Seasonal Forest\n")

cat("\nRecommended manuscript caption:\n")
cat(
  "Fig. 5 Principal component analysis of annual land-use and land-cover ",
  "composition from 1985 to 2022 for (a) the Cachoeira River basin, ",
  "(b) Dense Ombrophilous Forest, (c) Semideciduous Seasonal Forest, ",
  "and (d) Deciduous Seasonal Forest. Points represent individual years, ",
  "gray lines connect consecutive years, and arrows show the four land-use ",
  "and land-cover classes with the largest loading magnitudes in each ",
  "ordination. PC1 and PC2 percentages indicate the variance explained ",
  "within each spatial unit.\n",
  sep = ""
)

cat("\nFINAL v5 editorial checks:\n")
cat("  - Fig. 3 statistics are separated from panel titles and use plotmath.\n")
cat("  - Fig. 5 has final loading-label spacing and stronger repulsion.\n")
cat("  - DSF 1985-1986 annual change is audited above and not removed.\n")
cat("  - Inspect the PDF/vector files at 100% before submission.\n\n")
cat("============================================================\n")
