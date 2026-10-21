# Série temporal de cobertura florestal e PCA — Bacia do Rio Cachoeira (1985–2022)
# Entradas: data/raw/*.xlsx (percentuais anuais por classe e unidade espacial)
# Análises: OLS, Newey-West HAC, Durbin–Watson, Mann–Kendall, Sen, PCA
# Saídas: outputs/ (tabelas numéricas, loadings, Figuras 3–5, registros da sessão)
# A coluna "NA" dos XLSX é uma categoria residual, não um valor ausente.
# Ela permanece na PCA quando sua variância é maior que zero, como na análise original.



## Iniciar análise
rm(list = ls())

## Caminhos de trabalho
find_repo_root <- function() {
  candidates <- unique(c(getwd(), dirname(getwd()), dirname(dirname(getwd()))))
  hits <- candidates[file.exists(file.path(candidates, "data", "raw", "Bacia_percentual.xlsx"))]
  if (!length(hits)) stop(
    "Repository root not found. Open the project .Rproj or setwd() to its root before running."
  )
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}
repo_root <- find_repo_root()

## Dependências
pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr", "janitor",
  "lmtest", "sandwich", "Kendall", "trend",
  "ggplot2", "patchwork", "ggrepel", "openxlsx",
  "scales", "viridisLite", "grid", "stringr"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

## Dados de entrada
input_dir <- file.path(repo_root, "data", "raw")
arq_bacia <- file.path(input_dir, "Bacia_percentual.xlsx")
arq_fod   <- file.path(input_dir, "fod_percentual.xlsx")
arq_fes   <- file.path(input_dir, "fes_percentual.xlsx")
arq_fed   <- file.path(input_dir, "fed_percentual.xlsx")

for (p in c(arq_bacia, arq_fod, arq_fes, arq_fed)) {
  if (!file.exists(p)) stop("Required source file not found: ", p)
}

dir_out <- file.path(repo_root, "outputs")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## Leitura e padronização
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

## Verificação dos dados de entrada
validate_input <- function(x, label) {
  stopifnot(nrow(x) == 38L, identical(as.integer(x$year), 1985:2022))
  classes <- x |>
    dplyr::select(-dplyr::any_of(c("ano", "year", "forest"))) |>
    dplyr::select(dplyr::where(is.numeric))
  if (anyNA(classes) || any(!is.finite(as.matrix(classes))))
    stop(label, ": non-finite class percentage(s)")
  sums <- rowSums(classes)
  if (max(abs(sums - 100)) > 0.0001)
    stop(label, ": annual land-cover percentages do not sum to 100")
  if (any(x$forest < 0 | x$forest > 100))
    stop(label, ": forest proportions outside the 0–100 range")
  invisible(TRUE)
}
invisible(purrr::imap(dat, validate_input))

## Tendências temporais: OLS, HAC, DW, MK e Sen
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

## Elementos gráficos
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

## Figura 3: cobertura florestal anual
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

## Figura 4: variação anual da cobertura florestal
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

## PCA das classes de uso e cobertura
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

## Ajuste da PCA padronizada
run_pca <- function(df, unit_name) {
  # "NA" e classes de baixa frequência são mantidas se houver variância.
  # Trata-se de uma classe de origem, não de valores ausentes.
  x <- df |>
    select(-any_of(c("ano", "year", "forest"))) |>
    select(where(is.numeric))

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

  # Valores NA do R (ausência de dado) são tratados separadamente.
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

years_to_label <- c(1985, 1990, 2000, 2010, 2020, 2022)

n_loadings <- 4

## Seleção dos vetores da PCA
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

## Figura 5: ordenações temporais
make_pca_panel <- function(pca_obj, panel_label, panel_name) {
  z <- prepare_pca_panel(pca_obj, n_vectors = n_loadings)
  scores <- z$scores
  load <- z$load

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
    geom_hline(
      yintercept = 0, linetype = "dotted",
      linewidth = 0.32, color = "grey70"
    ) +
    geom_vline(
      xintercept = 0, linetype = "dotted",
      linewidth = 0.32, color = "grey70"
    ) +

    geom_path(
      aes(group = 1),
      linewidth = 0.40,
      color = "grey55",
      alpha = 0.85
    ) +

    geom_point(
      aes(color = year),
      size = 1.65,
      alpha = 0.92
    ) +

    geom_point(
      data = endpoints,
      mapping = aes(x = PC1, y = PC2, shape = endpoint),
      inherit.aes = FALSE,
      size = 2.55,
      stroke = 0.65,
      fill = "white",
      color = "black"
    ) +

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

p5 <- (
  (p5a + p5b) /
  (p5c + p5d)
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

## Verificação da primeira variação anual em DSF
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

ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL.tif"),
  p5,
  width = 174, height = 172, units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL.pdf"),
  p5,
  width = 174, height = 172, units = "mm",
  device = cairo_pdf,
  bg = "white"
)

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

ggsave(
  file.path(dir_out, "Fig5_PCA_Multipanel_REEC_v5_FINAL_preview.png"),
  p5,
  width = 174, height = 172, units = "mm",
  dpi = 300,
  bg = "white"
)

## Exportação das tabelas e dos loadings
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

## Relatório de processamento
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

writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))
writeLines(capture.output(tools::md5sum(c(arq_bacia, arq_fod, arq_fes, arq_fed))),
           file.path(dir_out, "input_md5.txt"))
cat("Reproducibility audit written to outputs/sessionInfo.txt and outputs/input_md5.txt\n")
