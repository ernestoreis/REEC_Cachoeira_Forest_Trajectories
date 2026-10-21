

# ============================================================
#  SÉRIE TEMPORAL (1985–2022)
#  Testes: Durbin–Watson (DW) + Newey–West HAC + Mann–Kendall / Sen
#  Entrada: 4 arquivos separados (Bacia, FED, FES, FOD)
#  Saída: Tabela ABNT (Excel) com β (pp/ano), p HAC, DW, tau MK, Sen slope
# ============================================================
setwd("C:/...") # Indique aqui a pasta de trabalho

rm(list = ls())

# ----------------------------
# 0) Pacotes 
# ----------------------------
pkgs <- c("readxl","dplyr","purrr","janitor","lmtest","sandwich","Kendall","trend","openxlsx")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

# ----------------------------
# 1) Caminhos (AJUSTE se necessário)
# ----------------------------
arq_bacia <- "Bacia_percentual.xlsx"
arq_fed   <- "fed_percentual.xlsx"
arq_fes   <- "fes_percentual.xlsx"
arq_fod   <- "fod_percentual.xlsx"

# Pasta de saída
dir_out <- "saida_serie_temporal"
dir.create(dir_out, showWarnings = FALSE)
out_xlsx <- file.path(dir_out, "Tabela_ABNT_Tendencia_1985_2022.xlsx")
out_csv  <- file.path(dir_out, "Tabela_ABNT_Tendencia_1985_2022.csv")

# ----------------------------
# 2) Leitura (detecta colunas ano e percentual)
# ----------------------------
# Regras de detecção:
# - coluna de ano: contém "ano" (ou "year")
# - coluna de floresta %: tenta achar "floresta"/"forest"; se não, pega a 1ª numérica que NÃO seja ano
detect_year_col <- function(nms) {
  cand <- nms[grepl("^ano$|ano_|_ano$|year|^t$|^time$", nms)]
  if (length(cand) == 0) NA_character_ else cand[1]
}

detect_value_col <- function(df, year_col) {
  nms <- names(df)
  
  # 1) tenta por nome
  cand1 <- nms[grepl("florest|forest|perc_florest|pct_florest|cobertura_florest|forest_pct|floresta_pct", nms)]
  cand1 <- setdiff(cand1, year_col)
  if (length(cand1) > 0) return(cand1[1])
  
  # 2) fallback: primeira coluna numérica diferente do ano
  num_cols <- nms[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, year_col)
  if (length(num_cols) == 0) NA_character_ else num_cols[1]
}

read_series_file <- function(path, label) {
  shs <- readxl::excel_sheets(path)
  # regra: usa a 1ª aba; se tiver mais, ainda assim é plug-and-play
  df_raw <- readxl::read_xlsx(path, sheet = shs[1])
  df <- janitor::clean_names(df_raw)
  
  ycol <- detect_year_col(names(df))
  if (is.na(ycol)) {
    stop("Não encontrei coluna de ANO em: ", label, " (", path, ").\n",
         "Colunas disponíveis: ", paste(names(df), collapse = ", "))
  }
  
  # Coerção: ano e valor devem virar numéricos (aceita vírgula)
  df <- df %>%
    dplyr::mutate(
      .year = suppressWarnings(as.numeric(gsub(",", ".", as.character(.data[[ycol]]))))
    )
  
  vcol <- detect_value_col(df, year_col = ycol)
  if (is.na(vcol)) {
    stop("Não encontrei coluna de percentual/valor em: ", label, " (", path, ").\n",
         "Colunas disponíveis: ", paste(names(df), collapse = ", "))
  }
  
  df <- df %>%
    dplyr::mutate(
      .value = suppressWarnings(as.numeric(gsub(",", ".", as.character(.data[[vcol]]))))
    ) %>%
    dplyr::filter(is.finite(.year), is.finite(.value)) %>%
    dplyr::arrange(.year)
  
  if (nrow(df) < 8) {
    stop("Série muito curta em ", label, ": n = ", nrow(df), " (preciso de >= 8 aprox. para testes).")
  }
  
  # garante vetor numérico puro (corrige o seu erro do sens.slope)
  df$.value <- as.numeric(df$.value)
  df$.year  <- as.numeric(df$.year)
  
  # mantém metadados do arquivo
  attr(df, "label") <- label
  attr(df, "year_col") <- ycol
  attr(df, "value_col") <- vcol
  df
}

# ----------------------------
# 3) Função de análise (DW + HAC + MK/Sen)
# ----------------------------
analisar_serie <- function(df, label) {
  df <- df %>% dplyr::arrange(.year)
  
  # Modelo linear simples (tendência média de longo prazo)
  mod <- stats::lm(.value ~ .year, data = df)
  
  # β (pp/ano)
  beta <- unname(stats::coef(mod)[[".year"]])
  
  # p-valor com erros robustos HAC (Newey–West)
  # lag automático = floor(4*(n/100)^(2/9)) (convenção comum)
  n <- nrow(df)
  lag_nw <- max(1, floor(4 * (n/100)^(2/9)))
  ct_hac <- lmtest::coeftest(mod, vcov. = sandwich::NeweyWest(mod, lag = lag_nw, prewhite = FALSE, adjust = TRUE))
  p_hac  <- unname(ct_hac[".year", "Pr(>|t|)"])
  
  # Durbin–Watson (autocorrelação serial de 1ª ordem nos resíduos)
  # obs: dwtest assume ordenação temporal
  dw <- tryCatch(
    unname(lmtest::dwtest(mod)$statistic),
    error = function(e) NA_real_
  )
  
  # Mann–Kendall (tendência monotônica; robusto a não-normalidade)
  mk <- Kendall::MannKendall(df$.value)
  tau_mk <- unname(mk$tau)
  p_mk   <- unname(mk$sl)
  
  # Sen slope (magnitude da tendência monotônica)
  # trend::sens.slope requer x numérico -> garantido
  sen <- trend::sens.slope(df$.value)
  sen_slope <- unname(sen$estimates)  # por passo temporal (ano-a-ano) => pp/ano
  
  tibble::tibble(
    Unidade = label,
    Periodo = paste0(min(df$.year, na.rm = TRUE), "–", max(df$.year, na.rm = TRUE)),
    N = n,
    `β (pp/ano)` = beta,
    `p-valor HAC (NW)` = p_hac,
    `DW` = dw,
    `tau (MK)` = tau_mk,
    `p-valor (MK)` = p_mk,
    `Sen slope (pp/ano)` = sen_slope
  )
}

# ----------------------------
# 4) Rodar análises nas 4 planilhas/arquivos
# ----------------------------
dados_list <- list(
  "Bacia" = read_series_file(arq_bacia, "Bacia"),
  "FED"   = read_series_file(arq_fed,   "FED"),
  "FES"   = read_series_file(arq_fes,   "FES"),
  "FOD"   = read_series_file(arq_fod,   "FOD")
)

res <- purrr::imap_dfr(dados_list, ~ analisar_serie(.x, .y))

# Ordena de forma “didática”
res <- res %>%
  dplyr::mutate(Unidade = factor(Unidade, levels = c("Bacia","FED","FES","FOD"))) %>%
  dplyr::arrange(Unidade)

print(res)

# ----------------------------
# 5) Exportação ABNT (Excel + CSV)
# ----------------------------
# Formatação “ABNT”: sem firula, números com casas consistentes, cabeçalho destacado.
res_fmt <- res %>%
  dplyr::mutate(
    `β (pp/ano)` = round(`β (pp/ano)`, 3),
    `p-valor HAC (NW)` = signif(`p-valor HAC (NW)`, 3),
    `DW` = round(`DW`, 3),
    `tau (MK)` = round(`tau (MK)`, 3),
    `p-valor (MK)` = signif(`p-valor (MK)`, 3),
    `Sen slope (pp/ano)` = round(`Sen slope (pp/ano)`, 3)
  )

# CSV
utils::write.csv(res_fmt, out_csv, row.names = FALSE, fileEncoding = "UTF-8")

# Excel (openxlsx)
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Tabela_ABNT")

openxlsx::writeData(wb, "Tabela_ABNT", res_fmt, startRow = 1, startCol = 1)

# Estilos
header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)
body_style <- openxlsx::createStyle(halign = "center", valign = "center")

openxlsx::addStyle(wb, "Tabela_ABNT", header_style, rows = 1, cols = 1:ncol(res_fmt), gridExpand = TRUE)
openxlsx::addStyle(wb, "Tabela_ABNT", body_style, rows = 2:(nrow(res_fmt)+1), cols = 1:ncol(res_fmt), gridExpand = TRUE)

openxlsx::setColWidths(wb, "Tabela_ABNT", cols = 1:ncol(res_fmt), widths = "auto")
openxlsx::freezePane(wb, "Tabela_ABNT", firstRow = TRUE)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\n✅ OK!\n")
cat("Arquivo Excel (ABNT): ", out_xlsx, "\n", sep = "")
cat("Arquivo CSV: ", out_csv, "\n", sep = "")
cat("\nObs.: script usa apenas DW + HAC(Newey–West) + MK/Sen.\n")