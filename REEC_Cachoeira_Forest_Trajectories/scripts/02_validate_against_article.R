# Confere os resultados numéricos calculados em outputs/ com reference/.
# Requer a execução prévia de scripts/01_reproduce_article.R.

repo_root <- getwd()
if (!file.exists(file.path(repo_root, "reference", "Table1_Trend_REEC.csv"))) {
  if (file.exists(file.path(dirname(repo_root), "reference", "Table1_Trend_REEC.csv")))
    repo_root <- dirname(repo_root)
}
if (!file.exists(file.path(repo_root, "reference", "Table1_Trend_REEC.csv"))) {
  stop("Open the project root before validating results")
}
read_original <- function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
expected1 <- read_original(file.path(repo_root,"reference","Table1_Trend_REEC.csv"))
expected2 <- read_original(file.path(repo_root,"reference","Table2_PCA_REEC.csv"))
actual1 <- read_original(file.path(repo_root,"outputs","Table1_Trend_REEC.csv"))
actual2 <- read_original(file.path(repo_root,"outputs","Table2_PCA_REEC.csv"))

compare <- function(expected, actual, col, tolerance = 5.01e-4, relative = FALSE) {
  e <- as.numeric(expected[[col]])
  a <- as.numeric(actual[[col]])
  if (!length(e) || !length(a) || length(e) != length(a) || anyNA(e) || anyNA(a))
    stop("Missing/incompatible values for: ", col)
  dif <- abs(e - a)
  limit <- if (relative) pmax(abs(e) * 0.012, 1e-30) else tolerance
  if (any(dif > limit)) {
    bad <- which(dif > limit)
    stop("Mismatch in ", col, " at rows: ", paste(bad, collapse=","),
         " [expected=", paste(e[bad],collapse=","),
         "; actual=", paste(a[bad],collapse=","), "]")
  }
  cat(sprintf("OK %-45s max_diff=%g\n", col, max(dif)))
}

stopifnot(identical(expected1[["Spatial unit"]], actual1[["Spatial unit"]]),
          identical(expected2[["Spatial unit"]], actual2[["Spatial unit"]]))
stopifnot(all(actual1$n == 38), all(actual1$Period == "1985-2022"))

for (col in c("beta (percentage points yr-1)", "Durbin-Watson",
              "Mann-Kendall tau", "Sen slope (percentage points yr-1)", "R-squared")) {
  compare(expected1, actual1, col)
}
for (col in c("HAC p-value", "Mann-Kendall p-value")) {
  compare(expected1, actual1, col, relative=TRUE)
}
for (col in c("PC1 variance (%)", "PC2 variance (%)", "Cumulative variance (%)")) {
  compare(expected2, actual2, col, tolerance=0.05001)
}
figure_filenames <- c("Fig3_ForestCover_Trends_REEC_v5_FINAL.tif",
                      "Fig4_Annual_Forest_Change_REEC_v5_FINAL.tif",
                      "Fig5_PCA_Multipanel_REEC_v5_FINAL.tif")
for (f in figure_filenames) {
  p <- file.path(repo_root,"outputs", f)
  if (!file.exists(p) || file.info(p)$size < 20000) stop("Figure absent or too small: ", p)
  cat("OK figure ", f, "\n", sep="")
}
cat("\nALL NUMERIC VALIDATIONS PASSED.\n")
writeLines(c(paste("Validation date:", Sys.time()), capture.output(sessionInfo()),
             "ALL NUMERIC VALIDATIONS PASSED"),
           file.path(repo_root,"outputs","validation_PASS.txt"))
