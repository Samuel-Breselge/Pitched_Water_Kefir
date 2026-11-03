# 1. Settings and data -----------------------------------------------------
packages <- c("vegan", "readr", "here", "openxlsx")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))
input_file <- here::here("data", "Data_metabolites.csv")
output_file <- here::here("output", "RDA_summary.xlsx")
factors <- c("Yeast", "LAB", "AAB", "Bif", "other")
metabolites <- c("Sucrose", "Fructose", "Glucose", "Ethanol", "Glycerol",
                 "Mannitol", "Lactic_Acid", "Acetic_Acid")
# Keep TRUE only if NA means not detected in the following three columns.
replace_nd_with_zero <- TRUE
nd_columns <- c("Mannitol", "Lactic_Acid", "Acetic_Acid")

dat <- readr::read_csv(input_file, show_col_types = FALSE)
absent <- setdiff(c("sample", factors, metabolites), names(dat))
if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
if (!all(vapply(dat[metabolites], is.numeric, logical(1)))) {
  stop("Metabolite columns must be numeric; check ND/text entries.")
}
dat$sample <- trimws(as.character(dat$sample))
if (anyNA(dat$sample) || any(dat$sample == "") || anyDuplicated(dat$sample)) {
  stop("Sample IDs must be present and unique.")
}
for (f in factors) dat[[f]] <- trimws(as.character(dat[[f]]))
control <- tolower(dat$Yeast) %in% c("media", "wk grains")
if ("short_name" %in% names(dat)) {
  control <- control | tolower(trimws(as.character(dat$short_name))) %in% c("media", "wk grains")
}
dat <- dat[!control, , drop = FALSE]
if (!nrow(dat) || anyNA(dat[factors]) || any(as.matrix(dat[factors]) == "")) {
  stop("No defined communities, or missing biological-factor labels.")
}
if (any(tolower(as.matrix(dat[factors])) %in% c("media", "wk grains"))) {
  stop("Control labels remain among defined-community factors.")
}
nd_replacements <- 0L
for (m in metabolites) {
  replace <- replace_nd_with_zero & m %in% nd_columns & is.na(dat[[m]])
  nd_replacements <- nd_replacements + sum(replace)
  dat[[m]][replace] <- 0
  dat[[m]] <- dat[[m]] / 1000  # Raw micrograms/mL -> g/L.
}
complete <- apply(as.matrix(dat[metabolites]), 1, function(x) all(is.finite(x)))
excluded_ids <- dat$sample[!complete]
dat <- dat[complete, , drop = FALSE]
for (f in factors) dat[[f]] <- factor(dat[[f]])  # 'none' is a valid absence level.
if (nrow(dat) < 3 || any(vapply(dat[factors], nlevels, integer(1)) < 2)) {
  stop("Too few complete samples or factor levels.")
}
Y <- as.matrix(dat[metabolites])
if (any(Y < 0) || any(apply(Y, 2, sd) == 0)) {
  stop("Negative concentrations or constant metabolite columns; inspect data.")
}
Y_scaled <- scale(Y)  # Centre and scale each metabolite to unit SD; no log transform.
X <- model.matrix(~ Yeast + LAB + AAB + Bif + other, data = dat)
rank_X <- qr(X)$rank
if (rank_X < ncol(X)) stop("Biological factors are inseparable (rank-deficient design).")
if (nrow(X) - rank_X < 2) stop("Insufficient residual degrees of freedom.")

# 2. Full model and unique drop-one contributions --------------------------
fit <- vegan::rda(Y_scaled ~ Yeast + LAB + AAB + Bif + other, data = dat)
r2 <- vegan::RsquareAdj(fit)
labels <- c(Yeast = "Yeast identity", LAB = "LAB identity", AAB = "AAB identity",
            Bif = "Bifidobacterium identity", other = "Z. mobilis presence/absence")
summary_table <- data.frame(
  Term = "Full model", df = rank_X - 1L,
  Variance_percent = 100 * r2$r.squared, Adjusted_R2 = r2$adj.r.squared)
for (f in factors) {
  remaining <- setdiff(factors, f)
  reduced <- vegan::rda(as.formula(paste("Y_scaled ~", paste(remaining, collapse = " + "))), data = dat)
  reduced_rank <- qr(model.matrix(reformulate(remaining), data = dat))$rank
  summary_table <- rbind(summary_table, data.frame(
    Term = unname(labels[f]), df = rank_X - reduced_rank,
    Variance_percent = 100 * (r2$r.squared - vegan::RsquareAdj(reduced)$r.squared),
    Adjusted_R2 = NA_real_))
}

# 3. Export one compact supplementary table ---------------------------
notes <- c(
  "Full model: percentage of total standardized metabolic variation accounted for by all five factors.",
  "Factor rows: unique percentage of total variation, calculated as full-model R2 minus R2 after removing that factor. These percentages need not sum to full-model R2.",
  "Adjusted R2 is reported for the full model only; df = model degrees of freedom or degrees of freedom removed with a factor.",
  "Eight metabolites were centred and scaled to unit SD. Defined communities only; grain and media controls excluded. No interaction terms or set adjustment.",
  "Associations attributed to Bifidobacterium and Z. mobilis may include experimental-set variation.",
  paste0("NA-to-zero in mannitol/lactic acid/acetic acid: ", replace_nd_with_zero,
         "; values replaced before complete-case filtering: ", nd_replacements, "."),
  paste0("Samples excluded for incomplete/nonfinite metabolites: ", length(excluded_ids),
         if (length(excluded_ids)) paste0(" (", paste(excluded_ids, collapse = ", "), ")") else "."))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "RDA summary", gridLines = FALSE)
openxlsx::writeData(wb, 1, "RDA of defined water kefir communities", startRow = 1)
openxlsx::writeData(wb, 1, paste0("n = ", nrow(dat), " communities; ", length(metabolites),
                               " metabolites; residual df = ", nrow(X) - rank_X), startRow = 2)
openxlsx::writeDataTable(wb, 1, summary_table, startRow = 4,
                         tableStyle = "TableStyleMedium2", keepNA = FALSE)
openxlsx::addStyle(wb, 1, openxlsx::createStyle(numFmt = "0.0"),
                   rows = 5:(4 + nrow(summary_table)), cols = 3, gridExpand = TRUE)
openxlsx::addStyle(wb, 1, openxlsx::createStyle(numFmt = "0.000"),
                   rows = 5:(4 + nrow(summary_table)), cols = 4, gridExpand = TRUE)
openxlsx::setColWidths(wb, 1, cols = 1:4, widths = c(38, 9, 23, 18))
openxlsx::addStyle(wb, 1, openxlsx::createStyle(textDecoration = "bold", fontSize = 13), rows = 1, cols = 1)
for (i in seq_along(notes)) {
  row <- 5 + nrow(summary_table) + i
  openxlsx::mergeCells(wb, 1, cols = 1:4, rows = row)
  openxlsx::writeData(wb, 1, notes[i], startRow = row)
  openxlsx::addStyle(wb, 1, openxlsx::createStyle(wrapText = TRUE, valign = "top"), rows = row, cols = 1)
  openxlsx::setRowHeights(wb, 1, rows = row, heights = 44)
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
print(summary_table, row.names = FALSE)
message("Saved one summary table: ", output_file)

