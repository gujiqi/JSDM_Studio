#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else NA_character_
root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE) else getwd()
setwd(root)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(root, "output", paste0("phylogeny_taxonomy_contract_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

record <- function(check, status, detail = "", path = "") {
  data.frame(
    check = check,
    status = status,
    detail = detail,
    path = path,
    stringsAsFactors = FALSE
  )
}

rows <- list()
add_row <- function(check, status, detail = "", path = "") {
  rows[[length(rows) + 1L]] <<- record(check, status, detail, path)
}

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

drop_id_col <- function(dat) {
  if (!ncol(dat)) return(dat)
  first <- tolower(names(dat)[1])
  first_chr <- as.character(dat[[1]])
  first_num <- suppressWarnings(as.numeric(first_chr))
  sequence_index <- all(!is.na(first_num)) && identical(as.integer(first_num), seq_len(length(first_num)))
  text_id <- !all(!is.na(first_num)) && !anyDuplicated(first_chr)
  id_names <- c("", "x", "x.1", "...1", "row.names", "rowname", "row_id", "id", "site", "site_id", "sample", "sample_id", "species", "species_id", "taxon", "taxon_id")
  if (first %in% id_names || sequence_index || text_id) dat <- dat[-1]
  dat
}

species_from_y <- function(path) {
  y <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(drop_id_col(y))
}

safe_source_adapter <- function() {
  source(file.path(root, "R", "hmschpc_adapter.R"), local = FALSE)
}

tryCatch({
  parse(file.path(root, "app.R"))
  add_row("app_parse", "pass", "app.R parsed successfully.", file.path(root, "app.R"))
}, error = function(e) {
  add_row("app_parse", "fail", conditionMessage(e), file.path(root, "app.R"))
})

tryCatch({
  safe_source_adapter()
  add_row("hmschpc_adapter_source", "pass", "R/hmschpc_adapter.R sourced successfully.", file.path(root, "R", "hmschpc_adapter.R"))
}, error = function(e) {
  add_row("hmschpc_adapter_source", "fail", conditionMessage(e), file.path(root, "R", "hmschpc_adapter.R"))
})

app_txt <- read_text(file.path(root, "app.R"))
add_row(
  "no_hmschpc_phylo_cov_upload_id",
  if (grepl("hmschpc_phylo_cov_file", app_txt, fixed = TRUE)) "fail" else "pass",
  "Hmsc-HPC GUI should expose taxonomy.csv and/or phylo_tree.nwk, not phylo_cov.csv."
)
add_row(
  "no_hmsc_c_file_gui_input",
  if (grepl("input\\$hmsc_C_file", app_txt)) "fail" else "pass",
  "Hmsc GUI no longer reads a user-supplied C matrix input."
)

upload_roots <- c(
  file.path(root, "examples", "00_ACTUAL_UPLOAD_DATA"),
  file.path(root, "examples", "00_UPLOAD_READY_book_chapter_cases")
)
for (upload_root in upload_roots) {
  old_files <- if (dir.exists(upload_root)) list.files(upload_root, pattern = "^phylo_cov\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  add_row(
    paste0("no_user_phylo_cov_files_", basename(upload_root)),
    if (length(old_files)) "fail" else "pass",
    paste0("Found ", length(old_files), " user-upload phylo_cov.csv files."),
    upload_root
  )
}

taxonomy_files <- unique(c(
  list.files(file.path(root, "examples", "00_ACTUAL_UPLOAD_DATA"), pattern = "^taxonomy\\.csv$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(root, "examples", "00_UPLOAD_READY_book_chapter_cases"), pattern = "^taxonomy\\.csv$", recursive = TRUE, full.names = TRUE)
))

if (!length(taxonomy_files)) {
  add_row("taxonomy_files_present", "fail", "No taxonomy.csv files found in example upload data.")
} else {
  add_row("taxonomy_files_present", "pass", paste0("Found ", length(taxonomy_files), " taxonomy.csv files."))
}

for (tax_path in taxonomy_files) {
  dir_path <- dirname(tax_path)
  y_path <- file.path(dir_path, "Y.csv")
  if (!file.exists(y_path)) {
    add_row("taxonomy_has_matching_y", "skip", "No Y.csv beside taxonomy.csv; this may be a common-input folder.", tax_path)
    next
  }
  tryCatch({
    sp <- species_from_y(y_path)
    tax <- utils::read.csv(tax_path, check.names = FALSE, stringsAsFactors = FALSE)
    tax_cols <- tolower(names(tax))
    rank_ok <- all(c("kingdom", "phylum", "class", "order", "family", "genus", "species") %in% tax_cols)
    if (!rank_ok) stop("taxonomy.csv must contain kingdom/phylum/class/order/family/genus/species columns.")
    C <- hmschpc_taxonomy_to_correlation(tax, sp)
    ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    add_row(
      "taxonomy_to_correlation",
      if (nrow(C) == length(sp) && min(ev) > 0) "pass" else "fail",
      paste0("species=", length(sp), "; matrix=", nrow(C), "x", ncol(C), "; min_eigen=", signif(min(ev), 4)),
      tax_path
    )
  }, error = function(e) {
    add_row("taxonomy_to_correlation", "fail", conditionMessage(e), tax_path)
  })
}

tree_files <- unique(c(
  list.files(file.path(root, "examples", "00_ACTUAL_UPLOAD_DATA"), pattern = "^(phylogeny|phylo_tree)\\.nwk$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(root, "examples", "00_UPLOAD_READY_book_chapter_cases"), pattern = "^(phylogeny|phylo_tree)\\.nwk$", recursive = TRUE, full.names = TRUE)
))

for (tree_path in tree_files) {
  dir_path <- dirname(tree_path)
  y_path <- file.path(dir_path, "Y.csv")
  if (!file.exists(y_path)) next
  tryCatch({
    sp <- species_from_y(y_path)
    C <- hmschpc_newick_to_correlation(tree_path, sp)
    ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    add_row(
      "newick_to_correlation",
      if (nrow(C) == length(sp) && min(ev) > 0) "pass" else "fail",
      paste0("species=", length(sp), "; matrix=", nrow(C), "x", ncol(C), "; min_eigen=", signif(min(ev), 4)),
      tree_path
    )
  }, error = function(e) {
    add_row("newick_to_correlation", "fail", conditionMessage(e), tree_path)
  })
}

tryCatch({
  sample_tax <- taxonomy_files[file.exists(file.path(dirname(taxonomy_files), "Y.csv"))][1]
  if (is.na(sample_tax)) stop("No taxonomy file with adjacent Y.csv found.")
  sp <- species_from_y(file.path(dirname(sample_tax), "Y.csv"))
  tax <- utils::read.csv(sample_tax, check.names = FALSE, stringsAsFactors = FALSE)
  C <- hmschpc_taxonomy_to_correlation(tax, sp)
  test_out <- file.path(out_dir, "derived_output_contract")
  hmschpc_write_derived_phylogeny(test_out, C = C, source = "taxonomy", input_file = sample_tax)
  expected <- file.path(test_out, c(
    "data/derived_phylogeny_correlation.csv",
    "tables/derived_phylogeny_correlation.csv",
    "standard/derived_phylogeny_correlation.csv",
    "tables/phylogeny_input_contract.csv"
  ))
  add_row(
    "derived_matrix_output_contract",
    if (all(file.exists(expected))) "pass" else "fail",
    paste0("Expected output files present: ", sum(file.exists(expected)), "/", length(expected)),
    test_out
  )
}, error = function(e) {
  add_row("derived_matrix_output_contract", "fail", conditionMessage(e), out_dir)
})

validation <- do.call(rbind, rows)
csv_path <- file.path(out_dir, "phylogeny_taxonomy_contract_validation.csv")
utils::write.csv(validation, csv_path, row.names = FALSE, fileEncoding = "UTF-8")

md_path <- file.path(out_dir, "phylogeny_taxonomy_contract_validation.md")
status_counts <- as.data.frame(table(validation$status), stringsAsFactors = FALSE)
names(status_counts) <- c("status", "n")
md <- c(
  "# Phylogeny and Taxonomy Input Contract Validation",
  "",
  paste0("Project: `", root, "`"),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Contract",
  "",
  "- User uploads either a real Newick tree or a taxonomy rank table.",
  "- User does not upload a species-by-species phylogenetic covariance/correlation matrix.",
  "- JSDM Studio derives the correlation matrix internally and writes it to output `data/`, `tables/` and `standard/`.",
  "",
  "## Status Counts",
  "",
  paste(capture.output(print(status_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "## Detailed CSV",
  "",
  paste0("`", csv_path, "`")
)
writeLines(md, md_path, useBytes = TRUE)

message("VALIDATION_CSV=", csv_path)
message("VALIDATION_MD=", md_path)
if (any(validation$status == "fail")) quit(status = 1)
