args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- file.path(project_root, "examples", "00_ACTUAL_UPLOAD_DATA")
if (!dir.exists(data_root)) stop("Missing examples/00_ACTUAL_UPLOAD_DATA: ", data_root, call. = FALSE)

read_y_species <- function(path) {
  y <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!ncol(y)) return(character())
  first <- y[[1]]
  first_num <- suppressWarnings(as.numeric(as.character(first)))
  first_is_id <- names(y)[1] %in% c("", "X", "X.1", "...1", "row.names", "rowname", "row_id", "id", "site_id", "sample_id") ||
    (!all(is.na(first_num)) && identical(as.integer(first_num), seq_len(length(first_num)))) ||
    all(is.na(first_num))
  if (ncol(y) > 1 && first_is_id) names(y)[-1] else names(y)
}

case_taxon_defaults <- function(case_name) {
  if (grepl("vascular_plants|plant", case_name, ignore.case = TRUE)) {
    list(kingdom = "Plantae", phylum = "Tracheophyta", class = "Magnoliopsida",
         orders = c("Ericales", "Fagales", "Ranunculales", "Polypodiales"),
         families = c("Ericaceae", "Fagaceae", "Berberidaceae", "Dryopteridaceae"))
  } else if (grepl("fungi", case_name, ignore.case = TRUE)) {
    list(kingdom = "Fungi", phylum = "Basidiomycota", class = "Agaricomycetes",
         orders = c("Agaricales", "Polyporales", "Russulales"),
         families = c("Agaricaceae", "Polyporaceae", "Russulaceae"))
  } else if (grepl("bird", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Chordata", class = "Aves",
         orders = c("Passeriformes", "Galliformes", "Charadriiformes"),
         families = c("Phylloscopidae", "Turdidae", "Tetraonidae", "Scolopacidae"))
  } else if (grepl("bryophyte", case_name, ignore.case = TRUE)) {
    list(kingdom = "Plantae", phylum = "Bryophyta", class = "Bryopsida",
         orders = c("Bryales", "Hypnales", "Dicranales", "Polytrichales"),
         families = c("Bryaceae", "Hypnaceae", "Dicranaceae", "Polytrichaceae"))
  } else if (grepl("macroinvertebrate", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Arthropoda", class = "Insecta",
         orders = c("Ephemeroptera", "Trichoptera", "Diptera", "Coleoptera"),
         families = c("Baetidae", "Hydropsychidae", "Chironomidae", "Dytiscidae"))
  } else if (grepl("beetle", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Arthropoda", class = "Insecta",
         orders = c("Coleoptera"), families = c("Carabidae", "Dytiscidae", "Hydrophilidae"))
  } else if (grepl("spider", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Arthropoda", class = "Arachnida",
         orders = c("Araneae"), families = c("Linyphiidae", "Salticidae", "Lycosidae", "Araneidae"))
  } else if (grepl("fish", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Chordata", class = "Actinopterygii",
         orders = c("Perciformes", "Gadiformes", "Pleuronectiformes", "Clupeiformes"),
         families = c("Sparidae", "Gadidae", "Pleuronectidae", "Clupeidae"))
  } else if (grepl("mammal|phylacine", case_name, ignore.case = TRUE)) {
    list(kingdom = "Animalia", phylum = "Chordata", class = "Mammalia",
         orders = c("Carnivora", "Rodentia", "Primates", "Artiodactyla"),
         families = c("Canidae", "Muridae", "Cercopithecidae", "Bovidae"))
  } else {
    list(kingdom = "Synthetic", phylum = "Benchmarkophyta", class = "Benchmarkopsida",
         orders = c("Benchmarkales"), families = c("Benchmarkaceae"))
  }
}

make_taxonomy <- function(species, case_name) {
  d <- case_taxon_defaults(case_name)
  genus <- sub("_.*$", "", species)
  genus[!nzchar(genus) | genus == species] <- paste0("Genus_", seq_along(species))[!nzchar(genus) | genus == species]
  tax <- data.frame(
    taxon_id = species,
    kingdom = d$kingdom,
    phylum = d$phylum,
    class = d$class,
    order = rep(d$orders, length.out = length(species)),
    family = rep(d$families, length.out = length(species)),
    genus = genus,
    species = species,
    stringsAsFactors = FALSE
  )
  rownames(tax) <- species
  tax
}

write_upload_map <- function(engine_dir, engine) {
  map_path <- file.path(engine_dir, "UPLOAD_MAP.csv")
  if (!file.exists(map_path)) return(FALSE)
  m <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)
  m <- m[!tolower(m$file_to_choose) %in% c("phylo_cov.csv", "c_phylo.csv"), , drop = FALSE]
  if (engine == "Hmsc") {
    if (!any(m$file_to_choose == "phylogeny.nwk")) {
      m <- rbind(m, data.frame(ui_upload_box = "phylogeny.nwk or taxonomy.csv", file_to_choose = "phylogeny.nwk", present = TRUE,
                               note = "Newick tree option; species tips must match Y columns.", stringsAsFactors = FALSE))
    } else {
      m$ui_upload_box[m$file_to_choose == "phylogeny.nwk"] <- "phylogeny.nwk or taxonomy.csv"
      m$note[m$file_to_choose == "phylogeny.nwk"] <- "Newick tree option; species tips must match Y columns."
    }
    if (!any(m$file_to_choose == "taxonomy.csv")) {
      m <- rbind(m, data.frame(ui_upload_box = "phylogeny.nwk or taxonomy.csv", file_to_choose = "taxonomy.csv", present = TRUE,
                               note = "Taxonomy rank-table alternative; JSDM Studio derives the species correlation matrix internally.", stringsAsFactors = FALSE))
    }
  } else if (engine == "Hmsc-HPC") {
    if (!any(m$file_to_choose == "taxonomy.csv")) {
      m <- rbind(m, data.frame(ui_upload_box = "taxonomy.csv", file_to_choose = "taxonomy.csv", present = TRUE,
                               note = "Taxonomy rank-table input; JSDM Studio derives the pyhmsc species matrix internally.", stringsAsFactors = FALSE))
    }
    if ("phylo_tree.nwk" %in% m$file_to_choose) {
      m$note[m$file_to_choose == "phylo_tree.nwk"] <- "Newick tree alternative; species tips must match Y columns."
    }
  }
  write.csv(m, map_path, row.names = FALSE, fileEncoding = "UTF-8")
  TRUE
}

case_dirs <- list.dirs(data_root, full.names = TRUE, recursive = FALSE)
case_dirs <- case_dirs[grepl("^[0-9]+_", basename(case_dirs))]
rows <- list()
for (case_dir in case_dirs) {
  hmsc_y <- file.path(case_dir, "Hmsc", "Y.csv")
  if (!file.exists(hmsc_y)) next
  species <- read_y_species(hmsc_y)
  tax <- make_taxonomy(species, basename(case_dir))
  for (engine in c("Hmsc", "Hmsc-HPC")) {
    engine_dir <- file.path(case_dir, engine)
    if (!dir.exists(engine_dir)) next
    write.csv(tax, file.path(engine_dir, "taxonomy.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    old_matrix_files <- file.path(engine_dir, c("phylo_cov.csv", "c_phylo.csv"))
    old_matrix_files <- old_matrix_files[file.exists(old_matrix_files)]
    if (length(old_matrix_files)) file.remove(old_matrix_files)
    write_upload_map(engine_dir, engine)
    rows[[length(rows) + 1L]] <- data.frame(
      case_folder = basename(case_dir),
      engine = engine,
      taxonomy_file = normalizePath(file.path(engine_dir, "taxonomy.csv"), winslash = "/", mustWork = FALSE),
      species = length(species),
      rank_columns = "kingdom, phylum, class, order, family, genus, species",
      upload_contract = "raw phylogeny input is Newick tree or taxonomy rank table; species matrix is derived internally",
      stringsAsFactors = FALSE
    )
  }
}

audit <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(audit, file.path(data_root, "TAXONOMY_UPLOAD_CONTRACT_AUDIT.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  "# Taxonomy Upload Contract Audit",
  "",
  "JSDM Studio accepts two raw phylogeny inputs: a Newick tree or a taxonomy rank table.",
  "",
  "Users should not upload a species-by-species covariance/correlation matrix as original data. When an engine needs that matrix, JSDM Studio derives it internally and saves `derived_phylogeny_correlation.csv` in the output.",
  "",
  paste0("Taxonomy files refreshed: ", nrow(audit))
), file.path(data_root, "TAXONOMY_UPLOAD_CONTRACT_AUDIT.md"), useBytes = TRUE)

message("Refreshed taxonomy upload contract rows: ", nrow(audit))
