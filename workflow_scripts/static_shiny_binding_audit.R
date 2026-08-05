txt <- paste(readLines("app.R", warn = FALSE), collapse = "\n")

rx <- function(pattern) {
  out <- regmatches(txt, gregexpr(pattern, txt, perl = TRUE))[[1]]
  if (length(out) == 1 && identical(out, "-1")) character() else out
}

get_id <- function(matches) {
  unique(sub(".*?[\"']([^\"']+)[\"'].*", "\\1", matches, perl = TRUE))
}

input_call_pattern <- paste0(
  "(fileInput|textInput|numericInput|selectInput|checkboxInput|sliderInput|",
  "actionButton|radioButtons|checkboxGroupInput|selectizeInput|textAreaInput|synced_numeric_slider)",
  "\\s*\\(\\s*[\"'][^\"']+[\"']"
)
output_call_pattern <- paste0(
  "(textOutput|verbatimTextOutput|DTOutput|plotOutput|uiOutput|downloadButton)",
  "\\s*\\(\\s*[\"'][^\"']+[\"']"
)

ui_ids <- get_id(rx(input_call_pattern))
output_ui_ids <- get_id(rx(output_call_pattern))
download_ui_ids <- get_id(rx("downloadButton\\s*\\(\\s*[\"'][^\"']+[\"']"))

input_refs <- unique(c(
  sub("input\\$([A-Za-z0-9_.]+).*", "\\1", rx("input\\$[A-Za-z0-9_.]+"), perl = TRUE),
  sub(".*input\\[\\[[\"']([^\"']+)[\"']\\]\\].*", "\\1", rx("input\\[\\[[\"'][^\"']+[\"']\\]\\]"), perl = TRUE)
))
output_handlers <- unique(sub("output\\$([A-Za-z0-9_.]+).*", "\\1", rx("output\\$[A-Za-z0-9_.]+\\s*<-"), perl = TRUE))
download_handlers <- unique(sub("output\\$([A-Za-z0-9_.]+).*", "\\1", rx("output\\$[A-Za-z0-9_.]+\\s*<-\\s*downloadHandler"), perl = TRUE))
render_handlers <- unique(sub("output\\$([A-Za-z0-9_.]+).*", "\\1", rx("output\\$[A-Za-z0-9_.]+\\s*<-\\s*render[A-Za-z]+"), perl = TRUE))

# file_preview_block() creates UI outputs with paste0(prefix, ...), and
# register_engine_file_preview() creates matching dynamic output handlers.
prefixes <- c("hmsc", "hmschpc", "jsdm", "gjam", "spocc", "sjsdm", "boral")
preview_input_ids <- paste0(prefixes, "_file_preview_choice")
preview_output_ids <- as.vector(outer(
  prefixes,
  c("_file_preview_selector", "_file_preview_meta", "_file_preview_table"),
  paste0
))
ui_ids <- unique(c(ui_ids, preview_input_ids))
input_refs <- unique(c(input_refs, preview_input_ids))
output_ui_ids <- unique(c(output_ui_ids, preview_output_ids))
output_handlers <- unique(c(output_handlers, preview_output_ids))
render_handlers <- unique(c(render_handlers, preview_output_ids))

refs_not_in_ui <- setdiff(input_refs, ui_ids)
ui_not_read <- setdiff(ui_ids, input_refs)
outputs_without_handler <- setdiff(output_ui_ids, output_handlers)
handlers_without_ui <- setdiff(output_handlers, output_ui_ids)
downloads_without_handler <- setdiff(download_ui_ids, download_handlers)
download_handlers_without_button <- setdiff(download_handlers, download_ui_ids)
non_download_without_render <- setdiff(setdiff(output_ui_ids, download_ui_ids), render_handlers)

cat("STRICT_SHINY_BINDING_AUDIT\n")
cat("ui_ids=", length(ui_ids), " input_refs=", length(input_refs),
    " refs_not_in_ui=", length(refs_not_in_ui),
    " ui_not_read=", length(ui_not_read), "\n", sep = "")
if (length(refs_not_in_ui)) cat("refs_not_in_ui:\n", paste(refs_not_in_ui, collapse = "\n"), "\n", sep = "")
if (length(ui_not_read)) cat("ui_not_read:\n", paste(ui_not_read, collapse = "\n"), "\n", sep = "")

cat("output_ui_ids=", length(output_ui_ids), " output_handlers=", length(output_handlers),
    " outputs_without_handler=", length(outputs_without_handler),
    " handlers_without_ui=", length(handlers_without_ui),
    " non_download_without_render=", length(non_download_without_render), "\n", sep = "")
if (length(outputs_without_handler)) cat("outputs_without_handler:\n", paste(outputs_without_handler, collapse = "\n"), "\n", sep = "")
if (length(handlers_without_ui)) cat("handlers_without_ui:\n", paste(handlers_without_ui, collapse = "\n"), "\n", sep = "")
if (length(non_download_without_render)) cat("non_download_without_render:\n", paste(non_download_without_render, collapse = "\n"), "\n", sep = "")

cat("download_buttons=", length(download_ui_ids), " download_handlers=", length(download_handlers),
    " downloads_without_handler=", length(downloads_without_handler),
    " handlers_without_button=", length(download_handlers_without_button), "\n", sep = "")
if (length(downloads_without_handler)) cat("downloads_without_handler:\n", paste(downloads_without_handler, collapse = "\n"), "\n", sep = "")
if (length(download_handlers_without_button)) cat("download_handlers_without_button:\n", paste(download_handlers_without_button, collapse = "\n"), "\n", sep = "")

bad_download_content <- character()
for (id in download_handlers) {
  pattern <- paste0(
    "output\\$", id, "\\s*<-\\s*downloadHandler[\\s\\S]*?",
    "(?=\\n\\s*output\\$|\\n\\s*observeEvent|\\n\\s*#|\\n\\})"
  )
  block <- rx(pattern)[1]
  if (!length(block) || is.na(block)) next
  ok <- grepl("copy_zip_to_download\\(", block) ||
    grepl("write\\.csv\\(", block) ||
    grepl("file\\.copy\\(", block)
  if (!ok) bad_download_content <- c(bad_download_content, id)
}

cat("download_content_static_bad=", length(bad_download_content), "\n", sep = "")
if (length(bad_download_content)) {
  cat("bad_download_content:\n", paste(bad_download_content, collapse = "\n"), "\n", sep = "")
}

failures <- length(refs_not_in_ui) +
  length(outputs_without_handler) +
  length(handlers_without_ui) +
  length(non_download_without_render) +
  length(downloads_without_handler) +
  length(download_handlers_without_button) +
  length(bad_download_content)

if (failures) stop("Strict Shiny binding audit failed with ", failures, " issue(s).", call. = FALSE)
cat("STRICT_AUDIT_OK\n")
