# Read the same Shiny label and help tags that are rendered beside each control.
studio_parameter_reference <- function(ui) {
  panes <- htmltools::tagQuery(ui)$find(".tab-pane")$selectedTags()
  rows <- list()

  has_class <- function(tag, value) {
    inherits(tag, "shiny.tag") && value %in% strsplit(tag$attribs$class %||% "", " +")[[1L]]
  }
  tag_text <- function(x) {
    if (inherits(x, "shiny.tag")) return(tag_text(x$children))
    if (is.list(x)) return(paste(vapply(x, tag_text, character(1)), collapse = " "))
    if (is.character(x)) return(paste(x, collapse = " "))
    ""
  }
  text_of <- function(x) trimws(gsub("[[:space:]]+", " ", tag_text(x)))

  walk <- function(x, workflow) {
    children <- if (inherits(x, "shiny.tag")) x$children else if (is.list(x)) x else return()
    for (i in seq_along(children)) {
      tag <- children[[i]]
      if (has_class(tag, "sync-pair") || has_class(tag, "shiny-input-container")) {
        q <- htmltools::tagQuery(tag)
        labels <- q$find("label")$selectedTags()
        id <- tag$attribs[["data-sync-id"]] %||% tag$attribs$id
        if (is.null(id) && length(labels)) id <- labels[[1L]]$attribs[["for"]]
        if (is.null(id)) {
          fields <- c(q$find("input")$selectedTags(), q$find("textarea")$selectedTags(),
                      q$find("select")$selectedTags())
          if (length(fields)) id <- fields[[1L]]$attribs$id
        }
        helps <- q$find(".param-help")$selectedTags()
        next_tag <- if (i < length(children)) children[[i + 1L]] else NULL
        explanation <- if (length(helps)) text_of(helps[[1L]]) else if (has_class(next_tag, "param-help")) text_of(next_tag) else ""
        label <- if (length(labels)) text_of(labels[[1L]]) else ""
        if (length(id) && nzchar(id) && nzchar(label) && nzchar(explanation)) {
          rows[[length(rows) + 1L]] <<- data.frame(workflow = workflow, input_id = id,
                                                  label = label, explanation = explanation)
        }
      } else {
        walk(tag, workflow)
      }
    }
  }

  for (pane in panes) {
    workflow <- pane$attribs[["data-value"]]
    if (!is.null(workflow)) walk(pane, sub(" Workflow$", "", workflow))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    workflow = "Hmsc", input_id = "hmsc_group_col_*", label = "Additional grouping level",
    explanation = "Select a distinct studyDesign.csv column for each additional nonspatial level."
  )
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$input_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}
