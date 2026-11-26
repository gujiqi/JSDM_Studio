# Presentation preferences are independent of analysis configs and output bundles.
studio_languages <- c(
  en = "English", "zh-CN" = "\u7b80\u4f53\u4e2d\u6587", fi = "Suomi",
  fr = "Fran\u00e7ais", ko = "\ud55c\uad6d\uc5b4", ja = "\u65e5\u672c\u8a9e",
  es = "Espa\u00f1ol", "pt-BR" = "Portugu\u00eas (Brasil)", de = "Deutsch",
  ru = "\u0420\u0443\u0441\u0441\u043a\u0438\u0439", ar = "\u0627\u0644\u0639\u0631\u0628\u064a\u0629",
  hi = "\u0939\u093f\u0928\u094d\u0926\u0940", bn = "\u09ac\u09be\u0982\u09b2\u09be",
  id = "Bahasa Indonesia", vi = "Ti\u1ebfng Vi\u1ec7t"
)

studio_preference_defaults <- list(locale = "en", textSize = "standard", inlineHelp = "translated")

studio_validate_preferences <- function(value) {
  if (!is.list(value)) return(studio_preference_defaults)
  choices <- list(locale = names(studio_languages), textSize = c("standard", "large", "largest"),
                  inlineHelp = c("translated", "english", "hidden"))
  setNames(lapply(names(choices), function(key) {
    item <- value[[key]]
    if (is.character(item) && length(item) == 1L && item %in% choices[[key]]) item
    else studio_preference_defaults[[key]]
  }), names(choices))
}

studio_preferences_file <- function() {
  override <- Sys.getenv("JSDM_STUDIO_PREFERENCES_FILE")
  if (nzchar(override)) override else file.path(tools::R_user_dir("JSDMStudio", "config"), "interface.json")
}

studio_read_preferences <- function(path = studio_preferences_file()) {
  if (!file.exists(path)) return(studio_preference_defaults)
  studio_validate_preferences(tryCatch(jsonlite::read_json(path), error = function(e) NULL))
}

studio_write_preferences <- function(value, path = studio_preferences_file()) {
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) return(FALSE)
  staged <- tempfile("studio-preferences-", tmpdir = directory)
  on.exit(unlink(staged), add = TRUE)
  tryCatch({
    jsonlite::write_json(studio_validate_preferences(value), staged, auto_unbox = TRUE, pretty = TRUE)
    suppressWarnings(file.rename(staged, path))
  }, error = function(e) FALSE)
}

studio_bind_preferences <- function(input, session) {
  observeEvent(input$studio_preference_request, {
    session$sendCustomMessage("studio.preferences.load", studio_read_preferences())
  })
  observeEvent(input$studio_preference_change, {
    event <- input$studio_preference_change
    if (!is.list(event) || !is.list(event$prefs)) return()
    revision <- event$revision
    if (!is.numeric(revision) || length(revision) != 1L || !is.finite(revision)) return()
    ok <- studio_write_preferences(event$prefs)
    session$sendCustomMessage("studio.preferences.saved", list(ok = ok, revision = revision))
  }, ignoreNULL = TRUE)
}

studio_language_options <- function() {
  lapply(names(studio_languages), function(locale)
    tags$option(value = locale, lang = locale, dir = if (locale == "ar") "rtl" else "ltr",
                unname(studio_languages[[locale]])))
}

studio_preferences_ui <- function() {
  tagList(
    div(id = "studio-toolbar", class = "studio-toolbar",
      tags$label(class = "studio-toolbar-label", `for` = "studio-language", "Language"),
      tags$select(id = "studio-language", class = "form-select form-select-sm",
        `aria-label` = "Interface language", translate = "no", studio_language_options()),
      tags$button(id = "studio-open-settings", class = "btn studio-settings-toggle", type = "button",
        title = "Settings", `aria-label` = "Settings", `aria-haspopup` = "dialog",
        studio_symbol("tuning", mode = "glyph"))),
    tags$dialog(id = "studio-settings", class = "studio-settings", `aria-labelledby` = "studio-settings-title",
      tags$header(
        h2(id = "studio-settings-title", "Settings"),
        tags$button(class = "btn studio-settings-close", type = "button", `aria-label` = "Close",
          `data-studio-dismiss` = NA, icon("times"))),
      div(class = "studio-settings-fields",
        div(
          tags$label(`for` = "studio-settings-language", "Interface language"),
          tags$select(id = "studio-settings-language", class = "form-select", translate = "no", studio_language_options())),
        div(
          tags$label(`for` = "studio-settings-size", "Text size"),
          tags$select(id = "studio-settings-size", class = "form-select",
            tags$option(value = "standard", "Standard"),
            tags$option(value = "large", "Large"),
            tags$option(value = "largest", "Extra large"))),
        div(
          tags$label(`for` = "studio-settings-help", "Inline parameter help"),
          tags$select(id = "studio-settings-help", class = "form-select",
            tags$option(value = "translated", "In the interface language"),
            tags$option(value = "english", "Original English"),
            tags$option(value = "hidden", "Hidden")))),
      p(class = "studio-settings-note", "Model formulas, input names, package logs and scientific exports keep their original language."),
      p(id = "studio-settings-status", class = "studio-settings-note", role = "status"),
      tags$footer(
        tags$button(id = "studio-reset-settings", class = "btn btn-outline-secondary", type = "button", "Reset to defaults"),
        tags$button(class = "btn btn-primary", type = "button", `data-studio-dismiss` = NA, "Done"))))
}
