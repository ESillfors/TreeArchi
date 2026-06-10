safe_filename <- function(x, max_chars = 120) {

  if (length(x) == 0 || is.na(x) || is.null(x)) {
    return("missing")
  }

  x <- as.character(x)

  x <- gsub("[^A-Za-z0-9_\\.-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  if (nchar(x) > max_chars) {
    x <- substr(x, 1, max_chars)
  }

  x
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    readr::write_csv(x, path),
    error = function(e) {
      alt <- file.path(dirname(path), paste0("ALT_", substr(basename(path), 1, 60)))
      readr::write_csv(x, alt)
      message("Could not write: ", path)
      message("Wrote alternative file instead: ", alt)
    }
  )
}

safe_write_lines <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    readr::write_lines(as.character(x), path),
    error = function(e) {
      alt <- file.path(dirname(path), paste0("ALT_", substr(basename(path), 1, 60)))
      readr::write_lines(as.character(x), alt)
      message("Could not write: ", path)
      message("Wrote alternative file instead: ", alt)
    }
  )
}

safe_ggsave <- function(filename, plot, width, height, dpi = 340, bg = "white") {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = bg),
    error = function(e) {
      alt <- file.path(dirname(filename), paste0("ALT_", substr(basename(filename), 1, 60)))
      ggplot2::ggsave(alt, plot, width = width, height = height, dpi = dpi, bg = bg)
      message("Could not save plot: ", filename)
      message("Saved alternative plot instead: ", alt)
    }
  )
}

read_messy_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  first_line <- readLines(path, n = 1, warn = FALSE)
  sep <- ifelse(
    length(strsplit(first_line, ";", fixed = TRUE)[[1]]) >
      length(strsplit(first_line, ",", fixed = TRUE)[[1]]),
    ";", ","
  )
  df <- read.table(
    path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM",
    check.names = FALSE
  )
  tibble::as_tibble(df)
}
