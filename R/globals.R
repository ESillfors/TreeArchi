set_analysis_paths <- function(data_path, out_base) {

  if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path)
  }

  if (!dir.exists(out_base)) {
    dir.create(out_base, recursive = TRUE)
  }

  list(
    data_path = data_path,
    out_base = out_base
  )
}
