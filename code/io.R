atomic_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  use_gzip <- grepl("[.]gz$", path, ignore.case = TRUE)
  temporary <- tempfile(
    tmpdir = dirname(path),
    fileext = if (use_gzip) ".tsv.gz" else ".tsv"
  )
  data.table::fwrite(
    x,
    temporary,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (use_gzip) "gzip" else "none"
  )
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
  invisible(path)
}
