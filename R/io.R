import_data <- function(file, verbose=FALSE) {
  if (verbose) cat(paste("\nReading data from", file, "\n"))
  d <- read_csv(file, col_types = cols(), na = c("", "NA", "n/a", "N/A", "None"), guess_max = 10000, progress = FALSE)
  if (verbose) cat(paste("   ", nrow(d), "rows,", ncol(d), "columns\n"))
  return(d)
}

read_excel_ph <- function(file, verbose=FALSE) {
  if (verbose) cat(paste("\nReading data from", file, "\n"))
  d <- readxl::read_excel(file, .name_repair = "minimal")

  coln <- colnames(d)
  dup <- which(duplicated(coln))
  if (length(dup) > 0) {
    stop(paste("Duplicated column names found:", paste(coln[dup], collapse = ", ")))
  }

  if (!("ID" %in% coln)) {
    stop("Column 'ID' not found in the input file. DDD identifiers must be stored in column called 'ID'.")
  }

  d <- d %>%
    rename(Name = ID)

  if (verbose) cat(paste("   ", nrow(d), "rows,", ncol(d), "columns\n"))
  d
}

