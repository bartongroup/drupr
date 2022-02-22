import_data <- function(file, verbose=FALSE) {
  if(verbose) cat(paste("\nReading data from", file, "\n"))
  d <- read_csv(file, col_types = cols(), na=c("", "NA", "n/a", "N/A", "None"), guess_max = 10000, progress = FALSE)
  if(verbose) cat(paste("   ", nrow(d), "rows,", ncol(d), "columns\n"))
  return(d)
}



write_res <- function(d, odir, file) {
  d %>%
    mutate_if(is.numeric, ~signif(.x, 4)) %>%
    write_csv(file.path(odir, file))
}

write_results <- function(pred, out_dir) {
  if(!dir.exists(out_dir)) dir.create(out_dir)

  write_res(pred$missing, out_dir, "missing.csv")
  write_res(pred$predictions, out_dir, "predictions.csv")
  write_res(pred$importance, out_dir, "importance.csv")
  write_res(pred$train_variables, out_dir, "train_variables.csv")
  write_res(pred$variable_comparison, out_dir, "variable_comparison.csv")
  write_res(pred$mismatched_levels, out_dir, "mismatched_levels.csv")
}

read_excel_ph <- function(file, verbose=FALSE) {
  if(verbose) cat(paste("\nReading data from", file, "\n"))
  d <- readxl::read_excel(file) %>%
    rename(Name = ID)
  if(verbose) cat(paste("   ", nrow(d), "rows,", ncol(d), "columns\n"))
  d
}

