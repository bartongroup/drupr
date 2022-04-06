#' @import tidyverse fitdistrplus randomForest glue readxl
#' @importFrom stats as.dist cor cutree hclust median na.omit predict quantile sd
#' @importFrom utils capture.output head setTxtProgressBar tail txtProgressBar


select <- dplyr::select

#' Drug properties prediction
#'
#' @param train_file Training data file. A CSV file with training data.
#' @param info_file Training data variable info file. A CSV file with two columns: 'variable' and 'type'. The first column contains variable names, the second 'id' for identifier column, 'response' for response variables or 'none' for variables not to be used.
#' @param test_file Test data file.
#' @param verbose Logical. If TRUE, progress information will be outputted to the console.
#' @param min_unique Minimum number of unique levels for a variable to be considered.
#' @param min_good Minimum number of non-missing values for a variable to be considered.
#' @param max_cat_levels Maximum number of levels (unique values) for a categorical variable to be considered.
#'
#' @return A list with tibbles containing results.
#' @export
ddu_prediction <- function(train_file, info_file, test_file, verbose = TRUE, min_unique = 2, min_good = 1500, max_cat_levels = 10) {

  train_raw <- import_data(train_file, verbose)
  test_raw <- import_data(test_file, verbose)
  case_mismatch(colnames(train_raw), colnames(test_raw))

  variable_info <- read_csv(info_file, show_col_types = FALSE, progress = FALSE)
  train_set <- process_training_data(train_raw, variable_info, verbose = verbose)
  test_set <- process_test_data(test_raw, train_set, verbose = verbose)

  predict_new_rf_exps(train_set, test_set, min_unique, min_good, max_cat_levels, verbose = verbose)
}



#' Merge pH7.4 and pH2 files.
#'
#' @param desc_file Main file with all data and pH7.4 values
#' @param ph2_file File with pH2 values.
#' @param moka_file File with Moka designations (optional)
#' @param verbose Logical. If TRUE, progress information will be outputted to the console.
#'
#' @return A tibble with merged values.
#' @export
merge_ph <- function(desc_file, ph2_file, moka_file = NULL, verbose = TRUE) {
  desc <- read_excel_ph(desc_file, verbose)
  ph2 <- read_excel_ph(ph2_file, verbose)
  if (!is.null(moka_file)) moka <- read_excel_ph(moka_file, verbose)

  ph2_sel <- ph2 %>%
    select("Name", where(is.numeric))
  ph_cols <- colnames(ph2_sel)[2:ncol(ph2_sel)]
  if (verbose) cat(paste("\npH columns found:", paste(ph_cols, collapse = ", "), "\n"))

  desc_cols <- desc %>%
    select(where(is.numeric)) %>%
    colnames()

  mtch <- ph_cols %in% desc_cols
  if (!all(mtch)) {
    cat(paste("\nError: pH columns not found in the descriptor file:", paste(ph_cols[!mtch], collapse = ","), "\n"))
    stop()
  }

  desc_sel <- desc %>%
    select("Name", all_of(ph_cols))
  desc_rest <- desc %>%
    select(-all_of(ph_cols))

  ph2_names <- glue::glue("G+_pH2_{ph_cols}")
  desc_names <- glue::glue("G+_pH7.4_{ph_cols}")

  if (verbose) cat(paste("\nCreating the following columns:\n  ", paste(ph2_names, collapse = ", "), "\n  ", paste(desc_names, collapse = ", "), "\nin the output file.\n"))

  colnames(ph2_sel) <- c("Name", ph2_names)
  colnames(desc_sel) <- c("Name", desc_names)

  res <- desc_rest %>%
    left_join(desc_sel, by = "Name") %>%
    left_join(ph2_sel, by = "Name")
  if (!is.null(moka_file)) {
    res <- res %>%
      left_join(moka, by = "Name")
  }

  res
}



#' Rename variables in a test file to match Lombardo training data.
#'
#' @param in_file Input file
#' @param rename_file Renaming file with two columns: lombardo_variable and current_variable.
#'
#' @return A tibble with data with renamed columns. Also, values in Moka_status/moka_ionState7.4 are converted to lowercase, as in Lombardo file.
#' @export
rename_for_vd <- function(in_file, rename_file, verbose = FALSE) {
  dat <- read_csv(in_file, show_col_types = FALSE, na = c("", "NA", "n/a", "N/A", "None"), guess_max = 10000)

  r <- read_csv(rename_file, show_col_types = FALSE)
  case_mismatch(r$current_variable, colnames(dat))
  r <- filter(r, current_variable %in% colnames(dat))

  if (verbose) {
    cat("\nRenaming the following columns:\n")
    print(as.data.frame(r))
  }

  ren <- set_names(r$current_variable, r$lombardo_variable)

  d <- dat %>%
    rename(all_of(ren))
  if (!is.null(d[["moka_ionState7.4"]])) d$moka_ionState7.4 <- tolower(d$moka_ionState7.4)
  d
}

