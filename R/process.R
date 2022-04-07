

# Main function to process raw data

process_training_data <- function(raw, info, cols_nolog = NULL,
                        missing_row_limit = 0.95, max_missing = 100, sigma_limit = 5,
                        cut_tree = 0.25, min_group_unique = 100,
                        verbose = FALSE) {

  if (verbose) cat("\nProcessing training data\n")
  if (verbose) cat(" ", nrow(raw), "compounds in training set\n")

  cols_id <- info %>% filter(type == "id") %>% pull(variable)
  cols_remove <- info %>% filter(type == "none") %>% pull(variable)

  dg <- raw

  # remove specified variables
  dg <- dg %>% select(-all_of(cols_remove))

  # properties of all variables
  props <- variable_properties(dg %>% select(-all_of(cols_id)), cols_nolog = cols_nolog) %>%
    left_join(info, by = c("original_variable" = "variable")) %>%
    mutate(response = type == "response") %>%
    mutate(response = replace_na(response, FALSE)) %>%
    select(-type) %>%
    mutate(predictor = !response)

  # find variables that lead to no more than missing_row_limit missing rows
  # and with no more than max_missing missing values
  # response variables are kept anyway
  min_good <- floor(nrow(raw) * missing_row_limit)
  mis_stats <- missing_stats(raw, cols_id)
  cols_mis <- mis_stats %>%
    filter(n_row_good < min_good & nbad > max_missing) %>%
    pull(variable)
  props <- props %>%
    mutate(mis = original_variable %in% cols_mis)

  # find useless "null" variables - all missing or only one unique value
  cols_null <- props %>% filter(null) %>% pull(original_variable)
  dg <- dg %>% select(-all_of(cols_null))

  cols_real <- props %>% filter(!null & real) %>% pull(original_variable)
  cols_int <- props %>% filter(!null & integer) %>% pull(original_variable)
  cols_numcat <- props %>% filter(!null & numcat) %>% pull(original_variable)
  cols_cat <- props %>% filter(!null & categorical) %>% pull(original_variable)

  # make tables for ids, numerical and categorical variables
  d_id <- dg[, cols_id]
  d_real <- dg[, cols_real]
  d_int <- dg[, cols_int]
  d_cat <- dg[, cols_cat] %>%
    mutate_all(convert_yes_no) %>%
    bind_cols(select(dg, all_of(cols_numcat))) %>%
    mutate_all(as_factor)

  # Convert some numerical variables into logarithms; change names
  d_numlog <- convert_log(d_real, props) %>%
    remove_outliers(sigma_limit) %>%
    bind_cols(d_int)

  var2log <- set_names(props$variable, props$original_variable)

  # detailed summaries of numerical and categorical variables
  num_props <- numeric_variable_summary(d_numlog)
  cat_props <- categorical_variable_summary(d_cat)

  # separate response and predictor variables
  cols_resp_num <- props %>% filter(!null & response & numeric) %>% pull(variable)
  cols_resp_cat <- props %>% filter(!null & response & categorical) %>% pull(variable)
  cols_pred_num <- props %>% filter(!null & predictor & numeric) %>% pull(variable)
  cols_pred_cat <- props %>% filter(!null & predictor & categorical) %>% pull(variable)
  d_resp <- bind_cols(d_cat[, cols_resp_cat], d_numlog[, cols_resp_num])
  d_pred <- bind_cols(d_cat[, cols_pred_cat], d_numlog[, cols_pred_num])

  # full table of data: id columns, response columns and predictor columns
  tab <- bind_cols(d_id, d_resp, d_pred)

  # all variable info
  vars <- props %>%
    left_join(bind_rows(cat_props, num_props) %>% select(variable, levels, top_count, zeroes), by = "variable")

  gr <- group_variables(tab, vars, cut_tree, min_group_unique)
  final_vars <- gr$variables %>% filter(variable %in% colnames(tab))
  n_sel <- final_vars %>% filter(selected & !mis & !null) %>% nrow()

  if (verbose) {
    cat("  Variables found:\n")
    cat(sprintf("    %3d not usable\n", length(cols_null)))
    cat(sprintf("    %3d real\n", length(cols_real)))
    cat(sprintf("    %3d of them converted to logarithm\n", sum(props$log_trans)))
    cat(sprintf("    %3d integer\n", length(cols_int)))
    cat(sprintf("    %3d numeric converted into categorical\n", length(cols_numcat)))
    cat(sprintf("    %3d categorical\n", length(cols_cat)))
    cat("  Grouping of similar variables perfomed:\n")
    cat(sprintf("    %3d variables selected for downstream analysis\n", n_sel))
  }

  list(
    all_variables = vars,
    variables = final_vars,
    cols_id = cols_id,
    missing_stats = mis_stats,
    hc = gr$hc,
    hc_groups = gr$hc_groups,
    tab = tab
  )
}



process_test_data <- function(raw, train, verbose = FALSE) {
  if (verbose) cat("\nProcessing test data\n")
  if (verbose) cat(" ", nrow(raw), "compounds in test set\n")

  dg <- raw %>% select(-all_of(train$cols_id))

  descriptor_vars <- train$variables %>%
    filter(predictor & selected & !mis & !null)

  test_vars <- tibble(original_variable = colnames(dg))

  # properties of all variables
  props <- descriptor_vars %>%
    filter(original_variable %in% test_vars$original_variable)

  varcomp <- select(descriptor_vars, original_variable) %>%
    add_column(in_train = TRUE) %>%
    full_join(test_vars %>% add_column(in_test = TRUE), by = "original_variable") %>%
    mutate(
      in_train = replace_na(in_train, FALSE),
      in_test = replace_na(in_test, FALSE)
    )
  not_in_test <- varcomp %>% filter(in_train & !in_test) %>% pull(original_variable)
  not_in_train <- varcomp %>% filter(!in_train & in_test) %>% pull(original_variable)

  cols_real <- props %>% filter(!null & real) %>% pull(original_variable)
  cols_int <- props %>% filter(!null & integer) %>% pull(original_variable)
  cols_numcat <- props %>% filter(!null & numcat) %>% pull(original_variable)
  cols_cat <- props %>% filter(!null & categorical) %>% pull(original_variable)

  # make tables for ids, numerical and categorical variables
  d_id <- raw[, train$cols_id]
  d_real <- dg[, cols_real]
  d_int <- dg[, cols_int]
  d_cat <- dg[, cols_cat]

  if (ncol(d_cat) > 0) {
    d_cat <- d_cat %>%
      mutate_all(convert_yes_no) %>%
      bind_cols(select(dg, all_of(cols_numcat))) %>%
      mutate_all(as_factor)
  }

  # Convert some numerical variables into logarithms; change names
  d_numlog <- convert_log(d_real, props) %>%
    bind_cols(d_int)

  var2log <- set_names(props$variable, props$original_variable)

  # detailed summaries of numerical and categorical variables
  num_props <- numeric_variable_summary(d_numlog)
  cat_props <- categorical_variable_summary(d_cat)

  # level mismatch
  cat_match <- cat_props %>%
    left_join(props, by = "variable")
  if (!is.null(cat_match$levels.x)) {
    mismatch <- cat_match %>%
      select(variable, test_levels = levels.x, train_levels = levels.y) %>%
      levels_mismatch() %>%
      filter(!matched) %>%
      select(-matched)
    cat_props <- cat_props %>% filter(!(variable %in% mismatch$variable))
    d_cat <- d_cat %>% select(-all_of(mismatch$variable))
  } else {
    mismatch <- tibble(variable = character(0))
  }

  # all variable info
  vars <- props %>%
    select(-c(levels, top_count, zeroes)) %>%
    left_join(bind_rows(cat_props, num_props) %>% select(variable, levels, top_count, zeroes), by = "variable")


  # full table of data: id columns, response columns and predictor columns
  tab <- bind_cols(d_id, d_numlog, d_cat)

  if (verbose) {
    cat("  Variables found:\n")
    cat(sprintf("    %3d variables in test set\n", nrow(test_vars)))
    cat(sprintf("    %3d selected variables in the training set\n",nrow(descriptor_vars)))
    cat(sprintf("    %3d test variables match training set\n", nrow(vars)))
    cat(sprintf("    %3d training variables not found in test set\n", length(not_in_test)))
    cat(sprintf("    %3d categorical variables with mismatched levels\n", nrow(mismatch)))
    if (nrow(mismatch) > 0) cat(paste("        ", paste(mismatch$variable, collapse = ", "), "\n"))
  }

  list(
    variables = vars,
    variable_comparison = varcomp,
    cols_id = train$cols_id,
    not_in_test = not_in_test,
    not_in_train = not_in_train,
    mismatched_levels = mismatch,
    tab = tab
  )
}

