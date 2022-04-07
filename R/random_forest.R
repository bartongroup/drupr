# Run random forest on the set 'set', for the response variable 'resp_var'
random_forest_model <- function(set, resp_var, min_unique, min_good, max_cat_levels = 10, sel = NULL, seed = 666, ...) {
  set.seed(seed)
  d <- select_predictors_for_models(set, resp_var, min_unique, min_good, max_cat_levels, sel = sel) %>%
    drop_na()
  rows <- d$Name
  d <- d %>% select(-Name)
  # random forest doesn't like non-standard variable names, need to translate
  vars <- colnames(d)
  x_var <- tibble(
    var_original = vars,
    var_replaced = paste0("X", seq_along(vars))
  ) %>%
    mutate(var_replaced = if_else(var_original == "response", "response", var_replaced))
  colnames(d) <- x_var$var_replaced
  rf <- randomForest::randomForest(response ~ ., data = d, ...)
  list(
    rows = rows,
    rf = rf,
    x_var = x_var
  )
}


# extract importance from a list of rfs
rf_importance <- function(rf) {
  vars <- names(rf)
  map_dfr(vars, function(v) {
    rf[[v]]$rf %>%
      randomForest::importance() %>%
      as_tibble(rownames = "var_replaced") %>%
      left_join(rf[[v]]$x_var, by = "var_replaced") %>%
      select(-var_replaced) %>%
      select(
        variable = var_original,
        inc_mse = `%IncMSE`,
        inc_node_purity = IncNodePurity
      ) %>%
      add_column(response_variable = v, .before = "variable")
  }) %>%
    mutate(response_variable = as_factor(response_variable))
}

plot_rf_importance_sel <- function(rimp_sel, what = "inc_mse", top, mx = as.numeric(NA), title = NULL) {
  rimp_sel %>%
    drop_na() %>%
    rename(y = !!what) %>%
    arrange(y) %>%
    mutate(variable = factor(variable, levels = variable)) %>%
    tail(n = top) %>%
    ggplot(aes(x = variable, y = y, shape = class)) +
    geom_segment(aes(xend = variable, yend = 0), colour = "grey80") +
    geom_point() +
    coord_flip() +
    theme_bw() +
    theme(plot.title = element_text(size = 8), panel.grid = element_blank(), legend.position = "none") +
    labs(y = what, x = NULL, title = title) +
    scale_shape_manual(values = c(1, 16), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, mx))

}

plot_rf_importance <- function(rimp, what = "inc_mse", top = 20) {
  resps <- levels(rimp$response)
  mx <- max(rimp[[what]])
  map(resps, function(resp) {
    rimp %>%
      drop_na() %>%
      filter(response == resp) %>%
      plot_rf_importance_sel(what, top, mx, title = resp)
  }) %>%
    plot_grid(plotlist = ., ncol = 2, align = "v")
}

plot_rf_importance_heatmap <- function(rimp, what = "inc_mse") {
  tab <- rimp %>%
    pivot_wider(id_cols = response, names_from = variable, values_from = !!sym(what)) %>%
    column_to_rownames("response")
  ggheatmap(tab, second.labels.y = rep("", nrow(tab)), second.labels.x = rep("", ncol(tab))) +
    scale_fill_viridis_c(option = "cividis")
}

# Random forest does not like weird variable names, so when run, it renames them
# into X1, X2, ..., Xn. The conversion table is stored in $x_var field of the
# returned list. When this is done for the training set, the test set must be
# renamed in the same way. This function renames tibble 'dat', following
# conversion table 'x_var'. It returns a list with 'data' - converted table with
# no IDs and 'names' - a vector of IDs.
make_rf_xdat <- function(dat, x_var) {
  xvar <- x_var[-1, ] %>%
    filter(var_original %in% colnames(dat))
  sdat <- dat[, c("Name", xvar$var_original)] %>% drop_na()
  rows <- sdat$Name
  sdat <- sdat[, -1]
  colnames(sdat) <- xvar$var_replaced
  list(
    data = sdat,
    names = rows
  )
}

# Predict experimental variable from new data using random forest
predict_new_rf_exp <- function(set, resp_var, newdat, min_unique = 2, min_good = 1500, max_cat_levels = 10, seed = 123, n_tree = 1000) {

  # match good variables from the training with variables in the new set
  train <- select_predictors_for_models(set, resp_var, min_unique, min_good, max_cat_levels) %>%
    select(-c(Name, response))
  common_vars <- intersect(colnames(train), colnames(newdat))

  mdl <- random_forest_model(set, resp_var, min_unique, min_good, max_cat_levels, sel = common_vars, seed = seed, ntree = n_tree, importance = TRUE)

  # need to rename variables for RF
  sdat <- make_rf_xdat(newdat, mdl$x_var)

  # prediction errors from all trees
  rf_pred <- predict(mdl$rf, newdata = sdat$data, predict.all = TRUE)
  pred_aggregate <- tibble(prediction = rf_pred$aggregate)
  pred_individual <- predict(mdl$rf, newdata = sdat$data, predict.all = TRUE)$individual %>%
    apply(1, function(x) {
      c(
        quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE),
        mean(x, na.rm = TRUE) + c(0, -1, 1) * sd(x, na.rm = TRUE)
      )
    }) %>%
    t()
  colnames(pred_individual) <- c("pred_med", "pred_lo", "pred_up", "pred_mean", "pred_sd_lo", "pred_sd_up")
  pred_individual <- as_tibble(pred_individual)

  mdl$prediction <- bind_cols(pred_aggregate, pred_individual) %>%
    add_column(Name = sdat$names, .before = 1)

  mdl$new_data <- sdat$data

  return(mdl)
}

missing_values <- function(newdat) {
  newdat %>%
    mutate(across(-Name, is.na)) %>%
    pivot_longer(-Name, values_to = "missing") %>%
    filter(missing) %>%
    mutate(value = "+") %>%
    pivot_wider(Name, values_fill = "-")
}


# returns models and predictions with errors for a set of experimental variables
predict_new_rf_exps <- function(train_set, test_set, min_unique = 2, min_good = 1500, max_cat_levels = 10, seed = 123, verbose = FALSE) {
  if (verbose) cat("\nPredicting responses:\n")
  resp_vars <- train_set$variables %>% filter(response) %>% pull(variable)
  mdls <- map(resp_vars, function(resp_var) {
    if (verbose) cat(paste("    ", resp_var, " "))
    pr <- predict_new_rf_exp(train_set, resp_var, test_set$tab, min_unique, min_good, max_cat_levels, seed)
    if (verbose) cat(paste(nrow(pr$prediction), "compounds predicted\n"))
    pr
  }) %>%
    set_names(resp_vars)


  pr <- map_dfr(resp_vars, function(resp_var) {
    mdls[[resp_var]]$prediction %>% add_column(response_variable = resp_var, .before = 1)
  }) %>%
    mutate(response_variable = factor(response_variable, levels = resp_vars))

  imp <- rf_importance(mdls)

  list(
    test_data = test_set$tab,
    train_variables = train_set$variables,
    variable_comparison = test_set$variable_comparison,
    models = mdls,
    predictions = pr,
    importance = imp,
    mismatched_levels = test_set$mismatched_levels,
    missing = missing_values(test_set$tab)
  )
}




