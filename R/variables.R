is_integer <- function(v) {
  if (!is.numeric(v)) return(FALSE)
  vn <- na.omit(v)
  sum(vn != as.integer(vn)) == 0
}

variable_properties <- function(d, cols_nolog = NULL) {
  props <- d %>%
    map_df(function(v) {
      tibble(
        class = class(v),
        good = sum(!is.na(v)),
        missing = sum(is.na(v)),
        unique = length(unique(na.omit(v))),
        integer = is_integer(v) & unique > 2,    # "integer" does not include 0/1 variables
      )
    }) %>%
    mutate(
      null = (good == 0 | unique == 1),  # useless "null" variables - all missing or only one unique value
      numeric = (class == "numeric" & unique > 2),          # proper numerical variables
      real = (class == "numeric" & unique > 2 & !integer),  # non-integer values with more than 2 levels
      numcat = (class == "numeric" & unique == 2),          # numerical variables with 2 levels are categorical
      categorical = (class == "character" & !null)          # proper categorical variables
    ) %>%
    add_column(original_variable = names(d), .before = "class") %>%
    mutate(class = as_factor(class))

  cols_real <- props %>% filter(real) %>% pull(original_variable)
  log_tab <- detect_log_distribution(d[, cols_real], cols_nolog)

  props %>%
    left_join(log_tab, by = "original_variable") %>%
    mutate(
      variable = if_else(is.na(variable), original_variable, variable),
      log_trans = replace_na(log_trans, FALSE)
    ) %>%
    relocate(variable)
}


numeric_variable_summary <- function(tab) {
  tab %>% map_dfr(function(v) {
    vn <- na.omit(v)
    tibble(
      min = min(vn),
      median = median(vn),
      mean = mean(vn),
      max = max(vn),
      good = sum(!is.na(v)),
      missing = sum(is.na(v)),
      zeroes = sum(v == 0, na.rm = TRUE),
      unique = length(unique(v))
    )
  }) %>%
    add_column(variable = names(tab), .before = "min")
}

categorical_variable_summary <- function(tab) {
  if (ncol(tab) > 0) {
    map_df(tab, function(v) {
      vn <- na.omit(v)
      tibble(
        good = sum(!is.na(v)),
        missing = sum(is.na(v)),
        unique = length(unique(vn)),
        levels = paste(levels(vn), collapse = ","),
        top_count = max(fct_count(vn)$n)
      )
    }) %>% add_column(variable = names(tab), .before = "good")
  } else {
    tibble(variable = character(0), good = integer(0), missing = integer(0), unique = integer(0), levels = character(0), top_count = integer(0))
  }
}

missing_stats <- function(d, id_cols, n_top = 30) {
  d <- d %>% select(-all_of(id_cols))
  mst <- d %>%
    map_df(function(v) tibble(ngood = sum(!is.na(v)), nbad = sum(is.na(v)))) %>%
    add_column(variable = names(d), .before = "ngood") %>%
    arrange(ngood) %>%
    head(n_top)
  nrowgood <- function(rownumber) {
    select(d, -c(mst$variable[1:rownumber])) %>%
      na.omit %>%
      nrow
  }
  nt <- NULL
  for (i in 1:nrow(mst)) {
    nt <- c(nt, nrowgood(i))
  }
  mst$n_row_good <- nt
  mst
}


convert_yes_no <- function(x) {
  map_chr(x, function(s) {
    if (is.na(s)) return(NA)
    if (str_detect(s, "^(Yes|No|High|Low|SupSat|Non_SupSat)")) {
      s <- str_remove(s, "\\s+\\(\\d+%\\)$")
    } else if (s == "NoSites") {
      s <- NA
    }
    s
  })
}

detect_log_distribution <- function(tab, cols_nolog = NULL, min_vgood = 30, max_vbad = 10) {
  colnames(tab) %>%
    map_dfr(function(nm) {
      v <- tab %>% pull(nm)
      ret <- FALSE
      nbad <- sum(v <= 0, na.rm = TRUE)
      vgood <- v[!is.na(v) & v > 0]
      if (length(vgood) >= min_vgood & nbad <= max_vbad & !(nm %in% cols_nolog)) {
        log.norm <- capture.output({
          fit.norm <- tryCatch(fitdistrplus::fitdist(vgood, "norm"), error = function(err) list(aic = 1e16))
        })
        log.lnorm <- capture.output({
          fit.lnorm <- tryCatch(fitdistrplus::fitdist(vgood, "lnorm"), error = function(err) list(aic = 1e16))
        })
        if (fit.norm$aic > fit.lnorm$aic) ret <- TRUE
      }
      tibble(original_variable = nm, log_trans = ret)
    }) %>%
    mutate(variable = if_else(log_trans, paste("log", original_variable), original_variable), .before = "log_trans")
}

convert_log <- function(tab, props) {
  logvars <- props %>% filter(log_trans) %>% pull(original_variable)
  for (v in logvars) {
    x <- tab[[v]]
    x[x <= 0] <- NA
    lg <- log10(x)
    tab[, v] <- lg
  }

  newnames <- tibble(original_variable = colnames(tab)) %>%
    left_join(props, by = "original_variable") %>%
    pull(variable)

  colnames(tab) <- newnames
  tab
}

remove_outliers <- function(tab, sigma_limit = 5) {
  map_dfc(names(tab), function(nm) {
    v <- tab[[nm]]
    outl <- which(abs(scale(v)) > sigma_limit)
    if (length(outl) > 0) v[outl] <- NA
    tibble(v) %>% set_names(nm)
  })
}


# try: find the smallest pka per compound
convert_pka <- function(v) {
  v %>%
    as.character() %>%
    str_split(";") %>%
    map_dbl(function(x) {
      if (is.na(x[[1]])) return(NA)
      x %>%
        as.numeric() %>%
        min()
    })
}

variable_correlation <- function(tab) {
  ct <- cor(tab, use = "pairwise.complete.obs") %>%
    as.data.frame()
  names <- rownames(ct)
  n <- length(names)
  colnames(ct) <- 1:n

  ct %>%
    as_tibble() %>%
    add_column(id = 1:n) %>%
    pivot_longer(-id) %>%
    set_names(c("i1", "i2", "correlation")) %>%
    mutate(i2 = as.integer(i2)) %>%
    filter(i1 != i2) %>%
    mutate(rep = i1 < i2) %>%
    mutate(var1 = names[i1], var2 = names[i2]) %>%
    select(var1, var2, correlation, rep)
}



# when NAs are dropped, some of the variables might end up having
# only one level, we need to reject them before modelling
reduce_for_levels <- function(tab, min_unique) {
  tab <- drop_na(tab)
  vars <- names(tab)
  bad_vars <- vars %>%
    map_dfr(function(vr) tibble(variable = vr, n = length(unique(tab[[vr]])))) %>%
    filter(n < min_unique) %>%
    pull(variable)
  tab %>% select(-all_of(bad_vars))
}

# For a given data set 'set', response variable 'resp_var', create a subset of the main set, with variables with at least 'min_unique' unique values, at least 'min_good' good (non-missing) values, at most 'max_cat_levels' levels in categorical variables. The subset is done of selected (cluster centroids) variables with no 'mis' flag.
select_predictors_for_models <- function(set, resp_var, min_unique, min_good, max_cat_levels, sel = NULL) {
  pred <- set$variables %>%
    filter(predictor & !null & !mis & selected & good >= min_good & (class == "numeric" | unique <= max_cat_levels))
  if (!is.null(sel)) pred <- filter(pred, variable %in% sel)
  pred_vars <- pull(pred, variable)
  reduce_for_levels(set$tab[, c("Name", resp_var, pred_vars)], min_unique) %>%
    rename(response = !!sym(resp_var)) %>%
    filter(!is.na(response))
}



test_min_good <- function(set, resp_var, min_unique = 2, max_cat_levels = 10, min_range = 10) {
  max_range <- nrow(set$tab)
  pb <- txtProgressBar(min = 1, max = max_range - min_range + 1, initial = 1, style = 3)
  tb <- map_dfr(min_range:max_range, function(mg) {
    x <- select_predictors_for_models(set, resp_var, min_unique = min_unique, min_good = mg, max_cat_levels = max_cat_levels)
    setTxtProgressBar(pb, mg - min_range + 1)
    c(min_good = mg, n_rows = nrow(x), n_cols = ncol(x) - 2)
  })
  close(pb)
  tb
}

# check if levels in test set are a subset of levels in the train set
levels_mismatch <- function(m) {
  m %>%
    mutate_at(vars(test_levels, train_levels), ~str_split(.x, ",")) %>%
    rowwise() %>%
    mutate(matched = all(test_levels %in% train_levels)) %>%
    mutate_at(vars(test_levels, train_levels), ~str_c(.x, collapse = ","))
}
