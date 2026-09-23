make_df_site_mc <- function(D, batch) {
  df <- tibble::as_tibble(D)
  df$site <- factor(batch)
  df
}

auc_ovr_macro <- function(truth, prob_mat) {
  lv <- levels(truth)
  stopifnot(all(lv %in% colnames(prob_mat)))
  prob_df <- as.data.frame(prob_mat)
  aucs <- vapply(lv, function(cl) {
    y <- ifelse(truth == cl, 1, 0)
    pred <- as.numeric(prob_df[[cl]])
    as.numeric(
      pROC::roc(
        response = y,
        predictor = pred,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )$auc
    )
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}

mc_sens_spec_macro <- function(truth, pred, lv = levels(truth)) {
  truth <- factor(truth, levels = lv)
  pred <- factor(pred, levels = lv)

  acc <- mean(pred == truth)
  sens_vec <- vapply(lv, function(cl) {
    TP <- sum(truth == cl & pred == cl)
    FN <- sum(truth == cl & pred != cl)
    if ((TP + FN) == 0) return(NA_real_)
    TP / (TP + FN)
  }, numeric(1))

  spec_vec <- vapply(lv, function(cl) {
    TN <- sum(truth != cl & pred != cl)
    FP <- sum(truth != cl & pred == cl)
    if ((TN + FP) == 0) return(NA_real_)
    TN / (TN + FP)
  }, numeric(1))

  list(
    Accuracy = acc,
    Sensitivity_Macro = mean(sens_vec, na.rm = TRUE),
    Specificity_Macro = mean(spec_vec, na.rm = TRUE)
  )
}

make_shared_site_splits <- function(batch, n_splits, base_seed) {
  batch <- factor(batch)
  sites <- levels(batch)
  set.seed(base_seed)
  purrr::map(seq_len(n_splits), function(s) {
    train_idx <- unlist(lapply(sites, function(site_level) {
      idx_site <- which(batch == site_level)
      sample(idx_site, size = floor(0.5 * length(idx_site)))
    }))
    list(train = train_idx, test = setdiff(seq_along(batch), train_idx))
  })
}

evaluate_methods_site_multiclass <- function(data_list, batch, p, r,
                                             n_splits = 20,
                                             trees = 500,
                                             mtry = NULL,
                                             base_seed = 7000 + r) {
  if (is.null(mtry)) mtry <- floor(sqrt(p))
  methods <- names(data_list)
  splits <- make_shared_site_splits(batch, n_splits, base_seed)

  rf_spec <- parsnip::rand_forest(mtry = mtry, trees = trees) |>
    parsnip::set_engine("ranger", probability = TRUE) |>
    parsnip::set_mode("classification")

  purrr::map_dfr(methods, function(m) {
    dfm <- make_df_site_mc(data_list[[m]], batch)
    rec <- recipes::recipe(site ~ ., data = dfm) |>
      recipes::step_zv(recipes::all_predictors()) |>
      recipes::step_normalize(recipes::all_predictors())

    wf <- workflows::workflow() |>
      workflows::add_model(rf_spec) |>
      workflows::add_recipe(rec)

    purrr::map_dfr(seq_along(splits), function(k) {
      tr <- splits[[k]]$train
      te <- splits[[k]]$test
      fit <- parsnip::fit(wf, data = dfm[tr, , drop = FALSE])
      pr <- stats::predict(fit, dfm[te, , drop = FALSE], type = "prob")

      lv <- levels(dfm$site)
      pred_cols <- paste0(".pred_", lv)
      prob_mat <- pr[, pred_cols, drop = FALSE]
      colnames(prob_mat) <- lv

      truth <- factor(dfm$site[te], levels = lv)
      auc_macro <- auc_ovr_macro(truth, prob_mat)
      pred_class <- factor(lv[max.col(as.matrix(prob_mat), ties.method = "first")], levels = lv)
      mets <- mc_sens_spec_macro(truth, pred_class, lv = lv)

      tibble::tibble(
        Rep = r,
        Split = k,
        Method = m,
        Task = "Site_MC",
        AUC = auc_macro,
        Accuracy = mets$Accuracy,
        Sensitivity_Macro = mets$Sensitivity_Macro,
        Specificity_Macro = mets$Specificity_Macro
      )
    })
  })
}

make_df_sex <- function(D, batch, sex) {
  df <- tibble::as_tibble(D)
  df$site <- factor(batch)
  df$sex <- factor(sex, levels = c("F", "M"))
  df
}

evaluate_methods_loso_sex <- function(data_list, batch, sex, p, r,
                                      trees = 100, mtry = NULL) {
  if (is.null(mtry)) mtry <- floor(sqrt(p))
  methods <- names(data_list)
  sites <- levels(factor(batch))

  rf_spec <- parsnip::rand_forest(mtry = mtry, trees = trees) |>
    parsnip::set_engine("ranger", probability = TRUE) |>
    parsnip::set_mode("classification")

  purrr::map_dfr(methods, function(m) {
    df_m <- make_df_sex(data_list[[m]], batch, sex)
    rec <- recipes::recipe(sex ~ ., data = df_m) |>
      recipes::step_rm(site) |>
      recipes::step_zv(recipes::all_predictors()) |>
      recipes::step_normalize(recipes::all_predictors())

    wf <- workflows::workflow() |>
      workflows::add_model(rf_spec) |>
      workflows::add_recipe(rec)

    purrr::map_dfr(sites, function(hold) {
      tr <- which(df_m$site != hold)
      te <- which(df_m$site == hold)
      fit <- parsnip::fit(wf, data = df_m[tr, , drop = FALSE])

      pr_df <- stats::predict(fit, df_m[te, , drop = FALSE], type = "prob")
      cls_df <- stats::predict(fit, df_m[te, , drop = FALSE], type = "class")

      pr_M <- as.numeric(pr_df$.pred_M)
      cls <- factor(cls_df$.pred_class, levels = c("F", "M"))
      truth <- factor(df_m$sex[te], levels = c("F", "M"))

      auc <- as.numeric(
        pROC::roc(
          response = truth,
          predictor = pr_M,
          levels = c("F", "M"),
          direction = "<",
          quiet = TRUE
        )$auc
      )

      tibble::tibble(
        Rep = r,
        Holdout = hold,
        Method = m,
        Task = "Sex_LOSO",
        AUC = auc,
        Accuracy = yardstick::accuracy_vec(truth, cls),
        Sensitivity = yardstick::sens_vec(truth, cls, event_level = "second"),
        Specificity = yardstick::spec_vec(truth, cls, event_level = "second")
      )
    })
  })
}
