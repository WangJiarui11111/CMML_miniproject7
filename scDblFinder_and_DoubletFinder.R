# Data path
input_rds <- "C:/Users/wangj/Desktop/hm-6k.rds"
scrublet_csv <- "C:/Users/wangj/Desktop/scrublet_result.csv"
outdir <- "C:/Users/wangj/Desktop/hm6k_project_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Set seed
set.seed(123)

# 1. Packages
options(Seurat.object.assay.version = "v3") # Use legacy Assay to ensure DoubletFinder compatibility under Seurat v5
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(patchwork)
library(pROC)
library(PRROC)
library(scDblFinder)
library(SingleCellExperiment)
library(DoubletFinder)
library(forcats)

# 2. Plot theme
theme_benchmark <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 1
      ),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.position = "right"
    )
}

method_levels <- c("scDblFinder", "DoubletFinder", "Scrublet")
method_cols <- c(
  "scDblFinder" = "#66C2A5",
  "DoubletFinder" = "#FC8D62",
  "Scrublet" = "#8DA0CB"
)

truth_cols <- c(
  "singlet" = "grey80",
  "doublet" = "#D73027"
)

# 3. Import hm-6k data
hm <- readRDS(input_rds)

hm_counts <- hm[[1]]
hm_label <- hm[[2]]

hm_label <- as.character(hm_label)
hm_label <- tolower(hm_label)

if (!all(hm_label %in% c("singlet", "doublet"))) {
  stop("hm_label must contain only 'singlet' and 'doublet'. Please check labels.")
}

if (!is.null(names(hm_label))) {
  if (all(colnames(hm_counts) %in% names(hm_label))) {
    hm_label <- hm_label[colnames(hm_counts)]
  }
}

if (length(hm_label) != ncol(hm_counts)) {
  stop("Length of hm_label does not match number of cells in hm_counts.")
}

hm_meta <- data.frame(
  truth = hm_label,
  is_doublet_truth = hm_label == "doublet"
)

rownames(hm_meta) <- colnames(hm_counts)

hm_df <- CreateSeuratObject(
  counts = hm_counts,
  meta.data = hm_meta,
  project = "hm6k"
)

hm_df$truth <- tolower(hm_df$truth)
hm_df$is_doublet_truth <- hm_df$truth == "doublet"

cat("Initial cells:", ncol(hm_df), "\n")
cat("Initial genes:", nrow(hm_df), "\n")
cat("Initial true doublet rate:", mean(hm_df$is_doublet_truth), "\n")

# 4. QC and preprocessing
hm_df[["percent.mt"]] <- PercentageFeatureSet(hm_df, pattern = "^MT-")

p_qc_before <- VlnPlot(
  hm_df,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) +
  plot_annotation(title = "QC metrics before filtering")

ggsave(
  filename = file.path(outdir, "Supplementary_QC_before_filtering.png"),
  plot = p_qc_before,
  width = 10,
  height = 4,
  dpi = 300
)

hm_df <- subset(
  hm_df,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7500 &
    percent.mt < 10
)

cat("Cells after QC:", ncol(hm_df), "\n")
cat("Genes after QC:", nrow(hm_df), "\n")
cat("True doublet rate after QC:", mean(hm_df$is_doublet_truth), "\n")

true_dbr <- mean(hm_df$is_doublet_truth)

hm_df <- NormalizeData(hm_df, verbose = FALSE)
hm_df <- FindVariableFeatures(
  hm_df,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
hm_df <- ScaleData(hm_df, verbose = FALSE)
hm_df <- RunPCA(hm_df, verbose = FALSE)
hm_df <- FindNeighbors(hm_df, dims = 1:20, verbose = FALSE)
hm_df <- FindClusters(hm_df, resolution = 0.5, verbose = FALSE)
hm_df <- RunUMAP(hm_df, dims = 1:20, verbose = FALSE)

hm_df$cluster_ref <- hm_df$seurat_clusters

p_umap_truth_initial <- DimPlot(
  hm_df,
  group.by = "truth",
  reduction = "umap",
  cols = truth_cols,
  pt.size = 0.4
) +
  ggtitle("Ground truth after QC") +
  theme_benchmark()

ggsave(
  filename = file.path(outdir, "UMAP_ground_truth_after_QC.png"),
  plot = p_umap_truth_initial,
  width = 5,
  height = 4,
  dpi = 300
)

# 5. scDblFinder
sce <- as.SingleCellExperiment(hm_df)

sce <- scDblFinder(
  sce,
  dbr = true_dbr
)

hm_df$scDblFinder_score <- colData(sce)$scDblFinder.score
hm_df$scDblFinder_class <- as.character(colData(sce)$scDblFinder.class)
hm_df$scDblFinder_class <- tolower(hm_df$scDblFinder_class)
hm_df$scDblFinder_pred <- hm_df$scDblFinder_class == "doublet"

cat("\nscDblFinder confusion table:\n")
print(table(hm_df$truth, hm_df$scDblFinder_class))

# 6. DoubletFinder 
sweep.res.list <- paramSweep(
  hm_df,
  PCs = 1:20,
  sct = FALSE
)

sweep.stats <- summarizeSweep(
  sweep.res.list,
  GT = FALSE
)

bcmvn <- find.pK(sweep.stats)

bcmvn$pK_numeric <- as.numeric(as.character(bcmvn$pK))

best_pK <- bcmvn$pK_numeric[which.max(bcmvn$BCmetric)]
best_pK <- as.numeric(best_pK)

cat("\nBest pK:", best_pK, "\n")

p_pk <- ggplot(
  bcmvn,
  aes(x = pK_numeric, y = BCmetric)
) +
  geom_line(color = "grey40") +
  geom_point(size = 2, color = "#FC8D62") +
  geom_vline(
    xintercept = best_pK,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "pK",
    y = "BCmetric",
    title = "DoubletFinder pK selection"
  ) +
  theme_benchmark()

ggsave(
  filename = file.path(outdir, "Supplementary_DoubletFinder_pK_selection.png"),
  plot = p_pk,
  width = 5,
  height = 4,
  dpi = 300
)

nExp_poi <- round(true_dbr * ncol(hm_df))
nExp_poi <- as.numeric(nExp_poi)

cat("Expected doublets for main DoubletFinder:", nExp_poi, "\n")

old_df_cols <- grep(
  "^pANN|^DF.classifications",
  colnames(hm_df@meta.data),
  value = TRUE
)

if (length(old_df_cols) > 0) {
  hm_df@meta.data[, old_df_cols] <- NULL
}

hm_df <- doubletFinder(
  hm_df,
  PCs = 1:20,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp_poi,
  reuse.pANN = FALSE,
  sct = FALSE
)

df_class_col <- grep(
  "^DF.classifications",
  colnames(hm_df@meta.data),
  value = TRUE
)

df_score_col <- grep(
  "^pANN",
  colnames(hm_df@meta.data),
  value = TRUE
)

df_class_col <- df_class_col[length(df_class_col)]
df_score_col <- df_score_col[length(df_score_col)]

hm_df$DoubletFinder_class <- as.character(hm_df@meta.data[[df_class_col]])
hm_df$DoubletFinder_score <- hm_df@meta.data[[df_score_col]]
hm_df$DoubletFinder_pred <- hm_df$DoubletFinder_class == "Doublet"

hm_df$DoubletFinder_class_plot <- ifelse(
  hm_df$DoubletFinder_pred,
  "doublet",
  "singlet"
)

cat("\nDoubletFinder confusion table:\n")
print(table(hm_df$truth, hm_df$DoubletFinder_class_plot))

# 7. Scrublet import

# Scrublet must be run on cells x genes matrix in Python.
# Seurat counts are genes x cells, so the matrix should be transposed in Python.
if (!file.exists(scrublet_csv)) {
  stop(
    "Scrublet result CSV not found. Please run Scrublet first and save scrublet_result.csv."
  )
}

scrublet_result <- read.csv(scrublet_csv)

required_scrublet_cols <- c("barcode", "Scrublet_score", "Scrublet_pred")

if (!all(required_scrublet_cols %in% colnames(scrublet_result))) {
  stop(
    "scrublet_result.csv must contain columns: barcode, Scrublet_score, Scrublet_pred"
  )
}

rownames(scrublet_result) <- scrublet_result$barcode

scrublet_result <- scrublet_result[colnames(hm_df), ]

if (!all(rownames(scrublet_result) == colnames(hm_df))) {
  stop("Scrublet barcodes are not aligned with Seurat object barcodes.")
}

hm_df$Scrublet_score <- scrublet_result$Scrublet_score
hm_df$Scrublet_pred <- as.logical(scrublet_result$Scrublet_pred)
hm_df$Scrublet_class <- ifelse(hm_df$Scrublet_pred, "doublet", "singlet")

cat("\nScrublet confusion table:\n")
print(table(hm_df$truth, hm_df$Scrublet_class))

# 8. Metric function with AUROC and AUPRC
calculate_metrics <- function(truth, pred, score = NULL, method = "method") {
  
  truth <- as.logical(truth)
  pred <- as.logical(pred)
  
  TP <- sum(truth == TRUE & pred == TRUE, na.rm = TRUE)
  FN <- sum(truth == TRUE & pred == FALSE, na.rm = TRUE)
  TN <- sum(truth == FALSE & pred == FALSE, na.rm = TRUE)
  FP <- sum(truth == FALSE & pred == TRUE, na.rm = TRUE)
  
  safe_divide <- function(x, y) {
    ifelse(y == 0, NA, x / y)
  }
  
  sensitivity <- safe_divide(TP, TP + FN)
  specificity <- safe_divide(TN, TN + FP)
  precision <- safe_divide(TP, TP + FP)
  
  f1 <- ifelse(
    is.na(precision) |
      is.na(sensitivity) |
      (precision + sensitivity) == 0,
    NA,
    2 * precision * sensitivity / (precision + sensitivity)
  )
  
  accuracy <- safe_divide(TP + TN, TP + TN + FP + FN)
  false_positive_rate <- safe_divide(FP, FP + TN)
  predicted_doublet_rate <- mean(pred, na.rm = TRUE)
  true_doublet_rate <- mean(truth, na.rm = TRUE)
  
  AUROC <- NA
  AUPRC <- NA
  
  if (!is.null(score)) {
    
    keep <- !is.na(truth) & !is.na(score)
    truth_use <- truth[keep]
    score_use <- score[keep]
    
    roc_obj <- tryCatch(
      roc(
        response = truth_use,
        predictor = score_use,
        levels = c(FALSE, TRUE),
        direction = "<",
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(roc_obj)) {
      AUROC <- as.numeric(auc(roc_obj))
    }
    
    pr_obj <- tryCatch(
      PRROC::pr.curve(
        scores.class0 = score_use[truth_use == TRUE],
        scores.class1 = score_use[truth_use == FALSE],
        curve = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(pr_obj)) {
      AUPRC <- pr_obj$auc.integral
    }
  }
  
  data.frame(
    method = method,
    TP = TP,
    FN = FN,
    TN = TN,
    FP = FP,
    sensitivity = sensitivity,
    recall = sensitivity,
    specificity = specificity,
    precision = precision,
    F1 = f1,
    accuracy = accuracy,
    FPR = false_positive_rate,
    predicted_doublet_rate = predicted_doublet_rate,
    true_doublet_rate = true_doublet_rate,
    AUROC = AUROC,
    AUPRC = AUPRC
  )
}

# 9. Doublet detection benchmark
metrics_scDblFinder <- calculate_metrics(
  truth = hm_df$is_doublet_truth,
  pred = hm_df$scDblFinder_pred,
  score = hm_df$scDblFinder_score,
  method = "scDblFinder"
)

metrics_DoubletFinder <- calculate_metrics(
  truth = hm_df$is_doublet_truth,
  pred = hm_df$DoubletFinder_pred,
  score = hm_df$DoubletFinder_score,
  method = "DoubletFinder"
)

metrics_Scrublet <- calculate_metrics(
  truth = hm_df$is_doublet_truth,
  pred = hm_df$Scrublet_pred,
  score = hm_df$Scrublet_score,
  method = "Scrublet"
)

benchmark_table <- bind_rows(
  metrics_scDblFinder,
  metrics_DoubletFinder,
  metrics_Scrublet
) %>%
  mutate(
    method = factor(method, levels = method_levels)
  )

benchmark_table_round <- benchmark_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(benchmark_table_round)

write.csv(
  benchmark_table,
  file.path(outdir, "benchmark_detection_metrics_full.csv"),
  row.names = FALSE
)

write.csv(
  benchmark_table_round,
  file.path(outdir, "benchmark_detection_metrics_rounded.csv"),
  row.names = FALSE
)

# 10. Detection benchmark plots

# 10.1 AUROC / AUPRC plot
auc_long <- benchmark_table %>%
  select(method, AUROC, AUPRC) %>%
  pivot_longer(
    cols = c(AUROC, AUPRC),
    names_to = "metric",
    values_to = "value"
  )

p_auc_auprc <- ggplot(
  auc_long,
  aes(x = method, y = value, fill = method)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(value, 3)),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~ metric, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = "Score",
    title = "AUROC and AUPRC"
  ) +
  theme_benchmark() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Figure_AUROC_AUPRC_comparison.png"),
  plot = p_auc_auprc,
  width = 7,
  height = 4,
  dpi = 300
)

# 10.2 Precision, Recall, F1, Specificity
metric_long <- benchmark_table %>%
  select(method, precision, recall, F1, specificity) %>%
  pivot_longer(
    cols = c(precision, recall, F1, specificity),
    names_to = "metric",
    values_to = "value"
  )

p_detection_metrics <- ggplot(
  metric_long,
  aes(x = method, y = value, fill = method)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(value, 3)),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~ metric, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = "Metric value",
    title = "Classification metrics"
  ) +
  theme_benchmark() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Figure_detection_metrics_comparison.png"),
  plot = p_detection_metrics,
  width = 9,
  height = 4,
  dpi = 300
)

# 10.3 Confusion matrix stacked plot, supplementary

confusion_long <- benchmark_table %>%
  select(method, TP, FP, FN, TN) %>%
  pivot_longer(
    cols = c(TP, FP, FN, TN),
    names_to = "confusion_type",
    values_to = "count"
  ) %>%
  mutate(
    confusion_type = factor(
      confusion_type,
      levels = c("TP", "FP", "FN", "TN")
    )
  )

p_confusion_stack <- ggplot(
  confusion_long,
  aes(x = method, y = count, fill = confusion_type)
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "TP" = "#D73027",
      "FP" = "#FC8D59",
      "FN" = "#91BFDB",
      "TN" = "#4575B4"
    )
  ) +
  labs(
    x = NULL,
    y = "Cell count",
    fill = "Type",
    title = "Confusion matrix components"
  ) +
  theme_benchmark() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = file.path(outdir, "Supplementary_confusion_matrix_stacked.png"),
  plot = p_confusion_stack,
  width = 6,
  height = 4,
  dpi = 300
)

# 10.4 Score boxplot
score_df <- data.frame(
  truth = hm_df$truth,
  scDblFinder = hm_df$scDblFinder_score,
  DoubletFinder = hm_df$DoubletFinder_score,
  Scrublet = hm_df$Scrublet_score
)

score_long <- score_df %>%
  pivot_longer(
    cols = c(scDblFinder, DoubletFinder, Scrublet),
    names_to = "method",
    values_to = "score"
  ) %>%
  mutate(
    method = factor(method, levels = method_levels)
  )

p_score_box <- ggplot(
  score_long,
  aes(x = truth, y = score, fill = truth)
) +
  geom_boxplot(
    outlier.size = 0.4,
    width = 0.65,
    alpha = 0.85
  ) +
  facet_wrap(~ method, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = truth_cols) +
  labs(
    x = "Ground truth",
    y = "Doublet score",
    title = "Doublet score by ground truth"
  ) +
  theme_benchmark() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Figure_doublet_score_boxplot.png"),
  plot = p_score_box,
  width = 8,
  height = 4,
  dpi = 300
)

# 10.5 UMAP truth vs prediction

p_truth <- DimPlot(
  hm_df,
  group.by = "truth",
  reduction = "umap",
  cols = truth_cols,
  pt.size = 0.35
) +
  ggtitle("Ground truth") +
  theme_benchmark()

p_scDblFinder <- DimPlot(
  hm_df,
  group.by = "scDblFinder_class",
  reduction = "umap",
  cols = truth_cols,
  pt.size = 0.35
) +
  ggtitle("scDblFinder") +
  theme_benchmark()

p_DoubletFinder <- DimPlot(
  hm_df,
  group.by = "DoubletFinder_class_plot",
  reduction = "umap",
  cols = truth_cols,
  pt.size = 0.35
) +
  ggtitle("DoubletFinder") +
  theme_benchmark()

p_Scrublet <- DimPlot(
  hm_df,
  group.by = "Scrublet_class",
  reduction = "umap",
  cols = truth_cols,
  pt.size = 0.35
) +
  ggtitle("Scrublet") +
  theme_benchmark()

fig_umap_truth_prediction <- 
  (p_truth | p_scDblFinder) /
  (p_DoubletFinder | p_Scrublet)

ggsave(
  filename = file.path(outdir, "Figure_UMAP_truth_vs_predictions.png"),
  plot = fig_umap_truth_prediction,
  width = 10,
  height = 8,
  dpi = 300
)

# 11. DoubletFinder sensitivity analysis
pANN_main_col <- df_score_col
homotypic.prop <- modelHomotypic(hm_df$seurat_clusters)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

cat("\nHomotypic proportion:", homotypic.prop, "\n")
cat("nExp main:", nExp_poi, "\n")
cat("nExp homotypic-adjusted:", nExp_poi.adj, "\n")

df_sensitivity_grid <- data.frame(
  setting = c(
    "0.5x_expected",
    "1.0x_expected_main",
    "1.5x_expected",
    "2.0x_expected",
    "homotypic_adjusted"
  ),
  nExp = c(
    round(0.5 * nExp_poi),
    nExp_poi,
    round(1.5 * nExp_poi),
    round(2.0 * nExp_poi),
    nExp_poi.adj
  )
)

df_sensitivity_grid$nExp <- as.numeric(df_sensitivity_grid$nExp)

df_sensitivity_metrics <- list()

df_sensitivity_metrics[["1.0x_expected_main"]] <- calculate_metrics(
  truth = hm_df$is_doublet_truth,
  pred = hm_df$DoubletFinder_pred,
  score = hm_df$DoubletFinder_score,
  method = "DoubletFinder_1.0x_expected_main"
) %>%
  mutate(
    setting = "1.0x_expected_main",
    nExp = nExp_poi,
    homotypic_prop = homotypic.prop
  )

for (i in seq_len(nrow(df_sensitivity_grid))) {
  
  setting_i <- df_sensitivity_grid$setting[i]
  nExp_i <- df_sensitivity_grid$nExp[i]
  
  if (setting_i == "1.0x_expected_main") {
    next
  }
  
  message("Running DoubletFinder sensitivity setting: ", setting_i)
  
  hm_df <- doubletFinder(
    hm_df,
    PCs = 1:20,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_i,
    reuse.pANN = pANN_main_col,
    sct = FALSE
  )
  
  df_class_cols_now <- grep(
    "^DF.classifications",
    colnames(hm_df@meta.data),
    value = TRUE
  )
  
  latest_class_col <- df_class_cols_now[length(df_class_cols_now)]
  
  clean_setting <- gsub("[^A-Za-z0-9]", "_", setting_i)
  
  class_col_new <- paste0("DoubletFinder_class_", clean_setting)
  pred_col_new <- paste0("DoubletFinder_pred_", clean_setting)
  
  hm_df[[class_col_new]] <- hm_df@meta.data[[latest_class_col]]
  hm_df[[pred_col_new]] <- hm_df@meta.data[[class_col_new]] == "Doublet"
  
  df_sensitivity_metrics[[setting_i]] <- calculate_metrics(
    truth = hm_df$is_doublet_truth,
    pred = hm_df@meta.data[[pred_col_new]],
    score = hm_df$DoubletFinder_score,
    method = paste0("DoubletFinder_", setting_i)
  ) %>%
    mutate(
      setting = setting_i,
      nExp = nExp_i,
      homotypic_prop = homotypic.prop
    )
}

df_sensitivity_table <- bind_rows(df_sensitivity_metrics) %>%
  mutate(
    expected_doublet_rate = nExp / ncol(hm_df),
    n_cells = ncol(hm_df)
  )

df_sensitivity_table_round <- df_sensitivity_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(df_sensitivity_table_round)

write.csv(
  df_sensitivity_table,
  file.path(outdir, "DoubletFinder_nExp_sensitivity_full.csv"),
  row.names = FALSE
)

write.csv(
  df_sensitivity_table_round,
  file.path(outdir, "DoubletFinder_nExp_sensitivity_rounded.csv"),
  row.names = FALSE
)

df_sensitivity_long <- df_sensitivity_table %>%
  select(setting, precision, recall, F1, specificity, predicted_doublet_rate) %>%
  pivot_longer(
    cols = -setting,
    names_to = "metric",
    values_to = "value"
  )

p_df_sensitivity <- ggplot(
  df_sensitivity_long,
  aes(x = setting, y = value, fill = setting)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(value, 3)),
    vjust = -0.3,
    size = 2.8
  ) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "DoubletFinder nExp setting",
    y = "Value",
    title = "DoubletFinder sensitivity to expected doublet number"
  ) +
  theme_benchmark() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Supplementary_DoubletFinder_nExp_sensitivity.png"),
  plot = p_df_sensitivity,
  width = 10,
  height = 6,
  dpi = 300
)

# 12. Downstream clustering after doublet removal
reprocess_seurat <- function(obj, dims_use = 1:20, resolution_use = 0.5) {
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution_use, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  
  return(obj)
}

hm_original_downstream <- hm_df

hm_scDblFinder_singlet <- subset(
  hm_df,
  subset = scDblFinder_pred == FALSE
)

hm_DoubletFinder_singlet <- subset(
  hm_df,
  subset = DoubletFinder_pred == FALSE
)

hm_Scrublet_singlet <- subset(
  hm_df,
  subset = Scrublet_pred == FALSE
)

set.seed(123)
hm_original_downstream <- reprocess_seurat(hm_original_downstream)

set.seed(123)
hm_scDblFinder_singlet <- reprocess_seurat(hm_scDblFinder_singlet)

set.seed(123)
hm_DoubletFinder_singlet <- reprocess_seurat(hm_DoubletFinder_singlet)

set.seed(123)
hm_Scrublet_singlet <- reprocess_seurat(hm_Scrublet_singlet)

# 13. Downstream UMAP clustering plots
p_cluster_original <- DimPlot(
  hm_original_downstream,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.35
) +
  ggtitle("Original") +
  theme_benchmark() +
  theme(legend.position = "none")

p_cluster_scDblFinder <- DimPlot(
  hm_scDblFinder_singlet,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.35
) +
  ggtitle("After scDblFinder removal") +
  theme_benchmark() +
  theme(legend.position = "none")

p_cluster_DoubletFinder <- DimPlot(
  hm_DoubletFinder_singlet,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.35
) +
  ggtitle("After DoubletFinder removal") +
  theme_benchmark() +
  theme(legend.position = "none")

p_cluster_Scrublet <- DimPlot(
  hm_Scrublet_singlet,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.35
) +
  ggtitle("After Scrublet removal") +
  theme_benchmark() +
  theme(legend.position = "none")

fig_downstream_umap_clusters <- 
  (p_cluster_original | p_cluster_scDblFinder) /
  (p_cluster_DoubletFinder | p_cluster_Scrublet) +
  plot_annotation(
    title = "UMAP clustering before and after doublet removal"
  )

ggsave(
  filename = file.path(outdir, "Figure_downstream_UMAP_clusters_after_removal.png"),
  plot = fig_downstream_umap_clusters,
  width = 10,
  height = 8,
  dpi = 300
)

# 14. Downstream clustering quality summary
calculate_cluster_quality <- function(obj, method_name, original_ncells) {
  
  meta <- obj@meta.data
  
  cluster_summary <- meta %>%
    group_by(seurat_clusters) %>%
    summarise(
      n_cells_cluster = n(),
      singlet_fraction = mean(is_doublet_truth == FALSE),
      doublet_fraction = mean(is_doublet_truth == TRUE),
      .groups = "drop"
    )
  
  data.frame(
    method = method_name,
    n_cells = ncol(obj),
    retained_cell_fraction = ncol(obj) / original_ncells,
    removed_cell_fraction = 1 - ncol(obj) / original_ncells,
    n_clusters = length(unique(obj$seurat_clusters)),
    mean_cluster_singlet_purity = weighted.mean(
      cluster_summary$singlet_fraction,
      cluster_summary$n_cells_cluster
    ),
    remaining_true_doublet_rate = mean(obj$is_doublet_truth)
  )
}

original_ncells <- ncol(hm_df)

cluster_quality_table <- bind_rows(
  calculate_cluster_quality(
    hm_original_downstream,
    "Original",
    original_ncells
  ),
  calculate_cluster_quality(
    hm_scDblFinder_singlet,
    "scDblFinder",
    original_ncells
  ),
  calculate_cluster_quality(
    hm_DoubletFinder_singlet,
    "DoubletFinder",
    original_ncells
  ),
  calculate_cluster_quality(
    hm_Scrublet_singlet,
    "Scrublet",
    original_ncells
  )
) %>%
  mutate(
    method = factor(
      method,
      levels = c("Original", "scDblFinder", "DoubletFinder", "Scrublet")
    )
  )

cluster_quality_table_round <- cluster_quality_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(cluster_quality_table_round)

write.csv(
  cluster_quality_table,
  file.path(outdir, "downstream_clustering_quality_full.csv"),
  row.names = FALSE
)

write.csv(
  cluster_quality_table_round,
  file.path(outdir, "downstream_clustering_quality_rounded.csv"),
  row.names = FALSE
)

cluster_quality_long <- cluster_quality_table %>%
  select(
    method,
    retained_cell_fraction,
    removed_cell_fraction,
    n_clusters,
    mean_cluster_singlet_purity,
    remaining_true_doublet_rate
  ) %>%
  pivot_longer(
    cols = -method,
    names_to = "metric",
    values_to = "value"
  )

downstream_cols <- c(
  "Original" = "grey60",
  "scDblFinder" = "#66C2A5",
  "DoubletFinder" = "#FC8D62",
  "Scrublet" = "#8DA0CB"
)

p_cluster_quality <- ggplot(
  cluster_quality_long,
  aes(x = method, y = value, fill = method)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = round(value, 3)),
    vjust = -0.3,
    size = 2.8
  ) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = downstream_cols) +
  labs(
    x = NULL,
    y = "Value",
    title = "Downstream clustering summary"
  ) +
  theme_benchmark() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "Figure_downstream_clustering_quality_summary.png"),
  plot = p_cluster_quality,
  width = 10,
  height = 6,
  dpi = 300
)

# 15. Save final object
saveRDS(
  hm_df,
  file.path(outdir, "hm6k_final_with_doublet_predictions.rds")
)

saveRDS(
  list(
    original = hm_original_downstream,
    scDblFinder_removed = hm_scDblFinder_singlet,
    DoubletFinder_removed = hm_DoubletFinder_singlet,
    Scrublet_removed = hm_Scrublet_singlet
  ),
  file.path(outdir, "hm6k_downstream_objects.rds")
)
