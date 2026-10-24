BiocManager::install("GEOquery")
BiocManager::install("Biobase")
BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db",
                       "STRINGdb", "biomaRt"), update = FALSE, ask = FALSE)
install.packages(c("caret", "randomForest", "pROC", "pheatmap",
                   "ggplot2", "ggrepel", "VennDiagram", "dplyr"))
# Load libraries
library(GEOquery)
library(Biobase)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(caret)
library(randomForest)
library(pROC)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(dplyr)

set.seed(24)
#1.Data Acquisition
gse = getGEO("GSE45827", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset = gse[[1]]

expr = exprs(eset)        # probes x samples matrix (already log2-normalized)
pheno = pData(eset)       # sample metadata

 # Inspect what group labels exist (subtype / tissue type is usually in
# 'characteristics_ch1' or 'title' / 'source_name_ch1')
head(pheno$characteristics_ch1)
table(pheno$characteristics_ch1)

# Map probe IDs to gene symbols using the platform annotation (fData)
feat <- fData(eset)
head(feat[, c("ID", "Gene symbol")])

#2.Pre processing and quality control
#Boxplot of expression distributions across samples
boxplot(expr, outline = FALSE, las = 2, main = "Expression distribution per sample",
        ylab = "log2 expression", col = "lightblue")
# PCA to check sample clustering / detect outliers
pca = prcomp(t(expr), scale. = TRUE)
pca_df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     group = pheno$characteristics_ch1)

ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA of GSE45827 samples", x = "PC1", y = "PC2")

#3.Subtype-Specific Group Definition
pheno$group <- case_when(
  grepl("Basal", pheno$characteristics_ch1.1, ignore.case = TRUE) ~ "Basal_TNBC",
  grepl("N/A", pheno$characteristics_ch1.1, ignore.case = TRUE) ~ "Normal",
  TRUE ~ "Other_subtype"
)
keep = pheno$group %in% c("Basal_TNBC", "Normal")
expr_sub = expr[, keep]
group_sub = factor(pheno$group[keep], levels = c("Normal", "Basal_TNBC"))
table(group_sub)# confirm sample sizes per group
t = as.vector(table(group_sub))
labels=c("Normal","Basal_TNBC")
bar = barplot(t, names.arg = labels,col = "steelblue", ylim = c(0,50))
text(x=bar, y=t,label = c(11,41), pos = 3,cex = 1.5)
#4.Differential expression analysis (limma)
design = model.matrix(~ group_sub)
fit = lmFit(expr_sub, design)
fit = eBayes(fit)

deg_results = topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
deg_results$Gene = feat[rownames(deg_results), "Gene symbol"]

# Define significance thresholds
sig_thresh_p = 0.05
sig_thresh_fc = 1  # log2FC, i.e., 2-fold change

deg_results$Significant = with(deg_results,
                                ifelse(adj.P.Val < sig_thresh_p & abs(logFC) > sig_thresh_fc,
                                       ifelse(logFC > 0, "Up", "Down"), "NS"))
table(deg_results$Significant)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  geom_vline(xintercept = c(-sig_thresh_fc, sig_thresh_fc), linetype = "dashed") +
  geom_hline(yintercept = -log10(sig_thresh_p), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: Basal-like/TNBC vs Normal (GSE45827)",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)")

#Heatmap of top DEGs
top_deg_ids <- rownames(deg_results)[deg_results$Significant != "NS"][1:50]
pheatmap(expr_sub[top_deg_ids, ],
         annotation_col = data.frame(Group = group_sub, row.names = colnames(expr_sub)),
         show_rownames = FALSE, scale = "row",
         main = "Top 50 DEGs: Basal-like/TNBC vs Normal")

#5.Biomarker Prioritization + Literature Cross-Check
deg_results_unique = deg_results %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Gene) %>%
  slice_min(adj.P.Val, n = 1) %>%
  ungroup()
top_candidates = deg_results_unique %>%
  filter(Significant != "NS") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20) # last line important in later stage, to be commented in current analysis

print(top_candidates[, c("Gene", "logFC", "adj.P.Val", "Significant")])

# Known breast cancer markers to check for overlap (positive control / sanity check)
known_markers = c("ESR1", "ERBB2", "MKI67", "GATA3", "BRCA1", "BRCA2", "PGR")
overlap_known = intersect(top_candidates$Gene, known_markers)
cat("Known markers recovered in top candidates:", overlap_known, "\n")

#6.Functional Enrichment Analysis (GO/KEGG)
sig_genes = unique(deg_results$Gene[deg_results$Significant != "NS"])
sig_genes = sig_genes[!is.na(sig_genes) & sig_genes != ""]

entrez_ids = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

go_enrich = enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, readable = TRUE)

kegg_enrich = enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa",
                          pvalueCutoff = 0.05)
dotplot(go_enrich, showCategory = 15, title = "GO Biological Process Enrichment")
dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")

#7.Classifier + Cross-Validated ROC/AUC evaluation
panel_genes = top_candidates$Gene[1:15]              # top 15-gene panel
panel_ids = rownames(deg_results)[deg_results$Gene %in% panel_genes]

model_data = as.data.frame(t(expr_sub[panel_ids, ]))
colnames(model_data) = make.names(colnames(model_data))
model_data$Class = group_sub

ctrl = trainControl(method = "cv", number = 5, classProbs = TRUE,
                     summaryFunction = twoClassSummary, savePredictions = "final")

rf_model = train(Class ~ ., data = model_data, method = "rf",
                  trControl = ctrl, metric = "ROC")

print(rf_model)
roc_obj = roc(rf_model$pred$obs, rf_model$pred$Basal_TNBC)
plot(roc_obj, main = paste0("ROC Curve — 15-Gene Panel (AUC = ",
                            round(auc(roc_obj), 3), ")"))
#5 fold cv
n_folds = 5
n_top_genes = 20   # number of DEGs to re-select in each outer fold

outer_folds = createFolds(group_sub, k = n_folds, list = TRUE)

# Storage for results
outer_results = data.frame(
  Fold = integer(),
  ROC = numeric(),
  Sens = numeric(),
  Spec = numeric(),
  Accuracy = numeric()
)

selected_genes_per_fold = list()   # to check how stable gene selection is
all_predictions = data.frame()     # to build one overall ROC curve

for (i in seq_along(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")
  
  test_idx  = outer_folds[[i]]
  train_idx = setdiff(seq_along(group_sub), test_idx)
  
  expr_train = expr_sub[, train_idx]
  expr_test  = expr_sub[, test_idx]
  group_train = group_sub[train_idx]
  group_test  = group_sub[test_idx]
  
  # ---- STEP 1: Gene selection using ONLY the training fold ----
  # This is the key fix: limma DEG selection happens fresh each fold,
  # never seeing the test samples.
  design_train = model.matrix(~ group_train)
  fit_train = lmFit(expr_train, design_train)
  fit_train = eBayes(fit_train)
  
  deg_train = topTable(fit_train, coef = 2, number = Inf, adjust.method = "BH")
  deg_train$Gene = feat[rownames(deg_train), "Gene symbol"]
  
  # Collapse to one probe per gene (lowest p-value) to avoid duplicate-probe
  # inflation, same fix as before
  deg_train_unique = deg_train %>%
    filter(!is.na(Gene), Gene != "") %>%
    group_by(Gene) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(adj.P.Val)
  
  fold_genes = deg_train_unique$Gene[1:n_top_genes]
  fold_probe_ids = deg_train_unique$ID[1:n_top_genes]  # adjust column name if needed
  # If deg_train_unique doesn't retain probe ID as "ID", use rownames instead:
  fold_probe_ids = rownames(expr_train)[match(fold_genes,
                                               feat[rownames(expr_train), "Gene symbol"])]
  
  selected_genes_per_fold[[i]] = fold_genes
  
  # ---- STEP 2: Build training/test sets using ONLY the selected genes ----
  train_data = as.data.frame(t(expr_train[fold_probe_ids, ]))
  colnames(train_data) = make.names(fold_genes)
  train_data$Class = group_train
  
  test_data = as.data.frame(t(expr_test[fold_probe_ids, ]))
  colnames(test_data) = make.names(fold_genes)
  test_data$Class = group_test
  
  # ---- STEP 3: Train classifier on training fold (with inner CV for tuning) ----
  inner_ctrl = trainControl(method = "cv", number = 5,
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary)
  
  rf_fold = train(Class ~ ., data = train_data, method = "rf",
                   trControl = inner_ctrl, metric = "ROC")
  
  # ---- STEP 4: Evaluate on the held-out OUTER test fold ----
  pred_class = predict(rf_fold, newdata = test_data)
  pred_prob  = predict(rf_fold, newdata = test_data, type = "prob")
  
  cm = confusionMatrix(pred_class, test_data$Class)
  roc_fold = roc(response = test_data$Class,
                  predictor = pred_prob[, levels(group_sub)[2]])
  
  outer_results = rbind(outer_results, data.frame(
    Fold = i,
    ROC = as.numeric(auc(roc_fold)),
    Sens = cm$byClass["Sensitivity"],
    Spec = cm$byClass["Specificity"],
    Accuracy = cm$overall["Accuracy"]
  ))
  
  all_predictions = rbind(all_predictions, data.frame(
    obs = test_data$Class,
    prob = pred_prob[, levels(group_sub)[2]]
  ))
}
print(outer_results)

#8.Independent validation on a 2nd dataset
gse2 = getGEO("GSE38959", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset2 = gse2[[1]]
expr2 = exprs(eset2)
pheno2 = pData(eset2)
feat2 = fData(eset2)

# Check available characteristics
grep("characteristics", colnames(pheno2), value = TRUE)
table(pheno2$characteristics_ch1)
head(pheno2$title, 40)   

pheno2$group = case_when(
  grepl("triple negative|TNBC|basal", pheno2$characteristics_ch1, ignore.case = TRUE) ~ "Basal_TNBC",
  grepl("normal", pheno2$characteristics_ch1, ignore.case = TRUE) ~ "Normal",
  TRUE ~ "Other"
)

table(pheno2$group)   # should show ~30 Basal_TNBC and ~13 Normal, no "Other"
group2 = factor(pheno2$group, levels = c("Normal", "Basal_TNBC"))
design2 = model.matrix(~ group2)
fit2 = eBayes(lmFit(expr2, design2))
deg2 = topTable(fit2, coef = 2, number = Inf, adjust.method = "BH")
deg2$Gene = feat2[rownames(deg2), "Gene symbol"]

total_genes_universe <- length(unique(intersect(
  na.omit(feat$`Gene symbol`),      # genes annotated in GSE45827's platform
  na.omit(feat2$`Gene symbol`)      # genes annotated in GSE38959's platform
)))
cat("Total shared gene universe (background for hypergeometric test):",
    total_genes_universe, "\n")
panel_genes      <- unique(top_candidates$Gene)                     # your curated ~20 genes
validated_genes  <- unique(deg2$Gene[deg2$adj.P.Val < 0.05 &        # significant in GSE38959
                                       abs(deg2$logFC) > 1])

venn_panel <- venn.diagram(
  x = list(
    "GSE45827_TopCandidates" = panel_genes,
    "GSE38959_SignificantDEGs" = validated_genes
  ),
  filename = NULL,
  fill = c("skyblue", "salmon"),
  main = "Overlap: Top Candidate Panel (GSE45827) vs Significant DEGs (GSE38959)",
  cat.pos = c(-20, 20),
  cex = 1.3, cat.cex = 0.8
)
grid::grid.newpage()
grid::grid.draw(venn_panel)

replicated_panel_genes <- intersect(panel_genes, validated_genes)
cat("\nCurated candidate genes that replicated in GSE38959 (",
    length(replicated_panel_genes), "of", length(panel_genes), "):\n")
print(replicated_panel_genes)
overlap_panel <- length(replicated_panel_genes)
p_panel <- phyper(
  overlap_panel - 1,
  length(panel_genes),
  total_genes_universe - length(panel_genes),
  length(validated_genes),
  lower.tail = FALSE
)
cat("Hypergeometric test p-value (panel overlap):", format(p_panel, scientific = TRUE), "\n")

sig_genes_ds1 <- unique(deg_results$Gene[deg_results$Significant != "NS"])   # GSE45827
sig_genes_ds2 <- unique(deg2$Gene[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 1]) # GSE38959

venn_full <- venn.diagram(
  x = list("GSE45827" = sig_genes_ds1, "GSE38959" = sig_genes_ds2),
  filename = NULL,
  fill = c("lightgreen", "salmon"),
  main = "Overlap of All Significant DEGs (GSE45827 vs GSE38959)",
  cex = 1.2, cat.cex = 1.1
)
grid::grid.newpage()
grid::grid.draw(venn_full)

overlap_full <- length(intersect(sig_genes_ds1, sig_genes_ds2))
cat("\nFull DEG overlap:", overlap_full, "genes\n")

p_full <- phyper(
  overlap_full - 1,
  length(sig_genes_ds1),
  total_genes_universe - length(sig_genes_ds1),
  length(sig_genes_ds2),
  lower.tail = FALSE
)
cat("Hypergeometric test p-value (full DEG overlap):", format(p_full, scientific = TRUE), "\n")

#b
common_genes = intersect(panel_genes, deg2$Gene)
val_ids = rownames(deg2)[deg2$Gene %in% common_genes]
val_data = as.data.frame(t(expr2[val_ids, ]))
colnames(val_data) <- make.names(feat2[val_ids, "Gene symbol"])
val_data$Class <- group2
val_data = val_data[, intersect(colnames(model_data), colnames(val_data))]
preds = predict(rf_model, newdata = val_data)

panel_ids_common = intersect(panel_ids, rownames(expr2))

# Validation matrix
val_data <- as.data.frame(t(expr2[panel_ids_common, ]))

# Make column names valid
colnames(val_data) <- make.names(colnames(val_data))

# Ensure the columns are in the same order as the model
train_vars <- rf_model$finalModel$xNames
val_data <- val_data[, train_vars]

deg_unique <- deg_results %>%
  filter(!is.na(Gene), Gene != "") %>%
  arrange(adj.P.Val) %>%
  distinct(Gene, .keep_all = TRUE)
panel <- deg_unique[1:15, ]

panel_genes  <- panel$Gene
panel_probes <- rownames(panel)

train_data <- as.data.frame(t(expr_sub[panel_probes, ]))
colnames(train_data) <- make.names(panel_genes, unique = TRUE)

# Add class labels
train_data$Class <- group_sub
library(caret)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_model <- train(
  Class ~ .,
  data = train_common,
  method = "rf",
  trControl = ctrl_smote,
  metric = "ROC"
)
expr2
feat2
# Find which panel genes exist in GSE38959
common_genes <- intersect(panel_genes, feat2[, "Gene symbol"])

# Map genes back to probes in GSE38959
probe_map2 <- data.frame(
  Probe = rownames(feat2),
  Gene = feat2[, "Gene symbol"],
  stringsAsFactors = FALSE
)

probe_map2 <- probe_map2 %>%
  filter(Gene %in% common_genes) %>%
  distinct(Gene, .keep_all = TRUE)

val_data <- as.data.frame(t(expr2[probe_map2$Probe, ]))
colnames(val_data) <- make.names(probe_map2$Gene, unique = TRUE)

# Reorder columns to match training
val_data <- val_data[, colnames(train_data)[-ncol(train_data)]]
train_common <- train_data[, c(common_genes, "Class")]
val_common <- val_data[, common_genes]

preds <- predict(rf_model, newdata = val_common)
confusionMatrix(preds, group2)

library(randomForest)

importance(rf_model$finalModel)

#Form better model
# Install devtools if you haven't already
install.packages("devtools")

# Install DMwR2 directly from GitHub
devtools::install_github("ltorgo/DMwR2")

# Load the library
library(DMwR2)
install.packages("themis")
library(themis)

################################################################################
# CLASS IMBALANCE CORRECTION — TRAINING DATA ONLY (GSE45827)
#
# PRINCIPLE: Resampling (upsampling / downsampling / SMOTE / class-weighting)
# is applied ONLY to training folds, via caret's trainControl(sampling = ...),
# which automatically excludes the held-out test fold from resampling in
# every CV iteration. This avoids leakage (synthetic/duplicated samples
# ending up in both train and test) and keeps evaluation honest.
#
# Requires: model_data (data frame of selected gene panel + Class column,
# built from GSE45827 as in the earlier classifier script)
################################################################################

library(caret)
library(randomForest)
library(pROC)

set.seed(42)

# Confirm the imbalance you're correcting for
table(model_data$Class)
# e.g., Normal: 11, Basal_TNBC: 41  -> ~3.7:1 imbalance

################################################################################
# 1. BASELINE — no correction (for comparison)
################################################################################
ctrl_original <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE,
                              summaryFunction = twoClassSummary,
                              savePredictions = "final")

rf_original <- train(Class ~ ., data = model_data, method = "rf",
                     trControl = ctrl_original, metric = "ROC")

################################################################################
# 2. UPSAMPLING — random oversampling of minority class (fold-safe)
################################################################################
ctrl_up <- trainControl(method = "cv", number = 5,
                        classProbs = TRUE,
                        summaryFunction = twoClassSummary,
                        savePredictions = "final",
                        sampling = "up")

rf_up <- train(Class ~ ., data = model_data, method = "rf",
               trControl = ctrl_up, metric = "ROC")

################################################################################
# 3. DOWNSAMPLING — random undersampling of majority class (fold-safe)
# CAUTION: with only ~11 minority samples, downsampling shrinks each
# training fold drastically (~22 total samples). Included for comparison,
# but likely to underperform given how little data remains.
################################################################################
ctrl_down <- trainControl(method = "cv", number = 5,
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = "final",
                          sampling = "down")

rf_down <- train(Class ~ ., data = model_data, method = "rf",
                 trControl = ctrl_down, metric = "ROC")

################################################################################
# 4. SMOTE — synthetic minority oversampling (fold-safe)
# Apply only on the reduced gene panel (already the case here, not the full
# transcriptome), since SMOTE's nearest-neighbor interpolation is unreliable
# in very high dimensions.
################################################################################
# install.packages("themis")  # if 'smote' sampling option requires it
ctrl_smote <- trainControl(method = "cv", number = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           savePredictions = "final",
                           sampling = "smote")

rf_smote <- train(Class ~ ., data = model_data, method = "rf",
                  trControl = ctrl_smote, metric = "ROC")

################################################################################
# 5. CLASS WEIGHTING — alternative to resampling (no synthetic/duplicated data)
################################################################################
class_counts <- table(model_data$Class)
class_weights <- c(
  Normal     = as.numeric(class_counts["Basal_TNBC"]) / sum(class_counts),
  Basal_TNBC = as.numeric(class_counts["Normal"])      / sum(class_counts)
)
# Inverse-frequency weights: minority class gets a higher weight

rf_weighted <- randomForest(Class ~ ., data = model_data,
                            classwt = class_weights, importance = TRUE)

# Cross-validated evaluation of the weighted model (manual CV since randomForest()
# doesn't plug directly into caret's classwt tuning)
folds <- createFolds(model_data$Class, k = 5)
weighted_cv_auc <- sapply(folds, function(idx) {
  train_fold <- model_data[-idx, ]
  test_fold  <- model_data[idx, ]
  m <- randomForest(Class ~ ., data = train_fold, classwt = class_weights)
  pred_prob <- predict(m, test_fold, type = "prob")[, "Basal_TNBC"]
  as.numeric(auc(roc(test_fold$Class, pred_prob)))
})
cat("Class-weighted RF — CV AUC per fold:", round(weighted_cv_auc, 3),
    "\nMean AUC:", round(mean(weighted_cv_auc), 3), "\n")

################################################################################
# 6. COMPARE ALL METHODS
################################################################################
comparison <- resamples(list(
  Original    = rf_original,
  Upsampled   = rf_up,
  Downsampled = rf_down,
  SMOTE       = rf_smote
))

summary(comparison)
bwplot(comparison, metric = "ROC",
       main = "Cross-Validated ROC by Class Imbalance Correction Method")
bwplot(comparison, metric = "Sens",
       main = "Cross-Validated Sensitivity by Correction Method")
bwplot(comparison, metric = "Spec",
       main = "Cross-Validated Specificity by Correction Method")

# Print a clean summary table including the weighted model
final_comparison <- data.frame(
  Method = c("Original", "Upsampled", "Downsampled", "SMOTE", "Class-Weighted"),
  Mean_ROC = c(
    mean(rf_original$resample$ROC),
    mean(rf_up$resample$ROC),
    mean(rf_down$resample$ROC),
    mean(rf_smote$resample$ROC),
    mean(weighted_cv_auc)
  )
)
print(final_comparison)

################################################################################
# 7. SELECT FINAL MODEL FOR INDEPENDENT VALIDATION
################################################################################
# Choose whichever method gives the best BALANCED performance (not just
# highest ROC) — check Sens and Spec jointly, not ROC alone, since ROC can
# look fine even when one class is poorly predicted.
#
# Example: if SMOTE gives the best balance, retrain on the FULL GSE45827
# panel data (not just CV folds) for use in independent validation:

final_model <- train(Class ~ ., data = model_data, method = "rf",
                     trControl = trainControl(method = "cv", number = 5,
                                              classProbs = TRUE,
                                              summaryFunction = twoClassSummary,
                                              sampling = "smote"),  # <- set to chosen method
                     metric = "ROC")

saveRDS(final_model, "final_model_for_validation.rds")
cat("\nFinal model saved. Use this object in the independent validation script.\n")

library(caret)
library(pROC)
library(dplyr)
gene_symbols_train <- feat[panel_ids, "Gene symbol"]

model_data <- as.data.frame(t(expr_sub[panel_ids, ]))
colnames(model_data) <- make.names(gene_symbols_train)   # <- gene symbols, not probe IDs
model_data$Class <- group_sub

anyDuplicated(colnames(model_data))
# See which genes are duplicated
dup_genes <- colnames(model_data)[duplicated(colnames(model_data)) |
                                    duplicated(colnames(model_data), fromLast = TRUE)]
dup_genes <- setdiff(unique(dup_genes), "Class")
print(dup_genes)
# Add probe ID as an explicit column BEFORE any dplyr operation
deg_results$ProbeID <- rownames(deg_results)

best_probes <- deg_results %>%
  filter(!is.na(Gene), Gene %in% gene_symbols_train) %>%
  group_by(Gene) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup()

# Now pull probe IDs from the explicit column, not from rownames()
panel_ids_clean <- best_probes$ProbeID[match(gene_symbols_train, best_probes$Gene)]
panel_ids_clean <- panel_ids_clean[!is.na(panel_ids_clean)]

# Sanity check BEFORE subsetting expr_sub — this would have caught the bug early
all(panel_ids_clean %in% rownames(expr_sub))   # should be TRUE
length(panel_ids_clean)                         # should roughly match your panel size

gene_symbols_clean <- feat[panel_ids_clean, "Gene symbol"]

model_data <- as.data.frame(t(expr_sub[panel_ids_clean, ]))
colnames(model_data) <- make.names(gene_symbols_clean)
model_data$Class <- group_sub

anyDuplicated(colnames(model_data))   # confirm 0  # confirm this is now 0
# Check for duplicate gene symbols mapping to multiple probes (common issue) —
# if duplicates exist, keep only the probe with the highest variance or
# lowest p-value per gene before this step
# ---- Reload the trained model if starting a fresh session ----
final_model <- readRDS("final_model_for_validation.rds")

################################################################################
# 1. CONFIRM NATURAL (UNCORRECTED) CLASS DISTRIBUTION OF VALIDATION SET
################################################################################
table(group2)
# Expected: Normal = 13, Basal_TNBC = 30 -- left as-is, NOT resampled

################################################################################
# 2. BUILD THE VALIDATION FEATURE MATRIX USING THE SAME GENE PANEL
################################################################################
# Match genes used in training; some may be missing/renamed in GSE38959's
# annotation, so intersect explicitly and report any dropped genes.
common_genes <- intersect(panel_genes, feat2[, "Gene symbol"])
missing_genes <- setdiff(panel_genes, common_genes)
if (length(missing_genes) > 0) {
  cat("Warning: the following panel genes were not found in GSE38959 and will",
      "be excluded from validation:\n", paste(missing_genes, collapse = ", "), "\n")
}

val_probe_ids <- rownames(feat2)[feat2[, "Gene symbol"] %in% common_genes]

val_data <- as.data.frame(t(expr2[val_probe_ids, ]))
colnames(val_data) <- make.names(feat2[val_probe_ids, "Gene symbol"])
val_data$Class <- group2

# Align column order/names with what the model was trained on
train_features <- setdiff(colnames(final_model$trainingData), ".outcome")
missing_in_val <- setdiff(train_features, colnames(val_data))
if (length(missing_in_val) > 0) {
  cat("Note: model expects features not present in validation data:\n",
      paste(missing_in_val, collapse = ", "),
      "\n-> Consider retraining the model using only genes common to both datasets.\n")
}
val_data <- val_data[, intersect(train_features, colnames(val_data)), drop = FALSE]

################################################################################
# 3. PREDICT ON THE INDEPENDENT (NATURAL-DISTRIBUTION) DATASET
################################################################################
pred_class <- predict(final_model, newdata = val_data)
pred_prob  <- predict(final_model, newdata = val_data, type = "prob")

################################################################################
# 4. PERFORMANCE METRICS — REPORT IMBALANCE-ROBUST METRICS, NOT JUST ACCURACY
################################################################################
# Raw accuracy is misleading on a 30:13 split (a model predicting "TNBC" for
# everyone would score ~70% accuracy without learning anything). Report
# sensitivity, specificity, balanced accuracy, and AUC together.

cm <- confusionMatrix(pred_class, group2, positive = "Basal_TNBC")
print(cm)

roc_val <- roc(response = group2, predictor = pred_prob[, "Basal_TNBC"])
auc_val <- auc(roc_val)

balanced_accuracy <- (cm$byClass["Sensitivity"] + cm$byClass["Specificity"]) / 2

performance_summary <- data.frame(
  Metric = c("Accuracy", "Sensitivity (TNBC detection)", "Specificity (Normal detection)",
             "Balanced Accuracy", "AUC"),
  Value = c(
    cm$overall["Accuracy"],
    cm$byClass["Sensitivity"],
    cm$byClass["Specificity"],
    balanced_accuracy,
    as.numeric(auc_val)
  )
)
print(performance_summary)

################################################################################
# 5. VISUALIZE
################################################################################
plot(roc_val, main = paste0("Independent Validation ROC — GSE38959 (AUC = ",
                            round(auc_val, 3), ")"))

# Confusion matrix heatmap
cm_table <- as.data.frame(cm$table)
library(ggplot2)
ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix — Independent Validation (GSE38959, natural distribution)")

################################################################################
# 6. COMPARE: DISCOVERY (NESTED CV) vs. INDEPENDENT VALIDATION PERFORMANCE
################################################################################
# Pull your nested CV summary numbers from the earlier script for direct
# side-by-side comparison in your report.
comparison_table <- data.frame(
  Dataset = c("GSE45827 (Nested CV)", "GSE38959 (Independent Validation)"),
  Sensitivity = c(NA, cm$byClass["Sensitivity"]),   # fill NA with your nested CV value
  Specificity = c(NA, cm$byClass["Specificity"]),   # fill NA with your nested CV value
  AUC = c(NA, as.numeric(auc_val))                  # fill NA with your nested CV value
)
print(comparison_table)