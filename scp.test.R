library(SCP)
options(reticulate.conda_binary = "/home/chenzh/miniconda3/condabin/conda",SCP_env_name = "R4.2")
data("pancreas_sub")
BPPARAM <- BiocParallel::MulticoreParam(workers = 10)


pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
pancreas_sub <- RunDynamicFeatures(srt = pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200,BPPARAM = BiocParallel::MulticoreParam(workers = 10))
pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
# ht <- DynamicHeatmap(
#   srt = pancreas_sub, lineages = c("Lineage1", "Lineage2"),
#   use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
#   species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
#   heatmap_palette = "viridis", cell_annotation = "SubCellType",
#   separate_annotation = list("SubCellType", c("Nnat", "Irx1")), separate_annotation_palette = c("Paired", "Set1"),
#   feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#   pseudotime_label = 25, pseudotime_label_color = "red",
#   height = 5, width = 2
# )
# print(ht$plot)
pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200,BPPARAM = BiocParallel::MulticoreParam(workers = 10))




ht1 <- DynamicHeatmap( srt = pancreas_sub,lineages = "Lineage1",n_split = 5,split_method = "kmeans-peaktime",cell_annotation = "SubCellType"
,row_names_side="right")
ht1$plot
panel_fix(ht1$plot, raster = TRUE, dpi = 50)

DynamicHeatmap( srt = pancreas_sub,lineages = "Lineage1",n_split = 5,split_method = "kmeans-peaktime",cell_annotation = "SubCellType",row_names_side="right",features_label=c("Birc5","Aurkb"))$plot


# default 
features_ordered <- ht1$feature_metadata %>% tbl_df() %>% pull(features)
cells_ordered <-  pancreas_sub@meta.data %>% tibble::rownames_to_column("cell") %>% arrange(Lineage1) %>% filter(!is.na(Lineage1)) %>% pull(cell)
fit_exp <- pancreas_sub@tools$DynamicFeatures_Lineage1$fitted_matrix %>% t()
fit_exp[features_ordered,cells_ordered] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,scale="none")

nlabel <- 20
index_from <- ceiling((length(features_ordered) / nlabel) / 2)
index_to <- length(features_ordered)
index <- unique(round(seq(from = index_from, to = index_to, length.out = nlabel)))




# ht2 <- DynamicHeatmap(
#   srt = pancreas_sub,
#   lineages = "Lineage1",
#   features = c("Sox9", "Neurod2", "Isl1", "Rbp4", "Pyy", "S_score", "G2M_score"),
#   cell_annotation = "SubCellType"
# )
# ht2$plot
# panel_fix(ht2$plot, height = 5, width = 5, raster = TRUE, dpi = 50)


ht3 <- DynamicHeatmap(srt = pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_split = 5, split_method = "kmeans", cluster_rows = TRUE, cell_annotation = "SubCellType",row_names_side="right")
ht3$plot

pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
ht4 <- DynamicHeatmap(
  srt = pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 6,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c("SubCellType", "Phase", "G2M_score"),
  cell_annotation_palette = c("Paired", "simspec", "Purples"),
  separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
  separate_annotation_palette = c("Paired", "Set1"),
  separate_annotation_params = list(height = unit(10, "mm")),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25,
  pseudotime_label_color = "red"
)
ht4$plot

ht5 <- DynamicHeatmap(
  srt = pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 6,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c("SubCellType", "Phase", "G2M_score"),
  cell_annotation_palette = c("Paired", "simspec", "Purples"),
  separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
  separate_annotation_palette = c("Paired", "Set1"),
  separate_annotation_params = list(width = unit(10, "mm")),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  flip = TRUE, column_title_rot = 45
)
ht5$plot

ht6 <- DynamicHeatmap(
  srt = pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  cell_annotation = "SubCellType",
  n_split = 5, split_method = "mfuzz",
  species = "Mus_musculus", db = "GO_BP",
  anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
)
ht6$plot



#' detailed explore
srt=pancreas_sub
lineages <- "Lineage1"
cell_bins <- 100
cell_density <- 1
cell_annotation <-  "SubCellType"
assay <- "RNA"
order_by = c("peaktime")
decreasing <- FALSE


cell_union <- unique(colnames(srt@assays[[1]])[apply(srt@meta.data[, lineages, drop = FALSE], 1, function(x) !all(is.na(x)))])
Pseudotime_assign <- rowMeans(srt@meta.data[cell_union, lineages, drop = FALSE], na.rm = TRUE)
cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union, cells = cell_union),
                                  Pseudotime_assign = Pseudotime_assign,
                                  srt@meta.data[cell_union, lineages, drop = FALSE])

bins <- cut(Pseudotime_assign, breaks = seq(min(Pseudotime_assign, na.rm = TRUE), max(Pseudotime_assign, na.rm = TRUE), length.out = cell_bins), include.lowest = TRUE)
ncells <- ceiling(max(table(bins), na.rm = TRUE) * cell_density)
message("ncell/bin=", ncells, "(", cell_bins, "bins)")

cell_keep <- unlist(sapply(levels(bins), function(x) {
  cells <- names(Pseudotime_assign)[bins == x]
  out <- sample(cells, size = min(length(cells), ncells))
  return(out)
}))
cell_metadata <- cell_metadata[cell_keep, , drop = FALSE]



cell_order_list <- list()
for (l in lineages) {
  # if (!is.null(reverse_ht) && (l %in% lineages[reverse_ht] || l %in% reverse_ht)) {
  cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
  cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], decreasing = TRUE), , drop = FALSE]
  # } else {
  #   cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
  #   cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], decreasing = FALSE), , drop = FALSE]
  # }
  cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), l)
}


cell_metadata <- cbind.data.frame(
  cell_metadata,
  cbind.data.frame(
    srt@meta.data[rownames(cell_metadata), c(intersect(cell_annotation, colnames(srt@meta.data))), drop = FALSE]
  )
)


DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
dynamic <- list()
DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > 20 & DynamicFeatures$r.sq > 0.2 & DynamicFeatures$dev.expl > 0.2 & DynamicFeatures$padjust < 0.05, , drop = FALSE]
dynamic[[l]] <- DynamicFeatures
features <- DynamicFeatures[["features"]]

all_calculated <- Reduce(intersect, lapply(lineages, function(l) rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]])))
features_tab <- table(features)
#features <- names(features_tab)[which(features_tab %in% (num_intersections %||% seq_along(lineages)))]


gene <- features[features %in% rownames(srt@assays[[assay]])]
meta <- features[features %in% colnames(srt@meta.data)]

feature_metadata <- data.frame(row.names = features, features = features)

for (l in lineages) {
  feature_metadata[rownames(dynamic[[l]]), paste0(l, order_by)] <- dynamic[[l]][, order_by]
  feature_metadata[rownames(dynamic[[l]]), paste0(l, "exp_ncells")] <- dynamic[[l]][, "exp_ncells"]
  feature_metadata[rownames(dynamic[[l]]), paste0(l, "r.sq")] <- dynamic[[l]][, "r.sq"]
  feature_metadata[rownames(dynamic[[l]]), paste0(l, "dev.expl")] <- dynamic[[l]][, "dev.expl"]
  feature_metadata[rownames(dynamic[[l]]), paste0(l, "padjust")] <- dynamic[[l]][, "padjust"]
}


feature_metadata[, order_by] <- apply(feature_metadata[, paste0(lineages, order_by), drop = FALSE], 1, max, na.rm = TRUE)
feature_metadata <- feature_metadata[order(feature_metadata[, order_by], decreasing = decreasing), , drop = FALSE]
feature_metadata <- feature_metadata[rownames(feature_metadata) %in% features, , drop = FALSE]
features <- rownames(feature_metadata)


mat_list <- list()
for (l in lineages) {
  fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, -1]
  rownames(fitted_matrix) <- paste0(rownames(fitted_matrix), l)
  mat_list[[l]] <- t(fitted_matrix[, features])
}
mat_raw <- do.call(cbind, mat_list)


mat <- FunPreheatmapNoLog(mat_raw)
mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)]), na.rm = TRUE) * ifelse(mat[is.infinite(mat)] > 0, 1, -1)
mat[is.na(mat)] <- mean(mat, na.rm = TRUE)
feature_split_by <- lineages
mat_split <- mat[, unlist(cell_order_list[feature_split_by]), drop = FALSE]



b <- ceiling(min(abs(quantile(mat, c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 2) / 2
colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))

flip <- FALSE
for (l in lineages) {
  ha_top_list[[l]] <- HeatmapAnnotation(
    Pseudotime = anno_simple(
      x = cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l],
      #col = pseudotime_col,
      which = ifelse(flip, "row", "column"),
      na_col = "transparent",
      border = TRUE
    ), which = ifelse(flip, "row", "column"), show_annotation_name = l == lineages[1], annotation_name_side = ifelse(flip, "top", "left")
  )
}
mat_split <- mat[, unlist(cell_order_list[feature_split_by]), drop = FALSE]
row_split_raw <- row_split <- feature_split <- NULL

n_split <- 5
feature_y <- feature_metadata[rownames(mat_split), order_by]
names(feature_y) <- rownames(mat_split)
km <- kmeans(feature_y, centers = n_split, iter.max = 1e4, nstart = 20)
row_split <- feature_split <- km$cluster

df <- data.frame(row_split = row_split, order_by = feature_metadata[names(row_split), order_by])
df_order <- aggregate(df, by = list(row_split), FUN = mean)
df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), , drop = FALSE]



split_levels <- c()
for (i in seq_len(nrow(df_order))) {
  raw_nm <- df_order[i, "row_split"]
  feature_split[feature_split == raw_nm] <- paste0("C", i)
  level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
  row_split[row_split == raw_nm] <- level
  split_levels <- c(split_levels, level)
}
row_split_raw <- row_split <- factor(row_split, levels = split_levels)
feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))

feature_metadata[["feature_split"]] <- feature_split

dend <- cluster_within_group(t(mat_split), row_split_raw)
cluster_rows <- dend
row_split <- length(unique(row_split_raw))

mat_cluster <- mat[, unlist(cell_order_list[cluster_features_by]), drop = FALSE]




#' extra the group information
library(ComplexHeatmap)
set.seed(372)
m = matrix(rnorm(120), nc = 12)
colnames(m) = 1:12
fa = rep(c("a", "b", "c"), times = c(2, 4, 6))
fa_col = c("a" = 2, "b" = 3, "c" = 4)
dend1 = cluster_between_groups(m, fa)
a=Heatmap(m, cluster_columns = dend1, column_split = 3,
          row_title = "cluster_between_groups",
          top_annotation = HeatmapAnnotation(foo = fa, col = list(foo = fa_col)))
column_order(a)


dend2 = cluster_within_group(m, fa)
Heatmap(m, cluster_columns = dend2, column_split = 3,
        row_title = "cluster_within_group",
        top_annotation = HeatmapAnnotation(foo = fa, col = list(foo = fa_col)))
#,t(srt@assays[[assay]]@data[intersect(cell_annotation, rownames(srt@assays[[assay]])) %||% integer(), rownames(cell_metadata), drop = FALSE])







##'  study RunDynamicFeatures 
#https://github.com/zhanghao-njmu/SCP/blob/main/R/SCP-analysis.R
library(SCP)
library(dplyr)
library(Seurat)
options(reticulate.conda_binary = "/home/chenzh/miniconda3/condabin/conda",SCP_env_name = "R4.2")
data("pancreas_sub")
BPPARAM <- BiocParallel::MulticoreParam(workers = 10)


pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")

srt <- pancreas_sub; lineages <- "Lineage1"; features = NULL; suffix = lineages;
n_candidates = 1000; minfreq = 5; family = NULL;slot = "counts"; assay = NULL; libsize = NULL;
BPPARAM <- BPPARAM; seed = 11


assay <- "RNA"



meta <- c()
gene <- c()


Y <- GetAssayData(srt, slot = slot, assay = assay)
if (is.null(libsize)) {
  status <- check_DataType(srt, assay = assay, slot = "counts")
  if (status != "raw_counts") {
    Y_libsize <- setNames(rep(1, ncol(srt)), colnames(srt))
  } else {
    Y_libsize <- colSums(GetAssayData(srt, slot = "counts", assay = assay))
  }
} 

if (length(meta) > 0) {
  Y <- rbind(Y, t(srt@meta.data[, meta, drop = FALSE]))
}

features_list <- c()
srt_sub_list <- list()
for (l in lineages) {
  srt_sub <- subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])])
  if (is.null(features)) {
    if (is.null(n_candidates)) {
      stop("'features' or 'n_candidates' must provided at least one.")
    }
    HVF <- VariableFeatures(FindVariableFeatures(srt_sub, nfeatures = n_candidates, assay = assay), assay = assay)
    HVF_counts <- srt_sub[[assay]]@counts[HVF, , drop = FALSE]
    HVF <- HVF[apply(HVF_counts, 1, function(x) {
      length(unique(x))
    }) >= minfreq]
    features_list[[l]] <- HVF
  } else {
    features_list[[l]] <- features
  }
  srt_sub_list[[l]] <- srt_sub
}
features <- unique(unlist(features_list))
gene <- features[features %in% rownames(srt[[assay]])]
meta <- features[features %in% colnames(srt@meta.data)]

if (slot == "counts") {
  gene_status <- status
}
gene_status <- status <- check_DataType(srt, assay = assay, slot = slot)
meta_status <- sapply(meta, function(x) {
  check_DataType(data = srt[[x]])
})
if (is.null(family)) {
  family <- rep("gaussian", length(features))
  names(family) <- features
  family[names(meta_status)[meta_status == "raw_counts"]] <- "nb"
  if (gene_status == "raw_counts") {
    family[gene] <- "nb"
  }
} 

for (i in seq_along(lineages)) {
  l <- lineages[i]
  srt_sub <- srt_sub_list[[l]]
  t <- srt_sub[[l, drop = TRUE]]
  t <- t[is.finite(t)]
  t_ordered <- t[order(t)]
  Y_ordered <- as_matrix(Y[features, names(t_ordered), drop = FALSE])
  l_libsize <- Y_libsize[names(t_ordered)]
  raw_matrix <- as_matrix(cbind(data.frame(pseudotime = t_ordered), t(Y_ordered)))
  
  gam_out <- list()
  for (n in seq_len(nrow(Y_ordered))) {
    print(n)
    #n=1
    feature_nm <- rownames(Y_ordered)[n]
    family_current <- family[feature_nm]
    if (min(Y_ordered[feature_nm, ]) < 0 && family_current %in% c("nb", "poisson", "binomial")) {
      warning("Negative values detected. Replace family with 'gaussian' for the feature: ", feature_nm, immediate. = TRUE)
      family_use <- "gaussian"
    } else {
      family_use <- family_current
    }
    
    l_libsize <- l_libsize
    
    sizefactror <- median(Y_libsize) / l_libsize
    mod <- mgcv::gam(y ~ s(x, bs = "cs") + offset(log(l_libsize)),family = family_use,data = data.frame(y = Y_ordered[feature_nm, ], x = t_ordered, l_libsize = l_libsize))
    
    pre <- predict(mod, type = "link", se.fit = TRUE)
    upr <- pre$fit + (2 * pre$se.fit)
    lwr <- pre$fit - (2 * pre$se.fit)
    upr <- mod$family$linkinv(upr)
    lwr <- mod$family$linkinv(lwr)
    res <- summary(mod)
    fitted <- fitted(mod)
    pvalue <- res$s.table[[4]]
    dev.expl <- res$dev.expl
    r.sq <- res$r.sq
    fitted.values <- fitted * sizefactror
    upr.values <- upr * sizefactror
    lwr.values <- lwr * sizefactror
    exp_ncells <- sum(Y_ordered[feature_nm, ] > min(Y_ordered[feature_nm, ]), na.rm = TRUE)
    peaktime <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
    valleytime <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])
    
    # ggplot(data = data.frame(
    #   x = t_ordered,
    #   raw = FetchData(srt_sub, vars = feature_nm, slot = "counts")[names(t_ordered), feature_nm, drop = TRUE],
    #   fitted = fitted.values,
    #   upr.values = upr.values,
    #   lwr.values = lwr.values,
    #   l_libsize = l_libsize
    # )) +
    #   geom_point(aes(x = x, y = raw), color = "black", size = 0.5) +
    #   geom_point(aes(x = x, y = fitted), color = "red", size = 0.5) +
    #   geom_path(aes(x = x, y = upr.values), color = "blue") +
    #   geom_path(aes(x = x, y = lwr.values), color = "green")
    # 
    # a <- data.frame(x = t_ordered, y = Y_ordered[feature_nm, ])
    # qplot(a$x, a$y)
    # length(unique(a$y) > 5)
    # a <- a[a$y > 0, ]
    # qplot(a$x, a$y)
    gam_out[[feature_nm]] <- list(
      features = feature_nm, exp_ncells = exp_ncells,
      r.sq = r.sq, dev.expl = dev.expl,
      peaktime = peaktime, valleytime = valleytime,
      pvalue = pvalue, fitted.values = fitted.values,
      upr.values = upr.values, lwr.values = lwr.values
    )
  }
  
  
  
  
  fitted_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["fitted.values"]]))
  colnames(fitted_matrix) <- rownames(Y_ordered)
  fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)
  
  upr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["upr.values"]]))
  colnames(upr_matrix) <- rownames(Y_ordered)
  upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)
  
  lwr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["lwr.values"]]))
  colnames(lwr_matrix) <- rownames(Y_ordered)
  lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)
  
  DynamicFeatures <- as.data.frame(do.call(rbind.data.frame, lapply(gam_out, function(x) x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")])))
  char_var <- c("features")
  numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
  DynamicFeatures[, char_var] <- lapply(DynamicFeatures[, char_var, drop = FALSE], as.character)
  DynamicFeatures[, numb_var] <- lapply(DynamicFeatures[, numb_var, drop = FALSE], as.numeric)
  rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
  DynamicFeatures[, "padjust"] <- p.adjust(DynamicFeatures[, "pvalue", drop = TRUE])
}
