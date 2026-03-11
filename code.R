

install.packages("Seurat")
install.packages("SeuratData")


library(Seurat)
library(SeuratData)


AvailableData()


InstallData("panc8")


data("panc8")

panc8

data("panc8")


panc8 <- UpdateSeuratObject(panc8)

panc8



dim(panc8)


table(panc8$celltype)


cell_counts <- table(panc8$celltype)
cell_counts


major_types <- names(cell_counts[cell_counts > 200])
major_types


panc_major <- subset(panc8, subset = celltype %in% major_types)

table(panc_major$celltype)

dim(panc_major)


panc_major[["percent.mt"]] <- PercentageFeatureSet(panc_major, pattern = "^MT-")


VlnPlot(
  panc_major,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)


FeatureScatter(panc_major, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(panc_major, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



panc_major <- subset(
  panc_major,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)


dim(panc_major)


panc_major <- NormalizeData(
  panc_major,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


panc_major <- FindVariableFeatures(
  panc_major,
  selection.method = "vst",
  nfeatures = 2000
)



head(VariableFeatures(panc_major))


VariableFeaturePlot(panc_major)


top10 <- head(VariableFeatures(panc_major), 10)

LabelPoints(
  plot = VariableFeaturePlot(panc_major),
  points = top10,
  repel = TRUE
)


panc_major <- ScaleData(panc_major)



panc_major <- RunPCA(
  panc_major,
  features = VariableFeatures(object = panc_major)
)

print(panc_major[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(panc_major, reduction = "pca")

ElbowPlot(panc_major)

dims = 1:20


DimHeatmap(panc_major, dims = 1:10, cells = 500, balanced = TRUE)

panc_major <- FindNeighbors(
  panc_major,
  dims = 1:20
)


panc_major <- FindClusters(
  panc_major,
  resolution = 0.5
)


panc_major <- RunUMAP(
  panc_major,
  dims = 1:20
)


DimPlot(
  panc_major,
  reduction = "umap",
  label = TRUE
)


DimPlot(
  panc_major,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE
)


saveRDS(panc_major, "panc8_processed.rds")

markers <- FindAllMarkers(panc_major)


markers <- FindAllMarkers(
  panc_major,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


head(markers)

library(dplyr)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)


top10



DoHeatmap(
  panc_major,
  features = top10$gene
)


write.csv(markers, "pancreas_cluster_markers.csv")




expr <- GetAssayData(panc_major, slot = "data")

expr <- as.data.frame(t(expr))

dim(expr)

labels <- panc_major$celltype

table(labels)


features <- unique(top10$gene)

expr_subset <- expr[, features]


dim(expr_subset)


ml_data <- expr_subset
ml_data$celltype <- labels


head(ml_data)



install.packages("caret")
library(caret)


expr <- GetAssayData(panc_major, slot = "data")


expr <- as.data.frame(t(expr))

dim(expr)


labels <- panc_major$celltype

table(labels)


features <- unique(top10$gene)

features <- unique(top10$gene)
expr_subset <- expr[, features]
dim(expr_subset)


ml_data <- expr_subset
ml_data$celltype <- labels

head(ml_data)

set.seed(123)

trainIndex <- createDataPartition(
  ml_data$celltype,
  p = 0.8,
  list = FALSE
)

train_data <- ml_data[trainIndex, ]
test_data <- ml_data[-trainIndex, ]


set.seed(123)

trainIndex <- createDataPartition(
  ml_data$celltype,
  p = 0.8,
  list = FALSE
)

train_data <- ml_data[trainIndex, ]
test_data <- ml_data[-trainIndex, ]

model_log <- train(
  celltype ~ .,
  data = train_data,
  method = "multinom"
)

pred_log <- predict(model_log, test_data)



ml_data$celltype <- as.factor(ml_data$celltype)



library(caret)

# remove near-zero variance genes
nzv <- nearZeroVar(train_data)

train_data <- train_data[, -nzv]
test_data <- test_data[, -nzv]

cor_matrix <- cor(train_data[, -ncol(train_data)])

high_cor <- findCorrelation(cor_matrix, cutoff = 0.9)

train_data <- train_data[, -high_cor]
test_data <- test_data[, -high_cor]

dim(train_data)

model_log <- train(
  celltype ~ .,
  data = train_data,
  method = "multinom",
  trace = FALSE
)

pred_log <- predict(model_log, test_data)

confusionMatrix(pred_log, test_data$celltype)

library(randomForest)

model_rf <- randomForest(
  celltype ~ .,
  data = train_data
)

pred_rf <- predict(model_rf, test_data)

confusionMatrix(pred_rf, test_data$celltype)


nzv <- nearZeroVar(train_data)

train_data <- train_data[, -nzv]

test_data <- test_data[, colnames(train_data)]

test_data <- test_data[, colnames(train_data)]

model_log <- train(
  celltype ~ .,
  data = train_data,
  method = "multinom",
  trace = FALSE
)

ml_data <- expr_subset
ml_data$celltype <- as.factor(labels)


str(ml_data$celltype)

set.seed(123)

trainIndex <- createDataPartition(
  ml_data$celltype,
  p = 0.8,
  list = FALSE
)

train_data <- ml_data[trainIndex, ]
test_data  <- ml_data[-trainIndex, ]


nzv <- nearZeroVar(train_data)

train_data <- train_data[, -nzv]
test_data  <- test_data[, colnames(train_data)]

model_log <- train(
  x = train_data[, -ncol(train_data)],
  y = train_data$celltype,
  method = "multinom",
  trace = FALSE
)
