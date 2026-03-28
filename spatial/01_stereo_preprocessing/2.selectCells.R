setwd("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/CalculateCells/SelectCells")
.libPaths(c("/cluster/apps/anaconda3/2021.05/envs/R-4.1.2/lib/R/library"))
library(fastmap)
# omics
library(dplyr)
library(Signac)
library(Seurat)
library(SummarizedExperiment)
library(GenomeInfoDb)
source("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/CalculateCells/SelectCells/xselectCells.R")


#--------------------------------------GW7--------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/GW7_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
#xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

# gw7-1
subdata1 <- subset(obj, cod_1>8000 & cod_1<12000 & cod_2>12000 & cod_2<17000)
FeaturePlot(subdata1,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(subdata1,reduction = "cod",feature = "ATOH1",raster=FALSE)
xSelectCells(subdata1, type = "GEM",feature = "nFeature_RNA")
cell1 = read.csv("./Barcode_selected_gw7-1.csv",header = F)
cells = as.character(cell1$V1)
obj1 = subset(obj,cells = cells)
saveRDS(obj1,file = "./RDS/SeuratObj_gw7-1.rds")

# gw7-2
subdata2 <- subset(obj, cod_1>9000 & cod_1<15000 & cod_2>4000 & cod_2<7000)
FeaturePlot(subdata2,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(subdata2,reduction = "cod",feature = "ATOH1",raster=FALSE)
xSelectCells(subdata2, type = "GEM",feature = "nFeature_RNA")
cell2 = read.csv("./Barcode_selected_gw7-2.csv",header = F)
cells = as.character(cell2$V1)
obj2 = subset(obj,cells = cells)
saveRDS(obj2,file = "./RDS/SeuratObj_gw7-2.rds")

# gw7-3
subdata3 <- subset(obj, cod_1>14000 & cod_1<18000 & cod_2>14000 & cod_2<20000)
FeaturePlot(subdata3,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(subdata3,reduction = "cod",feature = "ATOH1",raster=FALSE)
xSelectCells(subdata3, type = "GEM",feature = "nFeature_RNA")
cell3 = read.csv("./Barcode_selected_gw7-3.csv",header = F)
cells = as.character(cell3$V1)
obj3 = subset(obj,cells = cells)
saveRDS(obj3,file = "./RDS/SeuratObj_gw7-3.rds")


#--------------------------------------GW6--------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/GW6_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
#xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

# gw6
subdata <- subset(obj, cod_1>4000 & cod_1<6500 & cod_2>5000 & cod_2<7000)
FeaturePlot(subdata,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(subdata,reduction = "cod",feature = "ATOH1",raster=FALSE)
FeaturePlot(subdata,reduction = "cod",feature = "SOX2",raster=FALSE)
FeaturePlot(subdata,reduction = "cod",feature = "GDF10",raster=FALSE)
FeaturePlot(subdata,reduction = "cod",feature = "BNC2",raster=FALSE) #ex 
FeaturePlot(subdata,reduction = "cod",feature = "TUBB3",raster=FALSE) #ex
xSelectCells(subdata, type = "GEM",feature = "nFeature_RNA")
cell = read.csv("./Barcode_selected_gw6.csv",header = F)
cells = as.character(cell$V1)
obj = subset(obj,cells = cells)
saveRDS(obj,file = "./RDS/SeuratObj_gw6.rds")


#--------------------------------------GW9--------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/GW9_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
#xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

# gw9
subdata <- subset(obj, cod_1<17500 )
FeaturePlot(subdata,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
saveRDS(subdata,file = "./RDS/SeuratObj_gw9.rds")


#--------------------------------------GW10-------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/GW10_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
#xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

# gw10
subdata <- subset(obj, cod_1>10000 & cod_1<16000 & cod_2>8000 & cod_2<13000)
FeaturePlot(subdata,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
xSelectCells(subdata, type = "GEM",feature = "nFeature_RNA")
cell = read.csv("./Barcode_selected_gw10.csv",header = F)
cells = as.character(cell$V1)
obj = subset(obj,cells = cells)
saveRDS(obj,file = "./RDS/SeuratObj_gw10.rds")


#--------------------------------------Y31-------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/Y31_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

saveRDS(obj,file = "./RDS/SeuratObj_y31.rds")


#--------------------------------------Y32-------------------------------------
object = readRDS("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/bin30Data/Y32_bin30_seurat.rds")
counts = object@assays$Spatial@counts
#rownames(counts) = as.character(object@misc$raw_genename)
obj = CreateSeuratObject(counts = counts,
                         assay = "RNA",
                         meta.data = object@meta.data)
cod = obj@meta.data[,c("x","y")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod_"
)
obj@reductions[["cod"]] = drobj

# overview
FeaturePlot(obj,reduction = "cod",feature = "nFeature_RNA",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
#xSelectCells(obj,type = "GEM",feature = "nFeature_RNA")

saveRDS(obj,file = "./RDS/SeuratObj_y32.rds")


FeaturePlot(obj,reduction = "cod",feature = "BNC2",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "TUBB3",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SOX2",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "MNX1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "ISL1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "NTRK1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SLC5A7",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "FLT1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "CLDN5",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "ITGA1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "PDGFRB",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "COL1A1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "COL1A2",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "TTC6",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "PPP1R17",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "GDF7",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "GDF10",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "MSX1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "STMN2",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "ASCL1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "FGFR3",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "ATOH1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "FOXA2",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "SHH",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "ZIC1",raster=FALSE)
FeaturePlot(obj,reduction = "cod",feature = "MBP",raster=FALSE)





