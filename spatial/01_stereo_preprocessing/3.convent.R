setwd("/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/CalculateCells/SelectCells")
.libPaths(c("/cluster2/huanglab/jiamao/conda/envs/sceasy/lib/R/library"))
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')

# 指定文件夹路径
input_directory <- "/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/CalculateCells/SelectCells/RDS"
output_directory <- "/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/data/CalculateCells/SelectCells/H5AD"

# 获取该目录下所有RDS文件
files <- list.files(input_directory, pattern = "*.rds", full.names = TRUE)

# 遍历每个文件并进行转换
for (file in files) {
  # 获取文件名并修改扩展名
  output_file <- gsub("rds$", "h5ad", basename(file))
  # 设置输出文件路径
  output_path <- file.path(output_directory, output_file)
  # 如果输出文件已存在，则跳过
  if (file.exists(output_path)) {
    cat("Skip (already exists):", output_path, "\n")
    next
  }
  # 读取Seurat对象
  seurat_object <- readRDS(file)
  # 进行格式转换并保存
  sceasy::convertFormat(seurat_object, from = "seurat", to = "anndata", outFile = output_path)
  # 打印处理进度
  cat("Converted:", file, "->", output_path, "\n")
}