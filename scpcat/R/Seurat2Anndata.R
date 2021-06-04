library(Seurat)
library(SeuratDisk)

Seurat2AnnData<-function (filename=NULL, name = "seurat",outDir = "anndata")
{   
    if(!dir.exists(outDir)){
       	    dir.create(outDir)
    }
    message("loading rds file..")
    seurat<-readRDS(filename)
    metadata <- seurat@meta.data
    metadata$barcode <- Cells(seurat)
    write.table(metadata, file.path(outDir, "metadata.csv"),
        sep = ",", quote = FALSE, row.names = FALSE)
    message("INFO : SaveH5Seurat ...")
    dest <- file.path(outDir, paste0(name, ".h5Seurat"))
    SaveH5Seurat(seurat, filename = dest)
    message("INFO : Convert into h5ad ...")
    Convert(dest, dest = "h5ad")
    message("INFO : Done!")
}

