#' Smmit: single-cell multi-sample and multi-omics integration
#'
#' Smmit performs integration both across samples and modalities to produce a single UMAP space. It first uses harmony to integrate across samples and then uses Seurat weighted nearest neighbor function to integrate across modalities.
#'
#' @param obj A list of named Seurat objects. Each Seurat object contains single-cell multi-omics data from a single sample, processed using the standard Seurat or Signac pipeline. The name of each element is the sample name.
#' @param mode Either RNA_ATAC or RNA_ADT. RNA_ATAC performs integration for single-cell multi-omics data that jointly profile gene expression and chromatin accessibility. RNA_ADT performs integration for single-cell multi-omics data that jointly profile gene expression and protein abundances.
#' @param RNA.normalization.method Normalization method for RNA. See Seurat::NormalizeData()
#' @param RNA.nfeatures Number of variable features for RNA. See Seurat::FindVariableFeatures()
#' @param RNA.npcs Number of PCs for RNA. See Seurat::RunPCA()
#' @param RNA.integrate.dims Number of dimensions for RNA during the modality integration process. See Seurat::FindMultiModalNeighbors()
#' @param ATAC.nsvd Number of singular values for ATAC. See Signac::RunSVD()
#' @param ATAC.integrate.dims Number of dimensions for ATAC during the modality integration process. See Seurat::FindMultiModalNeighbors()
#' @param ADT.npcs Number of PCs for ADT. See Seurat::RunPCA()
#' @import Seurat Signac harmony 
#' @export
#' @return A Seurat object. The UMAP that integrates both across samples and modalities is stored in the wsnnumap reduction. The Seurat object can be directly used in downstream analysis such as cell clustering, cell type identification, and differential analysis.
#' @author Changxin Wan, Zhicheng Ji<zhicheng.ji@@duke.edu>

smmit <- function(obj,mode='RNA_ATAC',RNA.normalization.method='LogNormalize',RNA.nfeatures=2000,RNA.npcs=30,RNA.integrate.dims=30,ATAC.nsvd=50,ATAC.integrate.dims=30,ADT.npcs=NULL) {
  if (mode=='RNA_ATAC') {
    
    ### Adding sample names
    for (i in names(obj)) {
      obj[[i]]@meta.data$orig.ident <- i
    }
    
    ### Merge objects
    eval(parse(text=paste0('obj <- merge(obj[[1]],y = c(',paste0(paste0('obj[[',2:length(obj),']]'),collapse = ','),'), add.cell.ids = names(obj), project = "merge")')))
    
    ### RNA across sample integration
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj,normalization.method=RNA.normalization.method)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = RNA.nfeatures)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs=RNA.npcs)
    obj <- RunHarmony(object = obj, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'RNA', project.dim = FALSE)
    obj[["integrated_rna"]] <- CreateDimReducObject(embeddings = obj[['harmony']]@cell.embeddings, key = "integratedRNA_", assay = 'RNA')
    obj[['harmony']] <- NULL
    
    ### ATAC across sample integration
    DefaultAssay(obj) <- "ATAC"
    obj <- FindTopFeatures(obj)
    obj <- RunTFIDF(obj)
    obj <- RunSVD(obj,n=ATAC.nsvd)
    obj <- RunHarmony(object = obj, group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
    obj[["integrated_atac"]] <- CreateDimReducObject(embeddings = obj[['harmony']]@cell.embeddings, key = "integratedATAC_", assay = 'ATAC')
    obj[['harmony']] <- NULL
    
    ### Integrate across modalities
    obj <- FindMultiModalNeighbors(
      object = obj,
      reduction.list = list("integrated_rna", "integrated_atac"),
      dims.list = list(1:RNA.integrate.dims, 2:ATAC.integrate.dims),
      modality.weight.name = "RNA.weight",
      verbose = TRUE
    )
    
    ### Build a joint UMAP visualization
    obj <- RunUMAP(
      object = obj,
      nn.name = "weighted.nn",
      assay = "ATAC",
      verbose = TRUE,
      reduction.key = "wsnnUMAP_",
      reduction.name = "wsnnumap"
    )
    
  } else if (mode=='RNA_ADT') {
    
    ### ADT filtering and adding sample names
    adtlist <- NULL
    for (i in names(obj)) {
      adtlist <- c(adtlist,rownames(obj[[i]]$ADT@counts))
    }
    tabadtlist <- table(adtlist)
    selectadt <- names(tabadtlist)[tabadtlist==length(obj)]
    
    for (i in names(obj)) {
      obj[[i]][["ADT"]] <- CreateAssayObject(counts = obj[[i]]$ADT@counts[selectadt,])
      obj[[i]]@meta.data$orig.ident <- i
    }
    
    ### Merge objects
    eval(parse(text=paste0('obj <- merge(obj[[1]],y = c(',paste0(paste0('obj[[',2:length(obj),']]'),collapse = ','),'), add.cell.ids = names(obj), project = "merge")')))
    
    ### RNA across sample integration
    DefaultAssay(obj) <- 'RNA'
    obj <- NormalizeData(obj,normalization.method=RNA.normalization.method)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = RNA.nfeatures)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs=RNA.npcs)
    obj <- RunHarmony(object = obj, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'RNA', project.dim = FALSE)
    obj[["integrated_rna"]] <- CreateDimReducObject(embeddings = obj[['harmony']]@cell.embeddings, key = "integrateRNA_", assay = 'RNA')
    obj[['harmony']] <- NULL
    
    ### ADT across sample integration
    DefaultAssay(obj) <- 'ADT'
    VariableFeatures(obj) <- rownames(obj[["ADT"]])
    obj <- NormalizeData(obj,normalization.method = 'CLR', margin = 2)
    obj <- ScaleData(obj)
    if (is.null(ADT.npcs)) {
      obj <- RunPCA(obj, reduction.name = 'apca')  
    } else {
      obj <- RunPCA(obj, reduction.name = 'apca',npcs=ADT.npcs)
    }
    obj <- RunHarmony(object = obj, group.by.vars = 'orig.ident', reduction = 'apca', assay.use = 'ADT', project.dim = FALSE)
    obj[["integrated_adt"]] <- CreateDimReducObject(embeddings = obj[['harmony']]@cell.embeddings, key = "integrateADT_", assay = 'ADT')
    obj[['harmony']] <- NULL
    
    ### Integrate across modalities
    ncdim <- min(RNA.integrate.dims,ncol(obj[['integrated_adt']]@cell.embeddings))
    obj <- FindMultiModalNeighbors(
      object = obj,
      reduction.list = list("integrated_rna", "integrated_adt"),
      dims.list = list(1:ncdim, 1:ncdim),
      modality.weight.name = "weight",
      verbose = TRUE
    )
    
    ### Build a joint UMAP visualization
    obj <- RunUMAP(
      object = obj,
      nn.name = "weighted.nn",
      assay = "RNA",
      verbose = TRUE,
      reduction.key = "wsnnUMAP_",
      reduction.name = "wsnnumap"
    )    
  }
}
