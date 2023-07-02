library(Seurat)
library(magrittr)
# source("model.R")
#source("model_windows.R")

srt = SeuratDisk::LoadH5Seurat("pbmc3k.h5Seurat")
set.seed(111)
X = srt@assays$RNA@counts
X@x = log2(X@x + 1)
# for binary
# X@x = ifelse(X@x > 0, 1)
meta = srt@meta.data
if (length(grep("T", levels(Idents(srt)))) > 0) {
    meta$Cell_subtype = Idents(srt)
}
umap = srt@reductions$umap@cell.embeddings
meta = data.frame(meta, UMAP_1 = umap[, 1], UMAP_2 = umap[, 2])
#rm(srt)

# blacklist_genes = "^MT|^RP|^HSP|B2M|MALAT1"
model_paras = list(b0 = 0.01, g1 = 5, g2 = 1, g3 = 1/3, g4 = 2/3)
# if you have in-house housekeeping gene list, set to FALSE
# require 'Housekeeping_GenesHuman.csv' and change path in < model.R >
nHDP_init = initTree(X, 
                     num_topics = c(5, 4, 3), 
                     blacklist_genes = "^MT|^RP|^HSP|B2M|MALAT1", 
                     housekeeping_genes = TRUE)
# one-batch
nHDP_trained = minibatchInfer(X, 
                              nHDP_init, 
                              model_paras, 
                              max_est = 50, 
                              subtree_size = c(1, 20),
                              batch_size = ncol(X), 
                              n_epoch = 50, 
                              learning_rate = 0.01)
# mini-batch
nHDP_trained_mb = minibatchInfer(X, 
                                 nHDP_init, 
                                 model_paras, 
                                 max_est = 50, 
                                 subtree_size = c(1, 20),
                                 batch_size = 1500, 
                                 n_epoch = 50, 
                                 learning_rate = 0.01)
save(list = c("nHDP_init", "nHDP_trained", "nHDP_trained_mb", "X", "meta"), 
     file = "pbmc3k.RData")