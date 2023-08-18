# scGEM

scGEM is a nested tree-structured generative model that identifies subtype-specific and -shared gene coexpressing modules (GEMs) via scRNAseq data.

## Installation

Currently, scGEM works on MAC/Linux. For faster implementation on MAC, you need to switch your built-in BLAS library to Apple's vecLib. (https://mpopov.com/blog/2021/10/10/even-faster-matrix-math-in-r-on-macos-with-m1/). 

```R
devtools::install_github("hansolo-bioinfo/scGEM")
?initTree
?minibatchInfer
```

Note: if the `Error in Lazy database` pops up, try `.rs.restartR()`. For more information, please see (https://github.com/r-lib/devtools/issues/1980).


## Getting Started

To test the scGEM performance, we will use processed peripheral blood mononuclear cells (PBMC) 3K dataset as our toy case.

### 1) preprocessing

```R
library(Seurat)
library(SeuratData)
library(magrittr)
library(scGEM)

data("pbmc3k.final")
set.seed(111)
X <- pbmc3k.final@assays$RNA@counts
X@x <- log2(X@x + 1) # log-transform the count data

# for binary learning
# X@x <- ifelse(X@x > 0, 1)
```

### 2) initialization

```R
# assuming create a 5-4-3 tree structure
scGEM_init <- initTree(
    X,
    num_topics = c(5, 4, 3),
    blacklist_genes = "^MT|^RP|^HSP|B2M|MALAT1",
    rm.housekeeping_genes = TRUE
)
```

### 3) learning

For one-batch training, please set the `batch_size` equals to the sample size:

```R
# one-batch
scGEM_trained <- minibatchInfer(
    X,
    scGEM_init,
    model_paras = list(b0 = 0.01,
                       g1 = 5,
                       g2 = 1,
                       g3 = 1/3,
                       g4 = 2/3),
    max_est = 50,
    subtree_size = c(1, 20),
    batch_size = ncol(X),
    n_epoch = 50,
    learning_rate = 0.01
)
```

For mini-batch training (recommended), please do:

```R
# mini-batch
scGEM_trained_mb <- minibatchInfer(
    X,
    scGEM_init,
    model_paras = list(b0 = 0.01,
                       g1 = 5,
                       g2 = 1,
                       g3 = 1/3,
                       g4 = 2/3),
    max_est = 50,
    subtree_size = c(1, 20),
    batch_size = 1000,
    n_epoch = 50,
    learning_rate = 0.01
)
```

### 4) plot likelihood

Compare the total likelihood during each epoch in one-batch and mini-batch modes.

```R
# total likelihood
plot(scGEM_trained$likelihood, col = "blue", type = "b",
     ylab = "Likelihood", xlab = "Epoch")
points(scGEM_trained_mb$likelihood, col = "red", type = "b")
legend("bottomright", legend = c("One-batch", "Mini-batch"),
       lty = 1, col = c("blue", "red"))
```

<img width="479" alt="Screenshot 2023-08-18 at 10 57 25" src="https://github.com/hansolo-bioinfo/scGEM/assets/65295899/334c3181-cb3b-4e41-a492-43ad2ef57312">

### 5) extract distributions

To get the raw values:

```R
gene_over_gem <- scGEM_trained_mb$centroids       # 10640 x 85
gem_over_cell <- scGEM_trained_mb$count_matrix    # 85 x 2638
tree_relation <- scGEM_trained_mb$tree            # 86 x 2
gene_name     <- scGEM_trained_mb$gene
cell_name     <- scGEM_trained_mb$cell
gem_name      <- paste0("GEM", 1:85)

rownames(gene_over_gem) <- gene_name
colnames(gene_over_gem) <- gem_name
rownames(gem_over_cell) <- gem_name
colnames(gem_over_cell) <- cell_name
```

To get a relative proportion like in other topic models:

```R
gene_over_gem <- apply(gene_over_gem, 2, function(x) x/sum(x))
gem_over_cell <- apply(gem_over_cell, 2, function(x) x/sum(x))
```

Plot gene distribution for GEM 1

```R
barplot(sort(gene_over_gem[, 1], decreasing = T), ylab = "weight", xlab = "gene", main = "GEM 1")
```

<img width="643" alt="Screenshot 2023-08-18 at 12 28 38" src="https://github.com/hansolo-bioinfo/scGEM/assets/65295899/a191f3a2-ac96-4fbd-a91c-a5d9ace66aaa">

Plot GEM distribution for the 100th cell

```R
barplot(gem_over_cell[, 100], ylab = "weight", xlab = "GEM", 
        main = "100th cell: AAGATTACCGCCTT")
```

<img width="673" alt="Screenshot 2023-08-18 at 12 31 24" src="https://github.com/hansolo-bioinfo/scGEM/assets/65295899/adced07b-67a8-47bd-90db-46492ebfad6e">
