# scGEM

scGEM is a nested tree-structured generative model that identifies subtype-specific and -shared GEMs via scRNAseq data.

## Installation

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

For one-batch training, please set the `batch_size` equals to the sample size.

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

For mini-batch training, please do

# mini-batch
```R
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

```R
# total likelihood
plot(scGEM_trained$likelihood, col = "blue", type = "b",
     ylab = "Likelihood", xlab = "Epoch")
points(scGEM_trained_mb$likelihood, col = "red", type = "b")
legend("bottomright", legend = c("One-batch", "Mini-batch"),
       lty = 1, col = c("blue", "red"))
```
