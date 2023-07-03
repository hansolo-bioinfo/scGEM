# scGEM
Calculation of nested hierarchical dirichlet process for gene coexpressing modules (GEM) in single-cell transcriptome data 

# Usage
Please go through the <run_this.R> for the tutorial. For Window users, please source corresponding R codes prior to model building. More details please look into the source code.

#### required model input 
x: dgCMatrix of single cell count data (row name is gene symbol, column name is cell index).  
num_topics: vector of children topics per parent in each level.

#### model output 
centroids: gene by gem weight matrix (&beta;). Row name is gene symbol, column name is GEM name.  
count_mat: gem by cell count matrix (&theta;). Row name is GEM name, column name is cell idx.
