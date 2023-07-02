#' initialize the global tree based on nested mini-batch kmeans
#' 
#' @param x dgCMatrix of single cell count data.
#' (row name is gene name; column name is cell index)
#' @param num_topics vector of children topics per parent in each level.
#' (example: c(5, 4, 3))
#' @param min_cells minimum number of cells required for each topic, 
#' topics fail to meet this requirement will be removed.
#' (default is min(num_topics))
#' @param blacklist_genes an expression of set of genes to be removed.
#' (default is NULL, example: "^MT|^RPS")
#' @param housekeeping_genes whether or not consider housekeeping genes.
#' (default is FALSE, do not consider)
#' 
#' @importFrom mbkmeans mbkmeans
#' 
#' @return centroids and counts of topics and parent-children relationship
initTree <- function(
        x, 
        num_topics, 
        min_cells = NULL, 
        blacklist_genes = NULL,
        housekeeping_genes = FALSE
        ) {
    
    if (! is(x, "Matrix")) {
        stop("x is not in dgCMatrix format")
    }
    
    if (is.null(min_cells)) {
        min_cells <- min(num_topics)
    }
    
    if (! is.null(blacklist_genes)) {
        rm_genes <- grep(blacklist_genes, rownames(x))
        x <- x[-rm_genes, ]
    }
    
    if (! housekeeping_genes) {
        hk <- read.table("~/Lab/Data/Databases/Housekeeping_GenesHuman.csv", 
                         sep = ";", head = TRUE)
        hk <- hk$Gene.name
        rm_hk <- which(rownames(x) %in% hk)
        x <- x[-rm_hk, ]
    }
    
    # normalize
    x <- normData(x)
    N <- ncol(x)      
    L <- length(num_topics)
    max_topic <- sum(cumprod(num_topics))
    clusters <- rep(0, N)
    size <- rep(0, max_topic)
    centroids <- matrix(0, nr = nrow(x), nc = max_topic)
    
    tree <- data.frame(parent = 0, child = 0)
    
    for (l in 1:L) {
        cat("Building level", l, "...\n")
        
        K <- num_topics[l]    
        S <- sort(unique(clusters))
    
        for (s in S) {
            idx_s <- which(clusters == s)
            x_s <- x[, idx_s]
            
            if (length(idx_s) < min_cells) {
                next
            }
            
            #if (length(idx_s) >= 500) {
            #    batch <- 500
            #} else {
            #    batch <- length(idx_s)
            #}
            
            #mbk <- KMeans(x_s, 
            #              clusters    = K, 
            #              batch_size  = batch, 
            #              num_init    = 3,
            #              initializer = "random",
            #              max_iters   = 3,
            #              min_n       = min_cells
            #)
            
            if (K == 1) {
                mbk <- list()
                mbk$centroids <- Matrix::rowMeans(x_s)
                mbk$clusters <- rep(1, length(idx_s))
                mbk$K <- 1
            } else {
                mbk <- KMeans_L1(x_s, K, 3)
            }
            
            M <- max(clusters, na.rm = T)
            clusters[idx_s] <- mbk$clusters + M
            centroids[, (M + 1):(M + mbk$K)] <- mbk$centroids
            tbl <- table(mbk$clusters)
            size[(M + 1):(M + mbk$K)] <- as.integer(tbl)
            new_child <- as.integer(names(tbl)) + M
            
            # update tree
            tree_expand <- data.frame(parent = s, child = new_child)
            tree <- rbind(tree, tree_expand)
            
            # subtract the mean and change negative to zero and renormalize
            for (k in 1:mbk$K) {
                idx_k <- which(mbk$clusters == k)
                
                if (length(idx_k) > 1) {
                    x_k <- x_s[, idx_k]
                    
                    if (K == 1) {
                        centroid_k <- as.matrix(mbk$centroids)
                    } else {
                        centroid_k <- as.matrix(mbk$centroids[, k])
                    }
                    x_k@x <- x_k@x - centroid_k[x_k@i + 1]
                    x_k@x[x_k@x < 0] <- 0
                    x_k <- normData(x_k)
                    
                    start <- (x@p + 1)[idx_s[idx_k]]
                    end <- (x@p + 1)[idx_s[idx_k] + 1] - 1
                    x_idx <- unlist(mapply(":", start, end))
                    x@x[x_idx] <- x_k@x
                } else {
                    x[, idx_s[idx_k]] <- 0
                }
                
            }
        }
    }
    
    trim_topic <- which(size != 0)
    size <- size[trim_topic]
    centroids <- centroids[, trim_topic]
    
    # create the tree matrix 
    # row is the number of children under this node, col is the parent level
    n_topics <- length(size)
    tree_mat <- matrix(0, nr = n_topics, nc = n_topics)
    for (i in 1:n_topics) {
        idx <- i
        repeat {
            idx <- tree[tree$child == idx, "parent"]
            # if searching back to root, stop
            if (idx == 0) {
                break
            } else {
                tree_mat[idx, i] <- 1
            }
        }
    }
    
    # scale centroid
    result <- list(cell      = colnames(x),
                   gene      = rownames(x),
                   centroids = centroids * N,
                   size      = size,
                   tree      = tree,
                   tree_mat  = tree_mat)
    return(result)
}


#' Infer nHDP tree based on mini-batch variational inference
#' 
#' @param X dgCMatrix single cell count data to learn
#' @param nHDP_init initialized global nHDP tree.
#' @param batch_size sample size for one mini-batch training.
#' @param model_paras list of hyperparameters in the model.
#' (b0, g1, g2, g3, g4)
#' @param n_epoch number of epochs (default is 20).
#' @param subtree_size size of the search tree for each cell.
#' (example: c(5, 20), at least 5, at most 20)
#' @param max_est maximum number of estimation for topic percentages in 
#' each cell. (default is 25)
#' @param likelihood_eps
#' @param count_eps
#' 
#' @importFrom
#' 
#' @return
minibatchInfer <- function(
        X,
        nHDP_init, 
        model_paras,
        batch_size = 1000, 
        n_epoch = 20,
        subtree_size = c(1, 20),
        max_est = 50,
        likelihood_eps = .5*10^-3,
        count_eps = .5*10^-2,
        learning_rate = NULL,
        tol = 1e-4){
    
    miss_gene <- nHDP_init$gene[! nHDP_init$gene %in% rownames(X)]
    if (length(miss_gene) > 0) {
        X_add <- Matrix::Matrix(0, nrow = length(miss_gene), ncol = ncol(X), 
                                sparse = TRUE)
        rownames(X_add) <- miss_gene
        colnames(X_add) <- colnames(X)
        X <- Matrix::rbind2(X, X_add)
    }
    matched_genes <- match(nHDP_init$gene, rownames(X))
    #matched_genes <- match(nHDP_init$gene, rownames(X))
    X <- X[matched_genes, ]
    nHDP_tree <- nHDP_init
    nHDP_tree$cell <- colnames(X)
    N <- length(nHDP_tree$cell)
    nHDP_tree$batch_size <- batch_size
    
    # create a zero matrix for initialized centroid then update with parameters
    #X_centroid <- nHDP_init$centroids[matched_genes, ]
    #X_centroid[is.na(X_centroid)] <- 0
    #nHDP_tree$centroids <- X_centroid
    #nHDP_tree$gene <- rownames(X)
    
    # check hyperparameters
    for (para in c("b0", "g1", "g2", "g3", "g4")) {
        if (! para %in% names(model_paras)) {
            stop(paste("parameter", para, "not set yet"))
        }
    }
    
    likelihood <- rep(0, n_epoch)
    ite_num <- 1
    
    for (epoch in 1:n_epoch) {
        cat("***************************\n")
        cat("Epoch: ", epoch, "/", n_epoch, "\n", sep = "")
        if (batch_size == N) {
            batches <- list()
            batches[[1]] <- sample.int(N)
        } else {
            batches <- split(sample.int(N), ceiling(seq.int(N) / batch_size))
        }
        n_ite <- length(batches)
        count_matrix <- matrix(0, nr = length(nHDP_tree$size), nc = N) 
        
        for (ite in 1:n_ite) {
            cat(" [", ite, "/", n_ite, "]", sep = "")
            ite_num <- ite_num + 1
            
            if (is.null(learning_rate)) {
                step_size <- ite_num ^ (-2)
                
                # v0.1: NEW FEATURE
                # prevent the learning rate from too low
                if (step_size < 0.0001) {
                    step_size <- 0.0001
                }
            } else {
                step_size <- learning_rate
            }
            cat(", Step: ", step_size, "\n", sep = "")
            
            d_idx <- batches[[ite]]
            X_mb <- X[, d_idx]
            
            batch_updates <- updateTree(X_mb, 
                                    nHDP_tree, 
                                    model_paras, 
                                    subtree_size,
                                    max_est,
                                    likelihood_eps,
                                    count_eps,
                                    step_size)
            nHDP_tree <- batch_updates$new_tree
            count_matrix[, d_idx] <- batch_updates$count_matrix
            likelihood[epoch] <- likelihood[epoch] + batch_updates$likelihood
        }
        
        nHDP_tree$count_matrix <- count_matrix
        cat("Total Likelihood:", likelihood[epoch], "\n")
        
        
        if (epoch > 1) {
            add_ll <- (likelihood[epoch - 1] - likelihood[epoch]) / likelihood[epoch - 1]
            if (add_ll < tol) {
                cat("Converged at epoch ", epoch, "!\n", sep = "")
                likelihood <- likelihood[1:epoch]
                break
            }
        }
        
        invisible(gc(reset = TRUE))
    }
    
    nHDP_tree$likelihood <- likelihood
    #if (is.null(nHDP_tree$likelihood)) {
    #    nHDP_tree$likelihood <- likelihood
    #} else {
    #    nHDP_tree$likelihood <- c(nHDP_tree$likelihood, likelihood)
    #}
    
    return(nHDP_tree)
}

#' Update nHDP tree
updateTree <- function(
        x,
        nHDP_tree, 
        model_paras,
        subtree_size,
        max_est,
        likelihood_eps,
        count_eps,
        learning_rate) {
    
    b0 <- model_paras$b0
    g1 <- model_paras$g1
    g2 <- model_paras$g2
    g3 <- model_paras$g3
    g4 <- model_paras$g4
    
    n_genes <- length(nHDP_tree$gene)
    n_cells <- ncol(x)
    n_topics <- length(nHDP_tree$size)
    N <- length(nHDP_tree$cell)
    
    # process the global tree and get the priors
    priors <- getTreePrior(nHDP_tree, b0, g1)
    relation <- nHDP_tree$tree
    tree_mat <- nHDP_tree$tree_mat
    
    parent_level <- colSums(tree_mat)
    level_penalty <- digamma(g3) - digamma(g3 + g4) + 
        parent_level * (digamma(g4) - digamma(g3 + g4))
    
    idx <- rep(1:n_cells, times = diff(x@p))
    Xid <- split(x@i + 1, idx)
    Xcnt <- split(x@x, idx)
    
    # pbmcapply::pbmclapply
    batch_results <- pbmcapply::pbmclapply(
        1:n_cells, FUN = Estep2, Xid, Xcnt, n_topics, priors, g2, g3, g4, relation, 
        subtree_size, max_est, likelihood_eps, count_eps, tree_mat, level_penalty,
        mc.cores = parallel::detectCores() - 2
    )
    updates <- Mstep(Xid, n_topics, n_genes, nHDP_tree, batch_results)
    
    times <- N / n_cells
    if (learning_rate == 1) {
        nHDP_tree$size <- as.numeric(updates$W) * times
        nHDP_tree$centroids <- updates$B * times
    } else {
        for (k in 1:ncol(nHDP_tree$centroids)) {
            step_B <- times * updates$B[, k]
            # (1 - rho) * old + rho * ((1 - rho/10) * new + rho/10 * avg
            updates$B[, k] <- (1 - learning_rate) * nHDP_tree$centroids[, k] + learning_rate * ( (1 - learning_rate / 10) * step_B + learning_rate / 10 * mean(step_B))
        }
        
        nHDP_tree$size <- 
            (1 - learning_rate) * nHDP_tree$size + learning_rate * as.numeric(updates$W) * times
        
        nHDP_tree$centroids <- updates$B
    }
    
    return(list(new_tree = nHDP_tree, count_matrix = updates$M, likelihood = updates$L))
}


# ---- helper functions ----
Estep <- function(d, Xid, Xcnt, n_topics, priors, g2, g3, g4, relation,
                  subtree_size, max_est, likelihood_eps, count_eps, tree_mat, 
                  level_penalty) {
    d_id <- Xid[[d]]
    d_cnt <- Xcnt[[d]]
    
    this_ElnB <- priors$ElnB[d_id, ] # cell-specific ElnB
    this_ElnV <- digamma(1) - digamma(1 + g2) # local ElnV prior
    this_Eln1_V <- digamma(g2) - digamma(1 + g2) # local Eln(1-V) prior
    this_ElnP <- rep(-Inf, n_topics) # ElnP prior
    this_ElnP[1:5] <- this_ElnV + digamma(g3) - digamma(g3 + g4)
    
    # select subtree for this cell
    picked_nodes <- NULL
    lower_bound <- NULL
    this_weight <- array(0, dim = n_topics)
    keep_searching <- TRUE
    eps <- 2.2e-16
    # XX = sqrt(input$Xcnt[[i]])
    XX = log1p(1e4 * d_cnt / sum(d_cnt))
    
    while(keep_searching) {
        active_nodes <- which(this_ElnP > -Inf)
        
        #penalty <- t(t(this_ElnB[, active_nodes]) + this_ElnP[active_nodes])
        penalty <- this_ElnB[, active_nodes] + matrix(rep(this_ElnP[active_nodes], each = length(d_id)), nr = length(d_id), nc = length(active_nodes))
        C_act <- penalty
        penalty <- penalty * XX # gene-wise multiplication
        ElnPtop_act <- priors$ElnPtop[active_nodes]
        
        if (is.null(picked_nodes)) {
            score <- colSums(penalty) + ElnPtop_act
            picked_nodes <- active_nodes[which.max(score)]
            # add likelihood of best candidate
            lower_bound <- max(score) - ElnPtop_act[picked_nodes]
        } else {
            temp <- array(0, dim = n_topics)
            temp[active_nodes] <- 1:length(active_nodes)
            idx_clps <- temp[picked_nodes]
            n_active <- length(active_nodes)
            
            # gene-wise max scaled ElnB + ElnP
            if (length(idx_clps) == 1) {
                v <- penalty[, idx_clps]
            } else {
                v <- matrixStats::rowMaxs(penalty[, idx_clps])
            }
            
            # c(525, 1260, 1995, 2730, 3465, 4200, 4935, 5670, 6405)
            C_act <- C_act - v # gene-wise minus
            C_act <- exp(C_act)
            numerator <- C_act * penalty
            
            if (length(idx_clps) == 1) {
                numerator <- numerator + numerator[, idx_clps]
                denominator <- C_act + C_act[, idx_clps]
                v <- C_act[, idx_clps] * log(eps + C_act[, idx_clps])
            } else {
                numerator <- numerator + rowSums(numerator[, idx_clps])
                denominator <- C_act + rowSums(C_act[, idx_clps])
                v <- rowSums(C_act[, idx_clps] * log(eps + C_act[, idx_clps]))
            }
            
            h <- log(denominator) - (C_act * log(eps + C_act) + v) / denominator
            score <- colSums(numerator / denominator) + ElnPtop_act + colSums(h * XX)
            #score <- colSums(numerator / denominator) + ElnPtop_act + as.numeric(crossprod(h, XX))
            score[idx_clps] <- -Inf
            this_pick <- which.max(score)
            picked_nodes <- c(picked_nodes, active_nodes[which.max(score)])
            # add likelihood of best candidate
            lower_bound <- c(lower_bound, max(score) - ElnPtop_act[this_pick])
        }
        
        # update candidates according to new selected node
        last_pick <- picked_nodes[length(picked_nodes)]
        last_pick_parent <- relation[relation$child == last_pick, "parent"]
        sibling_nodes <- relation[relation$parent == last_pick_parent, "child"]
        left_sibling_nodes <- sibling_nodes[! sibling_nodes %in% intersect(sibling_nodes, picked_nodes)]
        
        # update weights of unselected siblings of last picked node
        this_weight[left_sibling_nodes] <- this_weight[left_sibling_nodes] + this_Eln1_V 
        this_ElnP[left_sibling_nodes] <- this_ElnP[left_sibling_nodes] + this_Eln1_V
        last_pick_children <- relation[relation$parent == last_pick, "child"]
        # update weights of children of last picked node
        this_weight[last_pick_children] <- this_ElnV 
        this_ElnP[last_pick_children] <- this_ElnV + level_penalty[last_pick_children]
        
        idx <- last_pick
        repeat {
            this_ElnP[last_pick_children] <- this_ElnP[last_pick_children] + this_weight[idx]
            idx <- relation[relation$child == idx, "parent"]
            # if searching back to root, stop
            if (idx == 0) {
                break
            }
        }
        
        # subtree size is between subtree_size[1] and subtree_size[2]
        # another condition is that the added likelihood is lower than 0.0005
        n_nodes <- length(lower_bound)
        if (n_nodes > subtree_size[1]) {
            if (n_nodes == subtree_size[2]) {
                keep_searching <- FALSE
            } else {
                added_likelihood <- lower_bound[n_nodes] - lower_bound[n_nodes - 1]
                if (abs(added_likelihood) / abs(lower_bound[n_nodes - 1]) < likelihood_eps) {
                    keep_searching <- FALSE
                }
            }
        }
    }
    # learn document parameters for subtree
    n_subtree <- length(picked_nodes)
    this_ElnB <- priors$ElnB[d_id, picked_nodes]
    this_ElnP <- this_ElnP[picked_nodes]
    cnt_old <- array(0, dim = n_subtree)
    keep_estimating <- TRUE
    estimate_time <- 0
    while(keep_estimating) {
        estimate_time <- estimate_time + 1
        #this_C <- t(t(this_ElnB) + as.numeric(this_ElnP))
        this_C <- this_ElnB + matrix(rep(this_ElnP, each = length(d_id)), nr = length(d_id), nc = length(this_ElnP))
        this_C <- this_C - matrixStats::rowMaxs(this_C)
        this_C <- exp(this_C)
        this_C <- this_C / rowSums(this_C)
        cnt <- colSums(this_C * XX)
        this_ElnP <- updateCellWeight(cnt, relation, picked_nodes, g2, g3, g4, tree_mat) 
        
        if (estimate_time > 1) {
            if (estimate_time == max_est) {
                keep_estimating <- FALSE
            } else {
                added_estimate <- sum(abs(cnt - cnt_old)) / sum(cnt)
                if (added_estimate < count_eps) {
                    keep_estimating <- FALSE    
                }
            }
        }
        cnt_old <- cnt
    }
    
    updates <- list(
        d_count_matrix = cnt,
        d_W = picked_nodes,
        d_B = this_C * XX
    )
    return(updates)    
}


Estep2 <- function(d, Xid, Xcnt, n_topics, priors, g2, g3, g4, relation,
                  subtree_size, max_est, likelihood_eps, count_eps, tree_mat, 
                  level_penalty) {
    d_id <- Xid[[d]]
    d_cnt <- Xcnt[[d]]
    
    this_ElnB <- priors$ElnB[d_id, ] # cell-specific ElnB
    this_ElnV <- digamma(1) - digamma(1 + g2) # local ElnV prior
    this_Eln1_V <- digamma(g2) - digamma(1 + g2) # local Eln(1-V) prior
    this_ElnP <- rep(-Inf, n_topics) # ElnP prior
    first_level_child <- relation[relation$parent == 0 & relation$child != 0, "child"]
    this_ElnP[1:max(first_level_child)] <- this_ElnV + digamma(g3) - digamma(g3 + g4)
    
    # select subtree for this cell
    picked_nodes <- NULL
    lower_bound <- NULL
    this_weight <- array(0, dim = n_topics)
    keep_searching <- TRUE
    eps <- 2.2e-16
    #XX = log2(d_cnt + 1)
    # XX = log1p(1e4 * d_cnt / sum(d_cnt))
    XX = d_cnt
    # XX = d_cnt / sum(d_cnt)
    
    while(keep_searching) {
        active_nodes <- which(this_ElnP > -Inf)
        penalty <- this_ElnB[, active_nodes][]
        penalty <- as.matrix(penalty)
        #penalty <- t(t(this_ElnB[, active_nodes]) + this_ElnP[active_nodes])
        
        collapse::setop(penalty, "+", this_ElnP[active_nodes], rowwise = TRUE)
        #penalty <- this_ElnB[, active_nodes] + matrix(rep(this_ElnP[active_nodes], each = length(d_id)), nr = length(d_id), nc = length(active_nodes))
        C_act <- penalty[]
        
        collapse::setop(penalty, "*", XX, rowwise = FALSE)
        #penalty <- penalty * XX # gene-wise multiplication
        
        ElnPtop_act <- priors$ElnPtop[active_nodes]
        
        if (is.null(picked_nodes)) {
            score <- colSums(penalty) + ElnPtop_act
            picked_nodes <- active_nodes[which.max(score)]
            # add likelihood of best candidate
            lower_bound <- max(score) - ElnPtop_act[picked_nodes]
        } else {
            temp <- array(0, dim = n_topics)
            temp[active_nodes] <- 1:length(active_nodes)
            idx_clps <- temp[picked_nodes]
            n_active <- length(active_nodes)
            
            # gene-wise max scaled ElnB + ElnP
            if (length(idx_clps) == 1) {
                v <- penalty[, idx_clps]
                collapse::setop(C_act, "-", v, rowwise = FALSE)
            } else {
                v <- matrixStats::rowMaxs(penalty[, idx_clps])
                collapse::setop(C_act, "-", v)
            }
            
            # c(525, 1260, 1995, 2730, 3465, 4200, 4935, 5670, 6405)
            #C_act <- C_act - v # gene-wise minus
            if (any(C_act > 700)) {
                C_act[C_act > 700] <- 700
            }
            
            C_act <- exp(C_act)
            
            numerator <- C_act[]
            denominator <- C_act[]
            
            collapse::setop(numerator, "*", penalty)
            #numerator <- C_act * penalty
            
            if (length(idx_clps) == 1) {
                #numerator <- numerator + numerator[, idx_clps]
                collapse::setop(numerator, "+", numerator[, idx_clps], rowwise = FALSE)
                # denominator <- C_act + C_act[, idx_clps]
                collapse::setop(denominator, "+", denominator[, idx_clps])
                v <- C_act[, idx_clps] * log(eps + C_act[, idx_clps])
            } else {
                # numerator <- numerator + rowSums(numerator[, idx_clps])
                collapse::setop(numerator, "+", rowSums(numerator[, idx_clps]), rowwise = FALSE)
                # denominator <- C_act + rowSums(C_act[, idx_clps])
                collapse::setop(denominator, "+", rowSums(denominator[, idx_clps]))
                v <- rowSums(C_act[, idx_clps] * log(eps + C_act[, idx_clps]))
            }
            
            #h <- log(denominator) - (C_act * log(eps + C_act) + v) / denominator
            h <- log(denominator)
            epslison <- C_act * log(eps + C_act) + v
            collapse::setop(epslison, "/", denominator)
            collapse::setop(h, "-", epslison)
            
            # score <- colSums(numerator / denominator) + ElnPtop_act + colSums(h * XX)
            collapse::setop(h, "*", XX, rowwise = FALSE)
            score <- collapse::setop(numerator, "/", denominator)
            score <- colSums(score) + ElnPtop_act + colSums(h)
            
            #score <- colSums(numerator / denominator) + ElnPtop_act + as.numeric(crossprod(h, XX))
            score[idx_clps] <- -Inf
            this_pick <- which.max(score)
            picked_nodes <- c(picked_nodes, active_nodes[which.max(score)])
            # add likelihood of best candidate
            lower_bound <- c(lower_bound, max(score) - ElnPtop_act[this_pick])
            
        }
        
        # update candidates according to new selected node
        last_pick <- picked_nodes[length(picked_nodes)]
        last_pick_parent <- relation[relation$child == last_pick, "parent"]
        sibling_nodes <- relation[relation$parent == last_pick_parent, "child"]
        left_sibling_nodes <- sibling_nodes[! sibling_nodes %in% intersect(sibling_nodes, picked_nodes)]
        
        # update weights of unselected siblings of last picked node
        this_weight[left_sibling_nodes] <- this_weight[left_sibling_nodes] + this_Eln1_V 
        this_ElnP[left_sibling_nodes] <- this_ElnP[left_sibling_nodes] + this_Eln1_V
        last_pick_children <- relation[relation$parent == last_pick, "child"]
        # update weights of children of last picked node
        this_weight[last_pick_children] <- this_ElnV 
        this_ElnP[last_pick_children] <- this_ElnV + level_penalty[last_pick_children]
        
        idx <- last_pick
        repeat {
            this_ElnP[last_pick_children] <- this_ElnP[last_pick_children] + this_weight[idx]
            idx <- relation[relation$child == idx, "parent"]
            # if searching back to root, stop
            if (idx == 0) {
                break
            }
        }
        
        # subtree size is between subtree_size[1] and subtree_size[2]
        # another condition is that the added likelihood is lower than 0.0005
        n_nodes <- length(lower_bound)
        if (n_nodes > subtree_size[1]) {
            if (n_nodes == subtree_size[2]) {
                keep_searching <- FALSE
            } else {
                added_likelihood <- lower_bound[n_nodes] - lower_bound[n_nodes - 1]
                if (abs(added_likelihood) / abs(lower_bound[n_nodes - 1]) < likelihood_eps) {
                    keep_searching <- FALSE
                }
            }
        }
    }
    # learn document parameters for subtree
    n_subtree <- length(picked_nodes)
    this_ElnB <- priors$ElnB[d_id, picked_nodes]
    this_ElnP <- this_ElnP[picked_nodes]
    cnt_old <- array(0, dim = n_subtree)
    keep_estimating <- TRUE
    estimate_time <- 0
    all_parents <- relation[match(picked_nodes, relation$child), "parent"]
    partition <- sort(unique(all_parents))
    while(keep_estimating) {
        estimate_time <- estimate_time + 1
        
        this_C <- this_ElnB[]
        collapse::setop(this_C, "+", this_ElnP, rowwise = TRUE)
        collapse::setop(this_C, "-", matrixStats::rowMaxs(this_C), rowwise = FALSE)
        this_C <- exp(this_C)
        collapse::setop(this_C, "/", rowSums(this_C), rowwise = FALSE)
        cnt_B <- this_C[]
        collapse::setop(this_C, "*", XX, rowwise = FALSE)
        cnt <- colSums(this_C)
        
        this_ElnP <- updateCellWeight(cnt, all_parents, partition, picked_nodes, g2, g3, g4, tree_mat) 
        
        if (estimate_time > 1) {
            if (estimate_time == max_est) {
                keep_estimating <- FALSE
            } else {
                added_estimate <- sum(abs(cnt - cnt_old)) / sum(cnt)
                if (added_estimate < count_eps) {
                    keep_estimating <- FALSE    
                }
            }
        }
        cnt_old <- cnt
    }
    
    updates <- list(
        d_count_matrix = cnt,
        d_W = picked_nodes,
        d_B = cnt_B * XX,
        d_LL = lower_bound[length(lower_bound)]
    )
    return(updates)    
}

Mstep <- function(Xid, n_topics, n_genes, nHDP_tree, batch_results) {
    new_count <- array(0, dim = n_topics)
    new_center <- matrix(0, nr = n_genes, nc = n_topics)
    count_matrix <- matrix(0, nr = n_topics, nc = length(Xid))
    batch_ll <- 0
    for (b in 1:length(Xid)) {
        picked_nodes <- batch_results[[b]]$d_W
        gene_idx <- Xid[[b]]
        new_count[picked_nodes] <- new_count[picked_nodes] + 1
        new_center[gene_idx, picked_nodes] <- new_center[gene_idx, picked_nodes] + batch_results[[b]]$d_B
        count_matrix[picked_nodes, b] <- batch_results[[b]]$d_count_matrix
        batch_ll <- batch_ll + batch_results[[b]]$d_LL
    }
    
    return(list(W = new_count, B = new_center, M = count_matrix, L = batch_ll))
}


#' Update the ElnP for cell
updateCellWeight <- function(cnt, all_parents, partition, picked_nodes, g2, g3, g4, tree_mat) {
    n_picked_nodes <- length(cnt)
    ElnP <- array(0, dim = n_picked_nodes)
    #all_parents <- relation[match(picked_nodes, relation$child), "parent"]
    
    bin_cnt1 <- cnt
    bin_cnt0 <- as.numeric(tree_mat[picked_nodes, picked_nodes] %*% cnt)
    Elnbin1 <- digamma(bin_cnt1 + g3) - digamma(bin_cnt1 + bin_cnt0 + g3 + g4)
    Elnbin0 <- digamma(bin_cnt0 + g4) - digamma(bin_cnt1 + bin_cnt0 + g3 + g4)
    
    # re-order weights
    stick_cnt <- bin_cnt1+bin_cnt0
    #partition <- sort(unique(all_parents))
    for (i in partition) {
        partition_idx <- which(all_parents == i)
        t1 <- stick_cnt[partition_idx]
        t1_sort <- sort(t1, decreasing = TRUE, index.return = TRUE)
        #t1_order <- order(t1, decreasing = TRUE)
        #t1 <- sort(t1, decreasing = TRUE)
        t1_order <- t1_sort$ix
        t1 <- t1_sort$x
        t3 <- t1 %>% rev() %>% cumsum() %>% rev()
        if (length(t3) > 1) {
            t4 <- c(t3[-1], 0)
            t5 <- c(0, digamma(t4[-length(t4)] + g2) - digamma(t1[-length(t1)] + t4[-length(t4)] + 1 + g2))
        } else {
            t4 <- 0
            t5 <- 0
        }
        weights <- digamma(t1 + 1) - digamma(t1 + t4 + 1 + g2) + cumsum(t5)
        ElnP[partition_idx[t1_order]] <- weights
    }
    ElnP <- ElnP + Elnbin1 + as.numeric(crossprod(tree_mat[picked_nodes, picked_nodes], (Elnbin0 + ElnP)))
    return(ElnP)
}

#' Calculate the log likelihood of the document given the parameters
#' 
#' @param x dgCMatrix of the used single cell data
#' @param Beta W x K matrix
#' @param Theta K x D matrix
#'
#' @return predictive log likelihood of all documents
calcLogLik <- function(X, Beta, Theta) {
    wordDist <- apply(Beta, 2, function(x) x/sum(x))
    wordDist[is.na(wordDist)] <- 0
    topicDist <- apply(Theta, 2, function(x) x/sum(x))
    dw <- crossprod(topicDist, t(wordDist))
    idx <- rep(1:ncol(X), times = diff(X@p))
    Xid <- split(X@i + 1, idx)
    Xcnt <- split(X@x, idx)
    Xcnt <- lapply(Xcnt, function(x) log2(x + 1))
    likelihood <- sapply(1:ncol(X), function(x) {
        w <- log(dw[x, Xid[[x]]] * Xcnt[[x]])
        w <- w[! is.infinite(w)]
        sum(w)
    })
    return(sum(likelihood))
}

#' Mini-batch kmeans for dgCMatrix data, clusters that are less than a 
#' required number of observations will be removed, the samples associated
#' with those clusters will be reassigned to the nearest cluster based on 
#' L1 distance. Final clusters are sorted based on sample size.
#' 
#' @param x dgCMatrix matrix
#' @param min_n minimum number of observations in each cluster
KMeans <- function(x, min_n, ...) {
    
    mbk <- suppressWarnings(
        mbkmeans::mbkmeans(x, ...)
    )
    centroids <- t(mbk$centroids)
    clusters <- mbk$Clusters
    
    # remove clusters that do not have enough sample size
    ck <- table(clusters) < min_n
    n_valid_cluster <- length(which(! ck))
    
    # make sure we can have more than two clusters
    if (n_valid_cluster <= 1) return(NULL)
    
    if (any(ck)) {
        rm_clusters <- which(ck)
        reassign_idx <- which(clusters %in% rm_clusters)
        
        reassign_clusters <- apply(
            as.matrix(x[, reassign_idx]), 
            MARGIN = 2, 
            FUN = function(g) {
                # distance <- colSums(abs(g - centroids))
                distance <- colSums((g - centroids) ^ 2)
                order(distance)[2]
            }
        )
        
        clusters[reassign_idx] <- reassign_clusters
    }
    
    # reorder cluster id based on sample size
    cnt <- sort(table(clusters), decreasing = TRUE)
    od <- names(cnt) %>% as.numeric()
    group_reorder <- match(clusters, od)
    
    output <- list(centroids = centroids[, od], 
                   clusters  = group_reorder,
                   K         = length(od))
    return(output)
}

#' KMeans clustering with L1 norm, the final group will be renamed based on 
#' decreasing order of the number of samples within the K clusters
#' 
#' @param mat dgCMatrix
#' @param K number of clusters
#' @param iter_max number of iteration (default is 3)
#' 
#' @import parallel
#' 
#' @return list of two attributes: centers (p x K) and groups (1 x n)
KMeans_L1 <- function(mat, K, iter_max = 3) {
    n <- ncol(mat)
    n_cores <- parallel::detectCores() - 1 # prevent overhead
    batches <- split(seq.int(n), ceiling(seq.int(n) / 20000))
    nb <- length(batches)
    
    if (n >= K) {
        # use the first K samples as the initial center
        b <- mat[, 1:K] 
    } else {
        # use the first n samples as the initial center, the rest be 0
        b <- mat[, 1:n]
        tmp <- matrix(0, nr = nrow(mat), nc = K - n) %>%
            Matrix::Matrix(sparse = TRUE)
        b <- cbind(b, tmp) %>% as.matrix() %>% unname()
        
        output <- list(centers = b,
                       groups  = 1:n)
        return(output)
    }
    
    # inner function for L1 calculation
    L1 <- function(i) {
        abs(b - mat[, i]) %>% 
            Matrix::colSums() %>%
            which.min()
    }
    
    # depracated
    L1_fast0 <- function() {
        l1_dist <- parallel::mclapply(b_lst, FUN = function(k) {
            Matrix::colSums(abs(mat - k))
        }, mc.cores = K)
        group <- do.call(rbind, l1_dist) %>% apply(2, which.min)
        return(group)
    }
    
    L1_fast <- function() {
        group <- rep(0, dim = n)
        for (b in 1:nb) {
            l1_dist <- parallel::mclapply(b_lst, FUN = function(k) {
                Matrix::colSums(abs(mat[, batches[[b]]] - k))
            }, mc.cores = K)
            group_b <- do.call(rbind, l1_dist) %>% apply(2, which.min)
            group[batches[[b]]] <- group_b
        }
        return(group)
    }
    
    # inner function for averaged center
    avgCenter <- function(idx) {
        if (length(idx) > 1) {
            Matrix::rowMeans(mat[, idx])
        } else if (length(idx) == 1) {
            mat[, idx]
        } else {
            matrix(0, nr = nrow(mat), nc = 1) %>%
                Matrix::Matrix(sparse = TRUE)
        }
    }
    
    for (i in 1:iter_max) {
        b_lst <- list()
        for (k in 1:K) b_lst[[k]] <- as.numeric(b[, k])
        
        # E-step: L1 KMeans for each sample and assign to closet center
        #group <- parallel::mclapply(1:n, L1, mc.cores = n_cores) %>%
        #    unlist()
        group <- L1_fast()
        
        # v0.1: NEW FEATURE
        # prevent the K in K-means decreases
        for (k in 1:K) {
            if (length(which(group == k)) == 0) {
                group[k] <- k
            }
        }
        
        group_lst <- lapply(1:K, function(x) which(group == x))
        
        # M-step: averaged center
        b <- parallel::mclapply(group_lst, avgCenter, mc.cores = n_cores) %>%
            do.call(cbind, .)
    }
    
    # reorder center by the number of samples within each cluster
    cnt <- sort(table(group), decreasing = TRUE)
    od <- names(cnt) %>% as.numeric()
    group_reorder <- match(group, od)
    b <- b[, od]
    
    output <- list(centroids = b, 
                   clusters  = group_reorder,
                   K         = length(od))
    return(output)
    # remove clusters that do not have enough sample size
    #ck <- table(group_reorder) < min_n
    #n_valid_cluster <- length(which(! ck))
    
    # make sure we can have more than two clusters
    #if (n_valid_cluster <= 1) return(NULL)
    
    #if (any(ck)) {
    #    rm_clusters <- which(ck)
    #    reassign_idx <- which(group_reorder %in% rm_clusters)
    
    #    reassign_clusters <- apply(
    #        as.matrix(x[, reassign_idx]), 
    #        MARGIN = 2, 
    #        FUN = function(g) {
    # distance <- colSums(abs(g - centroids))
    #            distance <- abs(g - b)
    #            order(distance)[2]
    #        }
    #    )
    
    #    group_reorder[reassign_idx] <- reassign_clusters
    #}
}

#' Normalize data as LogNormal in single cell (scale to 1e4 then take log)
#' 
#' @param x dgCMatrix format
#' 
#' @importFrom Matrix colSums
#'
#' @return log-normalized data in dgCMatrix
normData <- function(x) {
    if (all(x@x == 0)) {
        return(x)
    } else {
        # x@x <- log1p(x@x / rep.int(Matrix::colSums(x), diff(x@p)) * 1e4)
        x@x <- x@x / rep.int(Matrix::colSums(x), diff(x@p))
    }
    return(x)
}
#' Process the global tree and get global variables as prior
#' 
#' @param nHDP_tree nHDP tree.
#' @param beta0 the add-up when calculating the digamma function of count.
#' @param gamma1 top level Dirichlet Process concentration
#' 
#' @return list of global tree information
#' (ElnB: Elog(theta) (n_gene x n_topic)
#'  ElnPtop: ElnP (1 x n_topic))
getTreePrior <- function(
        nHDP_tree, 
        beta0, 
        gamma1) {
    n_topics <- length(nHDP_tree$size)
    tree <- nHDP_tree$tree
    ElnB <- apply(nHDP_tree$centroids, 2, function(x) {
        digamma(x + beta0) - digamma(sum(x + beta0))
    })
    
    # Elogp (global Elogpi), n_topics x 1 
    ElnPtop <- array(dim = n_topics) 
    # for each parent
    for (p in unique(tree$parent)) {
        child_nodes <- tree[tree$parent == p, "child"]
        child_cnt <- nHDP_tree$size[child_nodes] 
        tau1 <- child_cnt + 1
        tau2 <- rev(cumsum(rev(child_cnt[-1])))
        tau2 <- c(tau2, 0) + gamma1
        ElnV <- digamma(tau1) - digamma(tau1 + tau2)
        Eln1_V <- digamma(tau2) - digamma(tau1 + tau2)
        ElnPtop[child_nodes] <- ElnV + c(0, cumsum(Eln1_V[-length(Eln1_V)]))
    }
    
    result <- list(ElnB    = ElnB,
                   ElnPtop = ElnPtop)
    return(result)
}