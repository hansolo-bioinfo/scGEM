#' initialize the global tree based on nested mini-batch kmeans
#'
#' @param x dgCMatrix of single cell count data (gene by cell).
#' @param num_topics vector of children topics per parent in each level.
#' (example: c(5, 4, 3))
#' @param min_cells minimum number of cells required for each topic,
#' topics fail to meet this requirement will be removed.
#' (default is min(num_topics))
#' @param blacklist_genes an expression of set of genes to be removed.
#' (default is NULL, example: "^MT|^RPS")
#' @param rm.housekeeping_genes whether or not to remove housekeeping genes.
#' (default is TRUE)
#'
#' @import mbkmeans Matrix
#'
#' @return centroids and counts of topics and parent-children relationship
#' @export
#'
#' @example
#' scGEM_init <- initTree(
#'   x,
#'   num_topics = c(5, 4, 3),
#'   blacklist_genes = "^MT|^RP|^HSP|B2M|MALAT1",
#'   rm.housekeeping_genes = TRUE
#' )
initTree <- function(
        x,
        num_topics,
        min_cells = NULL,
        blacklist_genes = NULL,
        rm.housekeeping_genes = TRUE
) {

    if (! is(x, "dgCMatrix")) {
        stop("x is not in dgCMatrix format")
    }

    if (is.null(min_cells)) {
        min_cells <- min(num_topics)
    }

    if (! is.null(blacklist_genes)) {
        rm_genes <- grep(blacklist_genes, rownames(x))
        if (length(rm_genes) > 0) {
            x <- x[-rm_genes, ]
        }
    }

    if (rm.housekeeping_genes) {
        data("Housekeeping_GenesHuman")
        hk <- as.character(Housekeeping_GenesHuman$Gene.name)
        rm_hk <- which(rownames(x) %in% hk)
        if (length(rm_hk) > 0) {
            x <- x[-rm_hk, ]
        }
    }

    # normalize
    x <- .normData(x)
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

            if (length(idx_s) < min_cells | length(idx_s) < K) {
                next
            }

            if (K == 1) {
                mbk <- list()
                mbk$centroids <- Matrix::rowMeans(x_s)
                mbk$clusters <- rep(1, length(idx_s))
                mbk$K <- 1
            } else {
                mbk <- .KMeans_L1(x_s, K, 3)
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

            # break information and renormalize
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
                    x_k <- .normData(x_k)

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

#' Infer scGEM tree based on mini-batch variational inference
#'
#' @param X dgCMatrix single cell data to learn
#' @param scGEM_init initialized global scGEM tree.
#' @param batch_size sample size for one mini-batch training.
#' @param model_paras list of hyperparameters in the model.
#' b0: Dirichlet hyperparameter
#' g1: global Dirichlet process concentration
#' g2: local Dirichlet process concentration
#' g3: switching probability (stay)
#' g4: switching probability (searching down)
#' @param n_epoch number of epochs (default is 20).
#' @param subtree_size size of the search tree for each cell.
#' (example: c(5, 20), at least 5, at most 20)
#' @param max_est maximum number of estimation for topic percentages in
#' each cell. (default is 25)
#' @param likelihood_eps stopping criteria for added likelihood in lower bound
#' @param count_eps stopping criteria for theta estimation in posterior
#' @param learning_rate learning rate of each batch (use ite ^ (-2) if NULL)
#' @param tol stopping criteria for learning, (LL(n)-LL(n-1))/LL(n-1)
#'
#' @import Matrix pbmcapply collapse matrixStats
#'
#' @return GEM distribution over cells, Gene distribution over GEM
#' @example
#' scGEM_trained_mb <- minibatchInfer(
#'   x,
#'   scGEM_init,
#'   model_paras,
#'   max_est = 50,
#'   subtree_size = c(1, 20),
#'   batch_size = 1500,
#'   n_epoch = 50,
#'   learning_rate = 0.01
#' )
minibatchInfer <- function(
        X,
        scGEM_init,
        model_paras,
        batch_size = 1000,
        n_epoch = 20,
        subtree_size = c(1, 20),
        max_est = 50,
        likelihood_eps = .5*10^-3,
        count_eps = .5*10^-2,
        learning_rate = NULL,
        tol = 1e-4){

    miss_gene <- scGEM_init$gene[! scGEM_init$gene %in% rownames(X)]
    if (length(miss_gene) > 0) {
        X_add <- Matrix::Matrix(0, nrow = length(miss_gene), ncol = ncol(X),
                                sparse = TRUE)
        rownames(X_add) <- miss_gene
        colnames(X_add) <- colnames(X)
        X <- Matrix::rbind2(X, X_add)
    }
    matched_genes <- match(scGEM_init$gene, rownames(X))
    #matched_genes <- match(scGEM_init$gene, rownames(X))
    X <- X[matched_genes, ]
    scGEM_tree <- scGEM_init
    scGEM_tree$cell <- colnames(X)
    N <- length(scGEM_tree$cell)
    scGEM_tree$batch_size <- batch_size

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
        count_matrix <- matrix(0, nr = length(scGEM_tree$size), nc = N)

        for (ite in 1:n_ite) {
            cat(" [", ite, "/", n_ite, "]", sep = "")
            ite_num <- ite_num + 1

            if (is.null(learning_rate)) {
                step_size <- ite_num ^ (-2)

                if (step_size < 0.0001) {
                    step_size <- 0.0001
                }
            } else {
                step_size <- learning_rate
            }
            cat(", Step: ", step_size, "\n", sep = "")

            d_idx <- batches[[ite]]
            X_mb <- X[, d_idx]

            batch_updates <- .updateTree(X_mb,
                                        scGEM_tree,
                                        model_paras,
                                        subtree_size,
                                        max_est,
                                        likelihood_eps,
                                        count_eps,
                                        step_size)
            scGEM_tree <- batch_updates$new_tree
            count_matrix[, d_idx] <- batch_updates$count_matrix
            likelihood[epoch] <- likelihood[epoch] + batch_updates$likelihood
        }

        scGEM_tree$count_matrix <- count_matrix
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

    scGEM_tree$likelihood <- likelihood
    return(scGEM_tree)
}

#' Update scGEM tree
#'
#' @import pbmcapply collapse matrixStats
.updateTree <- function(
        x,
        scGEM_tree,
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

    n_genes <- length(scGEM_tree$gene)
    n_cells <- ncol(x)
    n_topics <- length(scGEM_tree$size)
    N <- length(scGEM_tree$cell)

    # process the global tree and get the priors
    priors <- .getTreePrior(scGEM_tree, b0, g1)
    relation <- scGEM_tree$tree
    tree_mat <- scGEM_tree$tree_mat

    parent_level <- colSums(tree_mat)
    level_penalty <- digamma(g3) - digamma(g3 + g4) +
        parent_level * (digamma(g4) - digamma(g3 + g4))

    idx <- rep(1:n_cells, times = diff(x@p))
    Xid <- split(x@i + 1, idx)
    Xcnt <- split(x@x, idx)

    # pbmcapply::pbmclapply
    batch_results <- pbmcapply::pbmclapply(
        1:n_cells, FUN = .Estep2, Xid, Xcnt, n_topics, priors, g2, g3, g4, relation,
        subtree_size, max_est, likelihood_eps, count_eps, tree_mat, level_penalty,
        mc.cores = parallel::detectCores() - 2
    )
    updates <- .Mstep(Xid, n_topics, n_genes, scGEM_tree, batch_results)

    times <- N / n_cells
    if (learning_rate == 1) {
        scGEM_tree$size <- as.numeric(updates$W) * times
        scGEM_tree$centroids <- updates$B * times
    } else {
        for (k in 1:ncol(scGEM_tree$centroids)) {
            step_B <- times * updates$B[, k]
            updates$B[, k] <- (1 - learning_rate) * scGEM_tree$centroids[, k] + learning_rate * ( (1 - learning_rate / 10) * step_B + learning_rate / 10 * mean(step_B))
        }

        scGEM_tree$size <-
            (1 - learning_rate) * scGEM_tree$size + learning_rate * as.numeric(updates$W) * times

        scGEM_tree$centroids <- updates$B
    }

    return(list(new_tree = scGEM_tree, count_matrix = updates$M, likelihood = updates$L))
}
