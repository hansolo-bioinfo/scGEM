#' Mini-batch kmeans for dgCMatrix data, clusters that are less than a
#' required number of observations will be removed, the samples associated
#' with those clusters will be reassigned to the nearest cluster based on
#' L1 distance. Final clusters are sorted based on sample size.
#'
#' @param x dgCMatrix matrix
#' @param min_n minimum number of observations in each cluster
#'
#' @import mbkmeans
#.KMeans <- function(x, min_n, ...) {
#
#    mbk <- suppressWarnings(
#        mbkmeans::mbkmeans(x, ...)
#    )
#    centroids <- t(mbk$centroids)
#    clusters <- mbk$Clusters
#
#    # remove clusters that do not have enough sample size
#    ck <- table(clusters) < min_n
#    n_valid_cluster <- length(which(! ck))
#
#    # make sure we can have more than two clusters
#    if (n_valid_cluster <= 1) return(NULL)

#    if (any(ck)) {
#        rm_clusters <- which(ck)
#        reassign_idx <- which(clusters %in% rm_clusters)
#
#        reassign_clusters <- apply(
#            as.matrix(x[, reassign_idx]),
#            MARGIN = 2,
#            FUN = function(g) {
#                # distance <- colSums(abs(g - centroids))
#                distance <- colSums((g - centroids) ^ 2)
#                order(distance)[2]
#            }
#        )
#
#        clusters[reassign_idx] <- reassign_clusters
#    }
#
#    # reorder cluster id based on sample size
#    cnt <- sort(table(clusters), decreasing = TRUE)
#    od <- names(cnt) %>% as.numeric()
#    group_reorder <- match(clusters, od)
#
#    output <- list(centroids = centroids[, od],
#                   clusters  = group_reorder,
#                   K         = length(od))
#    return(output)
#}

#' KMeans clustering with L1 norm, the final group will be renamed based on
#' decreasing order of the number of samples within the K clusters
#'
#' @param mat dgCMatrix
#' @param K number of clusters
#' @param iter_max number of iteration (default is 3)
#'
#' @import parallel Matrix
#'
#' @return list of two attributes: centers (p x K) and groups (1 x n)
.KMeans_L1 <- function(mat, K, iter_max = 3) {
    n <- ncol(mat)
    n_cores <- parallel::detectCores() - 1 # prevent overhead
    batches <- split(seq.int(n), ceiling(seq.int(n) / 2500))
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
    cnt <- sort(table(unlist(group)), decreasing = TRUE)
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
#' @import Matrix
#'
#' @return log-normalized data in dgCMatrix
.normData <- function(x) {
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
.getTreePrior <- function(
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
