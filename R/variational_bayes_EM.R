#' Expectation step for variational Bayesian inference in scGEM
#'
#' @import collapse matrixStats
.Estep2 <- function(d, Xid, Xcnt, n_topics, priors, g2, g3, g4, relation,
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
            if (any(C_act > 600)) {
                C_act[C_act > 600] <- 600
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

            if (all(is.infinite(score))) {
                break
            }

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
                if (added_likelihood / abs(lower_bound[n_nodes - 1]) < likelihood_eps) {
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

        this_ElnP <- .updateCellWeight(cnt, all_parents, partition, picked_nodes, g2, g3, g4, tree_mat)

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

.Mstep <- function(Xid, n_topics, n_genes, nHDP_tree, batch_results) {
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
.updateCellWeight <- function(cnt, all_parents, partition, picked_nodes, g2, g3, g4, tree_mat) {
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
