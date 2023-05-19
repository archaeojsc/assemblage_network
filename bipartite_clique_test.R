library(igraph)  
set.seed(8675309)  
g <- graph_from_edgelist(matrix(sample(LETTERS[1:10], 50, replace=T), ncol = 2), directed = FALSE)  
plot(g, edge.arrow.size=0.5)

cliques <- max_cliques(g)

cliqueBP <- matrix(c(rep(paste0("cl", seq_along(cliques)), sapply(cliques, length)), names(unlist(cliques))), ncol=2, )
bp <- graph_from_edgelist(cliqueBP, directed = F)
V(bp)$type <- grepl("cl", V(bp)$name)
plot(bp, layout=layout_as_bipartite)

bp.ind <- t(as_incidence_matrix(bp))
bp.adj <- bp.ind %*% t(bp.ind)

bp.adj.g <- graph_from_adjacency_matrix(bp.adj, mode = "undirected")
plot(simplify(bp.adj.g))
bp.adj.mis <- independent.vertex.sets(bp.adj.g)

sets <- lapply(bp.adj.mis, function(x) cliqueBP[cliqueBP[,1] %in% as_ids(x), 2])
sets[which(sapply(sets, length) == max(sapply(sets, length)))]
