# Some helper functions for working with the ape package

# In ape's world, there are edges that connect both tips and (internal-only)
# nodes.  The tips are numbered with the lower half plus one of the node/tip
# index values.  There is one more tip+node than edges.
# Every edge has a clade (a set of named tips) associated with it, and every
# entry in the vector of inferred sequences does too.  If we match up the clades
# we can match up edges and nodes.  There is one inferred sequence missing
# compared with the tree object: the node (and clade) at the root itself, before
# the sequence the tree is rooted on.

#' Get tips under each node
#'
#' Get a vector of tree tip labels for the clade associated with each tree node,
#' with all the vectors in a list.  The lower half of the list is for the tips
#' themselves, with the upper half the internal nodes in the tree.  This way the
#' positions in the list match up with that in the edge matrix for the tree.
#'
#' @param tree phylogenetic tree object from \code{\link[ape]{read.tree}}.
#' @return list of character vectors, one for each node (including tips themselves)
#' @export
tree_get_clades <- function(tree) {
  node_idxs <- 1:max(tree$edge)
  lapply(node_idxs, function(node) sort(tree$tip.label[tree_gather_tips(tree$edge, node)]))
}

#' Get tip indexes for a given node
#'
#' @return numeric vector of all tip index values linked from a particular node
#'   index
tree_gather_tips <- function(edges, node) {
  # what nodes does this node have an edge toward?
  edge_idx <- which(edges[, 1] == node)
  # if the given node links to any others, gather tips for those
  # otherwise return this node itself
  tips <- integer()
  if (length(edge_idx)) {
    for (next_node in edges[edge_idx, 2]) {
      tips <- c(tips, tree_gather_tips(edges, next_node))
    }
  } else {
    tips <- node
  }
  return(tips)
}

# Find all node index values that point towards the given tip
# this will always include node 1 since that's at the root.
tree_get_ancestors <- function(tree, tip_name) {
  which(sapply(tree_get_clades(tree), function(vec) tip_name %in% vec))
}
