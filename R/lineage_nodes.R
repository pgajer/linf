# A node key encodes fitted column indices, with zero reserved for a synthetic
# rare component. It never parses a feature ID or a display label.
linf.path.key <- function(path) paste(path, collapse = "/")

linf.unique.node.names <- function(raw, keys) {
  out <- raw
  collision <- duplicated(raw) | duplicated(raw, fromLast = TRUE)
  for (i in which(collision)) {
    candidate <- paste0(raw[[i]], " [node ", keys[[i]], "]")
    while (candidate %in% c(raw, out[-i])) candidate <- paste0(candidate, "'")
    out[[i]] <- candidate
  }
  out
}

linf.node.level <- function(paths, previous, feature.ids, feature.labels, rare.label, sep) {
  keys <- vapply(paths, function(path) {
    if (!length(path)) NA_character_ else linf.path.key(path)
  }, character(1))
  unique.keys <- sort(unique(keys[!is.na(keys)]))
  members <- match(unique.keys, keys)
  unique.paths <- paths[members]
  render <- function(path, features) paste(ifelse(path == 0L, rare.label,
                                                 features[pmax(path, 1L)]), collapse = sep)
  ids <- vapply(unique.paths, render, character(1), features = feature.ids)
  labels <- vapply(unique.paths, render, character(1), features = feature.labels)
  leaf <- vapply(unique.paths, function(path) utils::tail(path, 1L), integer(1))
  parent <- if (is.null(previous)) rep(NA_character_, length(members)) else unname(previous[members])
  nodes <- data.frame(
    node.id = unique.keys, parent.node.id = parent,
    feature.index = ifelse(leaf == 0L, NA_integer_, leaf),
    is.rare = leaf == 0L,
    terminal = !is.na(parent) & unique.keys == parent,
    lineage.id = linf.unique.node.names(ids, unique.keys),
    lineage.label = linf.unique.node.names(labels, unique.keys),
    stringsAsFactors = FALSE
  )
  # A carried node may also be reached from a newly refined parent in a stored
  # policy view. Keep every realized parent instead of selecting the first row.
  nodes$parent.node.ids <- I(lapply(unique.keys, function(key) {
    if (is.null(previous)) character() else unique(unname(previous[which(keys == key)]))
  }))
  nodes$parent.node.id <- vapply(nodes$parent.node.ids, function(x) {
    if (length(x) == 1L) x else NA_character_
  }, character(1))
  nodes$terminal <- vapply(seq_along(unique.keys), function(i) {
    parents <- nodes$parent.node.ids[[i]]
    length(parents) == 1L && !is.na(parents) && identical(parents, unique.keys[[i]])
  }, logical(1))
  nodes$path <- I(unname(unique.paths))
  list(keys = keys, nodes = nodes)
}

linf.store.node.level <- function(csts, view, paths, depth, sep = "__") {
  key.field <- paste0("node.ids.", view)
  node.field <- paste0("nodes.", view)
  previous <- if (depth == 1L) NULL else csts[[key.field]][[depth - 1L]]
  level <- linf.node.level(paths, previous, csts$feature.ids, csts$feature.labels,
                           csts$rare.label, sep)
  names(level$keys) <- names(csts$lineage.id)
  csts[[key.field]][[depth]] <- level$keys
  csts[[node.field]][[depth]] <- level$nodes
  idx <- match(level$keys, level$nodes$node.id)
  for (kind in c("id", "label")) {
    values <- level$nodes[[paste0("lineage.", kind)]][idx]
    names(values) <- names(level$keys)
    csts[[paste0("lineage.", kind, "s.", view)]][[depth]] <- values
    csts[[paste0("lineage.", kind, ".", view)]] <- values
  }
  csts
}

linf.activate.nodes <- function(csts) {
  view <- csts$low.freq.policy
  for (field in c("nodes", "node.ids", "lineage.ids", "lineage.labels", "lineage.id", "lineage.label")) {
    csts[[field]] <- csts[[paste0(field, ".", view)]]
  }
  csts
}

linf.node.paths <- function(csts, view, depth) {
  nodes <- csts[[paste0("nodes.", view)]][[depth]]
  keys <- csts[[paste0("node.ids.", view)]][[depth]]
  lapply(match(keys, nodes$node.id), function(i) {
    if (is.na(i)) integer() else nodes$path[[i]]
  })
}

# Import unambiguous pre-node fits on consumption. A serialized label shared by
# different paths is not recoverable safely and must not be silently merged.
linf.ensure.nodes <- function(csts) {
  fields <- c("nodes.pure", "nodes.absorb", "node.ids.pure", "node.ids.absorb")
  present <- fields %in% names(csts)
  if (any(present)) {
    if (!all(present)) stop("incomplete dCST node metadata; rebuild the fit", call. = FALSE)
    for (view in c("pure", "absorb")) {
      nodes <- csts[[paste0("nodes.", view)]]
      keys <- csts[[paste0("node.ids.", view)]]
      if (length(nodes) != csts$depth || length(keys) != csts$depth ||
          any(lengths(keys) != length(csts$lineage.id))) {
        stop("invalid dCST node hierarchy; rebuild the fit", call. = FALSE)
      }
      for (d in seq_len(csts$depth)) {
        node <- nodes[[d]]
        required <- c("node.id", "parent.node.id", "parent.node.ids", "feature.index",
                      "is.rare", "terminal", "lineage.id", "lineage.label", "path")
        if (!is.data.frame(node) || !all(required %in% names(node)) ||
            anyNA(node$node.id) || anyDuplicated(node$node.id) ||
            any(!is.na(keys[[d]]) & !keys[[d]] %in% node$node.id)) {
          stop("invalid dCST node metadata; rebuild the fit", call. = FALSE)
        }
        for (kind in c("id", "label")) {
          expected <- node[[paste0("lineage.", kind)]][match(keys[[d]], node$node.id)]
          actual <- csts[[paste0("lineage.", kind, "s.", view)]][[d]]
          if (!identical(unname(actual), unname(expected))) {
            stop("dCST assignments disagree with node metadata; rebuild the fit", call. = FALSE)
          }
        }
      }
    }
    return(csts)
  }
  fail <- function() stop("ambiguous legacy dCST hierarchy; rebuild it from its matrix with the current linf", call. = FALSE)
  if (csts$rare.label %in% c(csts$feature.ids, csts$feature.labels)) fail()
  sep <- csts$sep %||% "__"
  for (view in c("pure", "absorb")) {
    levels <- csts[[paste0("lineage.ids.", view)]]
    labels <- csts[[paste0("lineage.labels.", view)]]
    paths <- rep(list(integer()), length(csts$lineage.id))
    for (d in seq_len(csts$depth)) {
      for (i in seq_along(paths)) {
        child <- levels[[d]][[i]]
        if (is.na(child)) { paths[[i]] <- integer(); next }
        leaf <- child
        if (d > 1L) {
          parent <- levels[[d - 1L]][[i]]
          if (identical(child, parent)) next
          if (is.na(parent) || !startsWith(child, paste0(parent, sep))) fail()
          leaf <- substring(child, nchar(paste0(parent, sep)) + 1L)
        }
        index <- if (identical(leaf, csts$rare.label)) 0L else match(leaf, csts$feature.ids)
        if (is.na(index)) fail()
        paths[[i]] <- c(paths[[i]], as.integer(index))
      }
      keys <- vapply(paths, linf.path.key, character(1))
      for (values in list(levels[[d]], labels[[d]])) {
        groups <- split(keys, values)
        if (any(vapply(groups, function(x) length(unique(x)) > 1L, logical(1)))) fail()
      }
      csts <- linf.store.node.level(csts, view, paths, d, sep)
    }
  }
  linf.activate.nodes(csts)
}
