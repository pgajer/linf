#' Triple-Supported Graph of Unordered Pair States
#'
#' Construct a graph whose vertices are pure level-two unordered dominance
#' states. Edges join overlapping pairs, optionally requiring both pair
#' configurations to occur within their unique common leading triple.
#'
#' @param X Samples-by-features input accepted by [linf.udcsts()], with at
#'   least three features. Supply either `X` or both fitted assignments below.
#' @param pair.fit,triple.fit Results of [linf.udcsts()] at levels two and
#'   three, respectively, fitted to the same rows and dictionary. Pure IDs
#'   are used regardless of the fitted rare-state policy. Row identities,
#'   ordered ranks and dictionary agreement are checked. Anonymous fitted
#'   rows are matched by position; supply stable row names when fitting.
#' @param sample.ids Unique, nonempty sample IDs in input order. Defaults to
#'   fitted row names, or character row numbers for unnamed inputs.
#' @param coverage Optional target percentage in `(0, 100]`. Retain complete
#'   frequency tiers until at least this coverage is reached. Mutually
#'   exclusive with `states` and a nondefault `n0`.
#' @param n0 Minimum pair-state frequency, default one. Applied after counts
#'   on all input rows have been computed.
#' @param states Optional character vector of observed pair IDs to retain.
#'   `character(0)` selects an empty graph. Mutually exclusive with `coverage`
#'   and a nondefault `n0`.
#' @param adjacency `"observed_triple"` (default) or `"shared_feature"`.
#' @param min.support Positive integer minimum incidence on each side of an
#'   observed-triple edge. Also sets the threshold for witnessed faces.
#' @param feature.ids,feature.labels,backend Passed to [linf.udcsts()] when
#'   `X` is supplied. For existing fits, their feature metadata are used.
#'
#' @details
#' For pairs A and B sharing exactly one feature, the unique witness is
#' T = union(A, B). An observed-triple edge requires both incidences
#' I(A,T) and I(B,T) to be at least `min.support`. All marginals and
#' similarities refer to the full input sample set, including excluded pair
#' states. Selecting vertices never renormalizes these counts.
#'
#' All graph edges have `length = 1`. Jaccard and support columns are
#' evidence/display attributes, not path lengths. A witnessed face requires
#' all three pairs of one triple to be selected and pass the incidence
#' threshold. A triangle with three different witnesses is not a face.
#' Isolated vertices and disconnected components are preserved. No MST
#' repair or embedding is performed. Pure IDs use column indices and must
#' travel with the returned feature dictionary; realign columns to that
#' dictionary before fitting changed input column orders.
#'
#' @return A list of class `linf_udcst_graph` with `nodes` (pair ID, label,
#'   frequency, feature indices, component), `edges` (admitted edges),
#'   `candidate.edges` (all shared-feature edges on retained nodes), `faces`,
#'   full nonzero `incidence`, `pair.frequencies`, `triple.frequencies`,
#'   `membership` (sample ID and pure pair/triple IDs), `dictionary`, and
#'   `metadata`. Edges use integer `from`/`to` indices into `nodes` and retain
#'   pair IDs, witness ID, full marginal/incidence counts, minimum support,
#'   endpoint Jaccards, their mean/minimum/geometric mean, pair-conditional
#'   fractions and triple-allocation fractions. Faces use node indices `a`,
#'   `b`, `c` and retain their witness and minimum incidence. Metadata
#'   records settings, sample coverage, component sizes and isolate count.
#'
#' @examples
#' X <- rbind(c(A=6, B=3, C=1), c(A=6, B=1, C=3), c(A=1, B=6, C=3))
#' g <- linf.udcst.graph(X, coverage=60)
#' g$nodes
#' g$edges
#' g$faces
#' @export
linf.udcst.graph <- function(X = NULL, pair.fit = NULL, triple.fit = NULL,
                             sample.ids = NULL, coverage = NULL, n0 = 1L,
                             states = NULL,
                             adjacency = c("observed_triple", "shared_feature"),
                             min.support = 1L, feature.ids = NULL,
                             feature.labels = NULL, backend = c("auto", "dense", "sparse")) {
  adjacency <- match.arg(adjacency)
  backend <- match.arg(backend)
  fail <- function(message) stop(paste0("linf.udcst.graph: ", message), call. = FALSE)
  integer.scalar <- function(x) is.numeric(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x) && x >= 1 && x %% 1 == 0
  if (!integer.scalar(n0) || !integer.scalar(min.support))
    fail("n0 and min.support must be positive integers")
  if (sum(c(!is.null(coverage), !is.null(states), n0 != 1)) > 1L)
    fail("use only one of coverage, states, or nondefault n0")
  if (!is.null(coverage) && (!is.numeric(coverage) || is.complex(coverage) ||
      length(coverage) != 1L || !is.finite(coverage) || coverage <= 0 || coverage > 100))
    fail("coverage must be in (0, 100]")
  if (!is.null(X)) {
    if (!is.null(pair.fit) || !is.null(triple.fit)) fail("supply X or fitted assignments, not both")
    pair.fit <- linf.udcsts(X, 2L, 1L, "pure", feature.ids, feature.labels, backend)
    triple.fit <- linf.udcsts(X, 3L, 1L, "pure", feature.ids, feature.labels, backend)
  } else if (is.null(pair.fit) || is.null(triple.fit)) fail("supply X or both fitted assignments")
  valid.fit <- function(f, level) {
    if (!is.list(f) || !identical(f$level, as.integer(level)) ||
        !is.matrix(f$ordered.features) || ncol(f$ordered.features) != level ||
        nrow(f$ordered.features) < 1 || !is.numeric(f$ordered.features) ||
        any(!is.finite(f$ordered.features)) ||
        any(f$ordered.features %% 1 != 0) || any(f$ordered.features < 1) ||
        any(f$ordered.features > length(f$feature.ids)) ||
        anyDuplicated(f$feature.ids) || anyNA(f$feature.ids) ||
        length(f$feature.ids) < 3 || length(f$feature.labels) != length(f$feature.ids) ||
        anyNA(f$feature.labels)) fail("invalid fitted assignment structure")
    sets <- t(apply(f$ordered.features, 1, sort))
    if (any(apply(sets, 1, anyDuplicated) > 0) ||
        !identical(unname(sets), unname(f$feature.sets)) ||
        !identical(unname(apply(sets, 1, paste, collapse=",")), unname(f$pure.id)))
      fail("fitted pure IDs or feature sets disagree with ordered ranks")
  }
  valid.fit(pair.fit, 2L); valid.fit(triple.fit, 3L)
  if (!identical(pair.fit$feature.ids, triple.fit$feature.ids) ||
      !identical(pair.fit$feature.labels, triple.fit$feature.labels) ||
      !identical(pair.fit$input.dimnames, triple.fit$input.dimnames) ||
      !identical(unname(pair.fit$ordered.features), unname(triple.fit$ordered.features[,1:2,drop=FALSE])))
    fail("pair and triple fits must have identical row order, dictionary and leading ranks")
  N <- length(pair.fit$pure.id)
  if (is.null(sample.ids)) sample.ids <- pair.fit$input.dimnames[[1L]]
  if (is.null(sample.ids)) sample.ids <- as.character(seq_len(N))
  if (!is.character(sample.ids) || length(sample.ids) != N || anyNA(sample.ids) ||
      any(!nzchar(sample.ids)) || anyDuplicated(sample.ids)) fail("sample.ids must be unique nonempty character IDs")
  pf <- linf.udcst.frequencies(table(pair.fit$pure.id))
  tf <- linf.udcst.frequencies(table(triple.fit$pure.id))
  if (!is.null(coverage)) n0 <- max(pf$freq[pf$cumm.pct >= coverage])
  if (!is.null(states)) {
    if (!is.character(states) || anyNA(states) || anyDuplicated(states) ||
        any(!states %in% pf$ID)) fail("states must be distinct observed pair IDs")
    nodes <- pf[pf$ID %in% states,,drop=FALSE]
  } else nodes <- pf[pf$freq >= n0,,drop=FALSE]
  sets <- pair.fit$feature.sets[match(nodes$ID, pair.fit$pure.id),,drop=FALSE]
  nodes$node <- seq_len(nrow(nodes))
  nodes$label <- vapply(seq_len(nrow(nodes)), function(i)
    paste(pair.fit$feature.labels[sets[i,]], collapse=" + "), "")
  nodes$feature.a <- sets[,1]; nodes$feature.b <- sets[,2]
  # Sparse tabulation avoids a triples-by-pairs dense contingency matrix.
  keys <- paste(pair.fit$pure.id, triple.fit$pure.id, sep="|")
  counts <- table(keys); bits <- strsplit(names(counts), "|", fixed=TRUE)
  incidence <- data.frame(pair=vapply(bits, `[`, "", 1L),
    triple=vapply(bits, `[`, "", 2L), count=as.integer(counts), stringsAsFactors=FALSE)
  lookup <- stats::setNames(incidence$count, names(counts))
  pcounts <- stats::setNames(pf$freq,pf$ID); tcounts <- stats::setNames(tf$freq,tf$ID)
  getcount <- function(x, key) {v <- unname(x[key]); v[is.na(v)] <- 0; v}
  # Enumerate shared-feature candidates once, even for a missing triple.
  byfeature <- split(rep(seq_len(nrow(nodes)), each=2), as.vector(t(sets)))
  pairs <- lapply(byfeature, function(ii) if(length(ii)>1) t(utils::combn(ii,2)) else NULL)
  ij <- do.call(rbind, pairs)
  if (is.null(ij)) ij <- matrix(integer(), ncol=2)
  if (nrow(ij)) ij <- unique(t(apply(ij,1,sort)))
  edges <- data.frame(from=as.integer(ij[,1]),to=as.integer(ij[,2]))
  edges$pair.a <- nodes$ID[edges$from]; edges$pair.b <- nodes$ID[edges$to]
  edges$triple <- vapply(seq_len(nrow(edges)),function(i)
    paste(sort(unique(c(sets[edges$from[i],],sets[edges$to[i],]))),collapse=","), "")
  edges$n.a <- getcount(pcounts,edges$pair.a); edges$n.b <- getcount(pcounts,edges$pair.b)
  edges$n.triple <- getcount(tcounts,edges$triple)
  edges$i.a <- getcount(lookup,paste(edges$pair.a,edges$triple,sep="|"))
  edges$i.b <- getcount(lookup,paste(edges$pair.b,edges$triple,sep="|"))
  edges$min.support <- pmin(edges$i.a,edges$i.b)
  edges$jaccard.a <- edges$i.a/(edges$n.a+edges$n.triple-edges$i.a)
  edges$jaccard.b <- edges$i.b/(edges$n.b+edges$n.triple-edges$i.b)
  edges$jaccard.mean <- (edges$jaccard.a+edges$jaccard.b)/2
  edges$jaccard.min <- pmin(edges$jaccard.a,edges$jaccard.b)
  edges$jaccard.geometric <- sqrt(edges$jaccard.a*edges$jaccard.b)
  edges$pair.fraction.a <- edges$i.a/edges$n.a; edges$pair.fraction.b <- edges$i.b/edges$n.b
  edges$triple.fraction.a <- ifelse(edges$n.triple>0,edges$i.a/edges$n.triple,0)
  edges$triple.fraction.b <- ifelse(edges$n.triple>0,edges$i.b/edges$n.triple,0)
  edges$length <- rep(1,nrow(edges)); candidates <- edges
  if(adjacency == "observed_triple") edges <- edges[edges$min.support >= min.support,,drop=FALSE]
  faces <- data.frame(triple=character(),a=integer(),b=integer(),c=integer(),
                      min.support=integer(),n.triple=integer())
  for(T in unique(incidence$triple)) {
    z <- incidence[incidence$triple==T & incidence$count>=min.support,,drop=FALSE]
    ids <- match(z$pair,nodes$ID)
    if(nrow(z)==3 && !anyNA(ids)) faces[nrow(faces)+1,] <- list(T,ids[1],ids[2],ids[3],min(z$count),getcount(tcounts,T))
  }
  adj <- vector("list",nrow(nodes))
  for(i in seq_len(nrow(edges))) {a<-edges$from[i];b<-edges$to[i];adj[[a]]<-c(adj[[a]],b);adj[[b]]<-c(adj[[b]],a)}
  component <- integer(nrow(nodes)); nc <- 0L
  for(i in seq_len(nrow(nodes))) if(component[i]==0L) {
    nc <- nc+1L; todo <- i; component[i] <- nc; pos<-1L
    while(pos<=length(todo)) {u<-todo[pos];pos<-pos+1L;new<-adj[[u]][component[adj[[u]]]==0L];component[new]<-nc;todo<-c(todo,new)}
  }
  nodes$component <- component; rownames(nodes)<-NULL; rownames(edges)<-NULL
  structure(list(nodes=nodes,edges=edges,candidate.edges=candidates,faces=faces,
    incidence=incidence,pair.frequencies=pf,triple.frequencies=tf,
    membership=data.frame(sample.id=sample.ids,pair=pair.fit$pure.id,triple=triple.fit$pure.id),
    dictionary=data.frame(index=seq_along(pair.fit$feature.ids),id=pair.fit$feature.ids,label=pair.fit$feature.labels),
    metadata=list(reference.n=N,selected.n=sum(nodes$freq),coverage=100*sum(nodes$freq)/N,
      requested.coverage=coverage,n0=if(is.null(states))n0 else NULL,selected.states=nodes$ID,
      adjacency=adjacency,min.support=min.support,edge.length="unit",components=nc,
      component.sizes=as.integer(table(component)),isolates=sum(lengths(adj)==0),
      tie.rule="decreasing abundance then original column order",schema.version=1L)),
    class="linf_udcst_graph")
}
