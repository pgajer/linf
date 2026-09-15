linf.normalize.low.freq.policy <- function(low.freq.policy) {
  if (length(low.freq.policy) > 1L) {
    low.freq.policy <- low.freq.policy[[1L]]
  }

  match.arg(low.freq.policy, choices = c("pure", "absorb"))
}

linf.active.low.freq.view <- function(low.freq.policy) {
  if (identical(low.freq.policy, "pure")) {
    return("pure")
  }

  low.freq.policy
}
