.lavaan_sections <- matrix(c(
    "tau", "vector", "Thres_", "", "",
    "nu", "vector", "Intercept_", "", "",
    "lambda", "matrix", "", "_BY_", "",
    "theta", "matrix", "", "_WITH_", "",
    "alpha", "vector", "Intercept_", "", "",
    "beta", "matrix", "", "_ON_", "",
    "gamma", "matrix", "", "_ON_", "",
    "psi", "matrix", "", "_WITH_", "",
    "delta", "vector", "bla", "", ""
  ), ncol = 5, byrow = TRUE)

#' @method get_estimates lavaan
#' @export
get_estimates.lavaan <- function(x, ...){
  param_id <- lavInspect(x, "free")
  est_lavaan <- lavInspect(x, "est")
  Sigma <- lavTech(x, "vcov")
  keep <- lapply(param_id, `==`, 0)
  if(any(names(keep) %in% c("theta", "psi", "beta", "gamma"))){
    keep[na.omit(match(c("theta", "psi", "beta", "gamma"), names(keep)))] <-
      lapply(
      keep[na.omit(match(c("theta", "psi", "beta", "gamma"), names(keep)))],
      function(x){x[upper.tri(x)] <- TRUE
      x}
    )
  }

  estimate_labels <- lapply(names(param_id), function(x){
    label_section(sec = param_id[[x]],
                  pref = .lavaan_sections[.lavaan_sections[, 1] == x, 3],
                  mid = .lavaan_sections[.lavaan_sections[, 1] == x, 4],
                  post = .lavaan_sections[.lavaan_sections[, 1] == x, 5])
  })

  estimates <- data.frame(
    Estimate = unlist(est_lavaan),
    id = as.integer(unlist(param_id)),
    label = unlist(estimate_labels),
    keep = !unlist(keep),
    stringsAsFactors = FALSE)
  estimates <- estimates[estimates$keep, ]
  #estimates <- estimates[!estimates$id == 0, ]

  duplicates <- NULL
  if(any(duplicated(estimates$id))){
    dups <- estimates[duplicated(estimates$id), ]
    dups <- estimates[duplicated(estimates$id, fromLast = T), ]
    estimates <- estimates[!duplicated(estimates$id), ]
    duplicates <- lapply(unique(dups$id), function(o){dups$label[dups$id == o]})
    names(duplicates) <- estimates$label[match(dups$id, estimate$id)]
  }
  estimates <- estimates[order(estimates$id), ]



  coefs <- estimate$Estimate
  names(coefs) <- estimate$label

  Sigma <- Sigma[estimates$id, estimates$id]
  colnames(Sigma) <- estimates$label
  rownames(Sigma) <- estimates$label

  out <- list(estimate = coefs,
              Sigma = Sigma)
  if(!is.null(duplicates)) out[["label_synonyms"]] <- duplicates
  class(out) <- "gorica_estimate"
  out
}
