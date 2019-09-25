.lavaan_sections <- matrix(c(
    "tau", "vector", "___Thres___", "", "",
    "nu", "vector", "___Int___", "", "~1",
    "lambda", "matrix", "", "=~", "",
    "theta", "matrix", "", "~~", "",
    "alpha", "vector", "___Int___", "", "~1",
    "beta", "matrix", "", "~", "",
    "gamma", "matrix", "", "~", "",
    "psi", "matrix", "", "~~", "",
    "delta", "vector", "bla", "", ""
  ), ncol = 5, byrow = TRUE)

label_section <- function(sec, pref = "", mid = "", post = ""){
  #browser()
  if(mid == "~"){
    tmp <- matrix(paste0(pref, rep(rownames(sec), ncol(sec)), mid, rep(colnames(sec), each = nrow(sec)), post), nrow = dim(sec)[1], ncol = dim(sec)[2])
  } else {
    tmp <- matrix(paste0(pref, rep(colnames(sec), each = nrow(sec)), mid, rep(rownames(sec), ncol(sec)), post), nrow = dim(sec)[1], ncol = dim(sec)[2])
  }
  #tmp <- matrix(paste0(pref, rep(colnames(sec), each = nrow(sec)), mid, rep(rownames(sec), ncol(sec)), post), nrow = dim(sec)[1], ncol = dim(sec)[2])
  #if(mid == "_WITH_"){
  #  diag(tmp) <- gsub("^.+_WITH_", "Var_", diag(tmp))
  #}
  tmp
}

#' @method get_estimates lavaan
#' @export
get_estimates.lavaan <- function(x, ...){
  param_id <- lavInspect(x, "free")
  est_lavaan <- lavInspect(x, "est")
  Sigma <- lavTech(x, "vcov")

# From Mplus --------------------------------------------------------------
  ngroups <- lavInspect(x, what = "ngroups")
  nlevels <- lavInspect(x, what = "nlevels")
  if(ngroups > 1 | nlevels > 1){
    keep <- lapply(param_id, function(this_group){
      lapply(names(this_group), function(x){
        out <- this_group[[x]] == 0
        if(x %in% c("theta", "psi", "gamma")){
          out[upper.tri(out)] <- TRUE
      }
      out
      })
    })


    estimate_labels <- lapply(names(param_id), function(group_name) {
        lapply(names(param_id[[group_name]]), function(x) {
          label_section(
            sec = param_id[[group_name]][[x]],
            pref = .lavaan_sections[.lavaan_sections[, 1] == x, 3],
            mid = .lavaan_sections[.lavaan_sections[, 1] == x, 4],
            post = paste(.lavaan_sections[.lavaan_sections[, 1] == x, 5], group_name, sep = ".")
          )
        })
      })

  } else {
    keep <- lapply(names(param_id), function(x) {
      out <- param_id[[x]] == 0
      if (x %in% c("theta", "psi", "beta", "gamma")) {
        out[upper.tri(out)] <- TRUE
      }
      out
    })
    #origin_matrix <- unlist(lapply(1:length(param_id), function(x){
    #  rep(names(param_id)[x], length(param_id[[x]]))
    #}))
    origin_matrix <- unlist(lapply(1:length(param_id), function(x){
      rep(names(param_id)[x], length(param_id[[x]]))
    }))

    estimate_labels <- lapply(names(param_id), function(x){
      label_section(sec = param_id[[x]],
                    pref = .lavaan_sections[.lavaan_sections[, 1] == x, 3],
                    mid = .lavaan_sections[.lavaan_sections[, 1] == x, 4],
                    post = .lavaan_sections[.lavaan_sections[, 1] == x, 5])
    })

  }

# End from Mplus ----------------------------------------------------------





  estimate_labels <- gsub("^___(Thres|Int)___(threshold|intercept)", "", unlist(estimate_labels))
  estimates <- data.frame(
    Estimate = unlist(est_lavaan),
    id = as.integer(unlist(param_id)),
    label = estimate_labels,
    keep = !unlist(keep),
    stringsAsFactors = FALSE)
  estimates <- estimates[estimates$keep, ]

  duplicates <- NULL
  if(any(duplicated(estimates$id))){
    dups <- estimates[duplicated(estimates$id), ]
    dups <- estimates[duplicated(estimates$id, fromLast = T), ]
    estimates <- estimates[!duplicated(estimates$id), ]
    duplicates <- lapply(unique(dups$id), function(o){dups$label[dups$id == o]})
    names(duplicates) <- estimates$label[match(dups$id, estimates$id)]
  }
  estimates <- estimates[order(estimates$id), ]



  coefs <- estimates$Estimate
  names(coefs) <- estimates$label

  Sigma <- Sigma[estimates$id, estimates$id]
  colnames(Sigma) <- estimates$label
  rownames(Sigma) <- estimates$label

  out <- list(estimate = coefs,
              Sigma = Sigma)
  if(!is.null(duplicates)) out[["label_synonyms"]] <- duplicates
  class(out) <- "gorica_estimate"
  out
}
