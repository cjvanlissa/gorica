.mplus_tech1_sections <- c("tau" = "Thresholds", "nu" = "Intercepts", "lambda" = ".BY", "alpha" = "Means", "beta" = ".ON", "gamma" = ".ON", "psi" = ".WITH")

#' @method get_estimates mplus.model
#' @export
get_estimates.mplus.model <- function(x, ...){
  techs_missing <- c("tech1" = is.null(x[["tech1"]]), "tech3" = is.null(x[["tech3"]]))
  techs_missing <- sapply(names(techs_missing)[!techs_missing], function(i){!length(x[[i]]) > 0})
  if(any(techs_missing)) stop("Need the ", paste(names(techs_missing)[techs_missing], collapse = " and "), " output in order to derive the covariance matrix of the estimates from Mplus. Please re-run your analysis in Mplus, adding the following syntax:\nOUTPUT: tech1 tech3;")

  par_spec <- x$tech1$parameterSpecification

  if(all(is.na(match(names(.mplus_tech1_sections), names(par_spec))))){

  }
  #sections <- sections[match(names(sections), names(par_spec))]

  param_id <- do.call(rbind, lapply(names(.mplus_tech1_sections), function(x) {
    sec <- par_spec[[x]]
    if (length(rownames(sec) == 1) & rownames(sec)[1] == "1") {
      cbind(.mplus_tech1_sections[x], colnames(sec), as.vector(sec))
    } else {
      cbind(rep(paste0(rownames(sec), .mplus_tech1_sections[x]), each = ncol(sec)),
            rep(colnames(sec), nrow(sec)),
            as.vector(t(sec)))
    }
  }))

  param_id <- param_id[!(is.na(param_id[, 3]) | param_id[, 3] == "0"), ]

  param_id <- data.frame(label = paste(param_id[, 1], param_id[, 2], sep = "."),
                         id = as.integer(param_id[, 3]))


  estimate <- x$parameters$unstandardized
  estimate$label <- paste(estimate$paramHeader, estimate$param, sep = ".")
  estimate <- merge(param_id, estimate, by = "label", sort = FALSE)

  estimate$label <- gsub("[\\$\\.]", "_", estimate$label)

  coefs <- estimate$est
  names(coefs) <- estimate$label

  Sigma <- x$tech3$paramCov
  Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
  colnames(Sigma) <- estimate$label
  rownames(Sigma) <- estimate$label

  out <- list(estimate = coefs,
              Sigma = Sigma)
  class(out) <- "gorica_estimate"
  out
}

#' @method print gorica_estimate
#' @export
print.gorica_estimate <- function(x,
                         digits = 3,
                         na.print = "", ...){
  dat <- x$estimate
  dat <- formatC(dat, digits = digits, format = "f")
  print(dat, quote = FALSE)
}


