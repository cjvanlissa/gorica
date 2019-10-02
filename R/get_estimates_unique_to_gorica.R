# @export
#get_estimates <- function (x, ...)
#{
#  UseMethod("get_estimates", x)
#}

# @method get_estimates lmerMod
# @export
#' @method get_estimates lmerMod
#' @export
#' @import bain
#' @importFrom lme4 fixef
get_estimates.lmerMod <- function (x, ...)
{
  out <- list(estimate = fixef(x), Sigma = vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "lme4"
  out
}

#' @method get_estimates lavaan
#' @export
get_estimates.lavaan <- function(x, standardize = FALSE, ...){
  cl <- as.list(match.call()[-1])
  cl <- c(cl, list(retain_which = c("=~", "~", "~1", "~~", "=="),
                    split_multigroup_sigma = FALSE,
                    allow_between_constraints = TRUE))
  out <- do.call(lav_get_estimates, cl)
  names(out)[which(names(out) == "x")] <- "estimate"
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "lavaan"
  out
}
