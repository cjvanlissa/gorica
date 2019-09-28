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
