#' @method get_estimates lmerMod
#' @export
get_estimates.lmerMod <- function (x, ...)
{
  out <- list(estimate = fixef(x), Sigma = vcov(x))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "lme4"
  out
}
