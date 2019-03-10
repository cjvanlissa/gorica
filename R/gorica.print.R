#' @method print gorica
#' @export
print.gorica <- function(x,
                       digits = 3,
                       na.print = "", ...){
  dat <- as.matrix(x$fit)
  #fits <- x$fit
  #dat <- fits[, stats]
  #miss_val <- is.na(dat)
  dat <- formatC(dat, digits = digits, format = "f")
  #dat[miss_val] <- ""
  rownames(dat) <- paste0("H", c(1:(nrow(dat)-1), "u"))
  cat("Informative hypothesis test for an object of class ", class(x$model)[1], ":\n\n", sep = "")
  prmatrix(dat,
           quote = FALSE,
           na.print = na.print)

  cat("\nHypotheses:\n ", paste(rownames(dat)[-nrow(dat)], ": ", x$hypotheses, sep = "", collapse = "\n  "))

  if(!is.null(x[["warnings"]])){
    warning("Gorica analysis returned the following warnings:\n  ", paste(1:length(x$warnings), ". ", x$warnings, sep = "", collapse = "\n  "))
  }
}
