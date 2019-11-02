#' @title ormle
#' @description FUNCTION_DESCRIPTION
#' @param est PARAM_DESCRIPTION
#' @param covmtrx PARAM_DESCRIPTION
#' @param constr PARAM_DESCRIPTION
#' @param rhs PARAM_DESCRIPTION
#' @param nec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ormle
#' @keywords internal
ormle <-
  function(est,covmtrx,constr,rhs,nec){
  covmtrx = as.matrix(covmtrx)
  ginvcovmat <- ginv(covmtrx)
  Dmat = 2*ginvcovmat
  dvec = 2*(est%*% ginvcovmat)
  solveQP = solve.QP(Dmat, dvec = dvec, t(constr), rhs, meq = nec, factorized = FALSE)
  restrictedest = solveQP$solution
  names(restrictedest)=names(est)
  loglik = as.numeric( ( -length(est)/2*log(2*pi) )-( 0.5*log(det(covmtrx) ) )-( 0.5* t(est- restrictedest)%*%ginvcovmat%*% (est-restrictedest)) )

  out <- list(est=est, covmtrx=covmtrx, constr=constr, rhs=rhs, nec=nec, logLik=loglik,restrictedest=restrictedest)
  class(out) <- "ormle"
  out

   }

print.ormle <- function(x, digits = max(3, getOption("digits") - 3), ...){
cat("\n$est\n")
print(x$est)
cat("\n")
cat("\n$restrictedest\n")
print(x$restrictedest)
cat("\n")
invisible(x)
 }
