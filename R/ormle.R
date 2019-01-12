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
#' @export
ormle <-
  function(est,covmtrx,constr,rhs,nec){
  K=length(est)
  covmtrx = as.matrix(covmtrx)
  Dmat = 2*ginv(covmtrx)
  dvec = 2*(est%*% ginv(covmtrx))
  solveQP = solve.QP(Dmat, dvec = dvec, t(constr), rhs, meq = nec, factorized = FALSE)
  tildeQ = solveQP$solution
restrictedest=solveQP$solution
names(restrictedest)=names(est)
loglik =as.numeric( ( -K/2*log(2*pi) )-( 0.5*log(det(covmtrx) ) )-( 0.5* t(est- tildeQ)%*%ginv(covmtrx)%*% (est-tildeQ)) )

out <- list(est=est, covmtrx=covmtrx, constr=constr, rhs=rhs, nec=nec, logLik=loglik,restrictedest=restrictedest)
    class(out) <- "ormle"
    return(out)

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
