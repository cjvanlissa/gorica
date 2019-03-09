#' @title Evaluate Informative Hypotheses
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @param iter PARAM_DESCRIPTION, Default: 1e+05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' academic_awards <- within(academic_awards, {
#'   prog <- factor(prog, levels = 1:3, labels = c("General", "Academic",
#'                                                 "Vocational") ) } )
#'   zmath <- scale(academic_awards$math)
#'   #Model fitting
#'   model <- glm(num_awards ~ prog + zmath + prog * zmath, family = "poisson",
#'                data = academic_awards)
#'                new_gor <- gorica(x = model, hypothesis = "zmath=0&zmath+
#'                progAcademic:zmath=0&zmath+progVocational:zmath=0;
#'                progAcademic:zmath > progVocational:zmath&
#'                progVocational:zmath>0; progVocational:zmath<0")
#' }
#' @export
compare_hypotheses <-
  function(object, ..., iter = 100000){
    UseMethod("compare_hypotheses")
  }

#' @method compare_hypotheses ormle
#' @export
compare_hypotheses.ormle <-
  function(object, ..., iter = 100000){
    if (!inherits(object, "ormle") & !inherits(object, "list")) stop("object needs to be of class ormle or a list of ormle objects")
    if (iter < 1) stop("No of iterations < 1")
    if (inherits(object, "ormle")) objlist <- list(object, ...) else objlist <- object
    isorlm <- sapply(objlist, function(x) inherits(x, "ormle"))
    orlmlist <- objlist[isorlm]
    Call <- match.call()
    Call$iter <- NULL
    if (inherits(object, "ormle")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
    loglik <- -2*sapply(orlmlist, function(x) x$logLik)
    penalty <- 2*sapply(orlmlist, function(x) gorica_penalty3(x, iter = iter))
    gor <- loglik + penalty
    delta <- gor - min(gor)
    gorica_weights <- exp(-delta/2) / sum(exp(-delta/2))
    data.frame(misfit = loglik,complexity = penalty,gorica = gor,gorica_weights = round(gorica_weights,4))
  }

#' @method compare_hypotheses list
#' @export
compare_hypotheses.list <- function(object, ..., iter = 100000){
  if (all(sapply(object, class) == "ormle")) out <- compare_hypotheses.ormle(object, iter = iter)
  return(out)
}

#' @importFrom MASS ginv
#' @importFrom mvtnorm rmvnorm
#' @importFrom quadprog solve.QP
gorica_penalty3 <-
  function(object, iter = 100000, mc.cores = 1){
    if(any(object$constr != 0) | object$nec != 0){
      K <- length(object$est)
      Z <- rmvnorm(n = iter, mean = rep(0, K), sigma = object$covmtrx)
      ginvcovm <- ginv(object$covmtrx)
      Dmat2 = 2*ginvcovm

      dvec2 <- 2*(Z %*% t(ginvcovm))
      t_const <- t(object$constr)

      nact <- apply(dvec2, 1, function(z){
        solveQP2 = solve.QP(Dmat2,z,t_const,object$rhs,meq = object$nec,factorized = FALSE)
        if (solveQP2$iact[1] == 0){
          0
        } else {
          length(solveQP2$iact)
        }
      })
      sum((1:K)*sapply(1:K, function(x) sum(x == (K - nact)))/iter)
    } else {
      length(object$est)
    }

  }


#' @importFrom MASS ginv
#' @importFrom mvtnorm rmvnorm
#' @importFrom quadprog solve.QP
gorica_penalty2 <-
  function(object, iter = 100000, mc.cores = 1){
    if(any(object$constr != 0) | object$nec != 0){
      K <- length(object$est)
      Z <- rmvnorm(n = iter, mean = rep(0, K), sigma = object$covmtrx)
      ginvcovm <- ginv(object$covmtrx)
      Dmat2 = 2*ginvcovm
      nact <- apply(Z, 1, function(z){

        dvec2 = 2*(z %*% ginvcovm)
        solveQP2 = solve.QP(Dmat2, dvec2, t(object$constr), object$rhs, meq = object$nec, factorized = FALSE)
        if (solveQP2$iact[1] == 0){
          0
        } else {
          length(solveQP2$iact)
        }
      })

      dimsol <- K - nact
      LP <- sapply(1:K, function(x) sum(x == (dimsol)))/iter
      sum((1:K)*LP)

    } else {
      length(object$est)
    }

  }


gorica_penalty <-
  function(object, iter = 100000, mc.cores = 1){

    if (!(inherits(object, "ormle") )) stop("object needs to be of class ormle")
    if (all(object$constr == 0) & object$nec == 0){
      length(object$est)
    } else {
      if (iter < 1) stop("No of iterations < 1")

      est<-object$est
      K<-length(est)
      covmtrx <- object$covmtrx
      constr<-object$constr
      rhs = object$rhs
      nec = object$nec

      Z <- rmvnorm(n = iter, mean = rep(0, K), sigma = covmtrx)
      Dmat2 = 2*ginv(covmtrx)

      nact <- apply(Z, 1, function(z){

        dvec2 = 2*(z%*%ginv(covmtrx))
        solveQP2 = solve.QP(Dmat2,dvec2,t(constr),rhs,meq = nec,factorized = FALSE)
        if (solveQP2$iact[1] == 0){
          return(0)
        } else {
          return(length(solveQP2$iact))
        }
      })

      dimsol <- K - nact
      LP <- sapply(1:K, function(x) sum(x == (dimsol)))/iter
      sum((1:K)*LP)

    }
  }
