gorica <-
 function(object, ..., iter=100000){
   UseMethod("gorica")
 }


gorica_penalty <-
function(object, iter=100000, mc.cores=1){
     
     if (!(inherits(object, "ormle") )) stop("object needs to be of class ormle")
     if (all(object$constr == 0) & object$nec == 0){
     penalty <- length(object$est)
   } else {
     if (iter < 1) stop("No of iterations < 1")
      
      est<-object$est 
      K<-length(est)
      covmtrx <- object$covmtrx
      constr<-object$constr
      rhs=object$rhs
      nec=object$nec
    
      Z <- rmvnorm(n=iter, mean=rep(0, K), sigma=covmtrx)
      Dmat2=2*ginv(covmtrx)

      nact <- apply(Z, 1, function(z){ 

      dvec2=2*(z%*%ginv(covmtrx)) 
      solveQP2= solve.QP(Dmat2,dvec2,t(constr),rhs,meq =nec,factorized = FALSE)
      if (solveQP2$iact[1] == 0) return(0) else return(length(solveQP2$iact))
})

dimsol <- K - nact
LP <- sapply(1:K, function(x) sum(x == (dimsol)))/iter
penalty <- sum((1:K)*LP[])

}

return(penalty)

}



 gorica.ormle <-
 function(object, ..., iter=100000){
   if (!inherits(object, "ormle") & !inherits(object, "list")) stop("object needs to be of class ormle or a list of ormle objects")
   if (iter < 1) stop("No of iterations < 1")
   if (inherits(object, "ormle")) objlist <- list(object, ...) else objlist <- object
   isorlm <- sapply(objlist, function(x) inherits(x, "ormle"))
   orlmlist <- objlist[isorlm]  
   Call <- match.call()
   Call$iter <- NULL
   if (inherits(object, "ormle")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
   loglik <- -2*sapply(orlmlist, function(x) x$logLik)
   penalty <- 2*sapply(orlmlist, function(x) gorica_penalty(x, iter=iter))
   gorica <- loglik + penalty
delta <- gorica - min(gorica)
gorica_weights <- exp(-delta/2) / sum(exp(-delta/2))
data.frame(misfit=loglik,complexity=penalty,gorica=gorica,gorica_weights=round(gorica_weights,4))
 }


gorica.list <- function(object, ..., iter=100000){
   if (all(sapply(object, class) == "ormle")) out <- gorica.ormle(object, iter=iter)
   return(out)
 }