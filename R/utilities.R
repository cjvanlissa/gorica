with_env <- function(fun, envir=parent.frame(), ...) {
  stopifnot(is.function(fun))
  environment(fun) <- envir
  do.call(fun, list(...))
}

create_res <- function(){
  assign("Goricares", list(
    fit = NULL,
    call = cl,
    model = x,
    estimates = x,
    Sigma = Sigma,
    comparison = comparison
  ), envir = parent.env())
}


hypothesis_remove <- function(which_par){
  hypothesis$hyp_mat <<- lapply(hypothesis$hyp_mat, function(R){
    sweep(R[, -which_par, drop = FALSE], MARGIN = 1, as.vector(R[, which_par]))
  })
  # Drop parameters not in hypothesis
  x <<- x[-which_par]
  Sigma <<- Sigma[-which_par, -which_par]
}

hypothesis_adjust <- function(){
  browser()
  for(i in 1:length(hypothesis$hyp_mat)){
    hyp <- hypothesis$hyp_mat[[i]]
    all_zeroes <- rowSums(hyp[, -ncol(hyp)] == 0) == (ncol(hyp)-1)
    # If there are row(s) with only zeros...
    if(any(all_zeroes)){
      browser()
      positive_constant <- hypothesis[all_zeroes, ncol(hypothesis)] > 0
      if(any(positive_constant)){
        q_adjust <- which(apply(hypothesis[all_zeroes, -ncol(hypothesis), drop = FALSE], 2, function(y) any(y != 0)))

        q <- ldei(E= hypothesis[, -ncol(hypothesis), drop = FALSE],
                  F = t(hypothesis[, ncol(hypothesis), drop = TRUE]),
                  G = diag(1, n_pars)[zero_est, , drop = FALSE],
                  H = rep(0, length(zero_est)))$X

        qadj2=zero_est[which(zero_est!=q_adjust)]
        ##########################
        fff=matrix(0,nrow = M,ncol=nrow(BKcovx)) # Empty matrix; M = number of hypotheses, nrow(BKcovx) is number of estimates
        fff[j,c(q_adjust)]= q[q_adjust]
        fff[j,c(qadj2)]=0
        #########################



        q[q_adjust]=0


        bvec2=t(Amat%*%q)




      }
    }
  }
}
