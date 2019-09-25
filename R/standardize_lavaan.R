if(FALSE){


converted <- mplus2lavaan("c:/git_repositories/lav_test.inp", run = F)
converted$model <- gsub("\\*", "", converted$model)
#names(converted$data) <- toupper(names(converted$data))
lavobject <- sem(converted$model, data=converted$data,
                      group = "CLASS")

lav_standardize_parameters <- function (lavobject = NULL, partable = NULL, est = NULL, est.std = NULL,
          GLIST = NULL, cov.std = TRUE, ov.var = NULL, lavmodel = NULL,
          lavpartable = NULL, cov.x = NULL)
{
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
      }
      else {
        stop("lavaan ERROR: could not find `est' in lavpartable")
      }
    }
  }
  else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lavaan:::lav_object_inspect_est(lavobject)
    }
    if (lavmodel@conditional.x) {
      cov.x <- lavobject@SampleStats@cov.x
    }
  }
  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(GLIST)) {
    GLIST <- lavmodel@GLIST
  }
  if (is.null(est.std)) {
    est.std <- lavaan:::lav_standardize_lv(lavobject = lavobject,
                                  partable = partable, est = est, GLIST = GLIST, cov.std = cov.std,
                                  lavmodel = lavmodel, lavpartable = lavpartable)
  }
  out <- est.std
  N <- length(est.std)
  stopifnot(N == length(partable$lhs))
  VY <- computeVY(lavmodel = lavmodel, GLIST = GLIST, diagonal.only = TRUE)
  for (g in 1:lavmodel@nblocks) {
    ov.names <- vnames(lavpartable, "ov", block = g)
    lv.names <- vnames(lavpartable, "lv", block = g)
    if (is.null(ov.var)) {
      OV2 <- VY[[g]]
      zero.idx <- which(abs(OV2) < .Machine$double.eps)
      if (length(zero.idx) > 0L) {
        OV2[zero.idx] <- as.numeric(NA)
      }
      tmp.OV2 <- OV2
      neg.idx <- which(tmp.OV2 < 0)
      if (length(neg.idx) > 0L) {
        tmp.OV2[neg.idx] <- as.numeric(NA)
      }
      OV <- sqrt(tmp.OV2)
    }
    else {
      OV2 <- ov.var[[g]]
      OV <- sqrt(OV2)
    }
    if (lavmodel@conditional.x) {
      ov.names.x <- vnames(lavpartable, "ov.x", block = g)
      ov.names.nox <- vnames(lavpartable, "ov.nox", block = g)
      ov.names <- c(ov.names.nox, ov.names.x)
      OV2 <- c(OV2, diag(cov.x[[g]]))
      OV <- c(OV, sqrt(diag(cov.x[[g]])))
    }
    idx <- which(partable$op == "=~" & !(partable$rhs %in%
                                           lv.names) & partable$block == g)
    out[idx] <- out[idx]/OV[match(partable$rhs[idx], ov.names)]
    idx <- which((partable$op == "~" | partable$op == "<~") &
                   partable$lhs %in% ov.names & partable$block == g)
    out[idx] <- out[idx]/OV[match(partable$lhs[idx], ov.names)]
    idx <- which((partable$op == "~" | partable$op == "<~") &
                   partable$rhs %in% ov.names & partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$rhs[idx], ov.names)]
    rv.idx <- which(partable$op == "~~" & !(partable$lhs %in%
                                              lv.names) & partable$lhs == partable$rhs & partable$block ==
                      g)
    out[rv.idx] <- (out[rv.idx]/OV2[match(partable$lhs[rv.idx],
                                          ov.names)])
    if (cov.std) {
      if (!is.complex(est[rv.idx])) {
        RV <- sqrt(abs(est[rv.idx]))
      }
      else {
        RV <- sqrt(est[rv.idx])
      }
      rv.names <- partable$lhs[rv.idx]
    }
    idx.lhs <- which(partable$op == "~~" & !(partable$lhs %in%
                                               lv.names) & partable$lhs != partable$rhs & partable$block ==
                       g)
    if (length(idx.lhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.lhs] <- (out[idx.lhs]/OV[match(partable$lhs[idx.lhs],
                                               ov.names)])
      }
      else {
        out[idx.lhs] <- (out[idx.lhs]/RV[match(partable$lhs[idx.lhs],
                                               rv.names)])
      }
    }
    idx.rhs <- which(partable$op == "~~" & !(partable$rhs %in%
                                               lv.names) & partable$lhs != partable$rhs & partable$block ==
                       g)
    if (length(idx.rhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.rhs] <- (out[idx.rhs]/OV[match(partable$rhs[idx.rhs],
                                               ov.names)])
      }
      else {
        out[idx.rhs] <- (out[idx.rhs]/RV[match(partable$rhs[idx.rhs],
                                               rv.names)])
      }
    }
    idx <- which(partable$op == "~1" & !(partable$lhs %in%
                                           lv.names) & partable$block == g)
    out[idx] <- out[idx]/OV[match(partable$lhs[idx], ov.names)]
    idx <- which(partable$op == "|" & !(partable$lhs %in%
                                          lv.names) & partable$block == g)
    out[idx] <- out[idx]/OV[match(partable$lhs[idx], ov.names)]
    idx <- which(partable$op == "~*~" & !(partable$lhs %in%
                                            lv.names) & partable$block == g)
    out[idx] <- 1
  }
  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@def.function(x)
  }
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@ceq.function(x)
  }
  idx <- which((partable$op == "<" | partable$op == ">"))
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@cin.function(x)
  }
  out
}

}
