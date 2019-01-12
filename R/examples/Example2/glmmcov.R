 covarglmm <-
 function(formula,B,cluster,fact,std){
 UseMethod("covglmm")
 }




print.covglmm <- function(x, digits = max(3, getOption("digits") - 3), ...){
cat("\n$Call\n")
print(x$Call) 
cat("\n")
cat("\n$est\n")
print(x$est) 
cat("\n")
cat("\n$covmtrx\n")
print(x$covmtrx) 
cat("\n")
invisible(x)  
 }


covglmm<-
function(formula,B,cluster,fact,std){
Call = match.call()
cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula","data"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")


y <- cbind(model.response(mf, "numeric"))
x <- model.matrix(mt, mf, contrasts)  


if ( (any(cluster<=0)) | (length(cluster)!=1)) stop("There should be only one cluster with a positive index number")
data=eval(summary(formula)$call$data)
cluster=data[,cluster]
model=formula
K=length(summary(model)$coefficients[,1])
J=length(unique(cluster))
boot=matrix(c(NA),nrow=B,ncol=K)

for (b in 1:B) {
cat(paste("The index of the current bootstrap sample: ",b,"\n"))

c.data=cbind(data,cluster)

c <- sort(unique(c.data$cluster))

clust.group <- function(c) {
    c.data[c.data$cluster==c,]
}

clust.list <- lapply(c,clust.group)

n=sapply(clust.list, nrow)

index=sapply(n, sample,replace=TRUE)

for (i in 1:J) { clust.list[[i]]=clust.list[[i]][index[[i]],] }

DataB=do.call("rbind", clust.list)


if (all(fact>0)) {DataB[,fact]=factor(DataB[,fact])}

if (any(fact<0)) stop("The index number for the factor variables cannot be negative")


if (all(std>0)) {DataB[,std]=scale(DataB[,std])}

if (any(std<0)) stop("The index number for the standardized variables cannot be negative")


model=lmer(eval(summary(model)$call$formula),data=DataB,REML=FALSE)
boot[b,]=summary(model)$coefficients[,1]

}

est=apply(boot[,1:K],2,mean)
covmtrx=cov(boot[,])


out <- list(Call=Call,est=est,covmtrx=covmtrx,bootsample=boot)
class(out) <- "covglmm"
return(out)


}




