#Where to put data for automated tests with testthat? and use inst/testdata, then access the files with system.file("testdata",...,package="gorica"
#f <- list.files("tests/testthat/files/Examples", pattern = ".txt", recursive = T, full.names = T)
#file.rename(f, file.path(getwd(), "inst", "testdata", gsub("^.+?(\\d+)/(.+?)\\.txt$", "\\2_\\1\\.txt", f)))
#


# Yasin's code ------------------------------------------------------------

####Installing automatically the required packages (this code below can run after installing the package easypackages)


library(MASS)
library(quadprog)
library(limSolve)
library(matrixcalc)



old_goricacont <- function(inpfile, datfile){


  ####Reading the Input.txt file which should be in the same folder with the file "GoricaCont.R"

  Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,nrow=1))


  ####The number of cells in the data (for example, for a 2*2*2 contingency table D = 8).
  ####The dev/gorica_cont/data.txt file should be in the same folder with the file "GoricaCont.R"

  D=as.numeric(Input[1,1])


  ####The number of cells in the data must be a positive numeric value

  if (!is.numeric(D)) stop("The number of cells in the data must be a single numeric value")
  if (D<0 | D==0 | is.na(D)==TRUE) stop("The number of cells in the data must be a positive numeric value")



  ####The number of structural parameters(i.e.,cell probabilities with linear restrictions on cell probabilities
  ####and functions of cell probabilities with nonlinear restrictions on cell probabilities) must be a single numeric value


  K=as.numeric(Input[1,2])

  ####The number of structural parameters in the data must be a positive numeric value
  ####Note that D and K are equal to each other when we have restrictions on cell probabilities themselves but
  ####not functions of them

  if (!is.numeric(K)) stop("The number of structural parameters must be a single numeric value")
  if (K<0 | K==0 | is.na(K)==TRUE) stop("The number of structural parameters must be a positive numeric value")



  ####Providing the seed value in Input.txt file to duplicate the results in the paper

  seed=as.numeric(Input[1,3])
  if (!is.numeric(seed) | is.na(seed)==TRUE ) stop("The seed value must be a numeric value ")


  ####The number of bootstrap samples (often B = 1000)

  B=as.numeric(Input[1,4])

  if (!is.numeric(B)) stop("The number of bootstrap samples must be a single numeric value")
  if (B<0 | B==0 | is.na(B)==TRUE) stop("The number of bootstrap samples must be a positive numeric value")


  ####The number of iterations when calculating the penalty part of the GORICA (often T = 10000 or T = 100000)

  T=as.numeric(Input[1,5])

  if (!is.numeric(T)) stop("The number of iterations must be a single numeric value")
  if (T<0 | T==0 | is.na(T)==TRUE) stop("The number of iterations must be a positive numeric value")




  #------------------------------------------------------------------------------------------------
  options(warn=-1)



  ####Reading the data
  dataset=read.table(datfile,header=FALSE)

  ####The observations in the data
  Obs=dataset[,ncol(dataset)]

  ####Creating ObsF based on Obs. This will be used in the bootstrapping later below
  ####Note that the usual boot function in package boot does not work for contingency tables in the context of GORICA,
  ####since there is liear dependency between the cell probabilities and we may be interested in the covariance matrix of
  ####functions of cell probabilities.
  ObsF=rep(1:D,Obs)


  ####The number of observations in the data
  obstotal=length(ObsF)

  ####Assigning place for the estimates of cell probabilities in B bootstrap samples
  bootstrap=matrix(c(NA),nrow=B,ncol=D,byrow=TRUE)

  ####Assigning place for the data in bootstrap samples
  conttab=matrix(c(NA),nrow=1,ncol=D,byrow=TRUE)



  ####Assigning place for the estimates and x matrix will be used in bootstrapping below
  Q=matrix(c(NA), nrow=B, ncol=K)
  x=matrix(c(NA), nrow=B, ncol=K)


  ####Reread the Input.txt file to take into account the restrictions in the hypotheses under evaluation

  Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,nrow=1+K))


  ####If D!=K, then this means that we have linear restrictions on functions of cell probabilities
  ####Therefore we need to read the parameters used in evaluation
  ####When D = K we do not need to specify parameters used in evaluation, since we know that they are cell probabilities

  if (D!=K) {

    AA=eval(parse(text="Input[2:(1+K),1]")) }



  ####Performing bootstrapping

  for (b in 1:B) {

    index=c(1:obstotal)
    indexx=sample(1:obstotal,replace=TRUE)

    A=ObsF[indexx]

    for (i in 1:D) {
      conttab[,i]=length(A[A==i])
    }

    conttab=conttab/obstotal
    bootstrap[b,]=conttab


    for (i in 1:D) {
      x[i]=bootstrap[b,i]

    }

    if (D!=K) {

      for (i in 1:K) {
        Q[b,i]=eval(parse(text=AA[i]))
      }


    }

    #### If D=K the estimates in the b-th bootstrap sample are basically the estimated cell probabilities
    if (D==K) { Q[b,]=bootstrap[b,] }


    ###bootstrapping finishes here
  }

  if(any(is.na(Q))==TRUE) {cat("WARNING: The model parameters cannot be estimated for", length(which(is.na(Q)==TRUE)), "bootstrap sample(s).
Thus, the overall estimates and their covariance matrix are obtained based on", B-length(which(is.na(Q)==TRUE)), "bootstrap samples.\n")}


  #----------------------------------------------------------------------------------------

  ####Obtaining the overall estimates in bootstrapping procedure
  eta=rep(NA,K)
  for (i in 1:K) {
    eta[i]=mean(na.omit(Q[,i]))
  }

  etax=eta
  # etax = c(0.41462237659026, 0.392819543072976, 0.465424773344854, 0.502195103020689)


  ####Obtaining covariance matrix of bootstrap estimates
  ####Note that Q is a matrix containing all the estimates of (functions of) cell probabilities
  BKcovx=cov(na.omit(Q))
  # BKcovx <- eval(parse(text = 'structure(c(0.000163315057642907, -3.64422923510321e-06, -2.31938484805941e-05, -1.68125058063052e-06, -3.64422923510321e-06, 0.000351433075043721, 1.74670548454854e-05, -5.66514272792934e-05, -2.31938484805941e-05, 1.74670548454854e-05, 0.00246101758389758, -8.38347824524701e-05, -1.68125058063052e-06, -5.66514272792934e-05, -8.38347824524701e-05, 0.00482286906061497), .Dim = c(4L, 4L))'))

  ####If D!=K and any(eta=1), then we have restrictions on functions of cell probabilities (eta) and we estimated
  ####one of them as 1. This may happen when the original data contains zeros and the hypotheses under evaluation
  ####are functions of cell probabilities. In such a case we just give warning and say that the user needs to
  ####respecify the hypotheses under evaluation accordingly taking into account that the corresponding eta(s) are 1.

  if ((D!=K) & (any(eta==1)==TRUE)) {

    listeta2=as.matrix(list(etax)[[1]])
    listeta2=t(listeta2)
    A=paste("p",1:nrow(BKcovx),sep="")
    colnames(listeta2)=c(A)
    print(list(MLEs=listeta2))

    stop("Some of the eta parameters are estimated as one because of the empty cell(s). \n  The user needs to rewrite the hypotheses under evaluation.")  }


  ####Similarly if D!=K and some of etas are estimated as Infinity, we give warning and say that the user needs to
  ####respecify the hypotheses accordingly.


  if ((D!=K) & (any(eta==Inf)==TRUE)) {


    listeta2=as.matrix(list(etax)[[1]])
    listeta2=t(listeta2)
    A=paste("p",1:nrow(BKcovx),sep="")
    colnames(listeta2)=c(A)
    print(list(MLEs=listeta2))

    stop("Some of the eta parameters cannot be estimated because of the empty cell(s). \n  The user needs to rewrite the hypotheses under evaluation.")  }


  #if(length(c(which(x==0),which(x==Inf)))==length(etax)) stop("sdfsf")


  ####The code below are used to define whether the eta's are linearly dependent on each other in different conditions.



  wz=max(which(etax!=0))


  if (wz==-Inf) stop("The restrictions are solely on the eta parameters that are zero.
  The hypotheses under evaluation should be rewritten in the linear form.")



  if (ncol(Q)==1) {BKcov=var(Q)}
  if (ncol(Q)!=1) { BKcov=as.matrix(cov(Q[,]))}
  # If raw cell probabilities are used, drop the final element from the eta vector and Sigma
  if ( (D==K) & (etax[D]!=0) )
  {
    eta=eta[-c(D)]
    BKcov=as.matrix(BKcov[-c(D),-c(D)])
  }

  if ( (D==K) & (etax[D]==0) )
  {

    # If the final element is zero, remove the last nonzero element
    eta=eta[-c(wz)]
    BKcov=as.matrix(BKcov[-c(wz),-c(wz)])
  }


  Sz=0

  if (any(apply(BKcov, 1, function(y) all(y==0)))) {

    Sz=which(apply(BKcov, 1, function(y) all(y==0)))



    eta2=eta
    eta=eta[-c(Sz)]
    BKcov2=BKcov

    BKcov=as.matrix(BKcov[-c(Sz),-c(Sz)]) }


  if ((D!=K) & (is.positive.definite(BKcov)==FALSE) ) stop("The functions of cell probabilities (eta) are linearly dependent on each other. \n  The covariance matrix of the eta's is not positive definite.")



  #------------------------------------------------------------------------



  ####Reread the Input.txt file taking into account also the number of hypotheses under evaluation


  Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,nrow=2+K))


  ####When the hypotheses contain only linear restrictions on cell probabilities themselves (implied by K=D in the input file), one should
  ####delete the part "#parameters used in evaluation (eta)" in the input file

  if (D!=K) {M=as.numeric(Input[2+K,1])}
  if (D==K) {M=as.numeric(Input[2,1]) }

  if (D==K & is.na(M)==TRUE)stop("Hypotheses contain linear restrictions on cell probabilities, which is implied by K = D. \n  Specification of ''#Parameters used in evaluation (eta)'' in the input file should be deleted.")


  ####The number of hypotheses (M) must be a positive numeric value


  if (!is.numeric(M)) stop("The number of models must be a single numeric value.")
  if (M<0 | M==0 | is.na(M)==TRUE) stop("The number of models must be a positive numeric value.")


  ####If D!=K then we have linear restrictions on functions of cell probabilites. Therefore, the input file contain the part "#parameters used in evaluation (eta)"
  #### and thus the input file should be read accordingly



  if (D!=K) {


    Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,nrow=2+K+M))

    summ=sum(as.numeric(Input[(3+K):(2+M+K),1:2]))

    Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,skip=K+6,colClasses=rep("numeric", 30)))


  }



  ####If D=K then we have linear restrictions on cell probabilities themselves. The input file does not contain the part "#parameters used in evaluation (eta)"
  ####and thus the input file should be read accordingly



  if (D==K) {


    Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,nrow=2+M))

    summ=sum(as.numeric(Input[3:(2+M),1:2]))

    Input=as.matrix(read.table(inpfile,header=FALSE,fill=TRUE,skip=5,colClasses=rep("numeric", 30)))


  }

  #-----------------------------------------------------------------------


  ####Specifying the number of equality and inequality restrictions in each hypothesis under evaluation
  ####Always contain two columns. The first column represents equality restrictions, the second column represents inequality restrictions
  ####The rows denote the hypotheses. for example, 2 1 in the first row means that hypothesis H1 contains 2 equality restrictions and 1 inequality restriction
  ####If a row is 0 0, this means that this hypothesis is the unconstrained hypothesis which does not impose any restriction on (functions of) cell probabilities

  rest=matrix(c(as.numeric(Input[1:M,1:2])),nrow=M,ncol=2)




  ####Assigning place for the order-restricted estimates

  tildeeta=matrix(c(NA), nrow=M, ncol=nrow(BKcov))





  ####Calculating the part common to all hypotheses when calculating the order-restricted estimates
  ####(BKcov is the covariance matrix obtained with nonparametric bootstrapping )


  Dmat = 2*ginv(BKcov)

  sc = norm(Dmat,"O")

  dvec = 2*(eta%*%ginv(BKcov))



  ####Assigning places for log likelihood and penalty parts (M is the number of hypotheses under evaluation)

  logHm =matrix(c(NA), nrow=1, ncol=M)
  PTm =matrix(c(NA), nrow=1, ncol=M)

  ll=1
  fff=matrix(c(0),nrow=M,ncol=nrow(BKcovx))

  #################################################################M loop starts######################################################################
  for (j in 1:M) {


    ####If there is no equality and inequality restriction, the hypothesis is the unconstrained hypothesis
    ####The restriction matrix (Amat) and the right hand side values (bvec) are determined accordingly


    if (sum(rest[j,])==0) {

      Amat = matrix(c(rep(0,K)),nrow=1,ncol=K)
      bvec=0

    }

    ####If there is equality and/or inequality restriction in the hypothesis, it is not the unconstrained hypothesis
    ####The restriction matrix (Amat) and the right hand side values (bvec) are determined accordingly



    if (sum(rest[j,])!=0) {


      Amat=matrix(c(as.numeric(Input[(M+ll):(M+ll-1+sum(rest[j,])),1:K])),nrow=sum(rest[j,]),ncol=K)
      bvec=matrix(c(as.numeric(Input[(M+ll+nrow(Amat)):(M+ll-1+nrow(Amat)+1),1:nrow(Amat)])),nrow=1,ncol=sum(rest[j,]))

      AA=apply(Input, 2, function(y) all(is.na(y)))

      ###The number of structural parameters and the number of columns in the restriction matrix must be equal to each other

      if ((ncol(t(na.omit(t(Amat))))!=K) | ((length(AA[AA==FALSE])!=K) & K >1)) stop("The number of eta parameters (K) and the number of columns in restriction matrices are not equal to each other.")

    }



    ####The code below deal with determining the restriction matrices for different types of
    ####restrictions


    for (lst in 1:nrow(Amat)) {

      if (D==K & Amat[lst,K]!=0)

      {
        bvec[lst]=bvec[lst]+((-1)*Amat[lst,K])
        Amat[lst,]=Amat[lst,]+((-1)*Amat[lst,K]*rep(1,K))

      }

    }



    for (lst in 1:nrow(Amat)) {

      if (D==K & Amat[lst,K]==0)

      {
        bvec[lst]=bvec[lst]+((-1)*Amat[lst,K])
        bvec[lst]=bvec[lst]+((-1)*Amat[lst,wz])


        Amat[lst,]=Amat[lst,]+((-1)*Amat[lst,K]*rep(1,K))
        Amat[lst,1:(K-1)]=Amat[lst,1:(K-1)]+((-1)*Amat[lst,wz]*rep(1,K-1))

      }

    }


    ll=ll+1+sum(rest[j,])


    if ( (sum(rest[j,])==0) & (D==K) & any(as.numeric(Sz)>0)==TRUE ) {Amat2=matrix(c(as.numeric(Amat[,])), nrow=1,ncol=K-1-length(Sz)) }

    if ( (sum(rest[j,])==0) & (D==K) & any(as.numeric(Sz)>0)==FALSE ) {Amat2=matrix(c(as.numeric(Amat[,])), nrow=1,ncol=K-1) }

    if ( (sum(rest[j,])==0) & (D!=K) & any(as.numeric(Sz)>0)==TRUE )  {Amat2=matrix(c(as.numeric(Amat[,])), nrow=1,ncol=K-length(Sz)) }

    if ( (sum(rest[j,])==0) & (D!=K) & any(as.numeric(Sz)>0)==FALSE ) {Amat2=matrix(c(as.numeric(Amat[,])), nrow=1,ncol=K) }

    if ( (sum(rest[j,])!=0) & (D==K) & any(as.numeric(Sz)>0)==TRUE ) {Amat2=matrix(c(as.numeric(Amat[,-c(Sz)])), nrow=sum(rest[j,]),ncol=K-1-length(Sz)) }

    if ( (sum(rest[j,])!=0) & (D!=K) & any(as.numeric(Sz)>0)==TRUE ) {Amat2=matrix(c(as.numeric(Amat[,-c(Sz)])), nrow=sum(rest[j,]),ncol=K-length(Sz)) }

    if ( (sum(rest[j,])!=0) & (D==K) & any(as.numeric(Sz)>0)==FALSE ) {Amat2=matrix(c(as.numeric(Amat[,])), nrow=sum(rest[j,]),ncol=K-1) }

    if ( (sum(rest[j,])!=0) & (D!=K) & any(as.numeric(Sz)>0)==FALSE ) {Amat2=matrix(c(as.numeric(Amat[,])), nrow=sum(rest[j,]),ncol=K) }

    if (any(apply(Amat2, 1, function(y) all(y==0)))==FALSE) {bvec2=bvec}


    if (any(apply(Amat2, 1, function(y) all(y==0)))==TRUE)

    {


      Szxx=which(apply(Amat2, 1, function(y) all(y==0)))

      if (all(bvec[Szxx]<0 | bvec[Szxx]==0))

      {

        bvec2=bvec

      }

      if (any(bvec[Szxx]>0)) {

        if (length(Szxx)>1) { qadj=which(apply(Amat[Szxx,], 2, function(y) any(y!=0))) }
        if (length(Szxx)==1) { qadj=which(Amat[Szxx,]!=0)}

        G=diag(1,K)[c(Sz),]
        H=rep(0,length(Sz))
        q=ldei(E=Amat, F=t(bvec), G = G, H = H)$X


        qadj2=Sz[which(Sz!=qadj)]



        ##########################
        fff[j,c(qadj)]=q[qadj]
        fff[j,c(qadj2)]=0
        #########################



        q[qadj]=0


        bvec2=t(Amat%*%q)




      }


    }

    #--------------------------------------------------------------------------------------------------------


    ####Calculating the order-restricited estimates

    solveQP=solve.QP(Dmat/sc, dvec = dvec/sc, t(Amat2), bvec2, meq = rest[j,1], factorized = FALSE)
    tildeeta[j,]=solveQP$solution


    ####Calculating the log likelihoods

    if (D==K) {logHm[,j] =as.numeric( ( -(D-1)/2*log(2*pi) )
                                      -( 0.5*log(det(BKcov) ) )-( 0.5* t(eta- tildeeta[j,])%*%ginv(BKcov)%*%(eta-tildeeta[j,])) )  }

    if (D!=K) {logHm[,j] =as.numeric( ( -K/2*log(2*pi) )
                                      -( 0.5*log(det(BKcov) ) )-( 0.5* t(eta- tildeeta[j,])%*%ginv(BKcov)%*%(eta-tildeeta[j,])) )   }

    #----------------------------------------------------------------------------------


    #Calculating penalty parts

    bvecpen=bvec-bvec
    len = rep(NA, T)

    for (i in 1:T) {

      z=mvrnorm(1,rep(0,nrow(BKcov)),BKcov)


      Dmat2=2*ginv(BKcov)
      dvec2=2*(z[]%*%ginv(BKcov))

      sc = norm(Dmat2,"O")


      solveQP2= solve.QP(Dmat2/sc,dvec2/sc,t(Amat2),bvecpen,meq =rest[j,1] ,factorized = FALSE)
      Ksz=dim(BKcov)[1]


      len[i]=Ksz-length(solveQP2$iact[solveQP2$iact != 0])
    }


    levelprob = matrix(c(0),nrow=1,ncol=Ksz,byrow=TRUE)
    for (i in 1:Ksz) {
      levelprob[1, i]=sum(len[ ]==i)/T}

    PT = matrix(c(0), nrow = 1, ncol = Ksz, byrow = TRUE)
    for (i in 1:Ksz) {
      PT[1,i]=levelprob[1,i]*i }
    PTm[,j] = sum(PT[,])

    ######################################################M loop finishes###############################################
  }


  Log=c(logHm)
  Penalty=c(PTm)

  Goricam=rep(c(NA),M)
  Weightm=rep(c(NA),M)


  ####Calculating GORICA values
  Goricam=-2*Log+2*Penalty


  ####Calculating GORICA weights
  Weightm=(exp(-0.5*Goricam))/(sum(exp(-0.5*Goricam)))

  #---------------------------------------------------------------------------------


  ####Some notes are given in the output depending on different situations


  if ((length(eta)<length(etax))==TRUE & (D!=K)) cat("NOTE: Some of the eta parameters are estimated as zero.\nThe hypotheses under evaluation are rewritten due to empty cell(s).\n \n")


  if ((length(eta)<(length(etax)-1))==TRUE & (D==K)) cat("NOTE: Some of the eta parameters are estimated as zero.\nThe hypotheses under evaluation are rewritten due to empty cell(s).\n \n")


  ####The rest of the code for writting the output properly on R console

  cell=rep(1:D)

  if (Sz==0) {mm=max(cell)}
  if (any(Sz)>0) {mm=max(cell[-c(cell[which(etax==0)])])   }


  listeta2=as.matrix(list(etax)[[1]])
  listeta2=t(listeta2)
  A=paste("p",1:nrow(BKcovx),sep="")
  colnames(listeta2)=c(A)
  print(list(MLEs=listeta2))






  listBKcovx=as.matrix(list(BKcovx))[[1]]
  A=paste("p",1:nrow(BKcovx),sep="")
  colnames(listBKcovx)=c(A)
  rownames(listBKcovx)=c(A)
  print(list(Covariance_matrix=listBKcovx))

  tildeetax=matrix(c(NA),nrow=M,ncol=nrow(BKcovx))


  if (any(as.numeric(Sz)>0)==TRUE & (etax[K]!=0) ) {
    tildeetax[1:M,c(Sz)]=fff[1:M,c(Sz)]
    if (D==K) {tildeetax[1:M,-c(Sz,D)]=tildeeta[,]}
    if (D!=K) {tildeetax[1:M,-c(Sz)]=tildeeta[,]}



    if (D==K) {

      for (i in 1:M) {
        tildeetax[i,D]=1-sum(tildeetax[i,1:(D-1)])


      }




    }


    listtilde=as.matrix(list(tildeetax)[[1]])
    A=paste("H",1:M,sep="")
    B=paste("p",1:nrow(BKcovx),sep="")
    rownames(listtilde)=c(A)
    colnames(listtilde)=B
    print(list(Restricted_MLEs=listtilde))
  }


  ##############################################################################

  if (any(as.numeric(Sz)>0)==TRUE & (etax[K]==0) ) {
    tildeetax[1:M,c(which(etax==0))]=fff[1:M,c(which(etax==0))]
    if (D==K) {tildeetax[1:M,-c(Sz,D)]=tildeeta[,]}
    if (D!=K) {tildeetax[1:M,-c(Sz)]=tildeeta[,]}


    if (D==K) {

      for (i in 1:M) {
        tildeetax[i,c(mm)]=1-sum(na.omit(tildeetax[i,]))

      }

    }


    listtilde=as.matrix(list(tildeetax)[[1]])
    A=paste("H",1:M,sep="")
    B=paste("p",1:nrow(BKcovx),sep="")
    rownames(listtilde)=c(A)
    colnames(listtilde)=B
    print(list(Restricted_MLEs=listtilde))
  }
  ################################################################################


  if (any(as.numeric(Sz)>0)==FALSE) {
    if (D==K) {tildeetax[1:M,-c(D)]=tildeeta[,]}
    if (D!=K) {tildeetax[1:M,]=tildeeta[,]}


    if (D==K) {

      for (i in 1:M) {
        tildeetax[i,D]=1-sum(tildeetax[i,1:(D-1)])
      }

    }
    listtilde=as.matrix(list(tildeetax)[[1]])
    A=paste("H",1:M,sep="")
    B=paste("p",1:nrow(BKcovx),sep="")
    rownames(listtilde)=c(A)
    colnames(listtilde)=B
    print(list(Restricted_MLEs=listtilde))
  }
  #####################################################################################


  if (D == K) cat("NOTE: The",mm,". order-restricted MLE is defined as one minus the others, \n      because of the linear dependency among cell probabilities. \n \n")


  if (!all(listtilde>-1e-10)==TRUE) cat("WARNING: A negative value is encountered for the order-restricted MLEs. \nPlease check the order-restricted MLEs and inspect the hypothesis for which the order-restricted MLE(s) are negative.\nThe restrictions of the corresponding hypothesis might not be reasonable. \n \n")


  list1=as.matrix(list(Log)[[1]])
  list2=as.matrix(list(Penalty)[[1]])
  list3=as.matrix(list(Goricam)[[1]])
  list4=as.matrix(list(Weightm)[[1]])
  listt=cbind(list1,list2,list3,list4)
  A=paste("H",1:M,sep="")
  rownames(listt)=c(A)
  colnames(listt)=c("LogLikelihood","Penalty", "GORICA_value", "GORICA_weight")
  print(listt)
  #------------------------------------------------------------------------
  return(list(estimates= etax, gorica_weights = listt[, 4]))

}



# Tests -------------------------------------------------------------------

restore_data <- function(datfile){
  dataset=read.table(datfile,header=FALSE)
  out <- matrix(dataset$V3, nrow = max(dataset$V1), ncol = max(dataset$V2), byrow = TRUE)
  class(out) <- c("table", class(out))
  out
}

datf <- system.file("testdata", "data_1.txt", package="gorica")
inpf <- system.file("testdata", "input_1.txt", package="gorica")

# readLines(inpf)
doei <- capture_output({
  res_oldgor <- old_goricacont(inpfile = inpf,
                               datfile = datf)
})

tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "a:=x[2,1]/(x[1,1]+x[2,1]);b:=x[2,2]/(x[1,2]+x[2,2]);c:=x[2,3]/(x[1,3]+x[2,3]);d:=x[2,4]/(x[1,4]+x[2,4]);a > (b,c,d); a = b & c > d;a >b & b > c & c > d")

test_that("Estimates close", expect_equivalent(res$estimates, res_oldgor$estimates, tolerance = .01))

test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .05))



# Example 2 ---------------------------------------------------------------

datf <- system.file("testdata", "data_2.txt", package="gorica")
doei <- capture_output({
  res_oldgor <- old_goricacont(inpfile = system.file("testdata", "input_2.txt", package="gorica"),
                               datfile = datf)
})

tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[1,1]=x[1,2] & x[2,1]>x[2,2];x[1,1]>x[1,2] & x[2,1]>x[2,2]")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))


# Example 3 ---------------------------------------------------------------

datf <- system.file("testdata", "data_3.txt", package="gorica")
res_oldgor <- old_goricacont(inpfile = system.file("testdata", "input_3.txt", package="gorica"),
                             datfile = datf)


tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[1,1]=x[1,2] & x[2,1]>x[2,2];x[1,1]>x[1,2] & x[2,1]>x[2,2]")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .013))



# Example 4 ---------------------------------------------------------------

datf <- system.file("testdata", "data_4.txt", package="gorica")

expect_error(old_goricacont(inpfile = system.file("testdata", "input_4.txt", package="gorica"), datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=x[1,1]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);b:=x[1,2]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);c:=x[1,3]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);d:=x[1,4]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);a > (b,c,d); a = b & c > d;a >b & b > c & c > d"))


# Example 5 ---------------------------------------------------------------

datf <- system.file("testdata", "data_5.txt", package="gorica")
expect_error(old_goricacont(inpfile = system.file("testdata", "input_5.txt", package="gorica"), datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);b:=(x[1,2]*x[2,4])/(x[1,3]*x[2,2]);c:=(x[1,3]*x[2,4])/(x[1,4]*x[2,3]);a>b&b>c"))


# Example 6 ---------------------------------------------------------------

datf <- system.file("testdata", "data_6.txt", package="gorica")
expect_error(old_goricacont(inpfile = system.file("testdata", "input_6.txt", package="gorica"), datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=x[2,1]/(x[1,1]+x[2,1]);b:=x[2,2]/(x[1,2]+x[2,2]);c:=x[2,3]/(x[1,3]+x[2,3]);d:=x[2,4]/(x[1,4]+x[2,4]);a>(b,c,d)"))



# Example 9 ---------------------------------------------------------------

datf <- system.file("testdata", "data_9.txt", package="gorica")
inpf <- system.file("testdata", "input_9.txt", package="gorica")
res_oldgor <- old_goricacont(inpfile = system.file("testdata", "input_9.txt", package="gorica"),
                             datfile = datf)

#readLines(datf)
#readLines(inpf)

tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[2,1],x[2,2],x[2,3],x[2,4],x[1,1],x[1,2],x[1,3],x[1,4];x[2,1]>(x[2,2],x[2,3],x[2,4])&x[2,1]>0.3")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))


# Example 10 ---------------------------------------------------------------

datf <- system.file("testdata", "data_10.txt", package="gorica")
inpf <- system.file("testdata", "input_10.txt", package="gorica")
expect_error(old_goricacont(inpfile = inpf,
                            datfile = datf))

#readLines(datf)
#readLines(inpf)

tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);a<1;a=1;a>1"))


# Example 11 ---------------------------------------------------------------

datf <- system.file("testdata", "data_11.txt", package="gorica")
inpf <- system.file("testdata", "input_11.txt", package="gorica")
expect_error(old_goricacont(inpfile = inpf,
                            datfile = datf))

#readLines(datf)
#readLines(inpf)

tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);a<1;a=1;a>1"))


# Example 12 ---------------------------------------------------------------

datf <- system.file("testdata", "data_12.txt", package="gorica")
inpf <- system.file("testdata", "input_12.txt", package="gorica")
expect_error(old_goricacont(inpfile = inpf,
                            datfile = datf))

#readLines(datf)
#readLines(inpf)

tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[1,4])/(x[1,2]*x[1,3]);b:=(x[2,1]*x[2,4])/(x[2,2]*x[2,3]);a<1&b<1;a=1&b=1;a>1&b>1", comparison = "none"))

