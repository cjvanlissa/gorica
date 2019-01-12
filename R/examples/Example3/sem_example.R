#Installation of the Relavant Libraries

library(easypackages)
list.of.packages <- c("base", "blavaan", "boot", "coda", "foreign", "FRACTION", "lavaan", "lme4",
"MASS", "matrixcalc", "mvtnorm", "nlme", "quadprog", "R2OpenBUGS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
libraries(list.of.packages)


#Example 3: Structural Equation Modeling



#Estimation of Structural Parameters and Their Covariance Matrix


########################################################################
#Maximum Likelihood Estimation

#Specification of the data set
wechsler <- read.csv("wechsler.csv")



#Specification of the variables and the model
SEM.model <- '
Cry = ~ y1 + y2 + y3 + y4
Fld = ~ y2 + y3 + y5 + y6 + y7 + y8
Cry ~ edc + age
Fld ~ edc + age
'


#Model fitting
fit <- cfa(SEM.model, data = wechsler, std.lv = TRUE)


##############################################################


print("######Maximum Likelihood Estimation######")
print("###Factor Loadings###")

strest1 <- standardizedSolution(fit)[1:10,4]
strcovmtrx1 <- as.matrix(lavInspect(fit, "vcov.std.all")[1:10,1:10])
names(strest1)<-colnames(strcovmtrx1)


print(list(MLEs=strest1,Covariance_matrix_of_MLEs=strcovmtrx1))

#Obtaining order-restricted estimates using ormle

#Source code in ormle is reached
source("restrictedest.txt")

#Specification of the hypotheses containing factor loadings

# Hypothesis 1
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 10, byrow = TRUE)
rhs <- rep(0, 4)
nec <- 0
H1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)


# Hypothesis 2
constr <- matrix(c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 3, ncol = 10, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H2 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 7, ncol = 10, byrow = TRUE)
rhs <- rep(0, 7)
nec <- 0
H3 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 10)), nrow = 1, ncol = 10, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Evaluation of hypotheses containing factor loadings
set.seed(111)

source("Gorica.txt")

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu1, iter = 100000)))

#########################################################################################

print("######Maximum Likelihood Estimation######")
print("###Regression Coefficients###")

strest2 <- standardizedSolution(fit)[11:14,4]
strcovmtrx2 <- as.matrix(lavInspect(fit, "vcov.std.all")[11:14,11:14])
names(strest2)<-colnames(strcovmtrx2)


print(list(MLEs=strest2,Covariance_matrix_of_MLEs=strcovmtrx2))


#Specification of the hypotheses containing regression coefficients

#Hypothesis 4
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0), nrow = 3, ncol = 4, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H4 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#Hypothesis 5
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0,
0, 1, 0, 0,
0, 0, 0, -1), nrow = 5, ncol = 4, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H5 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#The unconstrained hypothesis
constr <- matrix(c(rep(0, 4)), nrow = 1, ncol = 4, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu2 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

# Evaluation of hypotheses containing regression coefficients
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H4, H5, Hu2, iter = 100000)))

#########################################################################################

print("######Maximum Likelihood Estimation######")
print("###Covariance between latent variables###")

strest3 <- standardizedSolution(fit)[25,4]
strcovmtrx3 <- as.matrix(lavInspect(fit, "vcov.std.all")[23,23])

names(strest3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
rownames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
colnames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]


print(list(MLEs=strest3,Covariance_matrix_of_MLEs=strcovmtrx3))


#Hypothesis 6
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 1
H6 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 7
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H7 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 8
constr <- matrix(c(-1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H8 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Evaluation of hypotheses containing correlation between latent variables
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H6, H7, H8, iter = 100000)))

###########################################################################################

#Nonparametric bootstrapping

set.seed(111)


#Performing nonparametric bootstrapping
boot <- bootstrapLavaan(fit, R = 1000, type = "nonparametric", FUN = function(x) {
standardizedSolution(x)$est }, verbose = TRUE, warn = TRUE)


print("######Nonparametric bootstrapping######")
print("###Factor Loadings###")

strest1 <- apply(boot, 2, mean)[1:10]
strcovmtrx1 <- cov(boot)[1:10,1:10]

names(strest1)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[1:10,1:10]))
rownames(strcovmtrx1)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[1:10,1:10]))
colnames(strcovmtrx1)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[1:10,1:10]))

print(list(Bootstrapped_estimates=strest1,Covariance_matrix_of_bootstrapped_estimates=strcovmtrx1))


#Specification of the hypotheses containing factor loadings

# Hypothesis 1
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 10, byrow = TRUE)
rhs <- rep(0, 4)
nec <- 0
H1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)


# Hypothesis 2
constr <- matrix(c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 3, ncol = 10, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H2 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 7, ncol = 10, byrow = TRUE)
rhs <- rep(0, 7)
nec <- 0
H3 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 10)), nrow = 1, ncol = 10, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Evaluation of hypotheses containing factor loadings
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu1, iter = 100000)))

###############################################################################################

print("######Nonparametric bootstrapping######")
print("###Regression coefficients###")

strest2 <- apply(boot, 2, mean)[11:14]
strcovmtrx2 <- cov(boot)[11:14,11:14]

names(strest2)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[11:14,11:14]))
rownames(strcovmtrx2)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[11:14,11:14]))
colnames(strcovmtrx2)<-colnames(as.matrix(lavInspect(fit, "vcov.std.all")[11:14,11:14]))


print(list(Bootstrapped_estimates=strest2,Covariance_matrix_of_bootstrapped_estimates=strcovmtrx2))


#Specification of the hypotheses containing regression coefficients

#Hypothesis 4
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0), nrow = 3, ncol = 4, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H4 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#Hypothesis 5
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0,
0, 1, 0, 0,
0, 0, 0, -1), nrow = 5, ncol = 4, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H5 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#The unconstrained hypothesis
constr <- matrix(c(rep(0, 4)), nrow = 1, ncol = 4, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu2 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#Evaluation of hypotheses containing regression coefficients
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H4, H5, Hu2, iter = 100000)))
###########################################################################################################

print("######Nonparametric bootstrapping######")
print("###Covariance between latent variables###")

strest3 <- apply(boot, 2, mean)[25]
strcovmtrx3 <- as.matrix(cov(boot)[23,23])

names(strest3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
rownames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
colnames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]

print(list(Bootstrapped_estimates=strest3,Covariance_matrix_of_bootstrapped_estimates=strcovmtrx3))

#Hypothesis 6
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 1
H6 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 7
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H7 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 8
constr <- matrix(c(-1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H8 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Evaluation of hypotheses containing correlation between latent variables
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H6, H7, H8, iter = 100000)))

###################################################################################################
#Gibbs sampling

set.seed(111)

#Performing gibbs sampling 
fit <- bsem(SEM.model, data = wechsler, n.chains = 3, burnin = 3000, sample = 27000,
std.lv = TRUE)

###################################################################################################
print("######Gibbs sampling######")
print("###Factor Loadings###")

strest1 <- standardizedSolution(fit)[1:10,4]
strcovmtrx1 <- lavInspect(fit, "vcov.std.all")[1:10, 1:10]

names(strest1) <- colnames(strcovmtrx1)


print(list(Gibbs_estimates=strest1,Covariance_matrix_of_gibbs_estimates=strcovmtrx1))

#Specification of the hypotheses containing factor loadings

# Hypothesis 1
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 10, byrow = TRUE)
rhs <- rep(0, 4)
nec <- 0
H1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)


# Hypothesis 2
constr <- matrix(c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 3, ncol = 10, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H2 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, -1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, -1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, -1), nrow = 7, ncol = 10, byrow = TRUE)
rhs <- rep(0, 7)
nec <- 0
H3 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 10)), nrow = 1, ncol = 10, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu1 <- ormle(est = strest1, covmtrx = strcovmtrx1, const = constr, nec = nec, rhs = rhs)

# Evaluation of hypotheses containing factor loadings
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu1, iter = 100000)))

####################################################################################################

print("######Gibbs sampling######")
print("###Regression coefficients###")

strest2 <- standardizedSolution(fit)[11:14,4]
strcovmtrx2 <- lavInspect(fit, "vcov.std.all")[11:14,11:14]
names(strest2)<-colnames(strcovmtrx2)

print(list(Gibbs_estimates=strest2, Covariance_matrix_of_gibbs_estimates=strcovmtrx2))

#Hypothesis 4
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0), nrow = 3, ncol = 4, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 0
H4 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#Hypothesis 5
constr <- matrix(c(1, 0, 0, 0,
0, 0, 1, 0,
1, 0, -1, 0,
0, 1, 0, 0,
0, 0, 0, -1), nrow = 5, ncol = 4, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H5 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

#The unconstrained hypothesis
constr <- matrix(c(rep(0, 4)), nrow = 1, ncol = 4, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu2 <- ormle(est = strest2, covmtrx = strcovmtrx2, const = constr, nec = nec, rhs = rhs)

# Evaluation of hypotheses containing regression coefficients
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H4, H5, Hu2, iter = 100000)))

######################################################################################################

print("######Gibbs sampling######")
print("###Covariance between latent variables###")

strest3 <- standardizedSolution(fit)[25,4]
strcovmtrx3 <- as.matrix(lavInspect(fit, "vcov.std.all")[23,23])

names(strest3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
rownames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]
colnames(strcovmtrx3)<-colnames(lavInspect(fit, "vcov.std.all"))[23]

print(list(Gibbs_estimates=strest3,Covariance_matrix_of_gibbs_estimates=strcovmtrx3))

#Hypothesis 6
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 1
H6 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 7
constr <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H7 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Hypothesis 8
constr <- matrix(c(-1), nrow = 1, ncol = 1, byrow = TRUE)
rhs <- rep(0,1)
nec <- 0
H8 <- ormle(est = strest3, covmtrx = strcovmtrx3, const = constr, nec = nec, rhs = rhs)

#Evaluation of hypotheses containing correlation between latent variables
set.seed(111)

print(list(Evaluation_of_the_set_of_hypotheses=gorica(H6, H7, H8, iter = 100000)))


