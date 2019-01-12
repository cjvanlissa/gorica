#Installation of the Relavant Libraries

list.of.packages <- c("base", "blavaan", "boot", "coda", "foreign", "FRACTION", "lavaan", "lme4",
"MASS", "matrixcalc", "mvtnorm", "nlme", "quadprog", "R2OpenBUGS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
library(base)
library(blavaan)
library(boot)
library(coda)
library(foreign)
library(FRACTION)
library(lavaan)
library(lme4)
library(MASS)
library(matrixcalc)
library(mvtnorm)
library(nlme)
library(quadprog)
library(R2OpenBUGS)


#Example 1: Poisson Regression Modeling

########################################################################
#Maximum Likelihood Estimation


#Specification of the data set

academic_awards <- read.csv("R/examples/Example1/academic_awards.csv")

#Specification of the variables in the model
academic_awards <- within(academic_awards, {
prog <- factor(prog, levels = 1:3, labels = c("General", "Academic", "Vocational") ) } )
zmath <- scale(academic_awards$math)

#Model fitting
model <- glm(num_awards ~ prog + zmath + prog * zmath, family = "poisson",
data = academic_awards)
class(model)

print("######Maximum Likelihood Estimation######")

strest <- model$coefficients[c(4,5,6)]
strest
strcovmtrx <- vcov(model)[c(4,5,6), c(4,5,6)]

print(list(MLEs=strest,Covariance_matrix_of_MLEs=strcovmtrx))



#Obtaining order-restricted estimates using ormle


#Source code in ormle is reached
#source("restrictedest.txt")


#Specification of the hypotheses under evaluation

# Hypothesis 1
constr <- matrix(c(1, 0, 0,
1, 1, 0,
1, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 3
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(0, 1, -1,
0, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
rhs <- rep(0, 2)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 0, -1), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 3)), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)
# Source code in gorica is reached, which is saved in the file "Gorica.txt"
#source("Gorica.txt")

gorica(H1, H2, H3, Hu, iter = 100000)

#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))

########################################################################################

#Nonparametric bootstrapping


set.seed(111)


#Specification of the model
boot.fn <- function(data, index) {
return(coef(glm(num_awards ~ prog + zmath + prog * zmath, family = c("poisson"),
data = academic_awards, subset = index) ) ) }


#Performing nonparametric bootstrapping
boot_sim <- boot(academic_awards, boot.fn, R = 1000, sim = "ordinary")
colnames(boot_sim$t) <- names(boot_sim$t0)


print("######Nonparametric bootstrapping######")


strest <- apply(boot_sim$t, 2, mean)[c(4,5,6)]
strcovmtrx <- cov(boot_sim$t)[c(4,5,6), c(4,5,6)]


print(list(Bootstrapped_estimates=strest,Covariance_matrix_of_bootstrapped_estimates=strcovmtrx))




#Obtaining order-restricted estimates using ormle
save(academic_awards, file="data/academic_awards.RData")

#Source code in ormle is reached
source("restrictedest.txt")


#Specification of the hypotheses under evaluation

# Hypothesis 1
constr <- matrix(c(1, 0, 0,
1, 1, 0,
1, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 3
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(0, 1, -1,
0, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
rhs <- rep(0, 2)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 0, -1), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 3)), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)
# Source code in gorica is reached, which is saved in the file "Gorica.txt"
source("Gorica.txt")



#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))

########################################################################################

#Gibbs Sampling

#Introducing the sample size N = 200 and the variables in the model file (see bugsmodelglm.txt)
N <- nrow(academic_awards)
prog <- academic_awards$prog
academic <- model.matrix(~ prog)[, 2]
vocational <- model.matrix(~ prog)[, 3]
zmath <- as.numeric(zmath)
num_awards <- as.numeric(academic_awards$num_awards)
data <- list("N", "academic", "vocational", "zmath", "num_awards")


#Initializing the model parameters
beta0 <- c(rnorm(1, 0, 1e-05))
beta <- c(rnorm(5, 0, 1e-05))
inits <- function(){
list(beta0 = beta0, beta = beta) }


#Performing gibbs sampling
gibbs.sim <- bugs(data, inits, model.file = "bugsmodelglm.txt", parameters = c("beta0", "beta"),
n.chains = 3, n.iter = 30000, n.burnin = 3000, debug = FALSE, codaPkg = TRUE)


#Estimates of the structural parameters
codaobject <- as.matrix(read.bugs(gibbs.sim) )

colnames(codaobject) <- c("progAcademic", "progVocational", "zmath", "progAcademic:zmath", "progVocational:zmath", "Intercept", "deviance")



print("######Gibbs sampling######")



strest <- apply(codaobject, 2, mean)[c(3,4,5)]
strcovmtrx <- cov(codaobject)[c(3,4,5), c(3,4,5)]

print(list(Gibbs_estimates=strest,Covariance_matrix_of_gibbs_estimates=strcovmtrx))





#Obtaining order-restricted estimates using ormle


#Source code in ormle is reached
source("restrictedest.txt")


#Specification of the hypotheses under evaluation

# Hypothesis 1
constr <- matrix(c(1, 0, 0,
1, 1, 0,
1, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
rhs <- rep(0, 3)
nec <- 3
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(0, 1, -1,
0, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
rhs <- rep(0, 2)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 3
constr <- matrix(c(0, 0, -1), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 3)), nrow = 1, ncol = 3, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)
# Source code in gorica is reached, which is saved in the file "Gorica.txt"
source("Gorica.txt")



#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))





