#Installation of the Relavant Libraries

library(easypackages)
list.of.packages <- c("base", "blavaan", "boot", "coda", "foreign", "FRACTION", "lavaan", "lme4",
"MASS", "matrixcalc", "mvtnorm", "nlme", "quadprog", "R2OpenBUGS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
libraries(list.of.packages)

#Example 2: Multilevel Regression Modeling

########################################################################
#Maximum Likelihood Estimation

set.seed(111)

#Specification of the data set
reading_ach <- read.csv("reading_achievement.csv")

#Specification of the variables in the model
zgeread <- scale(reading_ach$geread)
zage <- scale(reading_ach$age)
zgevocab <- scale(reading_ach$gevocab)
reading_ach$gender <- as.factor(reading_ach$gender)

#Model fitting
model <- lmer(zgeread ~ zage + zgevocab + gender + zage * zgevocab + zage * gender +
zgevocab * gender + (1 + zgevocab | school), data = reading_ach, REML = FALSE)

print("######Maximum Likelihood Estimation######")

strest <- summary(model)$coefficients[c(2,3,5,6,7),1]
strcovmtrx <- vcov(model)[c(2,3,5,6,7), c(2,3,5,6,7)]

print(list(MLEs=strest,Covariance_matrix_of_MLEs=strcovmtrx))

#Obtaining order-restricted estimates using ormle

#Source code in ormle is reached
source("restrictedest.txt")

#Specification of the hypotheses under evaluation

# Hypothesis 1
constr <- matrix(c(1, 0, 0, 0, 0,
1, 0, 0, 1, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 2
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(-1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
-1, 0, 0, -1, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)


# Hypothesis 3
constr <- matrix(c(0, 0, 1, 0, 0,
1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
-1, 0, 0, -1, 0,
0, -1, 0, 0, -1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0,5)
nec <- 1
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 5)), nrow = 1, ncol = 5, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)
# Source code in gorica is reached, which is saved in the file "Gorica.txt"
source("Gorica.txt")



#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))

#######################################################################################################

#Nonparametric bootstrapping

set.seed(111)

#Source code in covglmm is reached
source("glmmcov.txt")


#Specification of the model
covar <- covglmm(lmer(geread ~ age + gevocab + gender + age * gevocab + age * gender +
gevocab * gender + (1 + gevocab | school), data = reading_ach, REML = FALSE),
B = 1000, cluster = c(1), fact = c(2), std = c(3, 4, 5) )

colnames(covar$boot)<-names(summary(model)$coefficients[,1])



print("######Nonparametric bootstrapping######")

strest<- apply(covar$boot, 2, mean)[c(2,3,5,6,7)]
strcovmtrx<-cov(covar$boot)[c(2,3,5,6,7), c(2,3,5,6,7)]

print(list(Bootstrapped_estimates=strest,Covariance_matrix_of_bootstrapped_estimates=strcovmtrx))

#Specification of the hypotheses under evaluation

# Hypothesis 1
constr <- matrix(c(1, 0, 0, 0, 0,
1, 0, 0, 1, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 2
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(-1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
-1, 0, 0, -1, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)


# Hypothesis 3
constr <- matrix(c(0, 0, 1, 0, 0,
1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
-1, 0, 0, -1, 0,
0, -1, 0, 0, -1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0,5)
nec <- 1
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 5)), nrow = 1, ncol = 5, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)


#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))

#######################################################################

#Gibbs Sampling

set.seed(111)


#Introducing the sample size N = 10320 and the variables in the model file
N <- nrow(reading_ach)
school <- reading_ach[, 1]
gender <- model.matrix(~ factor(reading_ach$gender))[, 2]
zage <- as.numeric(scale(reading_ach$age) )
zgevocab <- as.numeric(scale(reading_ach$gevocab) )
zgeread <- as.numeric(scale(reading_ach$geread) )
data <- list("N", "school", "gender", "zage", "zgevocab", "zgeread")
#Initializing the model parameters
beta0 <- c(rnorm(1, 0, 1e-05) )
beta <- c(rnorm(6, 0, 1e-05) )
inits <- function(){
list(beta0 = beta0, beta = beta, prec.sigma2 = 1, prec.tau2 = 1) }

#Performing gibbs sampling 
gibbs.sim <- bugs(data, inits, model.file = "bugsmodelglmm.txt",
parameters = c("beta0", "beta", "prec.sigma2", "prec.tau2"),
n.iter = 10000, n.burnin = 1000, n.chains = 3, debug = FALSE, codaPkg = TRUE)


codaobject <- as.matrix(read.bugs(gibbs.sim) )

colnames(codaobject)<-c("zage","zgevocab","gender2", "zage:zgevocab","zage:gender2","zgevocab:gender2","Intercept","deviance","prec.sigma2","prec.tau2")


print("######Gibbs sampling######")


strest <- apply(codaobject, 2, mean)[c(1,2,4,5,6)]
strcovmtrx <- cov(codaobject)[c(1,2,4,5,6),c(1,2,4,5,6)]

print(list(Gibbs_estimates=strest,Covariance_matrix_of_gibbs_estimates=strcovmtrx))

#Specification of the hypotheses under evaluation

#Hypothesis 1
constr <- matrix(c(1, 0, 0, 0, 0,
1, 0, 0, 1, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 2
H1 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# Hypothesis 2
constr <- matrix(c(-1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
-1, 0, 0, -1, 0,
0, 1, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0, 5)
nec <- 0
H2 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)


# Hypothesis 3
constr <- matrix(c(0, 0, 1, 0, 0,
1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
-1, 0, 0, -1, 0,
0, -1, 0, 0, -1), nrow = 5, ncol = 5, byrow = TRUE)
rhs <- rep(0,5)
nec <- 1
H3 <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

# The unconstrained hypothesis
constr <- matrix(c(rep(0, 5)), nrow = 1, ncol = 5, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0
Hu <- ormle(est = strest, covmtrx = strcovmtrx, const = constr, nec = nec, rhs = rhs)

set.seed(111)


#Performing gorica to obtain the values of misfit, complexity, GORICA, and GORICA weights
print(list(Evaluation_of_the_set_of_hypotheses=gorica(H1, H2, H3, Hu, iter = 100000)))






































