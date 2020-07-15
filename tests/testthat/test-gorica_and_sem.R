if(FALSE){


library(gorica)
library(lavaan)

#  Specify and fit the confirmatory factor model
model1 <- '
    A =~ Ab + Al + Af + An + Ar + Ac
    B =~ Bb + Bl + Bf + Bn + Br + Bc
'

# Use the lavaan sem function to execute the confirmatory factor analysis
fit1 <- sem(model1, data = sesamesim, std.lv = TRUE)

# Inspect the parameter names (such that hypotheses can be specified using these)
coef(fit1)

# Formulate hypotheses
hypotheses1 <- "
A=~Ab > .6 & A=~Al > .6 & A=~Af > .6 & A=~An > .6 & A=~Ar > .6 & A=~Ac > .6 &
B=~Bb > .6 & B=~Bl > .6 & B=~Bf > .6 & B=~Bn > .6 & B=~Br > .6 & B=~Bc > .6
"

# Call gorica
set.seed(100)
out_c <- gorica(fit1, hypotheses1, comparison = "complement", standardize=TRUE) # TO DO doet laatste iets??
out_c
# TO DO Krijg nu een negatieve PT voor het complemet (zou iets tussen 0 en 12 moeten zijn, ws dichtbij 12)
# NB In restriktor gaat dit wel goed - zie regressie voorbeeld hieronder.
out_u <- gorica(fit1, hypotheses1, standardize=TRUE) # TO DO doet laatste iets??
out_u



library(restriktor)
#
# Generate data
n <- 10
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rnorm(n)
x5 <- rnorm(n)
x6 <- rnorm(n)
y <- 1 + x1 + x2 + x3 + x4 + x5 + x6 + rnorm(n)
#
# Fit regression model
fit.lm <- lm(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6)
#
#1.	Fitted unconstrained (lm or glm) object + character constraints
h1 <- "x1 > 0; x2 > 0; x3 > 0; x4 > 0; x5 > 0; x6 > 0"
h2 <- "x1 > 0.6; x2 > 0.6; x3 > 0.6; x4 > 0.6; x5 > 0.6; x6 > 0.6"
#
# goric
out_1_u <- goric(fit.lm, h1)
summary(out_1_u)
#
out_1_c <- goric(fit.lm, h1, comparison = "complement")
summary(out_1_c)
#
out_2_u <- goric(fit.lm, h2)
summary(out_2_u)
#
out_2_c <- goric(fit.lm, h2, comparison = "complement")
summary(out_2_c)
}
