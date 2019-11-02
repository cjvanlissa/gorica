# BELOW THE T.TEST INPUT FOR BAIN IS TESTED

# ===============================================================================================

# THE ONE SAMPLE T-TEST WITH A T.TEST OBJECT



x<-sesamesim$postnumb
ttest <- t_test(x)
set.seed(100)
z <- bain(ttest, "x=30; x>30; x<30")
z_gor <- gorica(ttest, "x=30; x>30; x<30")

x<-sesamesim$postnumb[which(sesamesim$sex==1)]
y<-sesamesim$postnumb[which(sesamesim$sex==2)]
ttest <- t_test(x,y,paired = FALSE, var.equal = FALSE)
set.seed(100)
z <- bain(ttest, "x=y; x>y; x<y")
z_gor <- gorica(ttest, "x=y; x>y; x<y")


# THE INDEPENDENT GROUPS T-TEST WITH A T.TEST OBJECT

x<-sesamesim$postnumb[which(sesamesim$sex==1)]
y<-sesamesim$postnumb[which(sesamesim$sex==2)]
ttest <- t_test(x,y,paired = FALSE, var.equal = TRUE)
set.seed(100)
z <- bain(ttest, "x=y; x>y; x<y")
z_gor <- gorica(ttest, "x=y; x>y; x<y")


sesamesim$sex<-as.factor(sesamesim$sex)
ttest <- bain:::t_test_old.formula(postnumb~sex,data=sesamesim,paired = FALSE, var.equal = TRUE)
ttest <- t_test(postnumb~sex,data=sesamesim,paired = FALSE, var.equal = TRUE)
set.seed(100)
zh<-bain(ttest, "group1=group2; group1>group2; group1<group2")
zh_gor <- gorica(ttest, "group1=group2; group1>group2; group1<group2")


test_that("t_test formula and normal interface same", {
  expect_equivalent(z_gor$fit$gorica_weights, zh_gor$fit$gorica_weights, tolerance = .1)
})

# =================================================================================================

# THE PAIRED SAMPLES T-TEST WITH A T.TEST OBJECT



x<-sesamesim$prenumb
y<-sesamesim$postnumb

ttest <- t_test(x,y,paired = TRUE)
set.seed(100)
z <- bain(ttest, "difference=0; difference>0; difference<0")
z_gor <- gorica(ttest, "difference=0; difference>0; difference<0")

test_that("paired t_test bain and gorica similar", {
  expect_equivalent(z$fit$PMPb, z_gor$fit$gorica_weights, tolerance = .1)
})

#==================================================================================================

# THE EQUIVALENCE TEST WITH A T.TEST OBJECT



x<-sesamesim$postnumb[which(sesamesim$sex==1)]
y<-sesamesim$postnumb[which(sesamesim$sex==2)]

ttest <- t_test(x,y,paired = FALSE, var.equal = TRUE)
set.seed(100)
z <- bain(ttest, "x - y > -1 & x - y < 1")
z_gor <- gorica(ttest, "x - y > -1 & x - y < 1")
z_gor

test_that("equal variance t_test gorica bain similar", {
  expect_true(z$fit$PMPb[1] > z$fit$PMPb[2] & z_gor$fit$gorica_weights[1] > z_gor$fit$gorica_weights[2])
})
