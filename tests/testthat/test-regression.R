

regr <- lm(postnumb ~ prenumb + funumb + peabody, sesamesim)
regr$call$formula
# UNSTANDARDIZED REGRESSION USING AN LM OBJECT
set.seed(100)
z<-bain(regr,"pre=fu=pea;pea > fu > pre; pre>fu>pea", standardize = FALSE)
z_gor<-gorica(regr,"pre=fu=pea;pea > fu > pre; pre>fu>pea", standardize = FALSE)

test_that("gorica and bain give similar results for regression", {
  expect_equivalent(z_gor$fit$gorica_weights, z$fit$PMPb, tolerance = .07)
})

# STANDARDIZED REGRESSION USING AN LM OBJECT


regr <- lm(postnumb ~ prenumb + funumb + peabody, sesamesim)
set.seed(100)
#sz<-bain(regr,"pre=fu=pea;pea > fu > pre; pre>fu>pea", standardize = TRUE)
test_that("Warning when lm has argument standardize", {
  expect_warning(gorica(regr,"pre=fu=pea;pea > fu > pre; pre>fu>pea", standardize = TRUE))})


# REGRESSION WITH THE INTERCEPT INCLUDED IN THE RESTRICTIONS



regr <- lm(postnumb ~ prenumb + peabody, sesamesim)
coef(regr)
set.seed(100)
sz<-bain(regr,"Int>5 & pre > pea", standardize = FALSE)
sz_gor<-gorica(regr,"Int>5 & pre > pea", standardize = FALSE)

test_that("gorica and bain give similar results for regression", {
  expect_equivalent(sz_gor$fit$gorica_weights, sz$fit$PMPb, tolerance = .01)
})
