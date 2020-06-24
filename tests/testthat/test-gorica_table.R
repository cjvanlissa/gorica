source("tests/testthat/files/old_goricacont.R")
datf <- "dev/examples/example 1/data.txt"
res_oldgor <- old_goricacont(inpfile = "dev/examples/example 1/input.txt",
                             datfile = datf)

restore_data <- function(datfile){
  dataset=read.table(datfile,header=FALSE)
  out <- matrix(dataset$V3, nrow = max(dataset$V1), ncol = max(dataset$V2), byrow = TRUE)
  class(out) <- c("table", class(out))
  out
}

tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "a:=x[2,1]/(x[1,1]+x[2,1]);b:=x[2,2]/(x[1,2]+x[2,2]);c:=x[2,3]/(x[1,3]+x[2,3]);d:=x[2,4]/(x[1,4]+x[2,4]);a > (b,c,d); a = b & c > d;a >b & b > c & c > d")



test_that("Estimates close", expect_equivalent(res$estimates, res_oldgor$estimates, tolerance = .01))

test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))



# Example 2 ---------------------------------------------------------------

datf <- "dev/examples/example 2/data.txt"
res_oldgor <- old_goricacont(inpfile = "dev/examples/example 2/input.txt",
                             datfile = datf)

tab <- t(restore_data(datf))
est <- get_estimates(tab)
est$Sigma

test_that("Ex 2 estimates equivalent", expect_equivalent(est$estimate, res_oldgor$estimates, tolerance = .01))


res <- gorica(tab, hypothesis = "x[1,1]>x[1,2] & x[2,1]>x[2,2]")

test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))

