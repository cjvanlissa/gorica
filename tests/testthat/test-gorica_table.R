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

test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .05))



# Example 2 ---------------------------------------------------------------

datf <- "dev/examples/example 2/data.txt"
res_oldgor <- old_goricacont(inpfile = "dev/examples/example 2/input.txt",
                             datfile = datf)


tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[1,1]=x[1,2] & x[2,1]>x[2,2];x[1,1]>x[1,2] & x[2,1]>x[2,2]")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))


# Example 3 ---------------------------------------------------------------

datf <- "dev/examples/example 3/data.txt"
res_oldgor <- old_goricacont(inpfile = "dev/examples/example 3/input.txt",
                             datfile = datf)


tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[1,1]=x[1,2] & x[2,1]>x[2,2];x[1,1]>x[1,2] & x[2,1]>x[2,2]")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .013))



# Example 4 ---------------------------------------------------------------

datf <- "dev/examples/example 4/data.txt"

expect_error(old_goricacont(inpfile = "dev/examples/example 4/input.txt", datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=x[1,1]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);b:=x[1,2]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);c:=x[1,3]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);d:=x[1,4]/(x[1,1]+x[1,2]+x[1,3]+x[1,4]);a > (b,c,d); a = b & c > d;a >b & b > c & c > d"))


# Example 5 ---------------------------------------------------------------

datf <- "dev/examples/example 5/data.txt"
expect_error(old_goricacont(inpfile = "dev/examples/example 5/input.txt", datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);b:=(x[1,2]*x[2,4])/(x[1,3]*x[2,2]);c:=(x[1,3]*x[2,4])/(x[1,4]*x[2,3]);a>b&b>c"))


# Example 6 ---------------------------------------------------------------

datf <- "dev/examples/example 6/data.txt"
expect_error(old_goricacont(inpfile = "dev/examples/example 6/input.txt", datfile = datf))


tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=x[2,1]/(x[1,1]+x[2,1]);b:=x[2,2]/(x[1,2]+x[2,2]);c:=x[2,3]/(x[1,3]+x[2,3]);d:=x[2,4]/(x[1,4]+x[2,4]);a>(b,c,d)"))



# Example 9 ---------------------------------------------------------------

# TO DO

datf <- "dev/examples/example 9/data.txt"
inpf <- "dev/examples/example 9/input.txt"
res_oldgor <- old_goricacont(inpfile = "dev/examples/example 9/input.txt",
                             datfile = datf)

readLines(datf)
readLines(inpf)

tab <- restore_data(datf)

res <- gorica(tab, hypothesis = "x[2,1]>(x[2,2],x[2,3],x[2,4])&x[2,1]>0.3")

res <- gorica(tab, hypothesis = "x[2,1],x[2,2],x[2,3],x[2,4],x[1,1],x[1,2],x[1,3],x[1,4];x[2,1]>(x[2,2],x[2,3],x[2,4])&x[2,1]>0.3")

test_that("Estimates close", expect_equivalent(matrix(res$estimates, nrow = 2), matrix(res_oldgor$estimates, nrow = 2, byrow = T), tolerance = .015))


test_that("Weights close", expect_equivalent(res$fit$gorica_weights, res_oldgor$gorica_weights, tolerance = .01))


# Example 10 ---------------------------------------------------------------

datf <- "dev/examples/example 10/data.txt"
inpf <- "dev/examples/example 10/input.txt"
expect_error(old_goricacont(inpfile = inpf,
                             datfile = datf))

readLines(datf)
readLines(inpf)

tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);a<1;a=1;a>1"))


# Example 11 ---------------------------------------------------------------

# TO DO

datf <- "dev/examples/example 10/data.txt"
inpf <- "dev/examples/example 10/input.txt"
expect_error(old_goricacont(inpfile = inpf,
                            datfile = datf))

readLines(datf)
readLines(inpf)

tab <- restore_data(datf)

expect_error(gorica(tab, hypothesis = "a:=(x[1,1]*x[2,2])/(x[1,2]*x[2,1]);a<1;a=1;a>1"))

