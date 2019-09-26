#' Bayes factors for informative hypotheses
#'
#' \code{gorica} is an acronym for "Bayesian informative hypotheses evaluation".
#' It uses the Bayes factor to evaluate hypotheses specified using equality and
#' inequality constraints among (linear combinations of) parameters in a wide
#' range of statistical models. A tutorial is provided by Hoijtink, Mulder,
#' van Lissa, and Gu (2018) retrievable from the Psychological Methods website
#' at \url{https://www.apa.org/pubs/journals/met/} or the gorica website at
#' \url{https://informative-hypotheses.sites.uu.nl/software/gorica/} Users are
#' advised to read this tutorial before using \code{gorica}.
#'
#' @param x An R object containing the outcome of a statistical analysis.
#' Currently, the following objects can be processed:
#' \itemize{
#' \item \code{lm()} objects (anova, ancova, multiple regression). In the
#' Example section it is elaborated which calls to \code{lm} can be processed by
#' \code{gorica} (other calls cannot be processed).
#' \item \code{t.test()} objects (Student's t-test, Welch's t-test, paired
#' samples t-test, one-group t-test, equivalence test).
#' \item A named vector containing the estimates resulting from a statistical
#' analysis.
#' Note that, named means that each estimate has to be labelled such that it can
#' be referred to in \code{hypotheses}.
#' }
#' @param hypothesis A character string containing the informative hypotheses to
#' evaluate (see Details).
#' @param ... Additional arguments (see Details).
#'
#' @details
#' \code{hypotheses} is a character string that specifies which informative
#' hypotheses have to be evaluated. A simple example is \code{hypotheses <- "a >
#' b > c; a = b = c;"} which specifies two hypotheses using three estimates with
#' names "a", "b", and "c", respectively.
#'
#' The hypotheses specified have to adhere to the following rules:
#' \enumerate{
#' \item Parameters are referred to using the names specified in \code{names()}.
#' \item Linear combinations of parameters must be specified adhering to the
#' following rules:
#'         \enumerate{ \item Each parameter name is used at most once.
#'                     \item Each parameter name may or may not be
#'                     pre-multiplied with a number.
#'                     \item A constant may be added or subtracted from each
#'                     parameter name.
#'                     \item A linear combination can also be a single number.}
#' Examples are: \code{3 * a + 5}; \code{a + 2 * b + 3 * c - 2}; \code{a - b};
#' and \code{5}.
#' \item (Linear combinations of) parameters can be constrained using <, >, and
#' =. For example, \code{a > 0} or
#' \code{a > b = 0} or \code{2 * a < b + c > 5}.
#' \item The ampersand & can be used to combine different parts of a hypothesis.
#' For example, \code{a > b & b > c} which is equivalent to \code{a > b > c} or
#' \code{a > 0 & b > 0 & c > 0}.
#' \item Sets of (linear combinations of) parameters subjected to the same
#' constraints can be specified using (). For
#' example, \code{a > (b,c)} which is equivalent to \code{a > b & a > c}.
#' \item The specification of a hypothesis is completed by typing ; For example,
#' \code{hypotheses <- "a > b > c; a = b = c;"}, specifies two hypotheses.
#' \item Hypotheses have to be compatible, non-redundant and possible. What
#' these terms mean will be elaborated below.
#' }
#'
#' \emph{The set of hypotheses has to be compatible}. For the statistical
#' background of this requirement see Gu, Mulder, Hoijtink (2018). Usually the
#' sets of hypotheses specified by researchers are compatible, and if not,
#' \code{gorica} will return an error message. The following steps can be used to
#' determine if a set of hypotheses is compatible:
#' \enumerate{
#' \item	Replace a range constraint, e.g., \code{1 < a1 < 3}, by an equality
#' constraint in which the parameter involved is equated to the midpoint of the
#' range, that is, \code{a1 = 2}.
#' \item Replace in each hypothesis the < and > by =. For example, \code{a1 = a2
#' > a3 > a4} becomes \code{a1 = a2 = a3 = a4}.
#' \item The hypotheses are compatible if there is at least one solution to the
#' resulting set of equations. For the two hypotheses considered under 1. and
#' 2., the solution is a1 = a2 = a3 = a4 = 2. An example of two non-compatible
#' hypotheses is \code{hypotheses <- "a = 0; a > 2;"} because there is no
#' solution to the equations \code{a=0} and \code{a=2}.
#' }
#'
#' \emph{Each hypothesis in a set of hypotheses has to be non-redundant.} A
#' hypothesis is redundant if it can also be specified with fewer constraints.
#' For example, \code{a = b & a > 0 & b > 0} is redundant because it can also be
#' specified as \code{a = b & a > 0}. \code{gorica} will work correctly if
#' hypotheses specified using only < and > are redundant. \code{gorica} will
#' return an error message if hypotheses specified using at least one = are
#' redundant.
#'
#' \emph{Each hypothesis in a set of hypotheses has to be possible.} An
#' hypothesis is impossible if estimates in agreement with the hypothesis do not
#' exist. For example: values for \code{a} in agreement with \code{a = 0 &
#' a > 2} do not exist. It is the responsibility of the user to ensure that the
#' hypotheses specified are possible. If not, \code{gorica} will either return an
#' error message or render an output table containing \code{Inf}'s.
#'
#' @return Bla
#' @author Caspar van Lissa, Yassin Altinisik, Rebecca Kuiper
#' @references All references are retrievable via or from
#' \url{https://informative-hypotheses.sites.uu.nl/software/gorica/}
#'
#' @seealso The gorica website at
#' @examples
#' print("bla")
#' @rdname gorica
#' @export
#' @importFrom stats as.formula coef complete.cases cov lm model.frame
#' model.matrix pt qt sd setNames summary.lm var vcov
#'
gorica <- function(x, hypothesis, comparison = "unconstrained", ...) {
  UseMethod("gorica", x)
}

#' @method gorica default
#' @export
gorica.default <- function(x,
                           hypothesis,
                           comparison = "unconstrained",
                           Sigma,
                           ...
)
{
  cl <- match.call()
  # Parse hypotheses --------------------------------------------------------
  #ren_estimate <- rename_estimate(x)
  if(!inherits(comparison, "character")|length(comparison) > 1){
    stop("Argument 'comparison' must be an atomic character string.")
  } else {
    comp_arg <- pmatch(comparison, c("unconstrained", "complement", "none"))
    if(is.na(comp_arg)) stop("Argument 'comparison' did not match one of the available options: 'unconstrained', 'complement', or 'none'.")
    comparison <- c("unconstrained", "complement", "none")[pmatch(comparison, c("unconstrained", "complement", "none"))]
  }
#browser()
  if(inherits(hypothesis, "character")){
    hyp_params <- params_in_hyp(hypothesis)
    coef_in_hyp <- charmatch(rename_function(hyp_params),
                             rename_function(names(x)))

    if(anyNA(coef_in_hyp)){
      stop("Some of the parameters referred to in the 'hypothesis' do not correspond to parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' did not match any parameters in 'x': ",
           paste(hyp_params[is.na(coef_in_hyp)], collapse = ", "),
           "\n  The parameters in object 'x' are named: ",
           paste(names(x), collapse = ", "))
    }
    if(any(coef_in_hyp == 0)){
      stop("Some of the parameters referred to in the 'hypothesis' matched multiple parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' matched multiple parameters in 'x': ",
           paste(hyp_params[coef_in_hyp == 0], collapse = ", "),
           "\n  The parameters in object 'x' are named: ",
           paste(names(x), collapse = ", "))
    }
    hypothesis <- parse_hypothesis(names(x), hypothesis)
  } else {
    if(inherits(hypothesis, "list") & !is.null(hypothesis[["hyp_mat"]]) & !is.null(hypothesis[["n_ec"]])){
      hypothesis$original_hypothesis <- matrix_to_hyp(hypothesis, names(x))
    } else {
      stop("Argument 'hypothesis' must either be a character string with inequality constraints, or a list with an element 'hyp_mat', consisting of a list of hypothesis matrices, and and element 'n_ec', consisting of an integer vector with the number of equality constraints for each hypothesis matrix in 'hyp_mat'.")
    }
  }
  #browser()
  hypotheses <- mapply(function(this_hyp, nec_num){
    ormle(x,
          Sigma,
          constr = this_hyp[, -ncol(this_hyp), drop = FALSE],
          nec = nec_num,
          this_hyp[, ncol(this_hyp)]
    )
  }, this_hyp = hypothesis$hyp_mat, nec_num = hypothesis$n_ec, SIMPLIFY = FALSE)

  hyp <- reverse_rename_function(hypothesis$original_hypothesis)

  if(comparison == "unconstrained"){
    hypotheses <- c(hypotheses,
                    list(ormle(est = x,
                               covmtrx = Sigma,
                               const = matrix(c(rep(0, length(x))), nrow = 1),
                               nec = 0,
                               rhs = 0)
                    ))
    hyp <- c(hyp, "Hu")
  }
  res <- compare_hypotheses(hypotheses)
  fit <- res$comparisons
  #browser()
  if(comparison == "complement"){
    complement <- do.call(comp, c(hypotheses[[1]], wt_bar = res[[1]][[2]]))
    fit <- rbind(fit, complement)
    hyp <- c(hyp, "Hc")
  }

  fit$gorica_weights <- compute_weights(fit$gorica)

  Goricares <- list(
    fit = fit,
    #BFmatrix = BFmatrix,
    #b = b,

    call = cl,
    model = x,
    hypotheses = hyp,
    #independent_restrictions = rank_hyp,
    estimates = x#,
    #n = n
  )
  class(Goricares) <- "gorica"
  Goricares
}



#' @method gorica htest
#' @export
gorica.htest <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {
    stop("To be able to run gorica on the results of an object returned by t.test(), you must first load the 'gorica' package, and then conduct your t.test. The standard t.test does not return group-specific variances and sample sizes, which are required by gorica. When you load the gorica package, the standard t.test is replaced by a version that does return this necessary information.")
  }

#' @method gorica t_test
#' @export
gorica.t_test <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {
    cl <- match.call()
    Args <- as.list(cl[-1])

    Args$x <- x$estimate
    Args$n <- x$n

    if(length(x$estimate) == 1){
      Args$Sigma <- x$v/x$n
      Args$group_parameters <- 0
      Args$joint_parameters <- 1
    } else {
      if (!x$method == " Two Sample t-test") {
        Args$Sigma <- lapply(x$v/x$n, as.matrix)
      } else {
        df <- sum(x$n) - 2
        v <- 0
        if (x$n[1] > 1)
          v <- v + (x$n[1] - 1) * x$v[1]
        if (x$n[2] > 1)
          v <- v + (x$n[2] - 1) * x$v[2]
        v <- v/df
        Args$Sigma <- lapply(v / x$n, as.matrix)
      }
      Args$group_parameters <- 1
      Args$joint_parameters <- 0
    }

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x
    class(Gorica_res) <- c("t_test", class(Gorica_res))
    Gorica_res
  }



#' @method gorica lm
#' @export
gorica.lm <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {

    cl <- match.call()
    Args <- as.list(cl[-1])

    Args$x <- coef(x)
    Args$Sigma <- vcov(x)

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x

    #if(!is.null(Warnings)){
    #  Gorica_res$Warnings <- Warnings
    #}
    class(Gorica_res) <- c("gorica_lm", class(Gorica_res))
    Gorica_res
  }

#' @method gorica mplus.model
#' @keywords internal
gorica.mplus.model <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {

    cl <- match.call()
    Args <- as.list(cl[-1])
    mplus_est <- get_estimates(x)
    Args$x <- mplus_est$estimate
    Args$Sigma <- mplus_est$Sigma

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x

    #if(!is.null(Warnings)){
    #  Gorica_res$Warnings <- Warnings
    #}
    class(Gorica_res) <- c("gorica_mplus", class(Gorica_res))
    Gorica_res
  }

#' @method gorica lavaan
#' @export
gorica.lavaan <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {

    cl <- match.call()
    Args <- as.list(cl[-1])
    mplus_est <- get_estimates(x)
    Args$x <- mplus_est$estimate
    Args$Sigma <- mplus_est$Sigma

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x

    #if(!is.null(Warnings)){
    #  Gorica_res$Warnings <- Warnings
    #}
    class(Gorica_res) <- c("gorica_lavaan", class(Gorica_res))
    Gorica_res
  }

#' @method gorica lmerMod
#' @export
gorica.lmerMod <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {

    cl <- match.call()
    Args <- as.list(cl[-1])

    Args$x <- fixef(x)
    Args$Sigma <- vcov(x)

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x

    #if(!is.null(Warnings)){
    #  Gorica_res$Warnings <- Warnings
    #}
    class(Gorica_res) <- c("gorica_lmerMod", class(Gorica_res))
    Gorica_res
  }

#' @method gorica model_estimates
#' @export
gorica.model_estimates <-
  function(x,
           hypothesis,
           comparison = "unconstrained",
           ...) {

    cl <- match.call()
    Args <- as.list(cl[-1])

    Args$x <- x$estimate
    Args$Sigma <- x$Sigma

    Gorica_res <- do.call(gorica, Args)
    Gorica_res$call <- cl
    Gorica_res$model <- x

    #if(!is.null(Warnings)){
    #  Gorica_res$Warnings <- Warnings
    #}
    class(Gorica_res) <- c("gorica_model_estimates", class(Gorica_res))
    Gorica_res
  }

