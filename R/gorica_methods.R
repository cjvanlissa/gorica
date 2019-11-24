#' Evaluate informative hypotheses using the GORICA
#'
#' GORICA is an acronym for "generalized order-restricted information criterion
#' approximation". It can be utilized to evaluate informative hypotheses, which
#' specify directional relationships between model parameters in terms of
#' (in)equality constraints.
#'
#' @param x An R object containing the outcome of a statistical analysis.
#' Currently, the following objects can be processed:
#' \itemize{
#' \item \code{lm()} objects (anova, ancova, multiple regression).
#' \item \code{t_test()} objects.
#' \item \code{lavaan} objects.
#' \item \code{lmerMod} objects.
#' \item A named vector containing the estimates resulting from a statistical
#' analysis, when the argument \code{Sigma} is also specified.
#' Note that, named means that each estimate has to be labelled such that it can
#' be referred to in \code{hypotheses}.
#' }
#' @param hypothesis A character string containing the informative hypotheses to
#' evaluate (see Details).
#' @param comparison A character string indicating what the \code{hypothesis}
#' should be compared to. Defaults to \code{comparison = "unconstrained"};
#' options include \code{c("unconstrained", "complement", "none")}.
#' @param ... Additional arguments passed to the internal function
#' \code{compare_hypotheses}.
#'
#' @details
#' The GORICA is applicable to not only normal linear models, but also applicable to generalized linear models (GLMs) (McCullagh & Nelder, 1989), generalized linear
#' mixed models (GLMMs) (McCullogh & Searle, 2001), and structural equation
#' models (SEMs) (Bollen, 1989). In addition, the GORICA can be utilized in the context of contingency tables for which (in)equality constrained hypotheses do not necessarily contain linear restrictions on cell probabilities, but instead often contain non-linear restrictions on cell probabilities.
#'
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
#' @return An object of class \code{gorica}, containing the following elements:
#' \itemize{
#' \item \code{fit}  A \code{data.frame} containing the loglikelihood, penalty
#' (for complexity), the GORICA value, and the GORICA weights. The GORICA
#' weights are calculated by taking into account the misfits and complexities of
#' the hypotheses under evaluation. These weights are used to quantify the
#' support in the data for each hypothesis under evaluation. By looking at the
#' pairwise ratios between the GORICA weights, one can determine the relative
#' importance of one hypothesis over another hypothesis.
#' \item \code{call}  The original function call.
#' \item \code{model}  The original model object (\code{x}).
#' \item \code{estimates}  The parameters extracted from the \code{model}.
#' \item \code{Sigma}  The asymptotic covariance matrix of the
#' \code{estimates}.
#' \item \code{comparison}  Which alternative hypothesis was used.
#' \item \code{hypotheses}  The hypotheses evaluated in \code{fit}.
#' }
#' @author Caspar van Lissa, Yasin Altinisik, Rebecca Kuiper
#' @references Altinisik, Y. (2018). Evaluation of Inequality Constrained
#' Hypotheses Using an Akaike-Type Information Criterion (Doctoral dissertation,
#' Utrecht University). ISBN: 978-90-393-6918-0.
#' \url{https://dspace.library.uu.nl/handle/1874/360604}
#'
#' Bollen, K. (1989). Structural equations with latent variables. New York, NY:
#' John Wiley and Sons.
#'
#' Kuiper, R. M., Hoijtink, H., & Silvapulle, M. J. (2011).
#' An Akaike-type information criterion for model selection under inequality
#' constraints. Biometrika, 98, 495-501.
#' \href{https://doi.org/10.1093/biomet/asr002}{doi:10.1093/biomet/asr002}
#'
#' Kuiper, R. M., Hoijtink, H., & Silvapulle, M. J. (2012).
#' Generalization of the order-restricted information criterion for multivariate
#' normal linear models. Journal of statistical planning and inference, 142(8),
#' 2454-2463. \href{https://doi.org/10.1016/j.jspi.2012.03.007}{
#' doi:10.1016/j.jspi.2012.03.007}
#'
#' McCullagh, P. & Nelder, J. (1989). Generalized linear models (2nd ed.). Boca
#' Raton, FL: Chapman & Hall / CRC.
#'
#' McCullogh, C. E., & Searle, S. R. (2001). Generalized linear and mixed
#' models. New York, NY: Wiley.
#' @examples
#' \dontshow{
#' # EXAMPLE 1. One-sample t test
#' ttest1 <- t_test(iris$Sepal.Length,mu=5)
#' gorica(ttest1,"x<5.8", iterations = 5)
#'
#' # EXAMPLE 2. ANOVA
#' aov1 <- aov(yield ~ block-1 + N * P + K, npk)
#' gorica(aov1,hypothesis="block1=block5;
#'    K1<0", iterations = 5)
#'
#' # EXAMPLE 3. gml
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fit <- glm(counts ~ outcome-1 + treatment, family = poisson())
#' gorica(fit, "outcome1 > (outcome2, outcome3)", iterations = 5)
#'
#' # EXAMPLE 4. ANOVA
#' res <- lm(Sepal.Length ~ Species-1, iris)
#' est <- get_estimates(res)
#' est
#' gor <- gorica(res, "Speciessetosa < (Speciesversicolor, Speciesvirginica)",
#' comparison = "complement", iterations = 5)
#' gor
#' }
#' \donttest{
#' # EXAMPLE 1. One-sample t test
#' ttest1 <- t_test(iris$Sepal.Length,mu=5)
#' gorica(ttest1,"x<5.8")
#'
#' # EXAMPLE 2. ANOVA
#' aov1 <- aov(yield ~ block-1 + N * P + K, npk)
#' gorica(aov1,hypothesis="block1=block5;
#'    K1<0")
#'
#' # EXAMPLE 3. gml
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fit <- glm(counts ~ outcome-1 + treatment, family = poisson())
#' gorica(fit, "outcome1 > (outcome2, outcome3)")
#'
#' # EXAMPLE 4. ANOVA
#' res <- lm(Sepal.Length ~ Species-1, iris)
#' est <- get_estimates(res)
#' est
#' gor <- gorica(res, "Speciessetosa < (Speciesversicolor, Speciesvirginica)",
#' comparison = "complement")
#' gor
#' }
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
  Goricares <- list(
    fit = NULL,
    call = cl,
    model = x,
    estimates = x,
    Sigma = Sigma,
    comparison = comparison
  )

  if(is.list(Sigma) & length(Sigma) == 1) Sigma <- Sigma[[1]]
  names(x) <- rename_function(names(x))
  colnames(Sigma) <- rownames(Sigma) <- names(x)
  # Parse hypotheses --------------------------------------------------------
  #ren_estimate <- rename_estimate(x)
  if(!inherits(comparison, "character")|length(comparison) > 1){
    stop("Argument 'comparison' must be an atomic character string.")
  } else {
    comp_arg <- pmatch(comparison, c("unconstrained", "complement", "none"))
    if(is.na(comp_arg)) stop("Argument 'comparison' did not match one of the available options: 'unconstrained', 'complement', or 'none'.")
    comparison <- c("unconstrained", "complement", "none")[pmatch(comparison, c("unconstrained", "complement", "none"))]
  }
  if(inherits(hypothesis, "character")){
    hypothesis <- rename_function(hypothesis)
    hyp_params <- params_in_hyp(hypothesis)
    coef_in_hyp <- sort(unique(charmatch(rename_function(hyp_params),
                             rename_function(names(x)))))
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
    # Drop parameters not in hypothesis
    x <- x[coef_in_hyp]
    Sigma <- Sigma[coef_in_hyp, coef_in_hyp]

    hypothesis <- parse_hypothesis(names(x), hypothesis)
  } else {
    if(inherits(hypothesis, "list") & !is.null(hypothesis[["hyp_mat"]]) & !is.null(hypothesis[["n_ec"]])){
      hypothesis$original_hypothesis <- matrix_to_hyp(hypothesis, names(x))
    } else {
      stop("Argument 'hypothesis' must either be a character string with inequality constraints, or a list with an element 'hyp_mat', consisting of a list of hypothesis matrices, and and element 'n_ec', consisting of an integer vector with the number of equality constraints for each hypothesis matrix in 'hyp_mat'.")
    }
  }
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
                               constr = matrix(c(rep(0, length(x))), nrow = 1),
                               nec = 0,
                               rhs = 0)
                    ))
    hyp <- c(hyp, "Hu")
  }
  res <- compare_hypotheses(hypotheses, ...)
  fit <- res$comparisons

  if(comparison == "complement"){
    complement <- do.call(comp, c(hypotheses[[1]], wt_bar = res[[1]][[2]]))
    fit <- rbind(fit, complement)
    hyp <- c(hyp, "Hc")
  }

  fit$gorica_weights <- compute_weights(fit$gorica)

  Goricares[c("fit", "hypotheses")] <- list(fit, hyp)
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
    #Args$n <- x$n

    if(length(x$estimate) == 1){
      Args$Sigma <- list(matrix(x$v/x$n))
    } else {
      if (!x$method == " Two Sample t-test") {
        Args$Sigma <- list(diag(x$v/x$n)) #lapply(x$v/x$n, as.matrix)
      } else {
        df <- sum(x$n) - 2
        v <- 0
        if (x$n[1] > 1)
          v <- v + (x$n[1] - 1) * x$v[1]
        if (x$n[2] > 1)
          v <- v + (x$n[2] - 1) * x$v[2]
        v <- v/df
        Args$Sigma <- list(diag(v / x$n)) #lapply(v / x$n, as.matrix)
      }
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
    if(!is.null(Args[["standardize"]])){
      if(Args[["standardize"]]){
        warning("Cannot standardize an object of class 'lm'. Using unstandardized coefficients.")
      }
    }
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
           standardize = FALSE,
           ...) {
    cl <- match.call()
    Args <- as.list(cl[-1])
    mplus_est <- get_estimates(x, standardize)
    Args$x <- mplus_est$estimate
    Args$Sigma <- mplus_est$Sigma
    Args$hypothesis <- force(hypothesis)
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

