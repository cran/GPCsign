#' Compute logLikelihood
#'
#' Computes and returns the  log-likelihood value, the covariance matrix of latent process and covariance structure of a Gaussian Process Classification (GPC) model.
#'
#' @param par vector contains the \code{coef.m} and the log of \code{coef.cov}.
#' @param f vector of binary observations (+/-1) corresponding to the class labels.
#' @param Xf a matrix representing the design of experiments.
#' @param covtype a character string specifying the covariance structure for the latent GP. Default is \code{matern_5_2}.
#' @param noise.var nugget effect. Default is 1e-6.
#' @param seed to fix the seed, default is \code{NULL}.
#' @param return.all an optional boolean. If \code{FALSE}, only the log-likelihood is returned; if \code{TRUE}, \code{K} and \code{cov.fun} are also returned. Default is \code{FALSE}.
#' @usage logLikFunc(par, f, Xf, covtype = "matern5_2", noise.var = 1e-6,
#'      seed = NULL, return.all = FALSE)
#'   
#' @return
#'    \item{logLik}{the log-likelihood.}
#'    \item{K}{the covariance matrix of latent process.}
#'    \item{cov.fun}{a DiceKriging object specifying the covariance structure.}
#' @name logLikFunc
#' @references Bachoc, F., Helbert, C. & Picheny, V. Gaussian process optimization with failures: classification and convergence proof. \emph{J Glob Optim} \bold{78}, 483–506 (2020). \doi{10.1007/s10898-020-00920-0}
#' @references Botev, Z.,  Belzile, L. TruncatedNormal: Truncated Multivariate Normal and Student Distributions. R package version 2.2.2 https://cran.r-project.org/package=TruncatedNormal
#' @references Botev, Z. I.  (2017), \emph{The normal law under linear restrictions:simulation and estimation via minimax tilting}, Journal of the Royal Statistical Society, Series B, \bold{79} (1), pp. 1-24.
#' @references Roustant, O., Ginsbourger, D. & Deville, Y. Contributors: Chevalier, C. , Richet, Y. DiceKriging: Kriging Methods for Computer Experiments. R package version 1.6.0. \url{https://CRAN.R-project.org/package=DiceKriging}.
#' 
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#' @export
#' 
#' @importFrom TruncatedNormal pmvnorm
#' @importFrom DiceKriging covStruct.create
#'
#' @examples
#' # ------------
#' # A 1D example
#' # ------------
#'
#' # Design of Experiments Xf and the corresponding signs f
#' Xf <- as.matrix(c(0.08, 0.27, 0.42, 0.65, 0.78, 0.84))
#' f <- c(1, -1, -1, 1, -1, -1)
#'
#' # loglikelihood and covariance matrix at Xf 
#' par <- c(coef.cov = 0.1, coef.m = 0)
#' result <- logLikFunc(par = par, f = f, Xf = Xf, return.all = TRUE)
#' K <- result$K
#' logLik <- result$logLik
#' print(logLik)
#'
logLikFunc <- function(par, f, Xf, covtype = "matern5_2", noise.var = 1e-6, seed = NULL, return.all = FALSE) {
  
  coef.m <- par[1]
  coef.cov <- exp(par[-1])
  
  d <- ncol(Xf)
  nobs <- nrow(Xf)
  
  cov.fun <- DiceKriging::covStruct.create(covtype = covtype, d = d, known.covparam = "All", var.names = colnames(Xf), coef.cov = coef.cov, coef.var = 1)
  K <- DiceKriging::covMatrix(object = cov.fun, X = Xf, noise.var = rep(noise.var, nobs))$C
  
  nobs <- length(f)
  lower <- rep(0, nobs)
  lower[f == -1] <- -Inf
  upper <- rep(Inf, nobs)
  upper[f == -1] <- 0
  
  set.seed(seed)
  logLikli <- TruncatedNormal::pmvnorm(lb = lower, ub = upper, mu= rep(coef.m,length(upper)), sigma = K, B=2000, type="qmc")
  logLik = log(logLikli[1])
  res <- list(logLik = logLik, K = K, cov.fun = cov.fun)
  
  if (!return.all) {
    output <- logLik
  } else {
    output <- res
  }
  
  return(output)
}