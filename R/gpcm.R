#' Fit and/or create a Gaussian Process Classification (GPC) model
#'
#' \code{gpcm} is used to fit GPC models. When parameters are known, the function creates a model using given parameters. Otherwise, they are estimated by Maximum Likelihood. In both cases, the result is a \code{gpcm} object.
#'
#' @param f a vector containing the binary observations (+/-1) corresponding to the class labels.
#' @param Xf a matrix representing the design of experiments.
#' @param covtype a character string specifying the covariance structure for the latent GP. Default is \code{matern_5_2}.
#' @param noise.var variance value standing for the homogeneous nugget effect. Default is 1e-6.
#' @param coef.cov an optional vector containing the values for covariance parameters for the latent GP. (See below).
#' @param coef.m an optional scalar corresponding to the mean value of the latent GP. If both \code{coef.cov} and \code{coef.m} are provided, no covariance parameter estimation is performed. If at least one of them is missing, both are estimated.
#' @param multistart an optional integer indicating the number of initial points from which running the \code{BFGS} for covariance parameter optimization. The multiple optimizations
#' will be performed in parallel provided that a parallel backend is registered (see package future)
#' @param seed to fix the seed, default is \code{NULL}.
#' @param lower (see below).
#' @param upper \code{lower}, \code{upper}: bounds for the covariance parameters (scalars or vectors), if \code{NULL} they are set to 0.2 and 3, respectively.
#' @param nsimu the number of samples of the latent process at observation points \code{Xf} to generate. Must be a non-null integer.
#' @param normalize a logical parameter indicating whether to normalize the input matrix \code{Xf}. If \code{TRUE}, the matrix will be normalized using \code{X.mean} and \code{X.std} values if given; otherwise, the mean and standard deviation are computed and used for normalization.
#' @param X.mean (see below).
#' @param X.std optional vectors containing mean and  standard deviation values for each column of the input matrix. If they are not provided, they are computed from the input matrix \code{Xf}.
#'
#' @usage gpcm(f, Xf, covtype = "matern5_2", noise.var = 1e-6,
#'      coef.cov = NULL, coef.m = NULL, multistart = 1,
#'      seed = NULL, lower = NULL, upper = NULL, nsimu = 100,
#'      normalize = TRUE, X.mean = NULL, X.std = NULL)
#'      
#' @return An object of class \code{gpcm}. See \code{\link{gpcm-class}}.
#' @references Bachoc, F., Helbert, C. & Picheny, V. Gaussian process optimization with failures: classification and convergence proof. \emph{J Glob Optim} \bold{78}, 483–506 (2020). \doi{10.1007/s10898-020-00920-0}.
#' @references Kotecha, J. H., Djuric, P. M. (1999). Gibbs Sampling Approach For Generation of Truncated Multivariate Gaussian Random Variables. \emph{IEEE Computer Society}, 1757–1760.
#' @references Wilhelm, S. tmvtnorm: Truncated Multivariate Normal and Student t Distribution. R package version 	1.6. \url{https://CRAN.R-project.org/package=tmvtnorm}.
#' @references Roustant, O., Ginsbourger, D. & Deville, Y. Contributors: Chevalier, C. , Richet, Y. DiceKriging: Kriging Methods for Computer Experiments. R package version 1.6.0. \url{https://CRAN.R-project.org/package=DiceKriging}.
#' @references Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. \emph{SIAM Journal on Scientific Computing}, \bold{16}, 1190–1208. \doi{10.1137/0916069}.
#'
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#'
#' @details
#' The generation of the matrix of samples of the latent process \code{Z_obs} is done using Gibbs sampling. See \code{rtmvnorm} function in \code{tmvtnorm} package.
#'
#' @import stats
#' @importFrom DiceKriging covStruct.create
#' @importMethodsFrom DiceKriging covMatrix
#' @import future.apply
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom future plan
#' @import methods
#' @export
#'
#' @examples
#' # ----------------------------------
#' # A 1D example - sinusoidal function
#' # ----------------------------------
#' sinusoidal_function <- function(x) {
#'   sin(4 * pi * x)
#' }
#'
#' # Design of Experiments Xf and the corresponding signs f
#' Xf <- as.matrix(c(0.07, 0.19, 0.42, 0.56, 0.81, 0.90))
#' f <- rep(1, length(Xf))
#' f[(sinusoidal_function(Xf) < 0)] <- -1
#'
#' # GPC model
#' GPCmodel <- gpcm(f, Xf, multistart = 3)
#'
#' # Graphics of predictions
#' x <- as.matrix(seq(0, 1, length.out = 101))
#' result <- predict(object = GPCmodel, newdata = x)
#' probabilities <- result$prob
#' index <- match(Xf, x)
#' plot(x, probabilities, pch = "-")
#' points(Xf[f == 1], probabilities[index[f == 1]], pch = 20, col = "blue")
#' points(Xf[f == -1], probabilities[index[f == -1]], pch = 20, col = "red")
#' abline(h = 0.5, lty = 2)
#' legend("topright",title = "DoE Xf",title.cex = 0.7, legend = c("+", "-"), 
#'      col = c("blue", "red"), pch = 20)
#'
#' \donttest{
#' # ----------------------------------
#' # A 2D example - Branin-Hoo function
#' # ----------------------------------
#'
#' # 30-points DoE, and the corresponding response
#' d <- 2
#' nb_PX <- 30
#' require(DiceDesign)
#' X <- lhsDesign(nb_PX, d, seed = 123)$design
#' Xopt <- maximinSA_LHS(X, T0 = 10, c = 0.99, it = 10000)
#' x <- Xopt$design
#' require(DiceKriging)
#' fx <- apply(x, 1, branin)
#' f <- ifelse(fx < 14, -1, 1)
#' Xf <- as.matrix(x)
#'
#' # Fit and create a GPC model without parallelisation
#' t0 <- proc.time()
#' GPCmodel <- gpcm(f, Xf, multistart = 3, seed = 123)
#' t1 = proc.time() - t0
#' cat(" time elapsed : ",t1[3])
#' print(GPCmodel)
#'
#' # Graphics - Predict probabilities
#' ngrid <- 50
#' x.grid <- seq(0, 1., length.out = ngrid)
#' grid <- as.matrix(expand.grid(x.grid, x.grid))
#' probabilities <- predict(GPCmodel, newdata = grid, light.return = TRUE)
#' filled.contour(x.grid, x.grid, matrix(probabilities, ngrid, ngrid),
#'                color.palette = function(n) hcl.colors(n, "RdYlBu", rev = FALSE),
#'                main = "probabilities map",
#'                plot.axes = {
#'                  axis(1)
#'                  axis(2)
#'                  points(Xf[f == 1, 1], Xf[f == 1, 2], col = "blue", pch = 21, bg = "blue")
#'                  points(Xf[f == -1, 1], Xf[f == -1, 2], col = "red", pch = 21, bg = "red")
#'                }
#' )
#'
#' # Fit and create a GPC model with parallelisation
#' ## Use multisession futures
#' require(future)
#' plan(multisession)
#' t0 = proc.time()
#' GPCmodel2 <-  gpcm(f,Xf, multistart = 3, seed = 123 )
#' t1 = proc.time() - t0
#' cat(" time elapsed : ",t1[3])
#' print(GPCmodel2)
#' ## Explicitly close multisession workers by switching plan
#' plan(sequential)
#' }
`gpcm` <- function(f, Xf, covtype = "matern5_2", noise.var = 1e-6, coef.cov = NULL, coef.m = NULL, multistart = 1, seed = NULL, lower = NULL, upper = NULL, nsimu = 100, normalize = TRUE, X.mean = NULL, X.std = NULL) {
  
  model <- new("gpcm")
  Xf <- as.matrix(data.frame(Xf))
  model@X <- Xf
  model@y <- as.matrix(f)
  model@d <- ncol(Xf)
  model@n <- nrow(Xf)
  model@call <- match.call()
  model@noise.flag <- (length(noise.var) != 0)
  model@noise.var <- as.numeric(noise.var)
  
  
  if (length(lower)!=model@d) lower <- rep(0.2, model@d)
  if (length(upper)!=model@d) upper <- rep(3., model@d)
  
  lower <- c(-2, log(lower))
  upper <- c(2, log(upper))
  
  model@lower <- lower
  model@upper <- upper
  
  if (normalize == TRUE) {
    if (!(is.null(X.mean) || is.null(X.std))) {
      Xf <- scale(Xf, center = X.mean, scale = X.std)
    }else {
      X.mean <- colMeans(Xf)
      X.std <- sqrt(apply(Xf, FUN = var, MARGIN = 2))
      Xf <- scale(Xf, center = X.mean, scale = X.std)
    }
  }
  else{
    X.mean <- NaN
    X.std <- NaN
  }
  
  model@X.mean <- X.mean
  model@X.std <- X.std
  
  validObject(model)
  
  isCovAndMean <- !(is.null(coef.cov) || is.null(coef.m))
  isMulti <- multistart > 1
  
  if (!(isCovAndMean) || (isCovAndMean & isMulti)) {
    model@param.estim = T
    control <- list(fnscale=-1)
    
    if (!isMulti) {
      if(isCovAndMean){
        par <- c(coef.m, log(coef.cov))
      }else{
        set.seed(seed)
        par <- lower + runif(model@d + 1) * (upper - lower)
      }
     
      res <- optim(
        fn = `logLikFunc`, control=control, par = par, lower = lower, upper = upper, method = "L-BFGS-B",
        f = f, Xf = Xf, covtype = covtype, noise.var = noise.var, seed = seed
      )
    }else {
      if(isCovAndMean){
        set.seed(seed)
        par <- matrix(rep(lower, each = multistart - 1), nrow = multistart - 1) + runif((model@d + 1) * (multistart - 1)) * matrix(rep(upper - lower, each = multistart - 1), nrow = multistart - 1)
        par <- rbind(par, c(coef.m, log(coef.cov)))
      }else{
        set.seed(seed)
        par <- matrix(rep(lower, each = multistart), nrow = multistart) + runif((model@d + 1) * multistart) * matrix(rep(upper - lower, each = multistart), nrow = multistart)
      }
      par <- lapply(apply(par, 1, list), unlist)
      
      
      lres <- future_lapply(
        X = par, FUN = optim, fn = `logLikFunc`, control=control, lower = lower, upper = upper, method = "L-BFGS-B",
        f = f, Xf = Xf, covtype = covtype, noise.var = noise.var, seed = seed, future.globals = TRUE, future.seed = TRUE, future.packages = c("DiceKriging", "mvtnorm")
      )
      
      
      err <- Reduce(cbind, lapply(lres, function(alist) is.atomic(alist)))
      lres <- lres[!err]
      err_null <- Reduce(cbind, lapply(lres, function(alist) is.null(alist$value)))
      if (any(err_null)) lres <- lres[!err_null]
      
      res.logcov <- Reduce(cbind, lapply(lres, function(alist) alist$par[-1]))
      test.logcov <- res.logcov < lower[-1]
      
      bestlogl <- Reduce(cbind, lapply(lres, function(alist) alist$value))
      if (!is.null(dim(test.logcov))) {
        bool.logcov <- apply(test.logcov, 2, any)
        bestlogl[bool.logcov] <- 1e6
      }
      
      res <- lres[[which.min(bestlogl)]]
    }
    coef.m <- res$par[1] # res$par[1]
    coef.cov <- exp(res$par[-1])
    logLik<- res$value
  }else {
    par <- c(coef.m, log(coef.cov))
    logLik <- logLikFunc(par, f, Xf)
    model@param.estim <- F
  }
  
  cov.fun <- DiceKriging::covStruct.create(covtype = covtype, d = model@d, known.covparam = "All", var.names = colnames(Xf), coef.cov = coef.cov, coef.var = 1)
  model@covariance <- cov.fun
  model@coef.cov <- coef.cov
  model@coef.m <- coef.m
  model@logLik <- logLik
  
  
  K <- DiceKriging::covMatrix(object = cov.fun, X = Xf, noise.var = rep(noise.var, model@n))$C
  chol.K <- NULL
  chol.K <- try(chol(K), TRUE)
  if (!is.numeric(chol.K)) {
    return(list(error = TRUE))
  }
  invK <- chol2inv(chol.K)
  model@invK <- invK
  model@K <- K
  
  
  l <- rep(-coef.m, nrow(Xf))
  u <- rep(Inf, nrow(Xf))
  l[f == -1] <- -Inf
  u[f == -1] <- -coef.m
  set.seed(seed)
  Z_obs <- t(tmvtnorm::rtmvnorm(lower = l, upper = u, sigma = K, n = nsimu, algorithm = "gibbs")) + coef.m
  model@Z_obs <- Z_obs
  model@l <- l
  model@u <- u
  validObject(model)
  return(model)
}
