## ----------------
## CLASS definition
## ----------------

## gpcm Class

#' Gaussian Process Classification (GPC) models class
#'
#' S4 class for GPC models.
#'
#' @section Objects from the Class: To create a \code{gpcm} object, use \code{\link{gpcm}}. See also this function for more details.
#'
#'
#' @slot d Object of class \code{"integer"}. The spatial dimension.
#' @slot n Object of class \code{"integer"}. The number of observations.
#' @slot X Object of class \code{"matrix"}. The design of experiments.
#' @slot y Object of class \code{"matrix"}. The vector of binary observations at design points (+/-1) corresponding to the class labels.
#' @slot X.std Object of class \code{"numeric"}. The vector of standard deviation values of design points.
#' @slot X.mean Object of class \code{"numeric"}. The vector of mean values of design points.
#' @slot call Object of class \code{"language"}. User call reminder.
#' @slot coef.m Object of class \code{"numeric"}. Mean coefficient of latent GP.
#' @slot coef.cov Object of class \code{"numeric"}. Covariance coefficients of latent GP.
#' @slot covariance Object of class \code{"covKernel"}. A DiceKriging object specifying the covariance structure.
#' @slot noise.flag Object of class \code{"logical"}. Are the observations noisy?
#' @slot noise.var Object of class \code{"numeric"}. Nugget effect.
#' @slot param.estim Object of class \code{"logical"}. \code{TRUE} if at least one parameter is estimated, \code{FALSE} otherwise.
#' @slot lower Object of class \code{"numeric"}. Lower bounds for covariance parameters estimation.
#' @slot upper Object of class \code{"numeric"}. Upper bounds for covariance parameters estimation.
#' @slot logLik Object of class \code{"numeric"}. Value of the log-Likelihood at its optimum.
#' @slot Z_obs Object of class \code{"matrix"}. A \code{nobs} * \code{nsimu} matrix of samples of the latent process at design points.
#' @slot l Object of class \code{"numeric"}.  Lower truncation points. Parameter to generate new \code{Z_obs}.
#' @slot u Object of class \code{"numeric"}. Upper truncation points. Parameter to generate new \code{Z_obs}.
#' @slot K Object of class \code{"matrix"}. Covariance matrix of design points. Parameter to generate new \code{Z_obs}
#' @slot invK Object of class \code{"matrix"}. The inverse of the matrix \code{K} whose Cholesky decomposition was given.
#'
#' @seealso { \code{\link{gpcm}} for more details about slots and to create a \code{gpcm} object. \code{{covStruct.create}} in \code{DiceKriging} to construct a covariance structure.}
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#' @exportClass gpcm
#'
#' @keywords classes
#'
setClass("gpcm",slots =
           list(
             d = "integer",
             n = "integer",
             ## data
             X = "matrix",
             y = "matrix",
             X.std ="numeric",
             X.mean = "numeric",
             call = "language",
             ## mean
             coef.m = "numeric",
             ## cov
             coef.cov = "numeric",
             ## covariance
             covariance = "covKernel",
             ## noisy observations
             noise.flag = "logical",
             noise.var = "numeric",
             ## optimization
             param.estim = "logical",
             lower = "numeric",
             upper = "numeric",
             logLik = "numeric",
             ## auxiliary variables
             Z_obs = "matrix",
             l = "numeric",
             u ="numeric",
             K = "matrix",
             invK = "matrix"
           )
)


validTrackObject <- function(object) {
  if (object@n <= object@d) {
    return("the number of experiments must be larger than the spatial dimension")
  }
  
  if (ncol(object@y) != 1) {
    return("the response must be 1-dimensional")
  }
  
  if (!identical(nrow(object@X), nrow(object@y))) {
    return("the number of observations is not equal to the number of experiments")
  }
  
  if (any(object@y!=1 & object@y!=-1)){
    return("at least one observation is not equal 1 or -1")
  }
  
  TRUE
}

setValidity("gpcm", validTrackObject)

##*****************************************************************************
##*****************************************************************************
##                        P R E D I C T  METHOD
##*****************************************************************************
##*****************************************************************************

#' Predict class probability at newdata for a Gaussian Process Classification (GPC) model
#'
#' Predicted probability of class 1. Optionally, conditional covariance based on a \code{gpcm} model and 95% quantiles of the probability of class 1 are returned.
#'
#' @param object an object of class \code{gpcm}.
#' @param newdata a vector, matrix of points to be predicted.
#' @param nsimu an optional integer indicating whether to resample latent GP at observation points and how many samples are required. If \code{NULL}, current samples are used. Default is \code{NULL}.
#' @param light.return an optional boolean. If \code{TRUE}, only \code{prob} is returned. Default is \code{FALSE}.
#' @param checkNames an optional boolean. If \code{TRUE}, a consistency test is performed between the names of newdata and the names of the experimental design (contained in \code{object@Xf}). Default is \code{FALSE}.
#' @param seed to fix the seed (used if \code{nsimu} is not \code{NULL}). Default is \code{NULL}.
#' @param ... no other argument for this method
#'
#' @return
#'   \item{prob}{the (averaged) probability of class 1 at \code{newdata}.}
#'   \item{lower95, upper95}{ 95% confidence bounds for the probability at \code{newdata}.}
#'   \item{probs}{a matrix of sample predicted probabilities.}
#'   \item{Zsimu_var, Zsimu_mean}{conditional variance vector and mean matrix of the latent GP at \code{newdata}.}
#'   \item{cov}{conditional covariance matrix at \code{newdata}.}
#'   \item{c}{an auxiliary matrix, containing all the covariances between \code{newdata} and design points \code{Xf}.}
#'   \item{lambda}{an auxiliary vector, product of the inverse covariance matrix \code{invK} returned by \code{object} and the unconditional covariance matrix \code{c} between \code{newdata} and design points \code{Xf}.}
#'   \item{kz}{an auxiliary matrix, corresponding to the unconditional covariance matrix at \code{newdata}.}
#'
#'
#' @references Bachoc, F., Helbert, C. & Picheny, V. Gaussian process optimization with failures: classification and convergence proof. \emph{J Glob Optim} \bold{78}, 483–506 (2020). \doi{10.1007/s10898-020-00920-0}.
#' @references Roustant, O., Ginsbourger, D. & Deville, Y. Contributors: Chevalier, C. , Richet, Y. DiceKriging: Kriging Methods for Computer Experiments. R package version 1.6.0. \url{https://CRAN.R-project.org/package=DiceKriging}.
#'
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#' @name predict
#' @aliases predict,gpcm-method predict.gpcm
#' @rdname predict.gpcm
#' @exportMethod predict
#' @exportS3Method predict gpcm
#' 
#' @usage \method{predict}{gpcm}(object, newdata, nsimu = NULL, 
#' light.return = FALSE, checkNames=FALSE, seed = NULL, ...)
#' 
#' @importFrom tmvtnorm rtmvnorm
#' @importMethodsFrom DiceKriging covMat1Mat2 covMatrix
#' 
#' @examples
#' # ----------------------------------
#' # A 2D example - Branin-Hoo function
#' # ----------------------------------
#'
#' # 30-points DoE, and the corresponding response
#' d <- 2
#' nb_PX <- 30
#' require(DiceDesign)
#' X <- lhsDesign(nb_PX, d, seed = 123)$design
#' Xopt <- maximinSA_LHS(X, T0 = 10, c = 0.99, it = 1000)
#' x <- Xopt$design
#' require(DiceKriging)
#' fx <- apply(x, 1, branin)
#' s <- ifelse(fx < 14, -1, 1)
#' f <- s
#' Xf <- as.matrix(x)
#'
#' # Bulding GPC model
#' GPCmodel <- gpcm(f = f, Xf = Xf, coef.m = -0.1, coef.cov=c(0.8,0.5))
#'
#' # Graphics - Predict probabilities
#' ngrid <- 50
#' x.grid <- seq(0, 1., length.out = ngrid)
#' grid <- as.matrix(expand.grid(x.grid, x.grid))
#' probabilities <- predict(object = GPCmodel, newdata = grid)$prob
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
predict.gpcm <- function(object, newdata, nsimu = NULL, light.return = FALSE, checkNames=FALSE, seed = NULL, ...) {
  
  Xf <- object@X
  
  if (checkNames) {
    newdata <- checkNames(X1 = Xf, X2 = newdata, X1.name = "the design", X2.name = "newdata")
  } else {
    newdata <- as.matrix(newdata)
    d.newdata <- ncol(newdata)
    if (!identical(d.newdata, object@d)) {
      stop("newdata must have the same numbers of columns than the experimental design")
    }
    if (!identical(colnames(newdata), colnames(Xf))) {
      ##  warning("column names mismatch between 'newdata' and the experimental design -
      ## the columns of 'newdata' are interpreted in the same order as the experimental design names")
      colnames(newdata) <- colnames(Xf)
    }
  }
  
  if (is.null(nsimu)) {
    Z_obs <- object@Z_obs
    nsimu <- ncol(Z_obs)
  } else {
    set.seed(seed)
    Z_obs <- t(tmvtnorm::rtmvnorm(lower = object@l, upper = object@u, sigma = object@K, n = nsimu, algorithm = "gibbs")) + object@coef.m
  }
  
  if (!any(is.nan(object@X.mean))) {
    newdata <- scale(newdata, center = object@X.mean, scale = object@X.std)
    Xf <- scale(object@X, center = object@X.mean, scale = object@X.std)
  }
  
  
  
  kstar <- DiceKriging::covMat1Mat2(object = object@covariance, X1 = Xf, X2 = newdata)
  kz <- DiceKriging::covMatrix(object = object@covariance, X = newdata)$C
  probs <- matrix(NA, nrow = nrow(newdata), ncol = nsimu)
  
  lambda <- object@invK %*% kstar
  Zsimu_c <- kz - t(kstar) %*% lambda
  Zsimu_var <- object@covariance@sd2 - colSums(kstar * lambda)
  Zsimu_var[Zsimu_var < 0] <- 0
  Zsimu_sd <- sqrt(Zsimu_var)
  
  Zsimu_mean <- object@coef.m + t(lambda) %*% (Z_obs - object@coef.m)
  Zbar <- (apply(Zsimu_mean, 2, function(x) x / Zsimu_sd))
  probs <- 1 - pnorm(-Zbar)
  if (nrow(newdata) == 1) probs <- matrix(probs, ncol = nsimu)
  
  
  q <- apply(probs, 1, quantile, probs = c(.05, .95))
  
  if (!light.return) {
    output <- list(prob = rowMeans(probs), lower95 = q[1, ], upper95 = q[2, ], probs = probs, Zsimu_var = Zsimu_var, Zsimu_mean = Zsimu_mean, cov = Zsimu_c, c = kstar, lambda = lambda, kz = kz)
  } else {
    output <- rowMeans(probs)
  }
  return(output)
}




if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object, ...) standardGeneric("predict")
  )
}

setMethod(f="predict", signature="gpcm",
          function(object, newdata, nsimu = NULL, light.return = FALSE, checkNames = FALSE, seed = NULL, ...) {
            predict.gpcm(object = object, newdata = newdata, nsimu = nsimu, light.return=light.return, checkNames=checkNames, seed = seed, ...)
          }
)



##*****************************************************************************
##****************************************************************************
##			                    U P D A T E  METHOD
##****************************************************************************
##*****************************************************************************

#' Update of a Gaussian Process Classification (GPC) model
#'
#' Update a \code{\link{gpcm}} object when one or many new observations are added.
#'
#' @param object an object of \code{\link{gpcm}} class.
#' @param newf a vector corresponding to the new binary observations (+/-1) at \code{newXf} locations. These locations can be new locations or existing ones.
#' @param newXf a matrix with \code{object@d} columns representing the locations to be updated. These locations can be new locations or existing ones.
#' @param newnoise.var optional scalar, nugget effect at new observations.
#' @param covandmean.reestim should the mean and covariance parameters be re-estimated? Default is \code{TRUE}.
#' @param multistart an optional integer indicating the number of initial points from which running the \code{BFGS} for covariance parameter optimization. Default is 1.
#' @param seed to fix the seed, default is \code{NULL}.
#' @param lower (see below).
#' @param upper \code{lower}, \code{upper}: bounds for the covariance parameters (scalars or vectors), if \code{NULL} they are set to 0.2 and 3, respectively.
#' @param nsimu an integer indicating the number of samples of the latent GP at observation points \code{Xf} to generate. Must be a non-null integer. Default is 100.
#' @param normalize a logical parameter indicating whether to normalize the input matrix \code{Xf}. If \code{TRUE}, the matrix will be normalized using \code{X.mean} and \code{X.std} values if given; otherwise, the mean and standard deviation are calculated and used for normalization.
#' @param newX.alreadyExist Boolean: indicate whether the locations \code{newXf} are all news or not. Default: \code{TRUE}, corresponding to existing locations in newX.
#' @param ... no other argument for this method
#' 
#' @return Updated \code{\link{gpcm}} object.
#'
#' @references Bachoc, F., Helbert, C. & Picheny, V. Gaussian process optimization with failures: classification and convergence proof. \emph{J Glob Optim} \bold{78}, 483–506 (2020). \doi{10.1007/s10898-020-00920-0}.
#' @references Kotecha, J. H., Djuric, P. M. (1999). Gibbs Sampling Approach For Generation of Truncated Multivariate Gaussian Random Variables. \emph{IEEE Computer Society}, 1757–1760.
#' @references Wilhelm, S. tmvtnorm: Truncated Multivariate Normal and Student t Distribution. R package version 	1.6. \url{https://CRAN.R-project.org/package=tmvtnorm}.
#' @references Roustant, O., Ginsbourger, D. & Deville, Y. Contributors: Chevalier, C. , Richet, Y. DiceKriging: Kriging Methods for Computer Experiments. R package version 1.6.0. \url{https://CRAN.R-project.org/package=DiceKriging}.
#' @references Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. \emph{SIAM Journal on Scientific Computing}, \bold{16}, 1190–1208. \doi{10.1137/0916069}.
#'
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#'
#' @seealso \code{\link[GPCsign]{gpcm}}
#' @name update
#' @aliases update,gpcm-method update.gpcm
#' @rdname update.gpcm
#' @exportMethod update
#' @exportS3Method update gpcm
#' 
#' @usage \method{update}{gpcm}(object, newf, newXf, newX.alreadyExist,
#'  newnoise.var, covandmean.reestim=TRUE, multistart = 1, seed = NULL,
#'   lower = NULL, upper = NULL, nsimu = 100, normalize = TRUE, ...)
#' @examples
#' # ----------------------------------
#' # A 1D example - sinusoidal function
#' # ----------------------------------
#'
#' # Test function
#' sinusoidal_function <- function(x) {
#'   sin(4 * pi * x)}
#'
#' # Desing of Experiment Xf and the corresponding sign f
#' Xf <- as.matrix(c(0.07, 0.19, 0.42, 0.56, 0.81, 0.90))
#' f <- rep(1,length(Xf)); f[(sinusoidal_function(Xf)<0)]<- -1
#'
#' # Builidng a GPC model
#' GPCmodel1 <- gpcm(f = f, Xf = Xf, coef.m=0, coef.cov=0.26)
#' print(GPCmodel1)
#'
#' # New points added to the gpcm object.
#' newXf <- as.matrix(c(0.1,0.5,0.7, 0.95))
#' newf <- rep(1,length(newXf)); newf[(sinusoidal_function(newXf)<0)]<- -1
#'
#' # Updating GPC model
#' NewGPCmodel <- update(object = GPCmodel1, newf = newf, newXf = newXf)
#' print(NewGPCmodel)

update.gpcm <- function(object, newf, newXf, newX.alreadyExist=TRUE, newnoise.var=NULL, covandmean.reestim=TRUE, multistart = 1, seed = NULL, lower = NULL, upper = NULL, nsimu = 100, normalize = TRUE, ...) {
  set.seed(seed)
  
  
  Xf <- rbind(object@X, as.matrix(newXf))
  f <- rbind(object@y, as.matrix(newf))
  if(newX.alreadyExist){
    duplicates_id <- duplicated(Xf)
    f <- f[!duplicates_id]
    Xf <- Xf[!duplicates_id,]
  }
  
  
  
  coef.cov <- object@coef.cov
  coef.m <- object@coef.m
  covtype <- object@covariance@name
  
  
  if (length(lower)!=object@d) lower <- rep(0.2, object@d)
  if (length(upper)!=object@d) upper <- rep(3., object@d)
  
  
  if (normalize == TRUE) {
    if (!any(is.nan(object@X.mean))) {
      X.mean <- object@X.mean
      X.std <- object@X.std
    }else {
      X.mean <- colMeans(Xf)
      X.std <- sqrt(apply(Xf, FUN = var, MARGIN = 2))
    }}else{
      X.mean <- NULL
      X.std <- NULL
    }
  
  if(!normalize && !any(is.nan(object@X.mean)))
  {
    coef.cov <- NULL
    coef.m <- NULL
  }
  
  if(covandmean.reestim && multistart==1)
  {
    coef.cov <- NULL
    coef.m <- NULL
  }
  
  if(is.null(newnoise.var)){
    noise.var <- object@noise.var
  }else{
    noise.var <- newnoise.var
  }
  GPCmod <- gpcm(f = f, Xf = Xf, covtype=covtype, coef.cov = coef.cov, coef.m = coef.m, noise.var = noise.var, seed=seed, nsimu = nsimu,
                 lower=lower, upper=upper, X.mean = X.mean, X.std = X.std, normalize = normalize,  multistart=multistart)
  return(GPCmod)
}



if(!isGeneric("update")) {
  setGeneric(name = "update",
             def = function(object, ...) standardGeneric("update")
  )
}

setMethod(f="update",signature="gpcm",
          function(object, newf,  newXf, newnoise.var=NULL, covandmean.reestim=TRUE, multistart = 1, seed = NULL,
                   lower = NULL, upper = NULL, nsimu = 100, normalize = TRUE, ...) {
            update.gpcm(object = object, newf = newf, newXf = newXf, newnoise.var=newnoise.var, covandmean.reestim=covandmean.reestim,
                        multistart = multistart,seed = seed, lower = lower, upper = upper, nsimu = nsimu,
                        normalize = normalize, ...)
          }
)

##*****************************************************************************
##*****************************************************************************
##                        S H O W  METHOD
##*****************************************************************************
##*****************************************************************************

#' Print values of a  Gaussian Process Classification (GPC) model
#'
#' Show method for \code{gpcm} object. Printing the main features of a GPC model.
#'
#' @param object an object of class \code{gpcm}. See \code{\link{gpcm}}.
#' @return
#' returns an invisible 'NULL'
#' @name show
#' @rdname show.gpcm
#' @seealso [gpcm()]
#' @author Morgane MENZ, Céline HELBERT, Victor PICHENY, François BACHOC. Contributors: Naoual SERRAJI.
#' @export
#' 
#' @examples 
#' ## 20-points DoE, and the corresponding response
#' d <- 2
#' nb_PX <- 20
#' require(DiceDesign)
#' x <- lhsDesign(nb_PX, d, seed = 123)$design
#' require(DiceKriging)
#' fx <- apply(x, 1, branin)
#' f <- ifelse(fx < 14, -1, 1)
#' Xf <- as.matrix(x)
#' 
#' ## GPC model 
#' model <- gpcm(f, Xf, coef.m=0, coef.cov=c(0.5,0.5))
#' 
#' ## print the result 
#' show(model)
`show.gpcm` <- function(object) {
  cat("\n")
  cat("Call:\n")
  print(object@call)
  
  cat("Mean  coeff.:\n")
  print(object@coef.m, quote=FALSE)
  
  show(object@covariance)
  
}



if(!isGeneric("show")) {
  setGeneric(name = "show",
             def = function(object) standardGeneric("show")
  )
}

setMethod(f="show", signature="gpcm",
          function(object){
            show.gpcm(object)
          }
)

