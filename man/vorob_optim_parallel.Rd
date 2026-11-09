\name{vorob_optim_parallel}
\Rdversion{1.4.1}
\alias{vorob_optim_parallel}

\title{ Parallel Vorob'ev criterion }
\description{
Evaluation of the parallel Vorob'ev criterion for some candidate points. To be used in optimization routines, like in \code{\link{max_vorob_parallel}}.
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise.
The criterion is the integral of the posterior Vorob'ev uncertainty.
}
\usage{
vorob_optim_parallel(x, integration.points, integration.weights = NULL,
intpoints.oldmean, intpoints.oldsd,
precalc.data, model, T,
new.noise.var = NULL, batchsize, alpha, current.vorob,
penalisation=NULL,typeEx=">")
}

\arguments{
  	\item{x}{Input vector of size batchsize*d at which one wants to evaluate the criterion. This argument is NOT a matrix.}
  	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points.}
  	\item{intpoints.oldmean}{Vector of size p corresponding to the kriging mean at the integration points before adding the batchsize points \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p corresponding to the kriging standard deviation at the integration points before adding  the batchsize points \code{x} to the design of experiments.}
  	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Target value (scalar). The criterion CANNOT be used with multiple thresholds.}
  	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
  	\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
    \item{alpha}{The Vorob'ev threshold.}
  	\item{current.vorob}{Current value of the vorob criterion (before adding new observations)}
  	\item{penalisation}{Optional penalization constant for type I errors. If equal to zero,  computes the Type II criterion.}
    \item{typeEx}{A character (">" or "<") identifying the type of excursion}

}
\details{
The first argument \code{x} has been chosen to be a vector of size batchsize*d (and not a matrix with batchsize rows and d columns) so that an optimizer like genoud can optimize it easily.
For example if d=2, batchsize=3 and \code{x=c(0.1,0.2,0.3,0.4,0.5,0.6)}, we will evaluate the parallel criterion at the three points (0.1,0.2),(0.3,0.4) and (0.5,0.6).
}

\value{Parallel vorob value}
\references{
Chevalier C., Ginsbourger D., Bect J., Molchanov I. (2013) \emph{Estimating and quantifying uncertainties on level sets using the Vorob'ev expectation and deviation with gaussian process models} mODa 10, Advances in Model-Oriented Design and Analysis, Contributions to Statistics, pp 35-43. \doi{10.1007/978-3-319-00218-7_5}

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern

Azzimonti, D., Ginsbourger, D., Chevalier, C., Bect, J., and Richet, Y. (2021). \emph{Adaptive design of experiments for conservative estimation of excursion sets}. Technometrics, 63(1):13-26. \doi{10.1080/00401706.2019.1693427}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Dario Azzimonti (IDSIA, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_vorob_parallel}} }

\examples{
#vorob_optim_parallel

set.seed(9)
N <- 20 #number of observations
T <- 80 #threshold
testfun <- branin

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design,
	response = response,covtype="matern3_2")

###we need to compute some additional arguments:
#integration points, and current kriging means and variances at these points
integcontrol <- list(n.points=50,distrib="vorob",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,
lower=c(0,0),upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
alpha <- obj$alpha
pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

batchsize <- 4
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
#one evaluation of the vorob_optim_parallel criterion
#we calculate the expectation of the future "vorob" uncertainty
#when 4 points are added to the doe
#the 4 points are (0.1,0.2) , (0.3,0.4), (0.5,0.6), (0.7,0.8)
vorob_optim_parallel(x=x,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,alpha=alpha,current.vorob=Inf)


#the function max_vorob_parallel will help to find the optimum:
#ie: the batch of 4 minimizing the expectation of the future uncertainty

}
