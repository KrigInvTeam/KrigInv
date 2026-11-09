\name{jn_optim_parallel}
\Rdversion{1.4.1}
\alias{jn_optim_parallel}

\title{ Parallel jn criterion }
\description{
Evaluation of the parallel jn criterion for some candidate points. To be used in optimization routines, like in \code{\link{max_sur_parallel}}. 
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the posterior expected \code{"jn"} uncertainty, which is the posterior expected variance of the excursion set's volume.
}
\usage{
jn_optim_parallel(x, integration.points, integration.weights = NULL, 
intpoints.oldmean, intpoints.oldsd, 
precalc.data, model, T, 
new.noise.var = NULL, batchsize, current.sur, ai_precalc)
}

\arguments{
  	\item{x}{Vector of size batchsize*d at which one wants to evaluate the criterion. This argument is NOT a matrix.}
  	\item{integration.points}{Matrix of points for numerical integration. Two cases are handled. If the \code{"jn"} importance sampling distribution has been used, this is a (2p)*d matrix containing p couples of integration points in D x D, where D is the integration domain. If instead a p*d matrix is given, then the function will perform the integration by using all \code{p^2} couples.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points.}
  	\item{intpoints.oldmean}{Vector of size p, or 2p, corresponding to the kriging mean at the integration points before adding the batchsize points \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p, or 2p, corresponding to the kriging standard deviation at the integration points before adding  the batchsize points \code{x} to the design of experiments.}
  	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Array containing one or several thresholds.}
  	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
  	\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
  	\item{current.sur}{Current value of the sur criterion (before adding new observations).}
  	\item{ai_precalc}{This is a matrix with ith row equal to \code{intpoints.oldmean-T[i]}. The argument DOES ALWAYS need to be filled, even if there is only one threshold T.}
}
\details{
The first argument \code{x} has been chosen to be a vector of size batchsize*d (and not a matrix with batchsize rows and d columns) so that an optimizer like genoud can optimize it easily.
For example if d=2, batchsize=3 and \code{x=c(0.1,0.2,0.3,0.4,0.5,0.6)}, we will evaluate the parallel criterion at the three points (0.1,0.2),(0.3,0.4) and (0.5,0.6).
}

\value{Parallel jn value}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}
\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_sur_parallel}} }

\examples{
#jn_optim_parallel

set.seed(9)
N <- 20 #number of observations
T <- c(80,100) #thresholds
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
n.points <- 200
integcontrol <- list(n.points=n.points,distrib="jn",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,
lower=c(0,0),upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points # (2n.points)*d matrix
integration.weights <- obj$integration.weights #vector of size n.points
pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)
nT <- 2 # number of thresholds
ai_precalc <- matrix(rep(intpoints.oldmean,times=nT),
    nrow=nT,ncol=length(intpoints.oldmean),byrow=TRUE)
ai_precalc <- ai_precalc - T  # substracts Ti to the ith row of ai_precalc

batchsize <- 4
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
#one evaluation of the sur_optim_parallel criterion
#we calculate the expectation of the future "sur" uncertainty 
#when 4 points are added to the doe
#the 4 points are (0.1,0.2) , (0.3,0.4), (0.5,0.6), (0.7,0.8)
jn_optim_parallel(x=x,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,current.sur=0,ai_precalc=ai_precalc)

# the criterion takes a negative value, which is normal.
# See the Technometrics paper in the references

#the function max_sur_parallel will help to find the optimum: 
#ie: the batch of 4 minimizing the expectation of the future uncertainty
}
