\name{timse_optim_parallel2}
\Rdversion{1.4.1}
\alias{timse_optim_parallel2}

\title{ Parallel timse criterion }

\description{
Evaluation of the parallel timse criterion for some candidate points, assuming that some other points are also going to be evaluated. 
To be used in optimization routines, like in \code{\link{max_timse_parallel}}. 
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the posterior timse uncertainty.
}
\usage{
timse_optim_parallel2(x, other.points, 
integration.points, integration.weights = NULL, 
intpoints.oldmean, intpoints.oldsd, precalc.data, 
model, T, new.noise.var = NULL,weight = NULL, 
batchsize, current.timse)
}

\arguments{
	\item{x}{Input vector of size d at which one wants to evaluate the criterion. This argument corresponds to only ONE point.}
	\item{other.points}{Vector giving the other \code{batchsize-1} points at which one wants to evaluate the criterion }
  	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points.}
  	\item{intpoints.oldmean}{Vector of size p corresponding to the kriging mean at the integration points before adding \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p corresponding to the kriging standard deviation at the integration points before adding \code{x} to the design of experiments.}
  	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Array containing one or several thresholds.}
  	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
	\item{weight}{
Vector of weight function (length must be equal to the number of lines of the matrix integration.points). 
If nothing is set, the imse criterion is used instead of timse. It corresponds to equal weights.
}
  	\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
  	\item{current.timse}{Current value of the timse criterion (before adding new observations)}

}
\details{
The first argument \code{x} has been chosen to be a vector of size d so that an optimizer like genoud can optimize it easily.
The second argument \code{other.points} is a vector of size (batchsize-1)*d corresponding to the batchsize-1 other points.
}

\value{Parallel timse value}
\references{
Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7). \doi{10.1115/1.4001873}

Picheny V. (2009) \emph{Improving accuracy and compensating for uncertainty in surrogate modeling}, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}
}

\author{
Victor Picheny (INRA, Toulouse, France)

Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_timse_parallel}} }

\examples{
#timse_optim_parallel2

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
integcontrol <- list(n.points=1000,distrib="timse",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,lower=c(0,0),upper=c(1,1),
model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

#we also need to compute weights. Otherwise the (more simple) 
#imse criterion will be evaluated
weight0 <- 1/sqrt( 2*pi*(intpoints.oldsd^2) )
weight <- 0
for(i in 1:length(T)){
  Ti <- T[i]
  weight <- weight + weight0 * exp(-0.5*((intpoints.oldmean-Ti)/sqrt(intpoints.oldsd^2))^2)
}

batchsize <- 4
other.points <- c(0.7,0.5,0.5,0.9,0.9,0.8)
x <- c(0.1,0.2)
#one evaluation of the timse_optim_parallel criterion2
#we calculate the expectation of the future "timse" uncertainty 
#when 1+3 points are added to the doe
#the 1+3 points are (0.1,0.2) and (0.7,0.5), (0.5,0.9), (0.9,0.8)
timse_optim_parallel2(x=x,other.points,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,weight=weight,
          batchsize=batchsize,current.timse=Inf)

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
timse_parallel.grid <- apply(X=x,FUN=timse_optim_parallel2,MARGIN=1,other.points,
          integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,weight=weight,
          batchsize=batchsize,current.timse=Inf)
z.grid <- matrix(timse_parallel.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(matrix(other.points,ncol=2,byrow=TRUE), col="red", pch=17, lwd=4,cex=2)

i.best <- which.min(timse_parallel.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of timse_parallel criterion (black) and of f(x)=T (blue)")
}
