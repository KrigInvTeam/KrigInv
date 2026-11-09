\name{vorob_optim_parallel2}
\Rdversion{1.4.1}
\alias{vorob_optim_parallel2}

\title{ Parallel Vorob'ev criterion }
\description{
Evaluation of the Vorob'ev criterion for some candidate points, assuming that some other points are also going to be evaluated. To be used in optimization routines, like in \code{\link{max_vorob_parallel}}.
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise.
The criterion is the integral of the posterior Vorob'ev uncertainty.
}
\usage{
vorob_optim_parallel2(x, other.points,
integration.points, integration.weights = NULL,
intpoints.oldmean, intpoints.oldsd, precalc.data,
model, T, new.noise.var = NULL,
batchsize, alpha, current.vorob,
penalisation = NULL, typeEx = ">")
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
	\item{T}{Target value (scalar). The criterion CANNOT be used with multiple thresholds.}
  	\item{new.noise.var}{Optional scalar value of the noise variance of the new observations.}
  	\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
  	\item{alpha}{The Vorob'ev threshold.}
    \item{current.vorob}{Current value of the vorob criterion (before adding new observations)}
    \item{penalisation}{Optional penalization constant for type I errors. If equal to zero,  computes the Type II criterion.}
    \item{typeEx}{A character (">" or "<") identifying the type of excursion}
}
\details{
The first argument \code{x} has been chosen to be a vector of size d so that an optimizer like genoud can optimize it easily.
The second argument \code{other.points} is a vector of size (batchsize-1)*d corresponding to the batchsize-1 other points.
}

\value{Parallel Vorob'ev value}
\references{
Chevalier C., Ginsbourger D., Bect J., Molchanov I. (2013) \emph{Estimating and quantifying uncertainties on level sets using the Vorob'ev expectation and deviation with gaussian process models} mODa 10, Advances in Model-Oriented Design and Analysis, Contributions to Statistics, pp 35-43. \doi{10.1007/978-3-319-00218-7_5}

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Dario Azzimonti (IDSIA, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_vorob_parallel}} }

\examples{
#vorob_optim_parallel2

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
other.points <- c(0.7,0.5,0.5,0.9,0.9,0.8)
x <- c(0.1,0.2)
#one evaluation of the vorob_optim_parallel criterion2
#we calculate the expectation of the future "vorob" uncertainty when
#1+3 points are added to the doe
#the 1+3 points are (0.1,0.2) and (0.7,0.5), (0.5,0.9), (0.9,0.8)
vorob_optim_parallel2(x=x,other.points,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,alpha=alpha,current.vorob=Inf)


n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
vorob_parallel.grid <- apply(X=x,FUN=vorob_optim_parallel2,MARGIN=1,other.points,
          integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,alpha=alpha,current.vorob=Inf)
z.grid <- matrix(vorob_parallel.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(matrix(other.points,ncol=2,byrow=TRUE), col="red", pch=17, lwd=4,cex=2)

i.best <- which.min(vorob_parallel.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of vorob_parallel criterion (black) and of f(x)=T (blue)")
}
