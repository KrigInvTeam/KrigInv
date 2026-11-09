\name{sur_optim_parallel2}
\Rdversion{1.4.1}
\alias{sur_optim_parallel2}

\title{ Parallel sur criterion }
\description{
Evaluation of the parallel sur criterion for some candidate points, assuming that some other points are also going to be evaluated. 
To be used in optimization routines, like in \code{\link{max_sur_parallel}}. 
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the expected future sur uncertainty when the candidate points are added.
}
\usage{
sur_optim_parallel2(x, other.points, 
integration.points, integration.weights = NULL, 
intpoints.oldmean, intpoints.oldsd, precalc.data, 
model, T, new.noise.var = NULL, 
batchsize, current.sur,ai_precalc = NULL)
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
\item{new.noise.var}{Optional scalar value of the noise variance of the new observations.}
\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
\item{current.sur}{Current value of the sur criterion (before adding new observations)}
\item{ai_precalc}{When multiple thresholds are used (i.e. when T is a vector), this is an nT*p matrix with ith row equal to \code{intpoints.oldmean-T[i]}. The argument does not need to be filled if only one threshold is used.}
}
\details{
The first argument \code{x} has been chosen to be a vector of size d so that an optimizer like genoud can optimize it easily.
The second argument \code{other.points} is a vector of size (batchsize-1)*d corresponding to the batchsize-1 other points.
}

\value{Parallel sur value}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_sur_parallel}} }

\examples{
#sur_optim_parallel2

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
integcontrol <- list(n.points=50,distrib="sur",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,lower=c(0,0),upper=c(1,1),
model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
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
other.points <- c(0.7,0.5,0.5,0.9,0.9,0.8)
x <- c(0.1,0.2)
#one evaluation of the sur_optim_parallel criterion2
#we calculate the expectation of the future "sur" uncertainty when 
#1+3 points are added to the doe
#the 1+3 points are (0.1,0.2) and (0.7,0.5), (0.5,0.9), (0.9,0.8)
sur_optim_parallel2(x=x,other.points,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,current.sur=Inf,ai_precalc=ai_precalc)

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
sur_parallel.grid <- apply(X=x,FUN=sur_optim_parallel2,MARGIN=1,other.points,
          integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,current.sur=Inf,ai_precalc=ai_precalc)
z.grid <- matrix(sur_parallel.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(matrix(other.points,ncol=2,byrow=TRUE), col="red", pch=17, lwd=4,cex=2)

i.best <- which.min(sur_parallel.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of sur_parallel criterion (black) and of f(x)=T (blue)")
}
