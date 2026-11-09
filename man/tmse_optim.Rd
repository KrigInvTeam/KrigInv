\name{tmse_optim}
\Rdversion{1.4.1}
\alias{tmse_optim}

\title{ Targeted MSE criterion }

\description{
Evaluation of the Targeted MSE criterion. To be used in optimization routines, 
like in \code{\link{max_infill_criterion}}
}

\usage{
tmse_optim(x, model, T, method.param = NULL)
}

\arguments{
  \item{x}{
Input vector at which one wants to evaluate the criterion. 
This argument can be either a vector of size d (for an evaluation at a single point) or a p*d matrix (for p simultaneous evaluations of the criterion at p different points).
}
  \item{model}{
An object of class \code{\link[DiceKriging]{km}} (Kriging model).
}
  \item{T}{
Array containing one or several thresholds.
}
  \item{method.param}{
Scalar tolerance around the targets T.
}
}

\value{targeted MSE value.
When the argument \code{x} is a vector the function returns a scalar.
When the argument \code{x} is a p*d matrix the function returns a vector of size p.
}

\references{
Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7). \doi{10.1115/1.4001873}

Picheny V. (2009) \emph{Improving accuracy and compensating for uncertainty in surrogate modeling}, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne
}

\author{ 
Victor Picheny (INRA, Toulouse, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)

Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

\examples{
#tmse_optim

set.seed(9)
N <- 20 #number of observations
T <- c(40,80) #thresholds
testfun <- branin

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

x <- c(0.5,0.4)#one evaluation of the tmse criterion
tmse_optim(x=x,T=T,model=model)

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
tmse.grid <- tmse_optim(x=x,T=T,model=model)
z.grid <- matrix(tmse.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,25,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)

i.best <- which.max(tmse.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of tmse criterion (black) and of f(x)=T (blue)")
}
