\name{bichon_optim}
\Rdversion{1.4.1}
\alias{bichon_optim}

\title{
Bichon et al.'s Expected Feasibility criterion
}

\description{
Evaluation of Bichon's Expected Feasibility criterion. To be used in optimization routines, like in \code{\link{max_infill_criterion}}.
}

\usage{
bichon_optim(x, model, T, method.param = 1)
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
Target value (scalar).
}
  \item{method.param}{
Scalar tolerance around the target T. Default value is 1.
}
}


\value{
Bichon EF criterion.
When the argument \code{x} is a vector, the function returns a scalar.
When the argument \code{x} is a p*d matrix, the function returns a vector of size p.
}

\references{
Bichon B.J., Eldred M.S., Swiler L.P., Mahadevan S., McFarland J.M. (2008) \emph{Efficient global reliability analysis for nonlinear implicit performance functions}, AIAA Journal 46(10), pp 2459-2468. \doi{10.2514/1.34321}
}
\author{ 
Victor Picheny (INRA, Toulouse, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)

Clement Chevalier (University of Neuchatel, Switzerland)
}


\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

\examples{
#bichon_optim

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

x <- c(0.5,0.4) #one evaluation of the bichon criterion
bichon_optim(x=x,T=T,model=model)

n.grid <- 20 # resolution. You may use a larger value.
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
bichon.grid <- bichon_optim(x=x,T=T,model=model)
z.grid <- matrix(bichon.grid, n.grid, n.grid)

#plots: contour of the criterion, DOE points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,25,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)

i.best <- which.max(bichon.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of Bichon criterion (black) and of f(x)=T (blue)")
}
