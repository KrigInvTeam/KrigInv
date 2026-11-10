\name{tsee_optim}
\Rdversion{1.4.1}
\alias{tsee_optim}

\title{ Two Sided Expected Exceedance criterion }

\description{
Evaluation of the Two-Sided Expected Exceedance criterion. To be used in optimization routines, like in \code{\link{max_infill_criterion}}.
}
\usage{
tsee_optim(x, model, T, method.param = NULL)
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
Not used. Accepted for compatibility with the other criteria called by \code{\link{max_infill_criterion}}.
}
}

\value{
tsee criterion. 
When the argument \code{x} is a vector the function returns a scalar.
When the argument \code{x} is a p*d matrix the function returns a vector of size p.
}


\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Yann Richet (IRSN, France)
}
\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

\examples{
#tsee_optim

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

x <- c(0.5,0.4)#one evaluation of the tsee criterion
tsee_optim(x=x,T=T,model=model)

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
tsee.grid <- tsee_optim(x=x,T=T,model=model)
z.grid <- matrix(tsee.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,25,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)

i.best <- which.max(tsee.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of tsee criterion (black) and of f(x)=T (blue)")
}
