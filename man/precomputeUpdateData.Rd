\name{precomputeUpdateData}
\Rdversion{1.4.1}
\alias{precomputeUpdateData}

\title{Useful precomputations to quickly update kriging mean and variance}

\description{
This function is used in combination with \code{\link{computeQuickKrigcov}} and computes an output list that serves as input in that function.
}

\usage{
precomputeUpdateData(model, integration.points)
}

\arguments{
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{integration.points}{
p*d matrix of points for numerical integration in the X space.
}
}

\value{A list with components:
\item{Kinv.c.olddata}{Matrix equal to K^(-1)*c where K is the non conditional covariance matrix at the design points and c is the non conditional covariances between the design points and the integration points.}
\item{Kinv.F}{Matrix equal to K^(-1)*F where F is a matrix with the values of the trend functions at the design points.}
\item{first.member}{Matrix with a complicated expression.}
}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{computeQuickKrigcov}}}

\examples{
#precomputeUpdateData

set.seed(9)
N <- 20 #number of observations
testfun <- branin

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

#the points where we want to compute prediction (if a point new.x is added to the doe)
n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
integration.points <- expand.grid(x.grid,y.grid)
integration.points <- as.matrix(integration.points)
precalc.data <- precomputeUpdateData(model=model,integration.points=integration.points)

#now we can compute quickly kriging covariances 
#between the integration.points and any other points
newdata <- matrix(c(0.6,0.6),ncol=2)
pred <- predict_nobias_km(object=model,newdata=newdata,type="UK",se.compute=TRUE)

kn <- computeQuickKrigcov(model=model,integration.points=integration.points,X.new=newdata,
                    precalc.data=precalc.data,F.newdata=pred$F.newdata,
                    c.newdata=pred$c)

z.grid <- matrix(kn, n.grid, n.grid)

#plots: contour of the covariances, DOE points and new point
#these covariances have been computed quickly with computeQuickKrigcov
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(newdata, col="red", pch=17, lwd=4,cex=3)
title("Kriging covariances with the point (0.6,0.6), in red")
}
