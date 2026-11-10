\name{predict_update_km_parallel}
\Rdversion{1.4.1}
\alias{predict_update_km_parallel}

\title{
Quick update of kriging means and variances when one or many new points are added to the DOE
}

\description{
The functions uses the kriging update formulas to quickly compute kriging mean and variances at points newdata, when \code{r} new points \code{newX} are added. 
}

\usage{
predict_update_km_parallel(newXmean, newXvar, newXvalue, 
Sigma.r, newdata.oldmean, newdata.oldsd, kn)
}

\arguments{
  \item{newXmean}{
Vector of size r: old kriging mean at points x_(n+1),...,x_(n+r).
}
  \item{newXvar}{
Vector of size r: kriging variance at points x_(n+1),...,x_(n+r).
}
  \item{newXvalue}{
Vector of size r: value of the objective function at x_(n+1),...,x_(n+r).
}
  \item{Sigma.r}{
An r*r matrix: kriging covariances between the points x_(n+1),...,x_(n+r).
}
  \item{newdata.oldmean}{
Vector: old kriging mean at the points \code{newdata} (before adding x_(n+1),...,x_(n+r))
}
  \item{newdata.oldsd}{
Vector: old kriging standard deviations at the points \code{newdata} (before adding x_(n+1),...,x_(n+r))
}
  \item{kn}{
Kriging covariances between the points \code{newdata} and the r points \code{newX}. These covariances can be computed using the function \code{\link{computeQuickKrigcov}}
}
}

\value{A list with the following fields:
\item{mean}{Updated kriging mean at points \code{newdata}}
\item{sd}{Updated kriging standard deviation at points \code{newdata}}
\item{lambda}{New kriging weight of x_(n+1),...,x_(n+r) for the prediction at points \code{newdata}}
}

\references{
Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{computeQuickKrigcov}}, \code{\link{precomputeUpdateData}}}

\examples{
#predict_update_km_parallel

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

#points where we want to compute prediction (if a point new.x is added to the doe)
n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
newdata <- expand.grid(x.grid,y.grid)
newdata <- as.matrix(newdata)
precalc.data <- precomputeUpdateData(model=model,integration.points=newdata)
pred2 <- predict_nobias_km(object=model,newdata=newdata,type="UK",se.compute=TRUE)
newdata.oldmean <- pred2$mean; newdata.oldsd <- pred2$sd

#the point that we are going to add
new.x <- matrix(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),ncol=2,byrow=TRUE)   
pred1 <- predict_nobias_km(object=model,newdata=new.x,type="UK",
se.compute=TRUE,cov.compute=TRUE)
newXmean <- pred1$mean; newXvar <- pred1$sd^2; newXvalue <- pred1$mean + 2*pred1$sd
Sigma.r <- pred1$cov

kn <- computeQuickKrigcov(model=model,integration.points=newdata,X.new=new.x,
                    precalc.data=precalc.data,F.newdata=pred1$F.newdata,
                    c.newdata=pred1$c)

updated.predictions <- predict_update_km_parallel(newXmean=newXmean,newXvar=newXvar,
                                         newXvalue=newXvalue,Sigma.r=Sigma.r,
                                         newdata.oldmean=newdata.oldmean,
                                         newdata.oldsd=newdata.oldsd,kn=kn)

#the new kriging variance is usually lower than the old one
updated.predictions$sd - newdata.oldsd 

z.grid1 <- matrix(newdata.oldsd, n.grid, n.grid)
z.grid2 <- matrix(updated.predictions$sd, n.grid, n.grid)

par(mfrow=c(1,2))

image(x=x.grid,y=y.grid,z=z.grid1,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid1,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
title("Kriging standard deviation")

image(x=x.grid,y=y.grid,z=z.grid2,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid2,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(new.x, col="red", pch=17, lwd=4,cex=2)
title("updated Kriging standard deviation")
}
