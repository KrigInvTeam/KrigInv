\Rdversion{1.4.1}
\name{vorobVol_optim_parallel2}
\alias{vorobVol_optim_parallel2}
\title{Compute volume criterion}
\usage{
vorobVol_optim_parallel2(x, other.points, integration.points,
  integration.weights = NULL, intpoints.oldmean, intpoints.oldsd,
  precalc.data, model, T, new.noise.var = NULL, batchsize, alpha,
  current.crit, typeEx = ">")
}
\arguments{
\item{x}{One point where to evaluate the criterion}

\item{other.points}{remaining points of the batch}

\item{integration.points}{\eqn{p*d} matrix with the integration points for evaluating numerically the integral of the criterion}

\item{integration.weights}{vector of size \eqn{p} containting the integration weights for the integral numerical evaluation}

\item{intpoints.oldmean}{Vector of size \eqn{p} corresponding to the kriging mean at the integration points before adding the batchsize points \code{x} to the design of experiments.}

\item{intpoints.oldsd}{Vector of size \eqn{p} corresponding to the kriging standard deviation at the integration points before adding the batchsize points \code{x} to the design of experiments.}

\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}

\item{model}{a km Model}

\item{T}{threshold}

\item{new.noise.var}{Optional scalar with the noise variance at the new observation}

\item{batchsize}{size of the batch of new points}

\item{alpha}{threshold on pn obtained with conservative estimates}

\item{current.crit}{starting value for criterion}

\item{typeEx}{a character (">" or "<") identifying the type of excursion}
}
\value{
The value of the criterion at \code{x}.
}
\description{
Compute the volume criterion. Useful for optimization routines.
}
\references{
Azzimonti, D. and Ginsbourger, D. (2018). \emph{Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation}. Journal of Computational and Graphical Statistics, 27(2), 255-267. \doi{10.1080/10618600.2017.1360781}

Azzimonti, D. (2016). \emph{Contributions to Bayesian set estimation relying on random field priors}. PhD thesis, University of Bern.

Azzimonti, D., Ginsbourger, D., Chevalier, C., Bect, J., and Richet, Y. (2021). \emph{Adaptive design of experiments for conservative estimation of excursion sets}. Technometrics, 63(1):13-26. \doi{10.1080/00401706.2019.1693427}

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}
}

\author{
Dario Azzimonti (IDSIA, Switzerland)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_futureVol_parallel}} }

\examples{
#vorobVol_optim_parallel2

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

# alpha, the pn threshold should be computed with conservativeEstimate
# Here it is fixed at 0.992364
alpha <- 0.992364

\donttest{
# You can compute it with the following code
CE_design=as.matrix (randtoolbox::sobol (n = 500*model@d,
                                      dim = model@d))
colnames(CE_design) <- colnames(model@X)

CE_pred = predict.km(object = model, newdata = CE_design,
                          type = "UK",cov.compute = TRUE)
CE_pred$cov <- CE_pred$cov +1e-7*diag(nrow = nrow(CE_pred$cov),ncol = ncol(CE_pred$cov))

Cestimate <- anMC::conservativeEstimate(alpha = 0.95, pred=CE_pred,
                           design=CE_design, threshold=T, pn = NULL,
                           type = ">", verb = 1,
                           lightReturn = TRUE, algo = "GANMC")
alpha <- Cestimate$lvs
}


pred <- predict_nobias_km(object=model,newdata=integration.points,
                          type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

batchsize <- 4
other.points <- c(0.7,0.5,0.5,0.9,0.9,0.8)
x <- c(0.1,0.2)
#one evaluation of the vorobVol_optim_parallel criterion2
#we calculate the expectation of the future volume vorobev uncertainty when
#1+3 points are added to the doe
#the 1+3 points are (0.1,0.2) and (0.7,0.5), (0.5,0.9), (0.9,0.8)
vorobVol_optim_parallel2(x=x,other.points=other.points,integration.points=integration.points,
                      integration.weights=integration.weights,
                      intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
                      precalc.data=precalc.data,T=T,model=model,
                      batchsize=batchsize,alpha=alpha)


n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
vorobVol_parallel.grid <- apply(X=x,FUN=vorobVol_optim_parallel2,MARGIN=1,other.points=other.points,
                             integration.points=integration.points,
                             integration.weights=integration.weights,
                             intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
                             precalc.data=precalc.data,T=T,model=model,
                             batchsize=batchsize,alpha=alpha)
z.grid <- matrix(vorobVol_parallel.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(matrix(other.points,ncol=2,byrow=TRUE), col="red", pch=17, lwd=4,cex=2)

# Note that we want to maximize this criterion.
i.best <- which.max(vorobVol_parallel.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of vorobVol_parallel criterion (black) and of f(x)=T (blue)")



}
