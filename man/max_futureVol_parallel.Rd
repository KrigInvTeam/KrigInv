\name{max_futureVol_parallel}
\Rdversion{1.4.1}
\alias{max_futureVol_parallel}
\title{Maximize parallel volume criterion}
\usage{
max_futureVol_parallel(lower, upper, optimcontrol = NULL, batchsize,
  integration.param, T, model, new.noise.var = 0, typeEx = ">")
}
\arguments{
\item{lower}{lower bounds of the domain}

\item{upper}{upper bounds of the domain}

\item{optimcontrol}{optional list of control parameters for optimization aspects, see \code{\link{max_vorob_parallel}} for details}

\item{batchsize}{size of the batch of new points}

\item{integration.param}{Optional list of control parameter for the computation of integrals, containing the fields \code{integration.points}: a p*d matrix corresponding to p integrations points and \code{integration.weights}: a vector of size p corresponding to the weights of these integration points.
If nothing is specified, default values are used (see: function \code{\link{integration_design}} for more details). }

\item{T}{threshold}

\item{model}{a km Model}

\item{new.noise.var}{Optional scalar with the noise variance at the new observation}

\item{typeEx}{a character (">" or "<") identifying the type of excursion}
}
\value{
A list containing \code{par}, the best set of parameters found, \code{value} the value of the criterion and \code{alpha}, the Vorob'ev quantile corresponding to the conservative estimate.
}
\description{
Maximizes the criterion \code{vorobVol_optim_parallel}.
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

\seealso{\code{\link{EGIparallel}},\code{\link{max_vorob_parallel}}}

\examples{
#max_futureVol_parallel

set.seed(9)
N <- 20 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design,
            response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=200,optim.option=2)
integcontrol <- list(distrib="timse",n.points=400,init.distrib="MC")
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                        lower=lower,upper=upper,model=model,
                                        T=T)

batchsize <- 5 #number of new points

\donttest{
  obj <- max_futureVol_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,
                            batchsize=batchsize,T=T,model=model,
                            integration.param=integration.param)
  #5 optims in dimension 2 !

  obj$par;obj$value #optimum in 5 new points
  new.model <- update(object=model,newX=obj$par,newy=apply(obj$par,1,testfun),
                      cov.reestim=TRUE)

  consLevel = 0.95; n_discrete_design=500*new.model@d
  CE_design=as.matrix (randtoolbox::sobol (n = n_discrete_design,
                              dim = new.model@d))
  colnames(CE_design) <- colnames(new.model@X)
  current.pred = predict.km(object = new.model,
                            newdata = CE_design,
                            type = "UK",cov.compute = TRUE)
  current.pred$cov <- current.pred$cov +1e-7*diag(nrow = nrow(current.pred$cov),
                                                  ncol = ncol(current.pred$cov))

  current.CE = anMC::conservativeEstimate(alpha = consLevel, pred=current.pred,
                                    design=CE_design, threshold=T, pn = NULL,
                                    type = ">", verb = 1,
                                    lightReturn = TRUE, algo = "GANMC")

  par(mfrow=c(1,2))
  print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,
                    cex.points=2.5,main="probability of excursion",consQuantile=obj$alpha)

  print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,
                    new.points=batchsize,col.points.end="red",cex.points=2.5,
                    main="updated probability of excursion",consQuantile=current.CE$lvs)

}
}
