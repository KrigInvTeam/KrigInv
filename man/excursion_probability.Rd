\name{excursion_probability}
\Rdversion{1.4.1}
\alias{excursion_probability}

\title{Excursion probability with one or many thresholds}

\description{
Probability that Gaussian random variables with some mean and variance are over a threshold \code{T}, or in an union of intervals.
If \code{T} is a vector of size p, T1,T2,...,Tp then the considered union of interval is (T1,T2) U ... U (Tp, +infty) if p is odd, and (T1,T2) U ... U (Tp-1, Tp) if p is even.
}

\usage{
excursion_probability(mn,sn,T)
}

\arguments{
 
  \item{mn}{Array of size k containing the expectations of the Gaussian random variables.}
	\item{sn}{Array of size k containing the standard deviations of the Gaussian random variables.}
 	\item{T}{Array containing one or several thresholds.}
}


\value{Array of size k containing the k excursion probabilities.}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{predict_nobias_km}}}


\examples{
#excursion_probability

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

some_points <- matrix(runif(20),ncol=2)
pred <- predict_nobias_km(object = model,newdata = some_points,
                type = "UK",se.compute = TRUE)
                
T <- c(60,80,100)
excursion_probability(mn = pred$mean,sn = pred$sd,T=T)
# probability to be in the interval [60,80] U [100, infty]
}