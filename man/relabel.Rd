\name{relabel}
\alias{relabel}
\alias{relabel.matrix}
\alias{relabel.default}
%% \alias{relabel.bcfa}
%% \alias{relabel.cfabugs}
%% \alias{as.mcmc.relabel}

\title{
  Relabel Factor Loadings from MCMC output
}

\description{
  Uses the relabeling method of Erosheva and Curtis to correct for sign
  invariance in MCMC draws from confirmatory factor analysis (CFA) models.
}

\usage{
relabel(obj, \dots)

\method{relabel}{default}(obj, \dots)

\method{relabel}{matrix}(obj, random = FALSE, max.iter = 100, \dots)

%% \method{relabel}{cfabugs}(obj, random = FALSE, max.iter=100, \dots)

%% \method{relabel}{bcfa}(obj, random = FALSE, max.iter=100, \dots)

%% \method{as.mcmc}{relabel}(x)
}

\arguments{
  \item{obj}{

    An object of the appropriate class. Currently, the method is
    implemented only for objects of S3 class \code{matrix}. The matrix
    should contain the MCMC draws from a Bayesian CFA model where no
    arbitrary restrictions are imposed on the loadings to enforce
    rotational invariance. For example, no loadings should have been
    arbitrarily constrained to be one. Columns of the matrix must have
    names of the form "Lam[1,1]", "Lam[1,2]", \ldots, where indices
    denote the structure of the loading matrix, i.e., missing idices
    indicate loadings that are constrained to be zero.

  }

  \item{random}{
    Logical value. If \code{TRUE}, random starting points will be
    used. If \code{FALSE}, starting values are chosen roughly to allow
    the largest loadings on each factor to have a positive sign or to
    allow the majority of loadings on a given factor to have a positive
    sign.
  }
  \item{max.iter}{
    Integer that specifies the maximum number of iterations of the algorithm before aborting procedure.
  }
  %% \item{x}{
  %%   an object of class \code{relabel}.
  %% }
  \item{\dots}{
    Not used.
  }
}

\details{
  See the references.
}

\value{
  An S3 object of class \code{relabel} with the following components:
  \item{Lam}{original (untransformed) draws of the loadings.}
  \item{nuLam}{transformed draws of the loadings.}
  \item{nu.init}{initial values for the sign-change parameters.}
  \item{nu}{final values of the sign change parameters used to compute \code{nuLam}.}
  \item{m}{final values of the mean parameters in the loss function.}
  \item{s}{final values of the standard deviation parameters in the loss function.}
  \item{factor.idx}{vector of integers indicating which columns of the posterior draws belong to which factor.}
  \item{iter}{number of iterations the algorithm took to converge.}
  \item{loss}{final value of the loss function after convergence.}
  \item{converged}{logical indicating whether the algorithm converged.}
}

\references{
  Erosheva, E. A. and Curtis, S. M. (2011)
  "A relabeling scheme for confirmatory factor analysis."
  Technical report #589. University of Washington, Dept. of Statistics.
}

\author{
  S. McKay Curtis and Elena Erosheva
}

\examples{

n <- 1000
p <- 8
set.seed(1)
mu <- as.numeric(t(cbind(
    matrix(rep(sample(c(-1, 1), size=n, replace=TRUE)*4, each=8), n, 8, byrow=TRUE),
    matrix(rep(sample(c(-1, 1), size=n, replace=TRUE)*4, each=8), n, 4, byrow=TRUE),
    matrix(rep(sample(c(-1, 1), size=n, replace=TRUE)*4, each=8), n, 4, byrow=TRUE))))
Lam <- matrix(rnorm(length(mu), mu, 1.0), n, 16, byrow=TRUE)
colnames(Lam) <- c(paste0("Lam[", 1:p, ",", 1, "]"),
                   paste0("Lam[", 1:4, ",", 2, "]"),
                   paste0("Lam[", 5:8, ",", 3, "]"))
\donttest{
par(mfrow=c(2, 8))
apply(Lam, 2, function(x) plot(density(x)))
}
## Relabeling removes the bimodality
out <- relabel(Lam)
\donttest{
par(mfrow=c(2, 8))
apply(out$nuLam, 2, function(x) plot(density(x)))
}

}

\keyword{ manip }
