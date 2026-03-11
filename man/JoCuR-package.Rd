\name{JoCuR-package}
\alias{JoCuR-package}
\alias{JoCuR}
\docType{package}
\title{
  Joint Cure Rate Modeling for Longitudinal Tumor Burden and Time-to-Event Data
}
\description{
  \code{JoCuR} provides tools to simulate data and fit a joint cure rate model
  with longitudinal tumor burden and time-to-event outcomes. The package
  includes Monte Carlo EM estimation with C++ acceleration and helper
  functions for marginal longitudinal trajectory summaries.
}
\details{
  The main user-facing workflow is:
  \itemize{
    \item simulate data with \code{CRsimulation()},
    \item estimate model parameters with \code{infer_CRJoint_MLE()},
    \item evaluate longitudinal expectations at a new time using \code{E_y()}.
  }

  The default estimation path uses the identifiable standardized implementation;
  legacy non-standard estimation remains available via the \code{non_std = 1}
  option in \code{infer_CRJoint_MLE()}.
}
\author{
Yixiang Qu

Maintainer: Yixiang Qu <yqu@unc.edu>
}
\references{
  Qu Y, et al. Change-Point Detection Using a Cure Rate Joint Model for
  Longitudinal Tumor Burden and Time-to-Event Data.
}
\keyword{ package }
\seealso{
  \code{\link{CRsimulation}}, \code{\link{infer_CRJoint_MLE}}, \code{\link{E_y}}
}
\examples{
  \dontrun{
    fmla.tte <- as.formula(survival::Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
    fmla.long <- as.formula(PCHG ~ 0 + Y0SCALE)
    sim_data <- CRsimulation(
      fmla.tte = fmla.tte,
      fmla.long = fmla.long,
      beta.tte = c(0.2), normal.tte = TRUE, sd.tte = 0.2,
      beta.y = c(-0.5), sd.y = 0.1, beta.y.cure = c(-0.2), sd.y.cured = 0.2,
      n = 100, censor.parameter = 0.5
    )
  }
}
