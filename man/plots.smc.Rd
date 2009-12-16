

\keyword{datasets}
\name{plots.smc}
\alias{plots.smc}
\non_function{}
\title{Example plots data set for CONIFERS growth model}
\usage{data(plots.smc)}

\description{
This a simple plots \code{\link{data.frame}} that contains basic plot attributes used in the CONIFERS growth model.

The \code{plots.smc} data frame has 4 rows and 6 columns.
}

\format{
This data frame contains the following columns:
\describe{
  \item{plot}{A numeric plot identifier.}
  \item{elevation}{Plot elevation, in feet.}
  \item{slope}{Plot slope, in percent.}
  \item{aspect}{Plot aspect, in degrees from North.}
  \item{whc}{Water holding capacity of the plot, in inches.}
  \item{map}{Mean annual precipitation, in inches.}
  \item{si30}{Site Index 30, in meters (optional). It must be present
  to project with variant=1 (SMC_VARIANT).}
}
}

\details{The example plot data represents plot level observations for
a field sample that can be used to generate future forest conditions
using the CONIFERS forest growth model. The data are in imperial
units. The required values are the \code{plot}, \code{elevation}, and
\code{slope}, and \code{aspect} measurements. Values that have not
been observed, or are missing (e.g. whc and map), are denoted with an
\code{NA} and can be imputed using the \code{\link{impute}} function,
which the exception of the \code{si30}. The current version of the
library does not include the functions to compute \code{si30} from
tree observations in the \code{\link{impute}} function.}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

Flewelling, J., R. Collier, B. Gonyea. Marshall D, and Turnblom
E. 2001. Height-Age Curves for Planted Stands of Douglas-fir, with
Adjustments for Density. (Working Paper 1). See 

\url{http://www.cfr.washington.edu/research.smc/working_papers/smc_working_paper_1.pdf}.


}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{impute}},
	     \code{\link{plants.smc}},
	     \code{\link{smc}}
}


\examples{
library( rconifers )
data(plots.smc)
print( plots.smc )
}
