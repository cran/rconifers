

\keyword{datasets}
\name{swo}
\alias{swo}
\non_function{}
\title{Species table for the SWO variant of the CONIFERS growth model}
\usage{data(swo)}
\description{
These are species codes used for the SWO variant of the CONIFERS growth model.

The \code{swo} data frame has 38 rows and 13 columns.
}
\format{
This data frame contains the following columns:
\describe{

\item{idx}{The index of the species code that is used within the
  CONIFERS source code.}
\item{code}{The alpha-numeric code that is used during data input,
  normally in the field data recorder. These are used in the
  summaries and outputs rather than the species names.}
\item{fsp}{The functional species code (found in the coefficients)
  that is used to project the plant record. Every species code is
  mapped to a single functional species code. See another manual page.}
\item{name}{The long name for the plant species.}
\item{organon}{The code that should be used to export the species to
  the ORGANON growth and yield model. See \url{http://www.cof.orst.edu/cof/fr/research/organon/index.htm} for more details.}
\item{cactos}{The code that should be used to export the species to
  the CACTOS growth and yield model. See \url{http://www.cnr.berkeley.edu/~wensel/cactos/cactoss.htm} for
  more details.}
\item{fvs}{The code that should be used to export the species to
  the FVS growth and yield model. See \url{http://www.fs.fed.us/fmsc/index.php} for more
  details.}
\item{endemic.mort}{ is a percent for the endemic mortality, where
(0.002) means that two tenths of a percent of the trees of that
functional species type will be "killed" for each year of projection.}
\item{max.sdi}{ is the maximum stand density index the functional
  species can obtain before competition induced mortality is applied.}
\item{browse.damage}{ is a percent for the browse damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of browse damage.}
\item{mechanical.damage}{ is a percent for the mechanical damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of mechanical damage.}
\item{genetic.worth.h}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}
\item{genetic.worth.d}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}
}
}

\details{
  The dataset is similar to the swo.txt file that is distirbuted
  with the GUI version of CONIFERS. A species lookup table is nothing
  more than a data.frame with specific column names. 
}



\references{
  
Ritchie, M. and J. Hamann. 2006. Modeling dynamics of competing
vegetation in young conifer plantations of northern California and
southern Oregon, USA. Canadian Journal of Forest Research 36(10):
2523-2532.

Ritchie, M. and J. Hamann. 2008. Individual-tree height-, diameter-
and crown-width increment equations for young Douglas-fir
plantations. New Forests 35(2):173-186.

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{\code{\link{plants.smc}},
	     \code{\link{plots.smc}}
	     \code{\link{plants.swo}},
	     \code{\link{plots.swo}},
	     \code{\link{set.variant}},
 	     \code{\link{smc}}, 
}

\examples{
library( rconifers )

# load the Southwest-Oregon species coefficients into R as a data.frame object
data( swo )
#view max sdi
swo$max.sdi
#change the max sdi for Douglas-fir to 600
swo$max.sdi[1]<-600

## set the variant to the SWO variant and use the new Douglas-fir max sdi
set.species.map(set.variant( 0 ))


}

\keyword{datasets}


