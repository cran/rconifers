

\keyword{misc}
\name{set.variant}
\alias{set.variant}

\title{Set the CONIFERS growth model variant}

\description{Assigns the variant for the CONIFERS forest growth model.}

\usage{
set.variant( var=0 )
}

\arguments{ \item{var}{is the CONIFERS variant. The only possible
  values for the current at the current time are 0
  (\code{SWO_CONIFERS}) and 1 (\code{SMC_CONIFERS}).}
}

\details{

  The \code{set.variant} function only accepts an argument of zero (0)
  or one (1). The default variant is the (\code{SWO\_CONIFERS})
  variant. This function returns a value of 0, 1 or an error message.

}

\value{
  none.
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

Vaughn, Nicholas. 2007. An individual-tree model to predict the annual
growth of young stands of Douglas-fir (Pseudotsuga menziesii (Mirbel)
Franco) in the Pacific northwest. M.S. Thesis, University of
Washington. 91p.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{\code{\link{rconifers}},
 	     \code{\link{smc}},
	     \code{\link{swo}}
}

\examples{
library( rconifers )

## set the variant to the SWO variant
## and calls species.map( swo )
set.variant( 0 )

## set the variant to the SMC variant
## and calls species.map( smc )
set.variant( 1 )

}

\keyword{models}


