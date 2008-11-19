
%%	$Id $	

\name{project}
\alias{project}

\title{Projects a sample.data object using the CONIFERS forest growth model}

\description{
  Projects a CONIFERS sample.data object.
}

\usage{
project( x,years=1,control=list(rand.err=0,
                        rand.seed=0,
                        endemic.mort=0,
                        sdi.mort=0) )
}
		   
\arguments{
  \item{x}{a sample.data object.}
  \item{years}{number of years to project the sample.data}
  \item{control}{A list of control parameters. See *Details*}
}

\details{
  The project function returns a sample.data object that contains the
  projected plant list. The plot data are simply copied from the input
  sample.data object, and the x0 member is updated if the sample has
  been projected past the point where stand growth becomes limited by
  competition. 

  By default this function projects the sample.data one year, with the
  default control parameters. The defaults project the sample with a
  random error component added to the height increment, with endemic
  mortality turned on and the stand density index induced mortality
  turned on as well. 

  The 'control' argument is a list that can supply any of the following
  components:

  \describe{
  \item{rand.err}{Non-negative integer. If 0, no random error will be
    applied in the height growth function. If 1, then a random error will
    be added to each plant, in each year to the height growth. This is
    where you include some explanation regarding the random error.}

  \item{rand.seed}{Non-negative integer. If 0, the random number
    generator will be seeded with the time. If >0, the random number
    generator will be seeded with this value. To use random error for each
    projection, you only need to initialize the random number generator
    once. If you wish to repeat the randomness, you need to reseed with
    the same seed, in the same sequence.}
  
  \item{endemic.mort}{Non-negative integer. If 0, no endemic mortality
    will be applied to the sample.data. If 1, the default mortality in the
    mortality settings table will be applied to the sample.data for each
    growth cycle.} 
  
  \item{sdi.mort}{Non-negative integer. If 0, no stand density index mortality
    will be applied to the sample.data. If 1, stand density index
    mortality will be applied to the sample.data for each growth
    cycle. This operation is applied independently of the endemic
    mortality.}    
 }

}

\value{returns a projected sample.data object.}


\references{

Hann, D.W. and C.H. Wang. 1990. Mortality equations for individual
      trees in southwest Oregon. Oregon State University, Forest
      Research Laboratory, Corvallis, Oregon. Research Bulletin
      67. 17p.
  
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

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{impute}},
	     \code{\link{plants}},	
	     \code{\link{plots}} 
	     \code{\link{project}},
	     \code{\link{rand.seed}},
	     \code{\link{rconifers}},
	     \code{\link{sample.data}},	     		  
	     \code{\link{set.species.map}},		  
	     \code{\link{set.variant}},		  
 	     \code{\link{smc}}, 
 	     \code{\link{summary.sample.data}},
	     \code{\link{swo}}, 
	     \code{\link{thin}}
}

\examples{
library( rconifers )

## set the variant to the SWO variant
set.variant( 0 )

# load the Southwest-Oregon species coefficients into R as a data.frame object
data( swo )

# set the species map 
sp.map <- list(idx=swo$idx,
                 fsp=swo$fsp,
                 code=as.character(swo$code),
                 em=swo$endemic.mort,
		 msdi=swo$max.sdi,
                 b=swo$browse.damage,
                 m=swo$mechanical.damage,
                 gwh=swo$genetic.worth.h,
                 gwd=swo$genetic.worth.d)
set.species.map( sp.map )

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots )
print( plots )

# load and display CONIFERS example plants
data( plants )
print( plants )

# randomly remove 10 crown.width observations to test
# the impute function
blanks <- sample( 1:nrow( plants ), 10, replace=FALSE )
plants[blanks,]$crown.width <- NA

# create the sample.data list object
sample.3 <- list( plots=plots, plants=plants, age=3, x0=0.0 )
class(sample.3)  <- "sample.data"

# fill in missing values
sample.3.imp <- impute( sample.3 )

# print the maximum stand density index for the current settings
print( calc.max.sdi( sample.3.imp ) )

# print a summary of the sample
print( sample.3.imp )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.23 <- project( sample.3.imp, 
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.23 )

## plot the diagnostics from the fit a linear dbh-tht model 
## before thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23$plants$tht ~ sample.23$plants$dbh ) )
par( opar )

## thin the stand to capture mortality, redistribute growth, 
## and possibly generate revenue
## Proportional thin for selected tree species, does not remove shrubs
sample.23.t1 <- thin( sample.23, 
	     control=list(type=1, target=50.0, target.sp="DF" ) )
print( sample.23.t1 )

## Proportional thin across all tree species
sample.23.t2 <- thin( sample.23, 
	     control=list(type=2, target=50.0 ) )
print( sample.23.t2 )

## Thin from below, by dbh, all species
sample.23.t3 <- thin( sample.23, 
	     control=list(type=3, target=50.0 ) )
print( sample.23.t3 )

## Thin from below, by dbh for species "PM"
sample.23.t4 <- thin( sample.23, 
	     control=list(type=4, target=50.0, target.sp="PM" ) )
print( sample.23.t4 )

## plot the diagnostics from the fit a linear dbh-tht model
## after proportional thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23.t2$plants$tht ~ sample.23.t2$plants$dbh ) )
par( opar )

## print the differences, by species
print( sp.sums( sample.23.t4 ) - sp.sums( sample.23 ) )

}

\keyword{models}


