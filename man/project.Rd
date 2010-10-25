
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
                        sdi.mort=0,
			genetic.gains=0) )
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

  \item{genetic.gains}{Non-negative integer. If 0, genetic gains will not be applied to 
    sample.data. If 1, genetic.gains will be applied to the sample.data for each growth
    cycle. This operation is applied independently of the other control variables. 
    This option only applies to Douglas-fir species in the SMC variant. }    

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
	     \code{\link{rand.seed}},
	     \code{\link{rconifers}},
	     \code{\link{sample.data}},	     		  
	     \code{\link{set.species.map}},		  
	     \code{\link{set.variant}},		  
 	     \code{\link{smc}}, 
	     \code{\link{swo}},
	     \code{\link{thin}}
}

\examples{

library( rconifers )

## set the variant to the SWO variant and set species map
set.variant( 0 )

## grow the data SWO variant test data
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"

# now, project the sample forward for 20 years
# with all of the options turned off
sample.swo.23 <- project( sample.swo.3, 20,
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))
print( sample.swo.23 )

## print the projected summaries
print( sample.swo.23 )

# now, project the sample forward for 20 years
# with random error on and duplicate run 100 times, finish by calculating
#summaries of DF basal area
reps<-100
dfba<-vector(mode="double", length=reps)
for ( i in 1:reps){
sample.swo.23wr <- project( sample.swo.3, 20,
	  control=list(rand.err=1,rand.seed=0,endemic.mort=1,sdi.mort=1))
	  
#do a sp.sums of each run
summi<-sp.sums(sample.swo.23wr)

# select the Douglas-fir basal area
dfba[i]<-summi["DF","ba"]
}
summary(dfba)



################################################################################
## now, reset the variant to SMC and grow with genetic gains option
set.variant(1)

## grow the data SMC variant test data
# load and display CONIFERS example plots
data( plots.smc )
print( plots.smc )

# load and display CONIFERS example plants
data( plants.smc )
print( plants.smc )

# create the sample.data list object
sample.smc.3 <- list( plots=plots.smc, plants=plants.smc, age=3, x0=0.0 )
class(sample.smc.3)  <- "sample.data"

# now, project the sample forward for 20 years
# with all of the options turned off
sample.smc.23.0 <- project( sample.smc.3, 20,
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0,genetic.gains=0))
print( sample.smc.23.0 )


## now, project the sample forward for 20 years
## with all of the options turned off
## except for the genetic gains
data( smc )
smc[smc$code == 'DF',]$genetic.worth.d <- 20
smc[smc$code == 'DF',]$genetic.worth.h <- 20
set.species.map( smc )

## project with genetic gains turned on.
sample.smc.23.1 <- project( sample.smc.3, 20,
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0,genetic.gains=1))
print( sample.smc.23.1 )

##get the differences for the species


}

\keyword{models}


