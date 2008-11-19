/********************************************************************************/
/*                                                                              */
/*  smc_model.c                                                                 */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/

/* 	$Id: smc_model.c 620 2009-01-05 17:25:02Z mritchie $	 */

//#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"


void get_age_cut(           
    unsigned long           *return_code,
    double                  h40,                 
    double                  si_30,
    unsigned long           *genetics_age_cut)
{

    *return_code=CONIFERS_SUCCESS;
    *genetics_age_cut = 0;

    if(si_30 <= 0.0 || h40 <= 0.0)
    {
        *return_code=CONIFERS_ERROR;
        return;
    }
    
    if(si_30 < 45.0)   /* is age > 10? (1) is age > 15? (2) */
    {
        if(h40 > 1.1) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 5.9) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 12.9)
        {
            *genetics_age_cut = 3;
        }
    }
    else if(si_30 < 55.0)
    {
        if(h40 > 2.3) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 7.7) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 17.1)
        {
            *genetics_age_cut = 3;
        }
    }
    else if(si_30 < 65.0)
    {
        if(h40 > 2.6) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 9.5) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 21.3)
        {
            *genetics_age_cut = 3;
        }
     }
     else if(si_30 < 75.0)
     {
        if(h40 > 3.0) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 11.9) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 26.1)
        {
            *genetics_age_cut = 3;
        }
     }
     else if(si_30 < 85.0)
     {
        if(h40 > 3.7) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 15.1) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 31.4)
        {
            *genetics_age_cut = 3;
        }
     }
     else
     {
        if(h40 > 4.9) 
        {
            *genetics_age_cut = 1;
        }
        if(h40 > 19.45) 
        {
            *genetics_age_cut = 2;
        }
        if(h40 > 37.2)
        {
            *genetics_age_cut = 3;
        }
     }

     return;
}

/********************************************************************************/
/*                  smc_calc_height_growth    D1                                */
/********************************************************************************/
/*  Description :   calc_height_growth                                          */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   August 28, 2007                                              */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  flew_site                   -   Flewellings site index                      */
/*  con_tpa                     -   conifer trees per acre                      */
/*  h40                         -   height of 40 tallest trees per acre         */
/*  cash                        -   crown area shrubs                           */
/*  random_norm_0_1             -   normal random variable with mean 0, var 1   */ 
/*  random_unif_0_1             -   uniform random var. with range 0, 1         */
/*  random_unif_0_1a            -   a second unif ran. var with range 0, 1      */
/*  ind_random                  -   indicator for random error (0 or 1) 1=yes   */
/*  prob_browse                 -   probability of browse, only if h< 4.5       */
/*  prob_top_damage             -   probability of top damage                   */
/* -----------------------------------------------------------------------------*/
/*  *height_growth              -   predicted height growth                     */
/*  vector<double> *coeffs_ptr  -   pointer to a vector of doubles that         */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : hg =  base * rh_mod * veg_mod * tpa_mod                           */
/*                                                                              */
/*  base=1.0/(b1 + b3*((si^-b2)/total_height) + b4*si^1-b2) + b5*total_height^b2-1)*/
/*  rh_mod  = 2/( 1 + exp(-(h0 + h1 *(pow(h40, 0.5))) * log(rel_ht)) )          */
/*  veg_mod = exp(-exp(v1+v2*h40)*sqrt(cash))                                   */
/*  tpa_mod = 1 + d1*(h40-d2)*(con_tpa-300)                                     */
/*                                                                              */
/*                                                                              */
/*  Source  : Vaughns Thesis, with modifications 11/2007                        */
/*  Coeffs  : HG                                                                */
/********************************************************************************/
void smc_calc_height_growth(
	unsigned long *return_code,	
    double  total_height, 
    double  basal_d,
    double  flew_site,
    double  con_tpa,
    double  h40,
    double  cash,
	double  catcon,
    double  cathw,
    double  catsh,
    double  random_norm_0_1,
    double  random_unif_0_1a,
    double  random_unif_0_1b,
    long    ind_random,
    double  prob_browse,
    double  prob_top_damage,

	//double  genetic_worth,
	
	double		    genetic_worth_h,
    double          genetic_worth_d,

	unsigned long use_genetic_gains,

	double  *height_growth,
    double  *coeffs_ptr,
	unsigned long plant_type,
    unsigned long genetics_age_cut)
{

/****************************  right now, not all are used     ******************************************/
    double  b0,  b1,  b2,  b3,  b4,  b5,  v1,  v2,  d1,  d2,
            h0,  h1,  b6, b12, b13, b14, b15;

    double  height_var;
    double  height_for_error;
    double  total_crown_area;
    int     broken;
    int     browsing;
    double  rel_ht;
    double  base;
    double  rh_mod;
    double  veg_mod;
    double  tpa_mod;
    double  hg;
    double  genetic_modifier;
    double  cat;
/**************************** initialize variables *******************************************************/
	*return_code=CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *height_growth = 0.0;
	    *return_code=CONIFERS_ERROR;
        return;
    }

    *height_growth = 0.0;                                 /* set to zero for default                     */
    height_var       = random_norm_0_1;                   /* height growth error by default is N(0,1)    */
    height_for_error = total_height;
    rel_ht           = total_height/h40;

    b0                  = 0.0;
    b1                  = 0.0;
    b2                  = 0.0;
    b3                  = 0.0;
    b4                  = 0.0;
    b5                  = 0.0;
    v1                  = 0.0;
    v2                  = 0.0;
    d1                  = 0.0;
    d2                  = 0.0;
    h0                  = 0.0;
    h1                  = 0.0;
    b6                  = 0.0;
    b12=0.0;
    b13=0.0;
    b14=0.0;
    b15=0.0;
    total_crown_area    = 0.0;
    broken              = 0;
    browsing            = 0;
    genetic_modifier    = 1.0;

    cat=(catcon+cathw+catsh)/SQ_FT_PER_ACRE;
/****************************** check for some errors *******************************************************/
    if( total_height < 0.0)
    {
        *height_growth    = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }

    if( cash < 0.0)
    {
        cash               = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }
/* convert cash to cover in shrubs as percent   */
     cash=100.0*cash/SQ_FT_PER_ACRE;

    if(total_height <= 0.5)                                       
    {
        *height_growth = 0.20;
		return;
    }

/****************************** set up the browse coefficient (indicator) ***********************************/
    if(prob_browse > random_unif_0_1a && total_height <= 4.5)
    {
        browsing = 1;
    }

/******************************  set up the top damage coefficient (indicator) ******************************/
    if(prob_top_damage > random_unif_0_1b)
    {
        broken   = 1;
    }

/******************************  turn off the random error component if ind_random=0 ************************/
    if(ind_random == 0) 
    {
	    height_var = 0.0;   
    }

/****************************** if broken then shut off random error ****************************************/
    if (broken) 
    {
	    height_var = 0.0;
    }
	
/****************************** if browsed shut off random error ********************************************/
    if (browsing) 
    {
	    height_var = 0.0;
    }

/*********************  the next stuff will keep random error acting ok. ***********************************/
/*********************	first we will make sure that we dont get a huge error term, then we  ***************/
/*********************	will limit the random normal component *********************************************/
    if(total_height > 15.0 )
    {
 	    height_for_error = 15.0;
    }

    if(height_var < -1.645)                                         /*  keep random normal within reason     */
    {
        height_var = 0.0;                                        /*  if large negative, set to -2.0       */
    }

    if(height_var > 1.645)
    {
        height_var = 0.0;                                         /*   if large positive, set to +2.0      */
    }
    
/* set the coefficients for the model  base b0-b5, veg v1 v2, tpa d1 d2, relht h0 and h1 */
/************************************************************************************************************/
    if ( plant_type == CONIFER || plant_type == HARDWOOD ) 
    {
/********************* this block for genetic gains stuff*******************************************/
        if(use_genetic_gains && plant_type == CONIFER)
        {
            /* set modifier from Gould et al.*/
            if(genetics_age_cut == 0 || genetic_worth_h < 0.0 || genetic_worth_h > 20.0)
            {
                genetic_modifier = 1.00 + 0.0062770*genetic_worth_h; /* stand is <5    */
            }
            else if(genetics_age_cut == 1)
            {
                genetic_modifier = 1.00 + 0.0062770*genetic_worth_h; /* stand is 5-10  */
            }
            else if( genetics_age_cut == 2)
            {
                genetic_modifier = 1.00 + 0.0031119*genetic_worth_h; /* stand is 10-15 */
            }
            else 
            {
                genetic_modifier = 1.00 + 0.0041740*genetic_worth_h; /* stand is 15+   */
            }
        }
/************************* end of genetic gains stuff **********************************************/
		b0  = coeffs_ptr[0];
		b1  = coeffs_ptr[1];
		b2  = coeffs_ptr[2];
		b3  = coeffs_ptr[3];
		b4  = coeffs_ptr[4];
		b5  = coeffs_ptr[5];
		v1  = coeffs_ptr[6];   /* made adjustment in coefficients*/
		v2  = coeffs_ptr[7];
		d1  = coeffs_ptr[8]/100000;
		d2  = coeffs_ptr[9];
		h0  = coeffs_ptr[10];
		h1  = coeffs_ptr[11];

		/* Nick Vaughns base height growth model for doug fir */
		/* base model check nov 15 2007 */
		/* rhmod checked on oct 22 2007 */
		/* veg mod checked on nov 15 2007 */
		base    = 1.0/(b1 + b3*(pow(flew_site, -b2)/total_height) + b4*pow(flew_site, 1-b2) + b5*pow(total_height, b2-1));

	    rh_mod  = 2.0/( 1 + exp(-(h0 + h1 *(pow(h40, 0.5))) * log(rel_ht)) );

	    veg_mod = exp(-exp(v1+v2*h40)*sqrt(cash));

	    tpa_mod = 1 + d1*(h40-d2)*(con_tpa-300);
    
		hg = base * rh_mod * veg_mod * tpa_mod * genetic_modifier;
    
    
	    b12 = coeffs_ptr[12] * (double)browsing;
		b13 = coeffs_ptr[13] * (double)broken ;                         
		b14 = coeffs_ptr[14] * height_var;                             
		b15 = coeffs_ptr[15] * (sqrt(height_for_error) * height_var)/2.0 ; 
		
        if((hg + b14 + b15) >=0.0)
        {
            *height_growth =   hg  
							+ (b12 + b13)   /* damage adjustments                    */
							+ (b14 + b15);  /* error term                            */ 
        }
        else
        {
            *height_growth = 0.0 + b12 + b13;
        }
    }
    else if (plant_type == SHRUB )   
    {
        b0 = coeffs_ptr[0];
        b1 = coeffs_ptr[1];
        b2 = coeffs_ptr[2];
	    b3 = coeffs_ptr[3];                                   /*          extra                        */
	    b4 = coeffs_ptr[4];
		
        *height_growth= b0/total_height 
                  + exp(b1 + b2*(log(total_height)) + b3*total_height*basal_d + b4*cat*cat);
    }
	else /* its non stocked or forb */
	{
		*height_growth=0.0;
	}
/************************************************************************************************************/
	if( total_height + *height_growth <= 0.5)                       /* if predicted height is <0.5 then        */
    {
        *height_growth = (-1.0)*(total_height-0.51 ) ;
    }
    return;
}


/********************************************************************************/
/*                  smc_calc_dbh_growth        D2                               */
/********************************************************************************/
/*  Description :   smc_calc_dbh_growth                                         */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   August 28,  2007                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code      - return code for calling function to check */
/*  double         total_height     - total tree height                         */
/*  double         h40              - predicted height increment (annual)       */
/*  double         crown_ratio      - crown ratio (currently not used)          */
/*  double         current_dbh      - initial dbh (inches)                      */
/*  double         current_d6       - initial d6 (inches)                       */
/*  double         d6_growth        - predicted d6 growth (inches)              */
/*  double         *pred_dbh_growth - predicted dbh growth (inches)             */
/*  vector<double> *coeffs_ptr      -   pointer to a vector of doubles that     */
/*                                      contain the coefficients for the        */
/*                                      functional species code                 */
/********************************************************************************/
/*  Formula : dbhgro =                                                          */ 
/*  Source  : Vaughns Thesis 08/2007                                            */
/*  Coeffs  : DG                                                                */
/********************************************************************************/
void smc_calc_dbh_growth(   
    unsigned long   *return_code,
    double          total_height,
    double          h40,
    double          flew_site,
    double          current_dbh,
    double          basal_area,
    double          contpa,
    
//	double          genetic_worth,

	double		    genetic_worth_h,
    double          genetic_worth_d,


    unsigned long   genetics_age_cut,
    unsigned long   use_genetic_gains,
    double          *pred_dbh_growth,
    double          *coeffs_ptr,
  	unsigned long   plant_type)
{

    double  temp_dbh;
    double  tnum;
    double  tden;
    double  dhrat;
    double  hrat;
    double  b0;
    double  b1;
    double  b2;
    double  b3;
    double  b4;
    double  b5;
    double  b6;
    double  b7;
    double  b8;
    double  b9;
    double  d1;
    double  d2;
    double  e1;
    double  genetic_modifier;
    unsigned long   big_i;  /* large tree indicator this will only effect growth for trees beyond the range of data */

	*return_code=CONIFERS_SUCCESS;

    big_i = 0;
    if(current_dbh > 9.0)
    {
        big_i = 1;
    }
    if( coeffs_ptr == NULL )
    {
        *pred_dbh_growth    = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }

    if( total_height < 4.5)
    {
        *pred_dbh_growth    =0.0;
        return; 
    }

    b0  = 0.0;
    b1  = 0.0;
    b2  = 0.0;
    b3  = 0.0;
    b4  = 0.0;
    b5  = 0.0;
    b6  = 0.0;
    b7  = 0.0;
    b8  = 0.0;
    b9  = 0.0;
    d1  = 0.0;
    d2  = 0.0;
    genetic_modifier = 1.0;

    if(use_genetic_gains && plant_type == CONIFER)
    {
       /* set modifier from Gould et al. */
       if( genetics_age_cut == 0 || genetic_worth_d < 0.0 || genetic_worth_d > 20.0)
       {
           genetic_modifier = 1.00 + 0.0101054*genetic_worth_d;   /* stand is below 5 */                  
       }
       else if(genetics_age_cut == 1)
       {
           genetic_modifier = 1.00 + 0.0101054*genetic_worth_d; /* stand is 5-10  */
       }
       else if( genetics_age_cut == 2)
       {
           genetic_modifier = 1.00 + 0.0033696*genetic_worth_d; /* stand is 10-15 */
       }
       else 
       {
           genetic_modifier = 1.00 + 0.0029435*genetic_worth_d; /* stand is 15+   */
       }
    }


    dhrat=current_dbh/total_height;
    hrat=total_height/h40;


    temp_dbh    = 0.0;
    
    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    b2  = coeffs_ptr[2];
    b3  = coeffs_ptr[3];
    b4  = coeffs_ptr[4];
    b5  = coeffs_ptr[5];
    b6  = coeffs_ptr[6];
    b7  = coeffs_ptr[7];
    b8  = coeffs_ptr[8];
    b9  = coeffs_ptr[9];
    d1  = coeffs_ptr[10]/100000;
    d2  = coeffs_ptr[11];
    e1  = coeffs_ptr[12];    /* this is a correction to force a peaking function beyond the range of dbh */
                             /*  to get this I basically graphically derived a term that would mimic approx. */
                             /*  the peak in Hann et al for dbh growth */

    tnum= b1*pow(current_dbh, b2)*exp(b3*current_dbh + b4*basal_area + b5*dhrat + b6*hrat + b7*h40 
                                    + e1*big_i*pow(current_dbh,2.5)) ; /*correction for big trees beyond data */
    tden= 1 + exp(b8 - b9*flew_site);
    temp_dbh= (tnum / tden) * (1 + d1*(h40-d2) * (contpa-300));

    if (temp_dbh <= 0.0)
    {
	    *pred_dbh_growth    = 0.0;
	    *return_code        =CONIFERS_ERROR;
    }

    *pred_dbh_growth        = (sqrt(current_dbh*current_dbh+temp_dbh) - current_dbh) * genetic_modifier;
    *return_code            = CONIFERS_SUCCESS;

//    if(current_dbh > 11.5)
//    {
//        *pred_dbh_growth = 0.2*exp(5.0 - 5.5 + 0.40*log(current_dbh-1.0) - 0.00044*current_dbh*current_dbh);
//        *return_code     = CONIFERS_SUCCESS;
//    }

    if( *pred_dbh_growth < 0 )
    {
        *pred_dbh_growth    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}





/********************************************************************************/
/*                  smc_calc_cr_growth          D5                              */
/********************************************************************************/

void smc_calc_cr_growth(
    unsigned long   *return_code,
    int             hcb_growth_on,
    double          total_height,
    double          height_growth,
    double          crown_ratio,
    double          conifer_ca,
    double          hardwood_ca,
    double          shrub_ca,
    double          uniform_0_1,
    double          *cr_growth,
    double          *coeffs_ptr )
{
    double      b0;              /* temporary coefficients                   */
    double      b1;
    double      b2;
    double      b3;
    double      b4; 
    double      b5;
    double      b6;
    double      b7;
    double      b8;
    double      b9;
    double      hcb_growth;      /* predicted hcb grow if crown receeds      */
    double      crown_length;    /* starting crown length                    */
    double      temp_exponent;   /* a temporary variable for logistic func   */ 
    double      prob_hcb;        /* probability hcb does not receed          */

/*   set the temporary variables     */

    hcb_growth  = 0.0;
    crown_length= crown_ratio * total_height;

/*    set the coefficients           */

    b0          = coeffs_ptr[0];
    b1          = coeffs_ptr[1];
    b2          = coeffs_ptr[2];
    b3          = coeffs_ptr[3];
    b4          = coeffs_ptr[4];
    b5          = coeffs_ptr[5];
    b6          = coeffs_ptr[6];
    b7          = coeffs_ptr[7];
    b8          = coeffs_ptr[8];
    b9          = coeffs_ptr[9];


    if( coeffs_ptr == NULL )
    {
        *cr_growth    = 0.0;
        *return_code  = INVALID_COEFF;
        return;
    }

    if( !hcb_growth_on )  /*  on - off switch because hcb change is biannual */
                          /*  off = no change in hcb, i.e. change in hcb=0   */
    {
        *cr_growth    = ( crown_length + height_growth ) / 
                        ( total_height + height_growth ) -
                        ( crown_length / total_height  );
        *return_code  = CONIFERS_SUCCESS;
        return;
    }


/* (1)   first predict change in height to crown base              */    

    hcb_growth  = b0 
                + b1 * crown_length 
                + b2 * 0.00001 * conifer_ca
                + b3 * 0.00001 * hardwood_ca
                + b4 * 0.00001 * shrub_ca;

    if( hcb_growth < 0.0 )          /*  this should NEVER happen!  */
    {
        hcb_growth   = 0.0;
        *cr_growth   = 0.0;
        *return_code = CONIFERS_ERROR;
	return;
    }
    if( hcb_growth > ( crown_length + height_growth ))
    {
    /*  MOD039   added the 0.1 fudge factor to keep cr non-zero & positive */
        hcb_growth   = crown_length + height_growth - 0.1; /*  may never happen! */
    }

/* (2)   next calculate the probability of change in crown base    */

    temp_exponent   = exp( b5 
                         + b6 * crown_ratio  
                         + b7 * conifer_ca  / 43560.0 
                         + b8 * hardwood_ca / 43560.0
                         + b9 * shrub_ca    / 43560.0);

/* (2a)  logistic function                                         */

    prob_hcb        = temp_exponent / ( 1.0 + temp_exponent );

/* (3)   if the uniform r.v. is greater then apply change in cr    */

    if( uniform_0_1 > prob_hcb )  
    {
        *cr_growth   = ( crown_length + height_growth   - hcb_growth ) / 
                       ( total_height + height_growth ) -
                       ( crown_length / total_height  );

        *return_code = CONIFERS_SUCCESS;
    }
    else
    {
        *cr_growth   = ( crown_length + height_growth ) / 
                       ( total_height + height_growth ) -
                       ( crown_length / total_height  );
        *return_code = CONIFERS_SUCCESS;
    }

}

/********************************************************************************/
/*                  smc_calc_cw_growth          D6                              */
/********************************************************************************/
void smc_calc_cw_growth(   
	unsigned long   *return_code,
	double          total_height,
	double		    height_growth,
	double		    crown_width,
	double		    ca_conifers,
	double		    ca_hardwoods,
	double		    ca_shrubs,
	double		    catcon,
	double		    unif_0_1,
	double          expf,
	double          ba,
	double          flew_site,
	double          h40,
	double		    con_tpa,
    double          *pred_cw_growth,
	double          *pred_mortality,
    double          *coeffs_ptr,
	unsigned long   plant_type)
{

	double  temp_cover_grow;
	double  temp_cover;
	double  vc;
	double  temp_cwg;
	double  c0;
	double  c1;
	double  c2;
	double  c3;
	double  c4;
	double  c5;
	double  c6;
	double  c7;
	double  c8;
	double  c9;
	double  c10;

/* this function currently does not have a random damage component but that is */
/*  what the extra parameters (c7, c8) are for and the uniform_0_1 variable    */

	*return_code	= CONIFERS_SUCCESS;

    if( coeffs_ptr == NULL )
    {
        *pred_cw_growth    = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }
    if(ca_conifers < 0.0 || ca_hardwoods < 0.0 || ca_shrubs < 0.0)
    {
 	    *pred_cw_growth	=0.0;
	    *return_code	= INVALID_INPUT_VAL;
	    return;
    }

    if(crown_width <= 0.0 )
    {
	    *pred_cw_growth	=0.0;
	    *return_code	= INVALID_INPUT_VAL;
	    return;
    }

    if( height_growth < 0.0 ) 
    {
        *pred_cw_growth = 0.0;
	    *return_code	= CONIFERS_SUCCESS;
	    return;
    }


    temp_cwg    = 0.0;
	ca_conifers = (ca_conifers/SQ_FT_PER_ACRE);
	ca_hardwoods= (ca_hardwoods/SQ_FT_PER_ACRE);
	ca_shrubs   = (ca_shrubs/SQ_FT_PER_ACRE);
	catcon      = catcon/SQ_FT_PER_ACRE;
    
	c0  = coeffs_ptr[0];
	c1  = coeffs_ptr[1];
	c2  = coeffs_ptr[2];
	c3	= coeffs_ptr[3];
	c4	= coeffs_ptr[4];
	c5	= coeffs_ptr[5];
	c6	= coeffs_ptr[6];
	c7	= coeffs_ptr[7];
	c8	= coeffs_ptr[8];
	c9  = coeffs_ptr[9];
	c10 = coeffs_ptr[10];

    if(plant_type == CONIFER || plant_type == HARDWOOD)  /* then use the swo version for trees */
    {
	  temp_cwg = pow(height_growth,c1)
		    *(c0 
			+ c2*sqrt(crown_width) 
			+ c3*ca_conifers*ca_conifers 
			+ c4*ca_hardwoods*ca_hardwoods 
			+ c5*ca_shrubs*ca_shrubs 
			+ c6*log(crown_width));

	  if( temp_cwg < 0.0 ) /* if pred growth is negative */
	  {
		 *pred_cw_growth    = 0.10;  /* Growth will be set to 0.10 */
		 *return_code        = CONIFERS_ERROR;
		 return;
	  }
	  *pred_cw_growth     = temp_cwg ;
	  *return_code        = CONIFERS_SUCCESS;
	  return;
    }
    else if (plant_type == SHRUB) /* if it is a shrub, do this*/
    {
	  vc = (MY_PI/(4.0*SQ_FT_PER_ACRE))*(crown_width*crown_width)*expf;
	  /* vc should be equal to ca_shrubs if one plant per plot */

	  temp_cover      = 100.0 * vc * exp(c0    + c1*ca_conifers + c2*ca_shrubs ); 
	  temp_cover_grow = temp_cover - 100.0 * vc;
	  
      if(expf > 0.0)
      {
	    temp_cwg = sqrt((temp_cover/100.0)*(SQ_FT_PER_ACRE*4.0/(MY_PI*expf))) - crown_width;
      }
	  if(temp_cwg+crown_width>8.0) /*keep shrubs from getting over 8 feet in diameter */
	  {
	    *pred_cw_growth = 0.0;
		*return_code     = CONIFERS_SUCCESS;
  	    return;
	  }

	  if(temp_cover_grow < 0.0) /* if cover is decreasing make a change in mortality */
	  {
		*pred_mortality =expf- (temp_cover/100.0)*(SQ_FT_PER_ACRE*4.0/(MY_PI*crown_width*crown_width));
        *pred_cw_growth =0.0;
		*return_code = CONIFERS_SUCCESS;
        return;
	  }

	  *pred_cw_growth  = temp_cwg;
	  *return_code     = CONIFERS_SUCCESS;
	  return;
    }
	else  /* then you don't have a valid plant type for smc version*/
	{
      *pred_cw_growth = 0.0;
	  *return_code=INVALID_PLANT_TYPE;
	  return;
	}
}


/********************************************************************************/
/*                  smc_calc_d6_growth          D3                              */
/********************************************************************************/
/*  Description :   calc_d6_growth for smc, basically a simple proportion       */    
/*                   of dbh or a function of ht depending on size               */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   Aug. 29, 2007                                               */
/*  Returns     :   void                                                        */
/*  Comments    :   refit by martin                                             */
/*                  diameter growth for tree species...                         */
/*  Arguments   :                                                               */
/*  unsigned long *return_code - return code for calling function to check      */
/*  double  height_growth    - predicted height increment (annual)              */
/*  double  total_height     - total height of the plant in feet                */
/*  double  dbh_growth       - predicted diameter growth in inches              */
/*  double  current_d6       - initial d6 (inches)                              */
/*  double  *pred_d6_growth  - predicted dbh growth (inches)                    */
/*  vector<double> *coeffs_ptr - pointer to a vector of doubles that            */
/*                              contain the coefficients for the                */
/*                              functional species code                         */
/********************************************************************************/
/*  Formula : See above                                                         */ 
/*  Source  :                                                                   */
/*  Coeffs  : D1, which is really the static d6-h, and d6-h-dbh coeffs borrowed */
/********************************************************************************/
void smc_calc_d6_growth(
    unsigned long   *return_code,
    double  height_growth,
    double  total_height,
    double  dbh,
    double  dbh_growth,
    double  d6,
    double  *pred_d6_growth,
    double  *c_v_d6_ht,     /* coeffs vector d6 as a func of ht */
	double  *c_v_d6_ht_dbh, /* coeffs vector d6 as a func of ht & dbh */
    unsigned long plant_type)
{

    double  dg;
    
    double  b0;
    double  b1;
    double  d0;
    double  d1;
    double  d2;


    /* initialize variables */
    *return_code=CONIFERS_SUCCESS;
    dg  = 0.0;
    
    b0  = 0.0;
    b1  = 0.0;
    d0  = 0.0;
    d1  = 0.0;
    d2  = 0.0;

    if( c_v_d6_ht == NULL || c_v_d6_ht_dbh == NULL )
    {
        *pred_d6_growth = 0.0;
	    *return_code=CONIFERS_ERROR;
        return;
    }
/******************************if predicted height growth is negative set diameter growth to zero************************/
    if (height_growth <= 0.0)
    {
        *pred_d6_growth = 0.0;
  	    return;
    }

    b0  = c_v_d6_ht[0]; /* these two are from the d6-ht function */
    b1  = c_v_d6_ht[1];
    d0  = c_v_d6_ht_dbh[0]; /* these three are from the d6-h-dbh fn  */
    d1  = c_v_d6_ht_dbh[1];
    d2  = c_v_d6_ht_dbh[2];

    if(total_height > 4.5 && dbh_growth >0.0 && plant_type == CONIFER)
    {
        //dg=exp(b0) * ( pow((total_height+height_growth-0.5), b1) - pow((total_height-0.5),b1) );
	    dg = (dbh + dbh_growth)/( d0 - d1 * exp ( -d2 * ( total_height+height_growth-4.5 ) ) ) - d6;
    }
    else
    {
	    if(height_growth > 0.0)
		{
	        dg=exp(b0) * ( pow((total_height+height_growth-0.5), b1) - pow((total_height-0.5),b1) );
		}
	    else
		{
	        dg=0.0;	
		}
    }
    
    if( dg < 0.0 )
    {
        dg = 0.0;
    }
    *pred_d6_growth = dg ;
}



/********************************************************************************/
/*                  smc_calc_d6_from_total_height       S4                      */
/********************************************************************************/
/*  Description :   Calculate a missing basal diameter from                     */
/*                  the total height.                                           */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   August 28, 2007                                             */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  *pred_d6                    -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula :   ln ( D6 ) = b0 * b1 * ln (H - 0.5)                              */
/*       Where   H=total height in feet                                         */
/*               D6=basal diameter at 6 inches above ground, in inches          */
/*                                                                              */
/*  Source: Ritchie, Static Height Diameter Equations S3                        */
/*  Coefficients: MS in coefficients file                                       */
/*  coefficients last fit in October 1999                                       */
/********************************************************************************/
void smc_calc_d6_from_total_height(       
    unsigned long   *return_code,
    double          total_height, 
    double          *pred_d6,
    double          *coeffs_ptr    )
{
    double  b0;
    double  b1;

    /* initialize variables */
    *return_code    = CONIFERS_ERROR;
    *pred_d6        = 0.0;
    b0              = 0.0;
    b1              = 0.0;

    if( coeffs_ptr == NULL )
    {
        *pred_d6        = 0.0;
        *return_code    = INVALID_COEFF;
        return;
    }
        
    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];

    if(total_height > 0.5)
    {
        *pred_d6 = exp( b0 + b1 * log( total_height - 0.5 ) );
        *return_code = CONIFERS_SUCCESS;
    }
    else
    {
        *pred_d6=0.1;
    }

    if( *pred_d6 < 0 )
    {
        *pred_d6    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}

/********************************************************************************/
/*                  smc_calc_d6_from_ht_and_dbh          S5                     */
/********************************************************************************/
/*  Description :   Calculate a missing basal diameter from                     */
/*                  the total height and DBH.                                   */
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   October 13, 1999                                            */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  *pred_d6                    -   predicted basal diameter from the function  */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula:  D6/dbh = 1.0/(b0 -b1 * exp ( -b2 ln(H - 4.5))                     */ 
/*       Where: D6=basal diameter at 6 inches measured in inches                */
/*              dbh= breast height diameter in inches                           */
/*              H = total tree or plant height of the largest stem in feet      */
/*  Source  : Ritchie Static Height Diameter Equations,   S5                    */
/*  Coeffs  : DH                                                                */
/********************************************************************************/
void smc_calc_d6_from_ht_and_dbh(       
    unsigned long   *return_code,
    double          total_height,
    double          dbh,
    double          *pred_d6,
    double          *coeffs_ptr    )
{

    double  b0;
    double  b1;
    double  b2;

    /* perform error check for the correct number of coeffs */
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        return;
    }

    if( total_height <= 4.5 )
    {
        *return_code = INVALID_INPUT_VAL;
        return;
    }

    b0  = coeffs_ptr[0];
    b1  = coeffs_ptr[1];
    b2  = coeffs_ptr[2];

    *pred_d6 = dbh / ( b0 - b1 * exp ( -b2 * ( total_height - 4.5 ) ) );
    *return_code = CONIFERS_SUCCESS;

    if( *pred_d6 < 0.0 )
    {
        *pred_d6    = 0.0;
        *return_code        = CONIFERS_ERROR;
    }
}

/********************************************************************************/
/*                  smc_calc_crown_width     S1                                 */
/********************************************************************************/
/*  Description :   smc_calc_crown_width                                        */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   November 30, 2007                                           */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -  return code for calling function to check    */
/*  d6_area                     -  basal area in inches of largest stem at base */
/*  total_height                -  total height of the plant                    */
/*  *pred_crown_area            -  predicted crown area                         */
/*  *pred_crown_width           -  predicted crown width in square feet         */
/*  vector<double> *coeffs_ptr  -  pointer to a vector of doubles that          */
/*                                   contain the coefficients for the           */
/*                                   functional species code                    */
/*                                                                              */
/********************************************************************************/
/*  Formula : log (ca) = b0 + b1 log (d6ba) + b2* log(h)   S1                   */ 
/*           where: ca      = crown area in square feet                         */
/*                  h       = height, feet                                      */
/*                  d6ba    = basal area of largest stem at 6 inches, in^2      */
/*  Source  : Uzoh and Ritchie 1996 PSW RP 227                                  */
/*  Coeffs  : CW   may need to add b3 due to minor LOF probs with df            */
/********************************************************************************/
void smc_calc_crown_width( 
    unsigned long   *return_code,
    double          d6_area,
    double          total_height,
    double          *pred_crown_width,
    double          *pred_crown_area,
    double          *coeffs_ptr,
	unsigned long   plant_type)
{

    double  b0;
    double  b1;
    double  b2;

        /* check for valid height */
        if(total_height <= 0.0)
        {
            *return_code      = INVALID_INPUT_VAL;
            *pred_crown_width = 0.0;
            *pred_crown_area  = 0.0;
            return;
        }

        /* check for valid coefficients */
        if( coeffs_ptr == NULL )
        {
            *return_code        = INVALID_COEFF;
            *pred_crown_width   = 0.0;
            *pred_crown_area    = 0.0;
            return;
        }

        b0 = coeffs_ptr[0];
        b1 = coeffs_ptr[1];
        b2 = coeffs_ptr[2];

        if( total_height < 0.51 )
        {
            *pred_crown_width = 0.25;
            *pred_crown_area  = 0.04908739;
        }
        else
        {
            *pred_crown_area = exp( b0 + 
                                    b1 * log( d6_area * 144.0 ) + 
                                    b2 * log( total_height ) );
            if(plant_type == SHRUB)
			{
				if( *pred_crown_area > 50.1 )
				{
					*pred_crown_area = 50.1 ; /* limit crown width to 8 feet for shrubs */
				}
			}
			else
			{
				if( *pred_crown_area > 2827.0 )
				{
			        *pred_crown_area = 2827.0 ; /* limit crown width to 60 feet */
				}
			}

			*pred_crown_width = sqrt( *pred_crown_area * ONE_OVER_PI * 4.0); 
        }
        *return_code = CONIFERS_SUCCESS;

        if( *pred_crown_width < 0.0 )
        {
            *pred_crown_width   = 0.0;
            *pred_crown_area    = 0.0;
            *return_code        = CONIFERS_ERROR;
			return;
        }
}

/********************************************************************************/
/*                  smc_calc_crown_ratio                                        */
/********************************************************************************/
/*  Description :   calc_crown_width                                            */   
/*  Author      :   Martin W. Ritchie                                           */
/*  Date        :   May 9, 2008                                                 */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :                                                               */
/*  return void                                                                 */
/*  unsigned long *return_code  -   return code for calling function to check   */
/*  total_height                -   total height of the subject tree            */
/*  d6                          -   basal diameter                              */
/*  double *pred_cr             -   pointer to predicted crown ratio            */
/*  vector<double> *coeffs_ptr       -   pointer to a vector of doubles that    */
/*                                  contain the coefficients for the            */
/*                                  functional species code                     */
/********************************************************************************/
/*  Formula : CL =b0*H**b1 *exp( b3 d6/height)                                  */ 
/*  Source  : Ritchie May 2008                                                  */
/*  Coeffs  : CR                                                                */
/********************************************************************************/
void smc_calc_crown_ratio(
    unsigned long *return_code,  
    double  total_height,
    double  d6,
    double  *pred_cr,
    double  *coeffs_ptr    )
{

    double b0, b1, b3;
    double ratio;
    double newCR;           
    double crown_length;

    *return_code = CONIFERS_SUCCESS;
    if( coeffs_ptr == NULL )
    {
        *return_code = INVALID_COEFF;
        *pred_cr = 0.0;
        return;
    }

    if( total_height <= 0.0 || d6 <= 0.0)
    {
        *return_code = INVALID_INPUT_VAL;
        *pred_cr = 0.0;
        return;
    }

    ratio=d6/total_height;

    b0 = coeffs_ptr[0];
    b1 = coeffs_ptr[1];
    b3 = coeffs_ptr[3];


    crown_length = b0*pow(total_height,b1) * exp(b3*ratio);

    if (total_height <= 1.0)
    {
        crown_length = b0*total_height;
    }

    newCR  = crown_length/total_height;
    
    if( newCR < 0.0 )
    {
        newCR = 0.0;
    }
    if( newCR > 1.0 )
    {
        newCR = 1.0;
    }

    *pred_cr = newCR;

}


void smc_calc_endemic_mortality(   
    unsigned long   *return_code,
    double          expansion_factor,
    double          *pred_mortality,
    double          *coeffs_ptr )
{

    double  b0;
    double  current_mort;
    /*  initialize parameters    */
    /* *pred_mortality = 0.0;    */
    b0              = coeffs_ptr[0];
    current_mort	= *pred_mortality;

    if( coeffs_ptr == NULL )
    {
        *pred_mortality     = 0.0;
        *return_code        = INVALID_COEFF;
        return;
    }

    /*   check that the starting expansion factor and mortality rate are in bounds  */
    if( expansion_factor < 0.0 || b0 < 0.0 || b0 > 0.10 || expansion_factor < current_mort)
    {
        *pred_mortality     = 0.0;
        *return_code        = INVALID_INPUT_VAL;
        return;
    }
    
    /*   this will calculate the number of trees to die using the endemic mortality and current mort if any */
    *pred_mortality     = (expansion_factor - current_mort)* b0;
    *return_code        = CONIFERS_SUCCESS;

    if( *pred_mortality < 0.0 )
    {
        *pred_mortality     = 0.0;
        *return_code        = CONIFERS_ERROR;
    }

}
