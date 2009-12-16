/********************************************************************************/
/*                                                                              */
/*  grow.c                                                                      */
/*  functions used to predict values for the CONIFERS growth model              */
/*                                                                              */
/********************************************************************************/

/* 	$Id: grow.c 644 2009-11-24 02:12:55Z hamannj $	 */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conifers.h"

/*
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
*/

/* local function declarations */
static void project_plot( 
   unsigned long           *return_code,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           endemic_mortality,      
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	        variant,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains);


static void swo_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,      
   int                     hbc_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err );

static void smc_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,      
   int                     hbc_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err, 
   struct SUMMARY_RECORD   *before_sums,
   unsigned long            use_genetic_gains,
   unsigned long            genetics_age_cut);



/********************************************************************************/
/* project the plant list                                                       */
/* this function will project each plot for one year                            */
/* to project the entire sample for more than one year, this function           */
/* needs to be called once for each year                                        */
/********************************************************************************/
void __stdcall project_plant_list( 
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   unsigned long           n_points,
   struct PLOT_RECORD      *plots_ptr,
   double		           *x0, 		   
   unsigned long           endemic_mortality,      
   unsigned long           sdi_mortality,          
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	        variant,
   unsigned long            use_genetic_gains    )
{
   unsigned long           i;
   struct  PLOT_RECORD     *plot_ptr;
   double                  max_sdi;
   struct SUMMARY_RECORD   before_sums;
   struct SUMMARY_RECORD   after_sums;
   double		           mortality_proportion;

   double  plants_removed;
   double  ba_removed;


   /* update the plot variabels before you procede */
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   calc_plot_stats_2(   return_code, 
   n_species,
   species_ptr,
   n_coeffs,
   coeffs_ptr,
   n_plants,
   plants_ptr,
   n_points,
   plots_ptr );

   if( *return_code != CONIFERS_SUCCESS )
   {
      return;
   }

   /* update the total summaries before the plant list is projected */
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   update_total_summaries(  return_code,
   n_points,
   n_plants,
   n_species,
   species_ptr,
   n_coeffs,
   coeffs_ptr,
   plants_ptr,
   &before_sums );

   if( *return_code != CONIFERS_SUCCESS )
   {
      return;
   }


   /* for each plot, project it forward one year */
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   plot_ptr = &plots_ptr[0];
   for( i = 0; i < n_points; i++, plot_ptr++ )
   {
      project_plot( return_code,
      n_plants,
      plants_ptr,
      plot_ptr,
      n_species,
      species_ptr,
      n_coeffs,
      coeffs_ptr,
      endemic_mortality,
      hcb_growth_on,          
      use_precip_in_hg,  
      use_rand_err,
      variant,    
      &before_sums,
      use_genetic_gains);      

      if( *return_code != CONIFERS_SUCCESS )
      {
	 return;
	 //continue;
      }

   }

    
   /* now, calculate the limiting sdi mortality value for  */
   /* all the trees in the plant list                      */
    
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   /* calcuate rd */
   /* update the total summaries before the plant list is projected */
   update_total_summaries(  return_code,
   n_points,
   n_plants,
   n_species,
   species_ptr,
   n_coeffs,
   coeffs_ptr,
   plants_ptr,
   &after_sums );

   /* need_error_trap_here */
   /* if max sdi limit switch is on, then      */
   /* MOD012 */
   //Rprintf( "if( sdi_mortality ) check == 1 @ %s, %d, then ()\n", __FILE__, __LINE__ );
   if( sdi_mortality )
   {
      calc_max_sdi( return_code,
      n_species,
      species_ptr,
      n_coeffs,
      coeffs_ptr,
      n_plants,
      plants_ptr,
      n_points,
      &max_sdi  );

      /* it is choking here. */
      if( *return_code != CONIFERS_SUCCESS )
      {
	 //Rprintf( "calc_max_sdi *return_code != CONIFERS_SUCCESS @ %s, %d\n", __FILE__, __LINE__ );
	 return;
      }


      /* this might have been changed since the code was ported */
      /* to the R interface, but a check needs to happen to	*/
      /* insure that x0 isn't not reset. A check will be placed */
      /* here to verify x0 has not been set.			*/      
      if( *x0 == 0.0 ) 
      {
	 //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
  	 //Rprintf( "*x0 == 0.0 -- calling calc_hann_wang_x0\n" );
  	 //Rprintf( "after_sums.sdi %f\n", after_sums.sdi );
  	 //Rprintf( "max_sdi = %lf\n", max_sdi );

	 /* this seems to be broken */
	 calc_hann_wang_x0( return_code,
	                    after_sums.sdi,
	                    after_sums.bh_expf,
	                    max_sdi,
	                    before_sums.rel_density,
	                    after_sums.rel_density,
	                    x0);
	 
	 if( *return_code != CONIFERS_SUCCESS )
	 {
	    //Rprintf( "*return_code != CONIFERS_SUCCESS, %s, %d\n", __FILE__, __LINE__ );
	    //Rprintf( "*x0 == 0.0 -- calling calc_hann_wang_x0\n" );
	    return;
	 }
      }
	
      if( *x0 > 0.0 )
      {
	 //Rprintf( "*x0 > 0.0 -- calling calc_sdi_mortality\n" );

	 calc_sdi_mortality( return_code,
	                    after_sums.qmd,
	                    after_sums.sdi,
	                    max_sdi,
	                    *x0,
	                    after_sums.bh_expf,
	                    &mortality_proportion);
	 if( *return_code != CONIFERS_SUCCESS )
	 {
	    return;
	 }

	 /*    this next section does the sdi mortality   */
	 plot_ptr = &plots_ptr[0];
	 for( i = 0; i < n_points; i++, plot_ptr++ )
	 {
	    thin_plot(  return_code,
	                n_plants,
	                plants_ptr,
	                plot_ptr,
	                n_species,
	                species_ptr,
	                n_coeffs,
	                coeffs_ptr,
	                //NULL,
	                0,
	                DO_SDI_MORT,
	                mortality_proportion, 
	                &plants_removed, 
	                &ba_removed );
            
	    if( *return_code != CONIFERS_SUCCESS )
	    {
	       continue;/* error trap here */
	    }
	 }
      }
   } /* end of the sdi_mortality switch */

}


/********************************************************************************/
/* project_plot                                                                 */
/********************************************************************************/
/*  Description :   projects all the plants on the plot forward one year        */
/*                  into the future                                             */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 12, 1999                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   NONE                                                        */
/*  Arguments   :   unsigned long *return_code  - pointer to a return code      */
/*                  unsigned long n_plants      - total number fo plants in the */
/*                      plants pointer array                                    */
/*                  struct PLANT_RECORD *plants_ptr - array of plants in the    */
/*                      sample to be projected                                  */
/*                  struct PLOT_RECORD  *plot_ptr   - pointer to the current    */
/*                      plot that is to be grown                                */
/*                  unsigned long n_species - size of the species_ptr           */
/*                  struct SPECIES_RECORD   *species_ptr - array of             */
/*                      SPECIES_RECORD's that hold species specific information */
/*                  unsigned long   n_coeffs - sizes of the coeffs_ptr array    */
/*                  struct COEFFS_RECORD *coeffs_ptr - array of coefficients    */
/*                      that are used to project the individual plants on the   */
/*                      plot.                                                   */
/********************************************************************************/
static void project_plot( 
   unsigned long           *return_code,
   unsigned long           n_plants,
   struct PLANT_RECORD     *plants_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   unsigned long           endemic_mortality,      
   int                     hcb_growth_on,          
   unsigned long           use_precip_in_hg,       
   unsigned long           use_rand_err,
   unsigned long	       variant,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains)       
{

   unsigned long   i;
   struct  PLANT_RECORD    *plant_ptr;
   double                  si_30;
   double                  h40;
   unsigned long           genetics_age_cut;

   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   genetics_age_cut=0;
   /* for each plot, project it forward one year */
   if(variant == CONIFERS_SMC)
   {
      h40   =  before_sums->height_40;
      si_30 =  plot_ptr->site_30;
      get_age_cut(return_code,
                  h40,
                  si_30,
                  &genetics_age_cut);
   }

   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );
   plant_ptr = &plants_ptr[0];
   for( i = 0; i < n_plants; i++, plant_ptr++ )
   {

      /* if you're not on the plot, don't grow it */
      /* you should really be using the function  */
      /* that returns the index values for the    */
      /* trees on the plot to get the thing done  */
      /* just iterate through the plant list and  */
      /* pitch the plant if it's not on the plot  */
      if( plant_ptr->plot != plot_ptr->plot )
      {
	    continue;
      }

      /* for each tree on the plot....    */
      switch( variant )
      {
	    case CONIFERS_SWO:
	      swo_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err  );  
		break;

	    case CONIFERS_SMC:
	      smc_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err,
			                before_sums,
                            use_genetic_gains,
                            genetics_age_cut);   
		break;

	    default:
	      swo_project_plant(return_code,
			                n_species, 
			                species_ptr, 
			                n_coeffs, 
			                coeffs_ptr,
			                plant_ptr, 
			                plot_ptr,
			                endemic_mortality,  
			                hcb_growth_on,      
			                use_precip_in_hg,    
			                use_rand_err  );     
		break;
      }

      if( *return_code != CONIFERS_SUCCESS )
      {
	 *return_code = FAILED_PROJECT_PLANT;  
	 //Rprintf( "FAILED_PROJECT_PLANT, %s, %d\n", __FILE__, __LINE__ );
	 return;
      }

      /* update the current tree values...        */
      plant_ptr->d6           += plant_ptr->d6_growth;
      plant_ptr->dbh          += plant_ptr->dbh_growth;
      plant_ptr->tht          += plant_ptr->tht_growth;
      plant_ptr->cr           += plant_ptr->cr_growth;
      plant_ptr->crown_width  += plant_ptr->cw_growth;
      plant_ptr->expf         -= plant_ptr->expf_change;
        
      /* update the extra data for the tree       */
      plant_ptr->basal_area    = plant_ptr->dbh * plant_ptr->dbh * FC_I;
      plant_ptr->d6_area       = plant_ptr->d6 * plant_ptr->d6 * FC_I;
      plant_ptr->crown_area    = plant_ptr->crown_width * plant_ptr->crown_width
				 * MY_PI / 4.0;
      plant_ptr->pct_cover     = 100.0 * plant_ptr->expf * plant_ptr->crown_area/SQ_FT_PER_ACRE;
   }

}


/********************************************************************************/
/* swo_project_plant                                                            */
/********************************************************************************/
/*  Description :   This function compute the plant growth for one year.        */
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   September 12, 1999                                          */
/*  Returns     :   void                                                        */
/*  Comments    :   updated values are height_growth, d6_growth, dbh_growth,    */
/*                  crown_ratio_growth, crown_width_growth, and max_crown_width */
/*  Arguments   :   unsigned long *return_code  - pointer to a return code      */
/*                  unsigned long n_plants      - total number fo plants in the */
/*                      plants pointer array                                    */
/*                  struct PLANT_RECORD *plants_ptr - array of plants in the    */
/*                      sample to be projected                                  */
/*                  struct PLOT_RECORD  *plot_ptr   - pointer to the current    */
/*                      plot that is to be grown                                */
/*                  unsigned long n_species - size of the species_ptr           */
/*                  struct SPECIES_RECORD   *species_ptr - array of             */
/*                      SPECIES_RECORD's that hold species specific information */
/*                  unsigned long   n_coeffs - sizes of the coeffs_ptr array    */
/*                  struct COEFFS_RECORD *coeffs_ptr - array of coefficients    */
/*                      that are used to project the individual plants on the   */
/*                      plot.                                                   */
/********************************************************************************/
static void swo_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,  
   int                     hcb_growth_on,      
   unsigned long           use_precip_in_hg,   
   unsigned long           use_rand_err )
{

   struct COEFFS_RECORD    *c_ptr;
   double                  awi;

   /* tree level stats */
   double  bat_total;
   double  bat_c;
   double  bat_h;
   double  bat_s;
   double  bat_c_h;
   double  cat_c;
   double  cat_h;
   double  cat_s;
   double  new_d6_area;

/*    unsigned long htidx = 0; */

   double  normal;
   double  browse_random_unif_0_1;
   double  top_dam_random_unif_0_1;

   double  bait[PLANT_TYPES];
   double  cait[PLANT_TYPES];
    
//   Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   /* get the supporting structures for the plant */
   c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

   /* don't go any further if you don't    */
   /* have a functional species code       */
   if( c_ptr == NULL )
   {
      *return_code = CONIFERS_ERROR;
      return;
   }

//   Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   /* get a uniform deviate for the browse */
   /* and one for the top damage           */
   normal                  = (double)gauss_dev();  
   browse_random_unif_0_1  = uniform_0_1();
   top_dam_random_unif_0_1 = uniform_0_1();

   /* calc the height growth...        */
   /* calc diam growth...              */
   /* calc crown recession...          */
   /* calc mortality...                */
   awi         =   plot_ptr->water_capacity * log( plot_ptr->mean_annual_precip );


    /* this function takes a vector [PLANT_TYPES][TALLER_PLANT_SIZE] */
    //get_taller_attribs( plant_ptr->tht, 
    //                    plot_ptr, 
    //                    bait, 
    //                    cait );
   //Rprintf( "%s, %d\n", __FILE__, __LINE__ );

    get_in_taller_attribs( plant_ptr, 
                            plot_ptr, 
                            bait, 
                            cait );

   bat_c       =   bait[CONIFER];
   bat_h       =   bait[HARDWOOD];
   bat_s       =   bait[SHRUB];

   cat_c       =   cait[CONIFER];
   cat_h       =   cait[HARDWOOD];
   cat_s       =   cait[SHRUB];

   bat_c_h     =   bat_c + bat_h;
   bat_total   =   bat_c + bat_h + bat_s;

   plant_ptr->tht_growth = 0.0;
   plant_ptr->d6_growth  = 0.0;
   plant_ptr->cw_growth  = 0.0;
   plant_ptr->cr_growth  = 0.0;
   plant_ptr->dbh_growth = 0.0;
   plant_ptr->expf_change= 0.0;

//Rprintf( "%s, %d\n", __FILE__, __LINE__ );
	if( is_tree( c_ptr ) || is_shrub( c_ptr))
	{
		calc_height_growth( return_code,
			                plant_ptr->tht, 
			                plant_ptr->cr,
			                plant_ptr->d6_area,
			                plot_ptr->mean_annual_precip, 
			                plot_ptr->water_capacity,
			                plot_ptr->slope,
			                plot_ptr->aspect,
						
			                plot_ptr->ca_c,
			                cat_c,

			                plot_ptr->ca_h,         
			                plot_ptr->ca_s,         
            
			                cat_h,
			                cat_s,

			                plant_ptr->d6,
			                plot_ptr->elevation,
			                normal,                   
			                browse_random_unif_0_1,   
			                top_dam_random_unif_0_1,  
			                use_rand_err,         

			                species_ptr[plant_ptr->sp_idx].browse_damage,
			                species_ptr[plant_ptr->sp_idx].mechanical_damage,

			                use_precip_in_hg,     
			                &plant_ptr->tht_growth,
			                c_ptr->ht_growth,
			                c_ptr->type); 
        if( *return_code != CONIFERS_SUCCESS )
        {

	        return;
        }

	}


	if( is_tree( c_ptr ) || is_shrub( c_ptr))
	{
		calc_d6_growth(	return_code,
						plant_ptr->tht_growth,
						plant_ptr->crown_width,
						plant_ptr->tht,
						plant_ptr->d6,
                        cat_s,
						cat_c,
						cat_h,
						plot_ptr->ca_s,
						plot_ptr->ca_c,
						plot_ptr->ca_h,
						plot_ptr->water_capacity,
						&plant_ptr->d6_growth,
						c_ptr->d6_growth,
						c_ptr->type);
         if( *return_code != CONIFERS_SUCCESS )
         {
	         return;
         }
	}
   /* predict the new dbh from     */
   /* the new d6 and total height  */
   if( is_tree( c_ptr ) )
   {
      calc_dbh_growth(  return_code,
						plant_ptr->tht,
						plant_ptr->tht_growth,
						plant_ptr->cr,
						plant_ptr->dbh,
						plant_ptr->d6,
						plant_ptr->d6_growth,
						&plant_ptr->dbh_growth,
						c_ptr->dbh_growth );
      if( *return_code != CONIFERS_SUCCESS )
      {
		return;
      }
            
      /* calculate the new crown ratio and    */
      /* calculate the difference             */
      calc_cr_growth(   return_code,
			            hcb_growth_on,                
			            plant_ptr->tht,
			            plant_ptr->tht_growth,
			            plant_ptr->cr,
			            plot_ptr->ca_c,
			            plot_ptr->ca_h,
			            plot_ptr->ca_s,
			            uniform_0_1(),
			            &plant_ptr->cr_growth,
			            c_ptr->cr_growth);

      if( *return_code != CONIFERS_SUCCESS )
      {
	    return;
      }

      calc_max_crown_width( return_code,
			                plant_ptr->dbh + plant_ptr->dbh_growth,
			                plant_ptr->tht + plant_ptr->tht_growth,
			                &plant_ptr->max_crown_width,
			                c_ptr->max_crown_width);

      if( *return_code != CONIFERS_SUCCESS )
      {
	    return;
      }



   }

   /* calculate the new crown area for     */
   /* the plant record...  need to calc    */
   /* a temp new d6 area, total height,    */
   /* and crown width                      */
   new_d6_area = plant_ptr->d6 * plant_ptr->d6 * FC_I;
   
   if( is_tree( c_ptr ) || is_shrub( c_ptr))
   {
		calc_cw_growth( return_code,
				        plant_ptr->tht,
				        plant_ptr->tht_growth,
				        plant_ptr->crown_width,    
	                    plot_ptr->ca_c,
				        plot_ptr->ca_h,
					    plot_ptr->ca_s,
					    cat_c,
				        uniform_0_1(),
				        &plant_ptr->cw_growth,
					    c_ptr->cw_growth,
                        c_ptr->type);
        if( *return_code != CONIFERS_SUCCESS )
        {
	       return;
        }
   }

//Rprintf( "before endemic_mortality == 1 check %s, %d\n", __FILE__, __LINE__ );
   if( endemic_mortality == 1 )    
   {

//Rprintf( "before swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );


      swo_calc_endemic_mortality(	return_code,
      plant_ptr->expf,
      &plant_ptr->expf_change,
      &species_ptr[plant_ptr->sp_idx].endemic_mortality);

//Rprintf( "after swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );
      
      if( *return_code != CONIFERS_SUCCESS )
      {
//Rprintf( "error in swo_calc_endemic_mortality %s, %d\n", __FILE__, __LINE__ );
	 return;
      }
   }
//Rprintf( "%s, %d\n", __FILE__, __LINE__ );

   *return_code = CONIFERS_SUCCESS;

}


/********************************************************************************/
/* smc_project_plant                                                            */
/********************************************************************************/
/*  Description :   computes plant growth components for SMC variant		*/
/*  Author      :   Jeff D. Hamann                                              */
/*  Date        :   August 28, 2007						*/
/*  Returns     :   void                                                        */
/*  Comments    :   updated values are height_growth, d6_growth, dbh_growth,    */
/*                  crown_ratio_growth, crown_width_growth, and max_crown_width */
/*  Arguments   :   unsigned long *return_code  - pointer to a return code      */
/*                  unsigned long n_plants      - total number fo plants in the */
/*                      plants pointer array                                    */
/*                  struct PLANT_RECORD *plants_ptr - array of plants in the    */
/*                      sample to be projected                                  */
/*                  struct PLOT_RECORD  *plot_ptr   - pointer to the current    */
/*                      plot that is to be grown                                */
/*                  unsigned long n_species - size of the species_ptr           */
/*                  struct SPECIES_RECORD   *species_ptr - array of             */
/*                      SPECIES_RECORD's that hold species specific information */
/*                  unsigned long   n_coeffs - sizes of the coeffs_ptr array    */
/*                  struct COEFFS_RECORD *coeffs_ptr - array of coefficients    */
/*                      that are used to project the individual plants on the   */
/*                      plot.                                                   */
/********************************************************************************/
static void smc_project_plant(  
   unsigned long           *return_code,
   unsigned long           n_species,
   struct SPECIES_RECORD   *species_ptr,
   unsigned long           n_coeffs,
   struct COEFFS_RECORD    *coeffs_ptr,
   struct PLANT_RECORD     *plant_ptr,
   struct PLOT_RECORD      *plot_ptr,
   unsigned long           endemic_mortality,  
   int                     hcb_growth_on,      
   unsigned long           use_precip_in_hg,   
   unsigned long           use_rand_err,
   struct SUMMARY_RECORD   *before_sums,
   unsigned long           use_genetic_gains,
   unsigned long           genetics_age_cut)
{

   struct COEFFS_RECORD    *c_ptr;

   /* tree level stats */
   double  bat_total;
   double  bat_c;
   double  bat_h;
   double  bat_s;
   double  bat_c_h;
   double  cat_c;
   double  new_d6_area;

   unsigned long htidx = 0;

   double  normal;
   double  browse_random_unif_0_1;
   double  top_dam_random_unif_0_1;
   double  basal_area;
   double  h40;
   double  tpa_con_stand;
   double  current_height;
   double  new_height;

   double  bait[PLANT_TYPES];
   double  cait[PLANT_TYPES];

    
   /* get the supporting structures for the plant */
   c_ptr = &coeffs_ptr[species_ptr[plant_ptr->sp_idx].fsp_idx];

   /* don't go any further if you don't    */
   /* have a functional species code       */
   if( c_ptr == NULL )
   {
      *return_code = CONIFERS_ERROR;
      return;
   }

   /* get a uniform deviate for the browse */
   /* and one for the top damage           */
   normal                  = (double)gauss_dev();
   browse_random_unif_0_1  = uniform_0_1();
   top_dam_random_unif_0_1 = uniform_0_1();


   tpa_con_stand  = before_sums->con_tpa;
   h40		      = before_sums->height_40;
   basal_area     = before_sums->basal_area;

    /* this function takes a vector [PLANT_TYPES][TALLER_PLANT_SIZE] */
    //get_taller_attribs( plant_ptr->tht, 
    //                    plot_ptr, 
    //                    bait, 
    //                    cait );

    get_in_taller_attribs(  plant_ptr,
                            plot_ptr,
                            bait,
                            cait );

   bat_c       =   bait[CONIFER];
   bat_h       =   bait[HARDWOOD];
   bat_s       =   bait[SHRUB];

   /* add for new model */
   //cat_c       =   plot_ptr->cait[CONIFER][htidx];
   cat_c       =   cait[CONIFER];


   bat_c_h     =   bat_c + bat_h;
   bat_total   =   bat_c + bat_h + bat_s;

   plant_ptr->tht_growth = 0.0;
   plant_ptr->d6_growth  = 0.0;
   plant_ptr->cw_growth  = 0.0;
   plant_ptr->cr_growth  = 0.0;    /* these initializations were added august 2008 by mwr */
   plant_ptr->dbh_growth = 0.0;
   plant_ptr->expf_change= 0.0;

   if( is_non_stocked(c_ptr) || plant_ptr->expf <= 0.0)
   {
       *return_code = CONIFERS_SUCCESS;
       return;
   }

   current_height = plant_ptr->tht; 
   new_height = 0.0;

   if( is_tree( c_ptr ) || is_shrub( c_ptr))
   {
    	smc_calc_height_growth( return_code,
                                plant_ptr->tht,
                                plant_ptr->d6,
                                plot_ptr->site_30,
	                            tpa_con_stand,
	                            h40,
	                            plot_ptr->ca_s,
		                        plot_ptr->cait[CONIFER][htidx],
                                plot_ptr->cait[HARDWOOD][htidx],
                                plot_ptr->cait[SHRUB][htidx],
	                            normal,                   
	                            browse_random_unif_0_1,   
	                            top_dam_random_unif_0_1,  
	                            use_rand_err,         
	                            species_ptr[plant_ptr->sp_idx].browse_damage,
	                            species_ptr[plant_ptr->sp_idx].mechanical_damage,
                                
								species_ptr[plant_ptr->sp_idx].genetic_worth_h,
								species_ptr[plant_ptr->sp_idx].genetic_worth_d,

                                use_genetic_gains,
                                &plant_ptr->tht_growth,
                                c_ptr->ht_growth,
		                        c_ptr->type,
                                genetics_age_cut) ;   
        if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
		     
   }

   
   if( is_tree( c_ptr ) )
   {
        if(plant_ptr->tht > 4.5 && plant_ptr->dbh <= 0.0)
        {
            if( plant_ptr->d6 <= 0.0 )
            {
                calc_d6_from_total_height(  return_code,
                                            plant_ptr->tht, 
                                            &plant_ptr->d6,
                                            c_ptr->d6_ht);
                if ( *return_code != CONIFERS_SUCCESS)
                {
                    return;
                }
            }
            calc_dbh_from_height_and_d6(return_code,
                                        plant_ptr->d6,
                                        plant_ptr->tht,
                                        &plant_ptr->dbh,
                                        c_ptr->d6_ht_dbh );
            if( *return_code != CONIFERS_SUCCESS )
            {
                return;
            }
        }

	/* plot_ptr->site_30 == 0.0, this will still grow! */
        smc_calc_dbh_growth(return_code,
			    plant_ptr->tht,
			    h40,
			    plot_ptr->site_30,
			    plant_ptr->dbh,
			    basal_area,
			    tpa_con_stand,
			    
			    species_ptr[plant_ptr->sp_idx].genetic_worth_h,
			    species_ptr[plant_ptr->sp_idx].genetic_worth_d,
			    
			    genetics_age_cut,
                            use_genetic_gains,
			    &plant_ptr->dbh_growth,
			    c_ptr->dbh_growth,
                            c_ptr->type);
	
	      
        if( *return_code != CONIFERS_SUCCESS )
        {
	        return;
        }
            
      /* calculate the new crown ratio and    */
      /* calculate the difference             */
                      
        smc_calc_cr_growth( return_code,
			                hcb_growth_on,                
			                plant_ptr->tht,
			                plant_ptr->tht_growth,
			                plant_ptr->cr,
			                plot_ptr->ca_c,
			                plot_ptr->ca_h,
			                plot_ptr->ca_s,
			                uniform_0_1(),
			                &plant_ptr->cr_growth,
			                c_ptr->cr_growth);

        if( *return_code != CONIFERS_SUCCESS )
        {
    	   return;
        }

   }

   new_d6_area = plant_ptr->d6 * plant_ptr->d6 * FC_I;

   if( is_tree( c_ptr ))   /* then call for d6 growth */
   {
     if( plant_ptr->d6 <= 0.0 )
     {
         calc_d6_from_total_height( return_code,
                                    plant_ptr->tht, 
                                    &plant_ptr->d6,
                                    c_ptr->d6_ht);
         if ( *return_code != CONIFERS_SUCCESS)
         {
             return;
         }
     }
     smc_calc_d6_growth(return_code,
			            plant_ptr->tht_growth,
			            plant_ptr->tht,
                        plant_ptr->dbh,
			            plant_ptr->dbh_growth,
			            plant_ptr->d6,
        	            &plant_ptr->d6_growth,
			            c_ptr->d6_ht,
			            c_ptr->d6_ht_dbh,
			            c_ptr->type);
     if ( *return_code != CONIFERS_SUCCESS)
     {
          return;
     }

   } 
   if( is_shrub( c_ptr)) /* set d6 as a static function */
   {
       new_height=current_height;
       if(plant_ptr->tht_growth>0.0)
       {
           new_height=current_height+plant_ptr->tht_growth;
           calc_d6_from_total_height(return_code,
                                     new_height, 
                                     &plant_ptr->d6,
                                     c_ptr->d6_ht);
           /* need_error_trap_here */
       }
       plant_ptr->d6_growth  = 0.0;
   }


   if( is_tree( c_ptr ) || is_shrub( c_ptr))
   {
	    smc_calc_cw_growth( return_code,
			                plant_ptr->tht,
			                plant_ptr->tht_growth,
			                plant_ptr->crown_width,
			                plot_ptr->ca_c,
			                plot_ptr->ca_h,
			                plot_ptr->ca_s,
			                cat_c,
			                uniform_0_1(),
			                plant_ptr->expf,
			                basal_area,
			                plot_ptr->site_30,
			                h40,
			                tpa_con_stand,
			                &plant_ptr->cw_growth,
			                &plant_ptr->expf_change,
			                c_ptr->cw_growth,
			                c_ptr->type);
        if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
   }

   if( endemic_mortality == 1 )    
   {
        smc_calc_endemic_mortality( return_code,
			                        plant_ptr->expf,
			                        &plant_ptr->expf_change,
 			                        &species_ptr[plant_ptr->sp_idx].endemic_mortality );
        if ( *return_code != CONIFERS_SUCCESS)
        {
            return;
        }
   }

   *return_code = CONIFERS_SUCCESS;

}


