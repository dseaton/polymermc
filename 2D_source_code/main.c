//Semi-Flexible Polymer

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "InputOutput.h"
#include "rand.h"
#include "model.h"
#include "metrop.h"
#include "wanglandau.h"

main(int argc, char *argv[])
{
 
	int i,j,m;
	int tmpbins; 
	int Change_EPSILON;
	FILE *ofp_run,*ofp_time;
	char s[512];
	double tmp_flat;
	
	//Error message if the number of arguments is incorrect
	if (argc != 5) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);

	//Get Random Number Seed
	shelltimeseed(atoi(argv[4]));
	//gettimeseed();
  
	//Initializes the System
	initialize();
  
	//Metropolis Sampling 
	if(MetropolisSampling == 1)
	{
		//Tscan(MTi,MTf,MdT,MDROP,MDROP,MSAMPS,MSEP);
		Seq(MTf,MDROP,MSAMPS,MSEP);
			
		//write_mol2(0);
		exit(1);	
	}
	
	//Initialize the Wang-Landau sampling parameters
	initWL();
	tmpbins = sampledBins();
	//fprintf(stderr,"MAIN: %d\n",tmpbins);
	
	ProductionRun = atoi(argv[3]);
	//starts the production run if it is turned on in the input file
	if( ProductionRun == 1)
	{
		production_run();
		exit(11);
	}
	else if( (ProductionRun != 0) && (ProductionRun != 1) )
	{
		fprintf(stderr,"ProductionRun is not equal to 1 or 0! Is set to: %d \n\n",ProductionRun);
		exit(11);
	};
	
	ofp_run=fopen("run.dat","a");
	fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));

	ofp_time=fopen("timeseries.dat","a");
	fprintf(ofp_time,"#N = %d : JA, Energy/N, Num unsampled bins, lnwlf\n",N);
	
	//Initialization
	numf=1;  //Keeps track of the number of interations  
	
	//for(j=0;j<D1BINS;j++)
	//	fprintf(stderr,"%d\t%g\t%g\n",j,j/invdWLD1+WLD1min,JA);
	
	//exit(100);
	
	//Main Wang-Landau sampling loop
	for(lnwlf=ModFactorInit;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
	{
		IterSweeps=0;	//Keeps track of the number of sweeps per iteration
		resetWL();		//Resets the WL histgram      
		reset_XYZ();	//Resets monomer positions by subtracting away the center of mass
		tmp_flat=0.0;	//Resets the flatness for the loop below
	  	Change_EPSILON=1;
		tmpbins = sampledBins();
		//if(tmpbins == 0)
		//	tmpbins = 10000000;//D1BINS*D2BINS;
		if(lnwlf == ModFactorInit)
			tmpbins=10000000;
			
		//fprintf(stderr,"MAIN: %d\n",tmpbins);
  
		while( tmpbins > NUM_UNSAMPLED_BINS)//!= 0 )
		{
			//fprintf(stderr,"MAIN: %d\n",tmpbins);
			
			//Performs a number () of WL sweeps
			sweepWL( D1BINS*D2BINS );
			IterSweeps+=1;		
			TotalSweeps+=1;
			//tmp_flat=flatWL();
			tmpbins = sampledBins();
	 		
			fprintf(ofp_time,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",TotalSweeps,JA,currEtot*invN,pairnum(),Rgyr2(),aveBAngle(),icos_core_count(),tmpbins,lnwlf);
			fflush(ofp_time);
			
			//Adjusts the mask to include a certain % of the sampled states, e.g. 0.999 after the 4 iteration, if 1.0, always adds new bins
			resize_mask(PERCENT_DOS_CUT,6);
			if(numf==6)
			{
				NUM_UNSAMPLED_BINS=0.1*D2BINS;
			}
			
			//Reduces shift in energy dependent moves, and stops adjusting after a certain iteration
			if( (numf%4)==0 && Change_EPSILON==1 && numf <= 12)
			{	
				Change_EPSILON=0;
				EPSILON = EPSILON*0.5;
				write_EMIN();
			}	
			else if(numf>12)	//Turns off the shifting in each E dependent move
				EPSILON==0.0;
				
			//mindiff=1.0;
			//for(i=0;i<D1BINS;i++)	
				//if(mindiff > diffvalues[i] && mask[i]==1)
					//mindiff=diffvalues[i];

			if( (IterSweeps % 10000 ) == 0)
			{
				write_DOS_H();
				write_diff();
			};
		
		};
		
		write_DOS_H();
		//write_diff();
		//write_EMIN();
		
		//thermoqs();
	
		fprintf(ofp_run,"%g\t%d\t%d\t%g\t%d\n",lnwlf,TotalSweeps,IterSweeps,tmp_flat,tmpbins);
		fflush(ofp_run);
		numf=numf+1;  //Used in write_DOS_H() to lable each DOS
	};
	fclose(ofp_run);
	fclose(ofp_time);

	//Prints out the final normalized density of states, thermodynamics, and order parameters
	//write_DOS_H();
	write_diff();
	thermoqs();

}


