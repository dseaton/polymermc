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
	int STARTED;  //Shows whether or not the restart has been applied
	FILE *ofp_run;
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
	
	//if JA greater than zero, code is currently not ready
	if(JA > 0.0)
	{
		fprintf(stderr,"You have to adjust the energy calculations if JA > 0.0\n");
		exit(1);
	};

	//Initialize the Wang-Landau sampling parameters
	initWL();
	nummoves=0;  //used in wlhybrid() to keep track of total number of moves

	ProductionRun = atoi(argv[3]);
	//starts the production run if it is turned on in the input file
	if( ProductionRun == 1)
	{
		production_run();
		exit(1);
	}
	else if( (ProductionRun != 0) && (ProductionRun != 1) )
	{
		fprintf(stderr,"ProductionRun is not equal to 1 or 0! Is set to: %d \n\n",ProductionRun);
		exit(1);
	};
	
	
	ofp_run=fopen("run.dat","a");
	fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));

	//Initialization
	numf=1;  //Keeps track of the number of interations  
	STARTED = 0;
    
	//Main Wang-Landau sampling loop
	for(lnwlf=ModFactorInit;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
    {
		IterSweeps=0;	//Keeps track of the number of sweeps per iteration
		resetWL();		//Resets the WL histgram      
		reset_XYZ();	//Resets monomer positions by subtracting away the center of mass
		tmp_flat=0.0;	//Resets the flatness for the loop below
	  	  
		//Decides if the code should be restarted or start from scratch
		if( (atoi(argv[2]) == 1) && (STARTED == 0) )
		{
			read_restart();  //Loads in the previously simulated data from the last checkpoint	  
			STARTED = 1;
			tmp_flat=flatWL();  //Checks to see if the previous restart file already had a flat histogram
			
			if( (lnwlf <= ModFactorFinal) && (tmp_flat >= Flatness) )
			{
				fprintf(stderr,"Run has already reached the final modification factor!\n\n");
				exit(1);
			};
		};
  
		//Writes out the restart file at the beginning of each iteration
		write_restart();
		//fprintf(stderr,"%g\t%d\t%f\t%f\n",lnwlf,numf,currEtot*invN,currnonbondE*invN);
		//This while loop is executed until the histogram is sampled
		while(  tmp_flat <= Flatness  )
		{
			//Performs a number () of WL sweeps
			sweepWL( D1BINS );
			IterSweeps+=1;		
			TotalSweeps+=1;
			tmp_flat=flatWL();
	 		
			if( (IterSweeps % 10000 ) == 0)
			{
				write_DOS_H();
				write_restart();
			};
		};
		
		write_DOS_H();
		write_restart();
		//thermoqs();
		//reportcounters(stderr);
		//resetcounters();
	
		fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",lnwlf,TotalSweeps*D1BINS,IterSweeps*D1BINS,tmp_flat,numbelow_flat);
		fflush(ofp_run);
		numf=numf+1;  //Used in write_DOS_H() to lable each DOS
	};
	fclose(ofp_run); 			

	//Prints out the final normalized density of states, thermodynamics, and order parameters
	//write_DOS_H();
	//write_restart();
	write_normDOS();
	thermoqs();

}


