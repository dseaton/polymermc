#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rand.h"
#include "model.h"
#include "metrop.h"
#include "wanglandau.h"

void production_run(void)
{	
	int i,j,tmpBins;
	int Prodsweeps=0;
	FILE *ofp_ptime;
	
	//Using input file, decide which iteration to import
	numf=0; //slightly different than in main, main counts an additional numf when none takes place (numf initial is 1 in main)
	for(lnwlf=ModFactorInit;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
		numf++;
	
	readg();
	numf=111;	//Stands for production run so that new DOS files are labeled
	//write_DOS_H();
	read_diff();
	//write_diff();
		
	//Find lowest and highest energy states
	for(j=0;j<D2BINS;j++)
	{
		LOWESTESTATE[j]=0;
		for(i=1;i<D1BINS;i++)	//Starts from i + 1
		{
			if(mask[i-1][j]==0 && mask[i][j]==1 )
				LOWESTESTATE[j]=i;
			if(mask[i-1][j]==1 && mask[i][j]==0 ) 
				HIGHESTESTATE[j]=i;
		}
	}
	
	//Production run initialization
	lnwlf = 0.0;
	EPSILON = 0.0;
	resetWL();
	EXPLORE=0;
	WLBINVISITS=ProductionBinSamps;	//This is done to make sampleBins() use the # moves per bin for production run
	initialize_prodrun();
	
	ofp_ptime=fopen("prodtimeseries.dat","w");
	//fprintf(ofp_ptime,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
	
	tmpBins = sampledBins();
	
	//fprintf(stderr,"%g\t%g\n",Rgyr2(),gyration_tensor());
	//exit(100);
	
	while( tmpBins > 0.5*D2BINS) //0.5*D2BINS )
	{
		sweepWL(D1BINS*D2BINS);
		Prodsweeps++;
		//if(Prodsweeps%50==0)
		//{
		//	write_DOS_H();
			//MUCAstep();
			//numf++;
			
		//}
		tmpBins = sampledBins();
		fprintf(ofp_ptime,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",JA,currEtot*invN,gyration_tensor()*invN,K1P,K2P,EEdist(),icos_core_count(),tmpBins);
		fflush(ofp_ptime);
	}
	
	fclose(ofp_ptime);
/*	
	for(j=0;j<D2BINS;j++)
	{
		for(i=0;i<D1BINS;i++)
		{
			if(prun_H[i][j] > 0 && mask[i][j]==1)
			{
				fprintf(stdout,"%g\t%g\t%18.10e\t%18.10e\t%18.10e\n",j/invdWLD2+WLD2min,i/invdWLD1+WLD1min,prun_Rgyr2[i][j]/prun_H[i][j],prun_pairnum[i][j]/prun_H[i][j],prun_H[i][j]);
			}
			else
			{
				fprintf(stdout,"%g\t%g\t%18.10e\t%18.10e\t%18.10e\n",j/invdWLD2+WLD2min,i/invdWLD1+WLD1min,0.0,0.0,0.0);
			}
		}
		fprintf(stdout,"\n");
	}
*/
 
	//write_DOS_H();
	MUCAstep();
	//numf=211;
	write_DOS_H();
	
	thermoqs();
	orderqs();
	
	exit(11);
}

void initialize_prodrun(void)
{
	int i,j,k;
	
	//allocate storage for an array of pointers
	prun_pairnum = malloc( D1BINS * sizeof(double *) );
	prun_EEd = malloc( D1BINS * sizeof(double *) );
	prun_Rgyr2 = malloc( D1BINS * sizeof(double *) );
	prun_K1P = malloc( D1BINS * sizeof(double *) );
	prun_K2P = malloc( D1BINS * sizeof(double *) );
	prun_H = malloc( D1BINS * sizeof(double *) );
	
	//Check to see if memory was allocated properly
	if ( (prun_pairnum == NULL) || (prun_EEd == NULL) || (prun_Rgyr2 == NULL) || (prun_K1P == NULL) || (prun_K2P == NULL) || (prun_H == NULL) )
    {
		fprintf(stderr,"\nFailure to allocate memory for 'Production run 2D Arrays'.  See 'initialize_prodrun()'.\n");
        exit(1);
    };
	
	//for each pointer, allocate storage for an array of doubles or ints
	for(k=0;k<D1BINS;k++)
	{
		prun_pairnum[k] = malloc( D2BINS * sizeof(double) );
		prun_EEd[k] = malloc( D2BINS * sizeof(double) );
		prun_Rgyr2[k] = malloc( D2BINS * sizeof(double) );
		prun_K1P[k] = malloc( D2BINS * sizeof(double) );
		prun_K2P[k] = malloc( D2BINS * sizeof(double) );
		prun_H[k] = malloc( D2BINS * sizeof(double) );
	};

	for(i=0;i<D1BINS;i++)
		for(j=0;j<D2BINS;j++)
		{
			prun_pairnum[i][j]=0.0;
			prun_EEd[i][j]=0.0;
			prun_Rgyr2[i][j]=0.0;
			prun_K1P[i][j]=0.0;
			prun_K2P[i][j]=0.0;
			prun_H[i][j]=0.0;
		};	
}

void prunbin( int Estate, int JAstate )
{
	//This function can only appear in the "accepted" portion of the WL function (otherwise, unphysical states might be measured)
	prun_H[Estate][JAstate] += 1.0;
	prun_pairnum[Estate][JAstate] += pairnum();
	prun_EEd[Estate][JAstate] += EEdist();
	//Gyration tensor - must calculate this first, in order to have K1P and K2P values
	prun_Rgyr2[Estate][JAstate] += gyration_tensor();
	prun_K1P[Estate][JAstate] += K1P;	//K1P calculated in gyration_tensor()
	prun_K2P[Estate][JAstate] += K2P;	//K2P "" ""
	
	//fprintf(stderr,"%d\t%d\t%18.10e\t%18.10e\t%18.10e\n",Estate,JAstate,prun_Rgyr2[Estate][JAstate]/prun_H[i][j],prun_pairnum[Estate][JAstate]/prun_H[i][j]);
}

void MUCAstep(void)
{
	int i,j;
	
	//This loop finds the number of unsampled bins
	for(i=0;i<D1BINS;++i)
		for(j=0;j<D2BINS;j++)
			if( (wlH[i][j]>0.0) && (mask[i][j]==1) )
				wllng[i][j] += log(wlH[i][j]);

	//resetWL();
}


//returns 1) wlH_min/wlH_avg and 2) the number of states below the flatness criteria
double flatWL(void)
{
	int i,j,k;
	double avg,min;
  
	k=0;  //number of sampled bins, used to calculate average
	min=1.0e300;  //the minimum sampled bin in the histogram wlH[][]
	avg=0.0;  //average height of the histogram wlH[][]
  
	//This loop finds the flatness and the number of unsampled bins
	for(i=0;i<D1BINS;++i)
		for(j=0;j<D2BINS;j++)
		{
			//minimum of the histogram
			if( (wlH[i][j]<min) )
			{
				min=wlH[i][j];
			};
			//average of the histogram
			if( (wlH[i][j]>0.0) )
			{
				avg+=wlH[i][j];
				k+=1;
			};	
		};
  
	avg/=1.0*k;  //Calculates the average height of the histogram (for  wlH > 0.0)
	if(k==0)
		avg=1.0;  //This keeps from dividing by zero
  
	//This loop counts the number of bins below the flatness criteria
  /*
    for(i=0;i<D1BINS;++i)
    for(j=0;j<D2BINS;++j)
    if( wlH[i][j]<Flatness*avg )
    numbelow_flat+=1.0;
    //Stores the percentage of states below the flatness criteria	
    numbelow_flat = (numbelow_flat)/(1.0*D1BINS*D2BINS);
  */
	//returns the current flatenss of the histogram		
	return min/avg;
  
}

//returns 1) wlH_min/wlH_avg and 2) the number of states below the flatness criteria
int sampledBins(void)
{
	int i=0,j=0,k=0;
		
	//This loop finds the number of unsampled bins
	for(i=0;i<D1BINS;++i)
		for(j=0;j<D2BINS;j++)
		{
			//Checking histogram in normal WL run
			if( (ProductionRun == 0) && (wlH[i][j]<=WLBINVISITS) && (mask[i][j]==1) )
				k+=1;
		
			//Checking histogram in normal WL run
			if( (ProductionRun == 1) && (prun_H[i][j]<=WLBINVISITS) && (mask[i][j]==1) )
				k+=1;
		};
	
	//fprintf(stderr,"%d\n",k);
	return k;
	
}

//reset mask (resets the mask by certain % of the sampled bins, and according to numf)
void resize_mask(double PERCENTAGE, int CURRNUMF)
{
	int i,j;
  
	if(EXPLORE == 1 && numf == CURRNUMF)
	{
		for(j=0;j<D2BINS;j++)
		{
			for(i=0;i<D1BINS;i++)
			{
				if(mask[i][j]==0 && i<HIGHESTESTATE[j])
				{
					LOWESTESTATE[j]=i;
					//fprintf(stderr,"%d\t%d\n",j,LOWESTESTATE[j]);
				}
			}
				
			//fprintf(stderr,"%d\tOLD: %d\t%g\t",j,LOWESTESTATE[j],LOWESTESTATE[j]/invdWLD1+WLD1min);
			//fprintf(stderr,"FACTOR:%d\t",(D1BINS - LOWESTESTATE[j]) - (int)( PERCENTAGE * (D1BINS - LOWESTESTATE[j]) ) );
			LOWESTESTATE[j] = LOWESTESTATE[j] + (D1BINS - LOWESTESTATE[j]) - (int)( PERCENTAGE * (D1BINS - LOWESTESTATE[j]) );
			//fprintf(stderr,"NEW: %d\t%g\n",LOWESTESTATE[j],LOWESTESTATE[j]/invdWLD1+WLD1min);
				
			//exit(100);
		}
		//exit(100);
		EXPLORE=0;

		//RESET Density of States
		for(j=0;j<D2BINS;j++)
			for(i=0;i<D1BINS;i++)
			{
				if(i>=LOWESTESTATE[j] && i<HIGHESTESTATE[j])
				{
					mask[i][j]=1;
				}
				else
				{
					mask[i][j]=0;
					wllng[i][j]=0.0;
				}
			}

	};
}


//reset accumulated histogram
void resetWL()
{
	int i,j;
  
	for(i=0;i<D1BINS;++i)
		for(j=0;j<D2BINS;j++)
			wlH[i][j]=0.0;
  
}


//reads in the g.dat file
void readg(void)
{
	int i,j;
	FILE *ifp;
	char s1[512];
	double tmp;
  
	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	//sprintf(s1,"DOS_H_iter%d.dat",numf);
	ifp=fopen(s1,"r");
	if ( ifp == NULL )
    {
		fprintf(stderr,"\nFailure to read in DOS values from file.  See 'readg()'.\n");
        exit(10);
    }
	
	//Reading in the DOS and histogram from the restart function
	fscanf(ifp,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&lnwlf,&TotalSweeps,&IterSweeps,&tmp);
	//fprintf(stderr,"#%d\t%lg\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,tmp);
	for(j=0;j<D2BINS;j++)
	{
		for(i=0;i<D1BINS;i++)
		{				
			//fscanf(ifp,"%lg\t%lg\t%lg\t%lg\n",&tmp,&tmp,&(wllng[i][j]),&(wlH[i][j]));
			fscanf(ifp,"%lg\t%lg\t%lg\t%lg\t%d\n",&tmp,&tmp,&(wllng[i][j]),&(wlH[i][j]),&(mask[i][j]));
			//fprintf(ifp,"%lg\t%lg\t%18.10e\t%18.10e\t%d\n",tmp,tmp,(wllng[i][j]),(wlH[i][j]),(mask[i][j]));
		}
	};
	fclose(ifp);

}

//writes out the DOS and the histogram
void write_DOS_H(void)
{
	int i,j;
	FILE *ofp;
	char s1[512];
	double maxg=-1.0e300;
	double ming=1.0e300;
  
	//uncomment the following line to output in log_10
	//tmp=log10(exp(1.0));
	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	ofp=fopen(s1,"w");
	//ofp=fopen("g.dat","w");
  
	//Find the minimum of the DOS
	for(i=0;i<D1BINS;++i)
		for(j=0;j<D2BINS;j++)
		{
			if( wllng[i][j]>maxg && mask[i][j]==1)
				maxg=wllng[i][j];
			//if( wllng[i]<ming )
				//ming=wllng[i];
		};
  	
	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#%d\t%g\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL());
	for(j=0;j<D2BINS;j++)
	{
		for(i=0;i<D1BINS;i++)
		{
			//Printing Out 1D Data
			//fprintf(ofp,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i][j]);
			//if(mask[i][j]==1)
				fprintf(ofp,"%g\t%g\t%18.10e\t%18.10e\t%d\n",j/invdWLD2+WLD2min,i/invdWLD1+WLD1min,(wllng[i][j]),wlH[i][j],mask[i][j]);
		}
		fprintf(ofp,"\n");
	};
	fflush(ofp);
	fclose(ofp);
}


void read_diff(void)
{
	int i,j;
	FILE *ifp;
	double tmp;
	
	ifp=fopen("diffvalues.dat","r");
	if ( ifp == NULL )
    {
		fprintf(stderr,"\nFailure to read in 'diff' values from file.  See 'read_diff()'.\n");
        exit(10);
    }
	
	for(j=0;j<D2BINS;j++)
	{
		for(i=0;i<D1BINS;i++)
		{
			//Printing Out 1D Data
			//fprintf(ofp,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i][j]);
			//if(mask[i][j]==1)
			fscanf(ifp,"%lg\t%lg\t%lg\t%lg\n",&tmp,&tmp,&(diffvalues[i][j]),&tmp);
			//fprintf(stderr,"%18.10e\n",diffvalues[i][j]);
		}
		fprintf(ifp,"\n");
	};
	
	fflush(ifp);
	fclose(ifp);
}


void write_diff(void)
{
	int i,j;
	FILE *ofp;
	
	ofp=fopen("diffvalues.dat","w");
	
	for(j=0;j<D2BINS;j++)
	{
		for(i=0;i<D1BINS;i++)
		{
			//Printing Out 1D Data
			//fprintf(ofp,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i][j]);
			//if(mask[i][j]==1)
			if(diffattempts[i][j]==0.0)
				diffattempts[i][j]=1.0;
			
			fprintf(ofp,"%g\t%g\t%18.10e\t%18.10e\n",j/invdWLD2+WLD2min,i/invdWLD1+WLD1min,diffvalues[i][j],JArate[i][j]/JAattempts[i][j]);
			//fprintf(ofp,"%g\t%18.10e\t%18.10e\t%18.10e\t%18.10e\t%18.10e\n",i/invdWLD1+WLD1min,diffvalues[i],diffvalues[i]/diffattempts[i],crankvalues[i],crankrate[i]/crankattempts[i],cutjoinrate[i]/cutjoinattempts[i]);
		
		}
		fprintf(ofp,"\n");
	};
	
	fflush(ofp);
	fclose(ofp);
}

void read_EMAX(void)
{
	int i,j;
	FILE *ifp;
	double tmp=0.0,Esigma=0.0;
	
	ifp=fopen("EMAX_METROPvalues.input","r");
	if ( ifp == NULL )
    {
		fprintf(stderr,"\nFailure to read in EMAX values from file.  See 'read_EMAX()'.\n");
        exit(10);
    }
	else
	{
		for(j=0;j<D2BINS;j++)
		{
			fscanf(ifp,"%lg\t%lg\t%lg\n",&(tmp),&(HIGHESTE[j]),&(Esigma));
			//fprintf(stderr,"%g\t%g\t%g\n",tmp,HIGHESTE[j],Esigma);
			HIGHESTESTATE[j] = (int) (( (HIGHESTE[j] + 3.0*Esigma)-WLD1min)*invdWLD1);
			if(HIGHESTESTATE[j] >=D1BINS)
				HIGHESTESTATE[j] = D1BINS;//(int) ((WLD1max-WLD1min)*invdWLD1);
		
			//fprintf(stderr,"%d\t%g\n",HIGHESTESTATE[j],HIGHESTE[j]);
		};
	}
	
	fflush(ifp);
	fclose(ifp);

}

void write_EMIN(void)
{
	int i,j;
	FILE *ofp;
	
	ofp=fopen("EMIN_WLvalues.input","w");
	
	for(j=0;j<D2BINS;j++)
	{
		fprintf(ofp,"%g\t%g\t%d\n",j/invdWLD2+WLD2min,LOWESTE[j],LOWESTESTATE[j]);
	};

	fflush(ofp);
	fclose(ofp);
}

//routine to initialize things for the Wang-Landau simulation
void initWL(void)
{
	int i=0,j=0,k=0,istate=0,iS=0;

	//make sure that the polymer is initialized
	initialize();
	TotalSweeps=0;		//Total number of sweeps
	IterSweeps=0;			//Number of sweeps per iteration
	MASKCUTLOWE=0;
	EPSILON = 0.001;	//Shift in energy dependent mc moves
	
	//Primary Binning Direction for WL Simulation
	D1BINS = (int)((WLD1max - WLD1min)/(dWLD1*invN));
	//Secondary Binning Direction for WL Simulation
	D2BINS = (int)((WLD2max - WLD2min)/(dWLD2));
	
	NUM_UNSAMPLED_BINS = 0.5*D2BINS;
	
	//Inverse Bin Width
	invdWLD1=1.0/(dWLD1*invN);
	invdWLD2=1.0/(dWLD2);
    
	//this value times Ebond-WLEbondmin, gives the bin number;
	//WLinvdEbond=(1.0*Lbond)/(WLEbondmax-WLEbondmin);
	//WLinvdEnonbond=(1.0*Lnonbond)/(WLEnonbondmax-WLEnonbondmin);

	//Dynamic Memory Allocation - this has to be done because an input file is used.
	//WL 2D Arrays - Histogram, Density of States, Mask, and Rawmask
	//allocate storage for an array of pointers
	wlH = malloc( D1BINS * sizeof(double *) );
	wllng = malloc( D1BINS * sizeof(double *) );
	mask = malloc( D1BINS * sizeof(int *) );
	LOWESTE = malloc( D2BINS * sizeof(double) );
	LOWESTESTATE = malloc( D2BINS * sizeof(int) );
	HIGHESTE = malloc( D2BINS * sizeof(double) );
	HIGHESTESTATE = malloc( D2BINS * sizeof(int) );
	
	//Check to see if memory was allocated properly
	if ( (wlH == NULL) || (wllng == NULL) || (mask == NULL) )
    {
		fprintf(stderr,"\nFailure to allocate memory for 'Wang-Landau 2D Arrays'.  See 'initWL()'.\n");
        exit(1);
    };
		
	diffvalues = malloc( D1BINS * sizeof(double *) );
	diffrate = malloc( D1BINS * sizeof(double *) );
	diffattempts = malloc( D1BINS * sizeof(double *) );
	
	/*
		crankvalues = malloc( D1BINS * sizeof(double *) );
		crankrate = malloc( D1BINS * sizeof(double *) );
		crankattempts = malloc( D1BINS * sizeof(double *) );
		cutjoinrate = malloc( D1BINS * sizeof(double *) );
		cutjoinattempts = malloc( D1BINS * sizeof(double *) );
	*/
		JArate = malloc( D1BINS * sizeof(double *) );
		JAattempts = malloc( D1BINS * sizeof(double *) );
	
		
	//for each pointer, allocate storage for an array of doubles or ints
	for(k=0;k<D1BINS;k++)
	{
		wlH[k] = malloc( D2BINS * sizeof(double) );
		wllng[k] = malloc( D2BINS * sizeof(double) );
		mask[k] = malloc( D2BINS * sizeof(int) );
		
		diffvalues[k] = malloc( D2BINS * sizeof(double) );
		diffrate[k] = malloc( D2BINS * sizeof(double) );
		diffattempts[k] = malloc( D2BINS * sizeof(double) );
		
		/*
			crankvalues[k] = malloc( D2BINS * sizeof(double) );
			crankrate[k] = malloc( D2BINS * sizeof(double) );
			crankattempts[k] = malloc( D2BINS * sizeof(double) );
			cutjoinrate[k] = malloc( D2BINS * sizeof(double) );
			cutjoinattempts[k] = malloc( D2BINS * sizeof(double) );
		*/
			JArate[k] = malloc( D2BINS * sizeof(double) );
			JAattempts[k] = malloc( D2BINS * sizeof(double) );
		
	};
		
	//Run the standard MC routine until the configuration has energy within the WL simulation energy range
	T=0.5;
	while( (currEtot*invN>WLD1max) || (currEtot*invN<WLD1min) )
    {
		fprintf(stderr,"Relaxing:  %g of %g \n",currEtot*invN,WLD1max);
		mchybrid();
    };

	//Initialize Arrays
	for(i=0;i<D1BINS;++i)
	{
		for(j=0;j<D2BINS;j++)
		{
			wllng[i][j]=0.0;
			wlH[i][j]=0.0;
			
			mask[i][j]=0;
			if(EXPLORE == 0)
				mask[i][j]=1;
			
			diffvalues[i][j]= DD;//(double)(((0.3 - 0.01)/D1BINS)*i) + 0.01;
			diffrate[i][j]= 0.0;	
			diffattempts[i][j]= 0.0;
			
			/*
				crankvalues[i][j]= 0.5*M_PI;//(double)(((0.3 - 0.01)/D1BINS)*i) + 0.01;
				crankrate[i][j]= 0.0;
				crankattempts[i][j]= 0.0;
				cutjoinrate[i][j]= 0.0;
				cutjoinattempts[i][j]= 0.0;
			*/
				JArate[i][j]= 0.0;	
				JAattempts[i][j]= 0.0;
			
		}
	};

	for(j=0;j<D2BINS;j++)
	{
		LOWESTE[j]=WLD1max;
		LOWESTESTATE[j]=0;
		HIGHESTE[j]=WLD1min;
		HIGHESTESTATE[j] = D1BINS;//(int) ((WLD1max-WLD1min)*invdWLD1);
		//HIGHESTESTATE[j]=0;
		//fprintf(stderr,"%d\t%d\t%g\n",j,LOWESTESTATE[j],LOWESTE[j]);
	};

	read_EMAX();

	//Initialize current energy state and mask
	//LOWESTE=WLD1max;//currEtot*invN;
	istate=(int) ((currEtot*invN-WLD1min)*invdWLD1);
	iS=(int) ((JA-WLD2min)*invdWLD2);
	//fprintf(stderr,"%g\t%d\n",JA,iS);
	//for(i=istate;i<D1BINS;i++)
		//for(j=jJA;j<D2BINS;j++)
			mask[istate][iS] = 1;
 
}


//attempts to run one WL hybrid move per bin
void sweepWL(double sweeps)
{
	int i;
      
	for(i=0;i<sweeps;++i)
    {
		wlhybrid();
	};
}


//excutes all of the combined MC moves, performs one step of each type of move, and one sweep of diffusion moves
void wlhybrid(void)
{
	int i;
  
	for(i=0;i<N;++i)
	{
		//wldiff();
		wldiffE();
		//wlsinglecrank();
	}
	
	wlJAchange();
	wlreptation();	
	wlcutjoin();
	
	//wlcrankshaft();	//Should be unaffected by non-zero JA (but not checked)
	//wlrandpivot();	//Needs to have bond-angle calculations added for non-zero JA
	
}

int WangLandau(double Ei, double Ef, double DBvalue, double JAi, double JAf)
{
	int iti,iS;//indices of initial config
	int fti,fS;//indices of final config
	int tmpindex;
	
	double lngi,lngf;  //values for the density of states, intial and final
	double R;
	
	//Primary direction index
	iti=(int) ((Ei*invN-WLD1min)*invdWLD1);
	fti=(int) ((Ef*invN-WLD1min)*invdWLD1);
	
	//Secondary direction index
	iS=(int) ((JAi-WLD2min)*invdWLD2);
	fS=(int) ((JAf-WLD2min)*invdWLD2);
	
	//fprintf(stderr,"%g\t%g\n",JA,currEtot*invN);
	
	//Keeping track of the number of moves
	
	//fprintf(stderr,"%d\t%d\t%d\t%d\t%g\t%g\t%d\n",iti,fti,iS,fS,Ei*invN,Ef*invN,D1BINS);
	
	//This statement simply prints out the lowest energy configuration
	if(Ef*invN<LOWESTE[fS] && ProductionRun==0)
	{
		writeinc(fS,Ef*invN);
		write_mol2(fS);
		LOWESTE[fS]=Ef*invN;
		//fprintf(stderr,"New Lowest E config %g\n\n",LOWESTE);
	};
	
	//  fprintf(stderr,"%d,%d >> %d,%d\n",iti,itj,fti,ftj);
		
	//inside of bounds
	//if((fti<0)||(fti>=D1BINS)||(iS<0)||(fS>=D2BINS))
	//if((fti<LOWESTESTATE[fS])||(fti>=D1BINS)||(iS<0)||(fS>=D2BINS))
	if((fti<LOWESTESTATE[fS])||(fti>=HIGHESTESTATE[fS])||(iS<0)||(fS>=D2BINS))
    {
		//reject, outside of bounds
		wlH[iti][iS]+=1.0;
		wllng[iti][iS]+=lnwlf;
		
		//if(ProductionRun==1)
		//	prunbin( iti, iS );
		
		return 0;
    }
	else
	{
		lngi=wllng[iti][iS];
		lngf=wllng[fti][fS];
		
		R=DBvalue*exp(lngi-lngf);
		
		if(randd1()<R)
		{
			//accept
			currEtot=Ef;
			//JA = JAf;
			
			if(mask[fti][fS]==0 && EXPLORE==1 )
			{
				mask[fti][fS]=1;
				diffvalues[fti][fS]=diffvalues[iti][iS];
			}
			wlH[fti][fS]+=1.0;
			wllng[fti][fS]+=lnwlf;
			
			if(ProductionRun==1)
				prunbin( fti, fS );
				
			return 1;
		}
		else
		{
			//reject
			wlH[iti][iS]+=1.0;
			wllng[iti][iS]+=lnwlf;
			
			//if(ProductionRun==1)
			//	prunbin( iti, iS );
				
			return 0;
		};
    };
	
}

//Changes the stiffness of the chain
void wlJAchange(void)
{
	int i=0,iES=0,iJAS=0;
    double Ei=0.0,Ef=0.0; //intial and final energy values(before and after)
	double Ebai=0.0,Ebaf=0.0;
	double JAi=0.0,JAf=0.0; //intial and final stiffness values(before and after)
	double delta=3.0;

	Ei = currEtot;
	JAi=JA;
	Ebai = bondangleEnergy();
	
	iES=(int) ((Ei*invN-WLD1min)*invdWLD1);
	iJAS=(int) ((JAi-WLD2min)*invdWLD2);
	
	//randomly choose to increase or decrease stiffness
	if(randd1() < 0.5)
	{
		JAf = JA - dWLD2*delta*randd1();//delta;
		//if(JAf < WLD2min)
		//	JAf = JAi;
	}
	else
	{
		JAf = JA + dWLD2*delta*randd1();//*delta;
		//if(JAf >=WLD2max)
		//	JAf = JAi;
	}
	
	if(JAf<WLD2min || JAf>WLD2max)
		JAf=JAi;
	
	JA = JAf;
	Ebaf = bondangleEnergy();
	Ef = Ei + (Ebaf - Ebai);  
	//Ef = Etot();  
	
	JAattempts[iES][iJAS] += 1.0;
	
	//fprintf(stderr,"MOVE: %g\t%g\n",JA,fabs(JAi-JAf));
	if(WangLandau(currEtot,Ef,1.0,JAi,JAf)==1)
	{
		JArate[iES][iJAS] += 1.0;
	}
	else
	{
		JA=JAi;
	}
}

//////////  Wang-Landau Energy Dependent Diffusion Move  //////////
//Displaces a random monomer by a small amount that depends on energy
void wldiffE(void)
{
	int i=0,istate=0,fstate=0,fJAS=0,moviestate=0;
	double ox=0.0,oy=0.0,oz=0.0; //old positions of the monomer
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double Wforward=0.0,Wbackward=0.0;
	//double epsilon = 0.001;
	double INFENERGY=1.0e10;
	
	//CONSTANT JA
	fJAS=(int) ((JA-WLD2min)*invdWLD2);
	
	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
	//Ei=Etot();
	Ei=currEtot;
	istate=(int) ((Ei*invN-WLD1min)*invdWLD1);
	Ebondi=localE(i);
	Enblocali=localnonbondEnergy(i);
	
	//diffuse the ith monomer
	ox=x[i];
	oy=y[i];
	oz=z[i];
	
	x[i]+=diffvalues[istate][fJAS]*(2.0*randd1()-1.0);
	y[i]+=diffvalues[istate][fJAS]*(2.0*randd1()-1.0);
	z[i]+=diffvalues[istate][fJAS]*(2.0*randd1()-1.0);
	
    //Calculate the final energy
	//Ef=Etot();
	Ebondf=localE(i);
	Enblocalf=localnonbondEnergy(i);
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
	fstate=(int) ((Ef*invN-WLD1min)*invdWLD1);
	
	//MUST reject the move in the following two scenarios
	//#1 - The final state is outside the energy boundaries
	//#2 - The new volume does not overlap with the old volume (detailed balance)
	
	//fprintf(stderr,"%18.10e\t%d\t%d\t%d\n",Ei*invN - Ef*invN,istate,fstate,fJAS);
	//fprintf(stderr,"%d\t%d\t%d\n",istate,fstate,fJAS);	
	//fprintf(stderr,"%g\n",fabs(Etot()-Ef));
	
	diffattempts[istate][fJAS] += 1.0;
	
	if( (fstate>=0) && (fstate<D1BINS) )
	{
		if( EPSILON != 0.0 )
		{
			if(Ef < Ei)
			{
				if(diffvalues[istate][fJAS] < 0.2)
				{
					diffvalues[istate][fJAS] = diffvalues[istate][fJAS]*(1.0 + 2.0*EPSILON);
				}
			}
			else if(diffvalues[istate][fJAS]>0.0001)
			{
				diffvalues[istate][fJAS] = diffvalues[istate][fJAS]*(1.0 - 1.0*EPSILON);
			};
		};
			
		if( (fabs(x[i]-ox) < diffvalues[fstate][fJAS]) && (fabs(y[i]-oy) < diffvalues[fstate][fJAS]) && (fabs(z[i]-oz) < diffvalues[fstate][fJAS]) )
		{
			Wforward = diffvalues[fstate][fJAS]*diffvalues[fstate][fJAS]*diffvalues[fstate][fJAS];
			Wbackward = diffvalues[istate][fJAS]*diffvalues[istate][fJAS]*diffvalues[istate][fJAS];
			
			//fprintf(stderr,"%g\n",Wforward/Wbackward);
			attemptdiff+=1;
			if(WangLandau(Ei,Ef,Wbackward/Wforward,JA,JA)==1)
			{
				acceptdiff+=1;
				diffrate[istate][fJAS] += 1.0;
			}
			else
			{
				//reject - return monomer to old position
				x[i]=ox;
				y[i]=oy;
				z[i]=oz;
			}
		}
		else
		{
			attemptdiff+=1;
			if(WangLandau(Ei,Ef,-1.0,JA,JA)==1)
			{
				acceptdiff+=1;
				diffrate[istate][fJAS] += 1.0;
			}
			else
			{
				//reject - return monomer to old position
				x[i]=ox;
				y[i]=oy;
				z[i]=oz;
			}
		}
	}
	else
	{
		attemptdiff+=1;
		if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
		{
			acceptdiff+=1;
			diffrate[istate][fJAS] += 1.0;
		}
		else
		{
			//reject - return monomer to old position
			x[i]=ox;
			y[i]=oy;
			z[i]=oz;
		}
	}	
}

void wlsinglecrank(void)
{
	int m=0;
    int istate=0,fstate=0,fJAS=0;
	double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double rx=0.0,ry=0.0,rz=0.0;
	double distr=0.0,rangle=0.0,c=0.0,s=0.0,u=0.0;
	double ox=0.0,oy=0.0,oz=0.0; //old positions of the monomer
	double ttx=0.0,tty=0.0,ttz=0.0;
	double Wforward=0.0,Wbackward=0.0;
	//double epsilon=0.001;
	
	//CONSTANT JA
	fJAS=(int) ((JA-WLD2min)*invdWLD2);
	
	//randomly choose two monomers
	m=(int) (randd1()*N);
	
	//we don't want m to be one of the ends
	while( (m==0) || (m==N-1))  
	{
		m=(int) (randd1()*N);
	};
	
	//set temporary positions
	ox=x[m];
	oy=y[m];
	oz=z[m];
	
	//Set initial energies
	Ei=currEtot;
	istate=(int) ((Ei*invN-WLD1min)*invdWLD1);
	
	//choose random rotation angle
	rangle=DC*M_PI*(2.0*randd1()-1.0);
	//rangle=crankvalues[istate][fJAS]*(2.0*randd1()-1.0);
	
	//Calculate the initial local bond and bond angle energy
	Ebondi = localE(m);
	//Calculate the initial local nonbonded energy
	Enblocali=localnonbondEnergy(m);
	
	//fprintf(stderr,"LIST:  \n%d\t%d\t%d\n",m-1,m,m+1);
	
	//////////Rotation
	//rotation axis is the vector between m1 and m2
	distr = sqrt( (x[m-1]-x[m+1])*(x[m-1]-x[m+1]) + (y[m-1]-y[m+1])*(y[m-1]-y[m+1]) + (z[m-1]-z[m+1])*(z[m-1]-z[m+1]) );
	
	//fprintf(stderr,"\nROTATION: %18.10e\t%18.10e\n%18.10e\t%18.10e\n%18.10e\t%18.10e\n%g\t%d\n\n",x[m-1],x[m+1],y[m-1],y[m+1],z[m-1],z[m+1],distr,m);
	
	rx = (x[m-1]-x[m+1])/distr;
	ry = (y[m-1]-y[m+1])/distr;
	rz = (z[m-1]-z[m+1])/distr;
	
	//fprintf(stderr,"ROTATION: %18.10e\t%18.10e\t%18.10e\t%18.10e\t%d\n",distr,rx,ry,rz,m);
	
	//angle constants
	c=cos(rangle);
	s=sin(rangle);
	u=1.0-c;
	
	/////Rotation
	//now we perform the rotation
	//translate the rotation segment to the origin
	x[m] -= x[m-1];
	y[m] -= y[m-1];
	z[m] -= z[m-1];
		
	//store old vector for the rotation matrix below
	ttx=x[m];
	tty=y[m];
	ttz=z[m];
		
	//Perform the rotation on the selected segment
	x[m] = (u*rx*rx + c)*ttx     + (u*ry*rx - s*rz)*tty  + (u*rz*rx + ry*s)*ttz;
	y[m] = (u*rx*ry + rz*s)*ttx  + (u*ry*ry + c)*tty     + (u*rz*ry - rx*s)*ttz;
	z[m] = (u*rx*rz - ry*s)*ttx  + (u*ry*rz + rx*s)*tty  + (u*rz*rz + c)*ttz;	
		
	//translate the rotation segment back to the original reference frame
	x[m] += x[m-1];
	y[m] += y[m-1];
	z[m] += z[m-1];
	
	//Calculate the initial local bond and bond angle energy
	Ebondf = localE(m);
	//Calculate the initial local nonbonded energy
	Enblocalf=localnonbondEnergy(m);
	
	//Final energy
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
	fstate=(int) ((Ef*invN-WLD1min)*invdWLD1);
	
	if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
	{
		//crankrate[istate][fJAS] += 1.0;
	}
	else
	{
		//reject - return monomer to old position
		x[m]=ox;
		y[m]=oy;
		z[m]=oz;
	}
	
	//fprintf(stderr,"%18.10e\t%d\t%d\t%d\n",Ei*invN - Ef*invN,istate,fstate,fJAS);
	
	/*
	crankattempts[istate][fJAS] += 1.0;
	
	if( (fstate>=0) && (fstate<D1BINS) )
	{
		if( EPSILON != 0.0 )
		{
			if(Ef < Ei)
			{
				if(crankvalues[istate][fJAS] < 0.5*M_PI)
				{
					crankvalues[istate][fJAS] = crankvalues[istate][fJAS]*(1.0 + 2.0*EPSILON);
					//crankratedown[istate] = crankratedown[istate]+1.0;
				}
			}
			else
			{
				crankvalues[istate][fJAS] = crankvalues[istate][fJAS]*(1.0 - 1.0*EPSILON);
				//crankrateup[istate] = crankrateup[istate]+1.0;
			};
		};
		
		if( rangle < crankvalues[fstate][fJAS] )  //Make sure move can take you back to old position (meaning fstate must overlap current angle)
		{
			Wforward = crankvalues[fstate][fJAS];
			Wbackward = crankvalues[istate][fJAS];
			
			//fprintf(stderr,"%g\n",Wforward/Wbackward);
			if(WangLandau(Ei,Ef,Wbackward/Wforward,JA,JA)==1)
			{
				crankrate[istate][fJAS] += 1.0;
			}
			else
			{
				//reject - return monomer to old position
				x[m]=ox;
				y[m]=oy;
				z[m]=oz;
			}
		}
		else
		{
			if(WangLandau(Ei,Ef,-1.0,JA,JA)==1)
			{
				crankrate[istate][fJAS] += 1.0;
			}
			else
			{
				//reject - return monomer to old position
				x[m]=ox;
				y[m]=oy;
				z[m]=oz;
			}
		}
	}
	else
	{
		if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
		{
			crankrate[istate][fJAS] += 1.0;
		}
		else
		{
			//reject - return monomer to old position
			x[m]=ox;
			y[m]=oy;
			z[m]=oz;
		}
	};
	*/
}


//////////  Wang-Landau Diffusion Move  //////////
//Displaces a random monomer by a small amount 
void wldiff(void)
{
	int i=0;
	double ox=0.0,oy=0.0,oz=0.0; //old positions of the monomer
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables

	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
	//Ei=Etot();
	Ei=currEtot;
	Ebondi=localE(i);
	Enblocali=localnonbondEnergy(i);

	//diffuse the ith monomer
	ox=x[i];
	oy=y[i];
	oz=z[i];

	x[i]+=DD*(2.0*randd1()-1.0);
	y[i]+=DD*(2.0*randd1()-1.0);
	z[i]+=DD*(2.0*randd1()-1.0);

    //Calculate the final energy
	//Ef=Etot();
	Ebondf=localE(i);
	Enblocalf=localnonbondEnergy(i);
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
		
	attemptdiff+=1;
  
	if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
    {
		acceptdiff+=1;
	}
	else
    {
		//reject - return monomer to old position
		x[i]=ox;
		y[i]=oy;
		z[i]=oz;
    };	
	
}

//////////  Wang-Landau Reptation Move  //////////
//performs one slithering snake MC attempt
void wlreptation(void)
{
	int i=0;
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double rx=0.0,ry=0.0,rz=0.0,bondr=0.0; //used for random direction reptated monomer
	double eta1=0.0,eta2=0.0,etasq=0.0;
		
	//store the temporary positions
	for(i=0;i<N;i++)
	{
		tx[i]=x[i];
		ty[i]=y[i];
		tz[i]=z[i];
	};
	
	//choose a random direction using the Marsaglia method
	//give random unit vectors for the relative position of the moved endpoint
	etasq = 2.0;  //simply set to a value greater than 1.0 (the while condition below)
    while( etasq > 1.0 )
	{
		eta1 = 1.0 - 2.0*randd1();
		eta2 = 1.0 - 2.0*randd1();
	
		etasq = eta1*eta1 + eta2*eta2;
	};

	//these are the new unit vectors
	rx = 2.0*eta1*sqrt(1.0-etasq);
	ry = 2.0*eta2*sqrt(1.0-etasq);
	rz = 1.0 - 2.0*etasq;
	
	//Set initial energies
	//Ei=Etot();
	Ei=currEtot;
	
	//choose whether to move the 0 monomer or N-1 monomer
	if(randd1()<0.5)
	{
		//move the 0th monomer
		//measure the bond length of 0-1 bond
		bondr=sqrt( (x[0]-x[1])*(x[0]-x[1]) + (y[0]-y[1])*(y[0]-y[1]) + (z[0]-z[1])*(z[0]-z[1]) );
	    
		//store local bond energy of the 0th monomer (even though there shouldn't be any change)
		Ebondi = localE(0);
		//store local non-bond energy 
		Enblocali=localnonbondEnergy(0);
		
		//move x[1] to x[0], etc...
		for(i=0;i<N-1;i++)
		{
			x[i]=tx[i+1];
			y[i]=ty[i+1];
			z[i]=tz[i+1];
		};
		
		//new tail has random direction, but same bond length
		x[N-1]=x[N-2]+bondr*rx;
        y[N-1]=y[N-2]+bondr*ry;
        z[N-1]=z[N-2]+bondr*rz;
	
		//Calculate new bond energy
		Ebondf=localE(N-1);
		//Calculate new non-bond energy
		Enblocalf=localnonbondEnergy(N-1);
	}
	else
	{
		//move the N-1 monomer	
		//measure the bond length of (N-2) to (N-1) bond
		bondr=sqrt( (x[N-2]-x[N-1])*(x[N-2]-x[N-1]) + (y[N-2]-y[N-1])*(y[N-2]-y[N-1]) + (z[N-2]-z[N-1])*(z[N-2]-z[N-1]) );
		
		//store local energy change of the N-1 monomer (even though it will not change)
		Ebondi = localE(N-1);
		//store local non-bond energy 
		Enblocali=localnonbondEnergy(N-1);
		
		//move x[1] to x[0], etc...
		for(i=1;i<N;i++)
		{
			x[i]=tx[i-1];
			y[i]=ty[i-1];
			z[i]=tz[i-1];
		};
		
		//new tail has random direction, but same bond length
		x[0]=x[1]+bondr*rx;
        y[0]=y[1]+bondr*ry;
        z[0]=z[1]+bondr*rz;
	
		//Calculate new change in bond energy
		Ebondf=localE(0);
		//Calculate new non-bond energy
		Enblocalf=localnonbondEnergy(0);
	}

    //Calculate the final energy (Note: the bond energy should not change, so it isn't calculated)
	//Ef=Etot();
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
	
    attemptsnake+=1;
	
	if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
    {
		acceptsnake+=1;
	}
	else
    {
		//reject - return monomer to old position
		for(i=0;i<N;i++)
		{
			x[i]=tx[i];
			y[i]=ty[i];
			z[i]=tz[i];
		};
	};	
	
}


//////////  Wabg-Landau Random Pivot  //////////
//Performs a random pivot by selecting a single monomer, and rotating the chain around a random direction
//Note: previously the left or right portion of the chain was chosen at random
//However, now the shortest chain length is chosen.
void wlrandpivot(void)
{
	int m=0,i=0,j=0,direction=0,start=0,end1=0,end2=0;
	double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double eta1=0.0,eta2=0.0,etasq=0.0,rx=0.0,ry=0.0,rz=0.0;
	double distr=0.0,rangle=0.0,c=0.0,s=0.0,u=0.0;
	double ttx=0.0,tty=0.0,ttz=0.0;

	//set temperorary variables
	for(j=0;j<N;++j)
    {
		tx[j]=x[j];
		ty[j]=y[j];
		tz[j]=z[j];
    };

	//choose random monomer
	m=(int) (randd1()*(N-2)+1.0);  //note, never choses monomers on the ends, as this would do nothing to the internal degrees of freedom
	
	//choose random rotation angle
	rangle=DC*2.0*M_PI*(2.0*randd1()-1.0);

	//select the shortest portion of the chain to the left or right of the random monomer
	if(m < N/2.0)
	{
		direction = -1;
		start=m;
		end1=0;
		end2=N-1;
	}
	else
	{
		direction = 1;
		start=m;
		end1=N-1;
		end2=0;
	};
		
	//choose a random direction using the Marsaglia method
	etasq = 2.0;  //simply set to a value greater than 1.0 (the while condition below)
    while( etasq > 1.0 )
	{
		eta1 = 1.0 - 2.0*randd1();
		eta2 = 1.0 - 2.0*randd1();
	
		etasq = eta1*eta1 + eta2*eta2;
	};
	
	//these are the new unit vectors
	rx = 2.0*eta1*sqrt(1.0-etasq);
	ry = 2.0*eta2*sqrt(1.0-etasq);
	rz = 1.0 - 2.0*etasq;
			
	//angle constants
	c=cos(rangle);
	s=sin(rangle);
	u=1.0-c;

	//set energies
	Ei=currEtot;
	
	//Calculate the initial local bond energy (NOTE: this cannot go in the rotation loop)
/*	i=start;
	while( i != end1 )
	{
		i=i+direction;
		j=i-direction;
		distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
		Ebondi = Ebondi + bondpot(distr);	
		//fprintf(stderr,"%d\t%d\t%d\t%f\t%f\n",m,i,i-direction,distr,bondpot(distr));
	};
*/
	//Local bond and bond angle energy
	Ebondi = localE(m);
	
/////Rotation
	//now we perform the rotation
	i=start;
	while( i != end1 )
	{
		i=i+direction;
		
		//Calculate the initial local nonbond energy
		for(j=end2;j!=start;j+=direction)
		{	
			distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
			Enblocali = Enblocali + nonbondpot(distr);
			//fprintf(stderr,"%d\t%d\n",i,j);
		};
		
		//translate the rotation segment to the origin
		x[i] -= x[start];
		y[i] -= y[start];
		z[i] -= z[start];
		
		//store old vector for the rotation matrix below
		ttx=x[i];
		tty=y[i];
		ttz=z[i];
		
		//Perform the rotation on the selected segment
		x[i] = (u*rx*rx + c)*ttx     + (u*ry*rx - s*rz)*tty  + (u*rz*rx + ry*s)*ttz;
		y[i] = (u*rx*ry + rz*s)*ttx  + (u*ry*ry + c)*tty     + (u*rz*ry - rx*s)*ttz;
		z[i] = (u*rx*rz - ry*s)*ttx  + (u*ry*rz + rx*s)*tty  + (u*rz*rz + c)*ttz;	
	
		//translate the rotation segment back to the original reference frame
		x[i] += x[start];
		y[i] += y[start];
		z[i] += z[start];
	
		//Final local nonbond energy
		for(j=end2;j!=start;j+=direction)
		{	
			distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
			Enblocalf = Enblocalf + nonbondpot(distr);
			//fprintf(stderr,"%d\t%d\n",i,j);
		};
	};

	//Calculate the final local bond energy (NOTE: this cannot go in the rotation loop)
/*	i=start;
	while( i != end1 )
	{
		i=i+direction;
		j=i-direction;
		distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
		Ebondf = Ebondf + bondpot(distr);	
		//fprintf(stderr,"%d\t%d\t%d\t%f\t%f\n",m,i,i-direction,distr,bondpot(distr));
	};
*/
	//Local bond and bond angle energy
	Ebondf = localE(m);
	
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);	
	//Ef=Etot();

	attemptpivot+=1;

	if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
	{
		acceptpivot+=1;
		if( (acceptpivot%1000) == 0 )  //this insures that errors do not accumulate over time
		{
			currEtot = Etot();
		};
	}
	else
    {
		//reject
		for(j=0;j<N;++j)
		{
			x[j]=tx[j];
			y[j]=ty[j];
			z[j]=tz[j];
		};
	};

}


//////////  Wang-Landau Crankshaft Move  //////////
//Selects two random monomers and then rotates the chain between them.  Rotation axis is the axis 
//between the two bonds, and rotation angle is random.
void wlcrankshaft(void)
{
	int i=0,j=0,m1=0,m2=0,start=0,end=0;
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double rx=0.0,ry=0.0,rz=0.0;
	double distr=0.0,rangle=0.0,c=0.0,s=0.0,u=0.0;
	double ttx=0.0,tty=0.0,ttz=0.0;
	
	//set temporary positions
	for(j=0;j<N;++j)
    {
		tx[j]=x[j];
		ty[j]=y[j];
		tz[j]=z[j];
    };

	//randomly choose two monomers
	m1=(int) (randd1()*N);
	m2=(int) (randd1()*N);
	
	//choose random rotation angle
	rangle=DC*2.0*M_PI*(2.0*randd1()-1.0);
	
	//below are the cases where the move wouldn't change anything
	while( (m1==m2) || (abs(m1-m2)==N-1) || (abs(m1-m2)==1) )  
	{
		m1=(int) (randd1()*N);
		m2=(int) (randd1()*N);
	};

	//set indices for energy calculation
	if( m1<m2 )
	{
		start=m1;
		end=m2;
	}
	else
	{
		start=m2;
		end=m1;
	};

	//Set initial energies
	Ei=currEtot;

//////////Rotation
	//rotation axis is the vector between m1 and m2
	distr = sqrt( (x[m1]-x[m2])*(x[m1]-x[m2]) + (y[m1]-y[m2])*(y[m1]-y[m2]) + (z[m1]-z[m2])*(z[m1]-z[m2]) );
	rx = (x[m1]-x[m2])/distr;
	ry = (y[m1]-y[m2])/distr;
	rz = (z[m1]-z[m2])/distr;
		
	//angle constants
	c=cos(rangle);
	s=sin(rangle);
	u=1.0-c;
	
	//Calculate the initial local bond energy
/*	for(i=start;i<end;i++)
	{
		distr = sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );
		Ebondi = Ebondi + bondpot(distr);
	};
*/
	//Calcs local bond and bond angle energies
	Ebondi = localE(m1) + localE(m2);
	
/////Rotation
	//now we perform the rotation
	for(i=start+1;i<end;i++)
	{
		//calculate initial local nonbonded energy
		for(j=0;j<N;j++)
			if( (j<start) || (j>end) )
			{
				distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
				Enblocali = Enblocali + nonbondpot(distr);
				//fprintf(stderr,"%d\t%d\t%20.16f\t%20.16f\t%20.16f\n",i,j,distr,nonbondpot(distr),Enblocali);
			};
		
		//translate the rotation segment to the origin
		x[i] -= x[start];
		y[i] -= y[start];
		z[i] -= z[start];
		
		//store old vector for the rotation matrix below
		ttx=x[i];
		tty=y[i];
		ttz=z[i];
		
		//Perform the rotation on the selected segment
		x[i] = (u*rx*rx + c)*ttx     + (u*ry*rx - s*rz)*tty  + (u*rz*rx + ry*s)*ttz;
		y[i] = (u*rx*ry + rz*s)*ttx  + (u*ry*ry + c)*tty     + (u*rz*ry - rx*s)*ttz;
		z[i] = (u*rx*rz - ry*s)*ttx  + (u*ry*rz + rx*s)*tty  + (u*rz*rz + c)*ttz;	
	
		//translate the rotation segment back to the original reference frame
		x[i] += x[start];
		y[i] += y[start];
		z[i] += z[start];
	
		//calculate final local nonbonded energy
		for(j=0;j<N;j++)
			if( (j<start) || (j>end) )
			{
				distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
				Enblocalf = Enblocalf + nonbondpot(distr);
				//fprintf(stderr,"%d\t%d\t%20.16f\t%20.16f\t%20.16f\n",i,j,distr,nonbondpot(distr),Enblocalf);
			};
	};
	
	//Calculate the initial local bond energy
/*	for(i=start;i<end;i++)
	{
		distr = sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );
		Ebondf = Ebondf + bondpot(distr);
	};
*/	
	//Calcs local bond and bond angle energies
	Ebondf = localE(m1) + localE(m2);
	
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);

	attemptside+=1;

	if(WangLandau(Ei,Ef,1.0,JA,JA)==1)
	{
		acceptside+=1;
		if( (acceptside%1000) == 0 )
		{
			currEtot = Etot();
		};
	}
	else
    {
		//reject
		for(j=0;j<N;++j)
		{
			x[j]=tx[j];
			y[j]=ty[j];
			z[j]=tz[j];
		};
    };
}

//Cut and join move designed to rearrange the bonds at high densities 
void wlcutjoin(void)
{
	int i=0,j=0,m=0,istate=0,fstate=0,fJAS=0,neighfrwd=0,neighbkwd=0;
	int end=0;
	double dx=0.0,dy=0.0,dz=0.0;
	double dist=0.0;
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double Wfrwd=0.0,Wbkwd=0.0;
	
	//CONSTANT JA
	fJAS=(int) ((JA-WLD2min)*invdWLD2);
	
	//set temporary positions
	for(j=0;j<N;++j)
    {
		tx[j]=x[j];
		ty[j]=y[j];
		tz[j]=z[j];
    };

	//These statements decide which end to connect with
	if( randd1() < 0.5 )
		end = N-1;
	else
		end = 0;

	//randomly choose a monomer (excluding ends)
//	m=(int) (randd1()*(N-2) + 1);
	neighfrwd=0;
	for(i=0;i<N;i++)
	{
		Nlist[i]=-1;
		if(fabs(i-end)!=1 && fabs(i-end)!=0)
		{
			dx=x[i]-x[end];
			dy=y[i]-y[end];
			dz=z[i]-z[end];
			dist = sqrt(dx*dx+dy*dy+dz*dz);
			if( dist < Rcutbond )
			{
				Nlist[neighfrwd]=i;	
				neighfrwd++;
				//fprintf(stderr,"%d\t%d\t%g\t%d\n",i,end,dist,neighfrwd);
			};
		};
	};
	
	//Statement decides if it is worth making the move (if zero neighbors, don't do anything)
	if(neighfrwd > 0)
	{
		m = Nlist[ (int)(neighfrwd*randd1()) ];
		//fprintf(stderr,"%d\t%d\n",m,neighfrwd);
	
		//Generate a neighbor list for selecting which monomer to perform the move with
	
		//Set initial energies
		//Ei=Etot();
		Ei=currEtot;
		istate=(int) ((Ei*invN-WLD1min)*invdWLD1); 
		Ebondi=localE(m);
		//Local nonbond energy calculation
		dist=sqrt( (x[m]-x[end])*(x[m]-x[end]) + (y[m]-y[end])*(y[m]-y[end]) + (z[m]-z[end])*(z[m]-z[end]) );
		Enblocali = nonbondpot(dist);

		//Relabeling of bonds:  Switches labels in the direction of the newly bonded end
		if(end == N-1)  //Bonds toward N-1
		{
			j=m;
			for(i=N-1;i>m;i--)
			{
				j++;  //NOTE: list increases for bonds to the right (they go towards N-1 end)
				x[i]=tx[j];
				y[i]=ty[j];
				z[i]=tz[j];
				//fprintf(stderr,"%d\t%d\t%d\t%g\n",mi,i,j,sqrt(dist));
			}
		}
		else if(end == 0)  //Bonds toward 0
		{
			j=m;
			for(i=0;i<m;i++)
			{
				j--;  //NOTE: list decreases for bonds to the left (they go toward 0 end)
				x[i]=tx[j];
				y[i]=ty[j];
				z[i]=tz[j];
				//fprintf(stderr,"%d\t%d\t%d\t%g\n",mi,i,j,sqrt(dist));
			}
		}

		//Calculate the final energy
		//Ef=Etot();
		Ebondf=localE(m);
		//Final local nonbond energy calculation
		dist=sqrt( (x[m]-x[end])*(x[m]-x[end]) + (y[m]-y[end])*(y[m]-y[end]) + (z[m]-z[end])*(z[m]-z[end]) );
		Enblocalf = nonbondpot(dist);
		Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
		fstate=(int) ((Ef*invN-WLD1min)*invdWLD1);
		//fprintf(stderr,"%g\t%g\t%g\t%g\n",Ei,Ef,Enonbondf,Ebondf);
	
		neighbkwd=0;
		for(i=0;i<N;i++)
		{
			Nlist[i]=-1;
			if(fabs(i-end)!=1 && fabs(i-end)!=0)
			{
				dx=x[i]-x[end];
				dy=y[i]-y[end];
				dz=z[i]-z[end];
				dist = sqrt(dx*dx+dy*dy+dz*dz);
				if( dist < Rcutbond )
				{
					Nlist[neighbkwd]=i;	
					neighbkwd++;
					//fprintf(stderr,"%d\t%d\t%g\t%d\n",i,end,dist,neighbkwd);
				};
			};
		};
	
		Wfrwd = (double) neighfrwd;
		Wbkwd = (double) neighbkwd;
		//fprintf(stderr,"%d\t%d\t%g\n",neighfrwd,neighbkwd,Wbkwd/Wfrwd);
		
	}
	else //if zero neighbors, don't do anything
	{
		Wfrwd = 1.0;
		Wbkwd = 1.0;
		//fprintf(stderr,"%d\t%d\t%g\n",neighfrwd,neighbkwd,Wbkwd/Wfrwd);
	
		Ei = currEtot;
		Ef = currEtot;
	}
		
	//cutjoinattempts[istate][fJAS] += 1.0;
	//attemptcutjoin+=1;
  
	if(WangLandau(Ei,Ef,Wfrwd/Wbkwd,JA,JA)==1)
	{
		//acceptcutjoin+=1;
		//cutjoinrate[istate][fJAS] += 1.0;
	}
	else
	{
		//reject
		for(j=0;j<N;++j)
		{
			x[j]=tx[j];
			y[j]=ty[j];
			z[j]=tz[j];
		};
	};	
	//fprintf(stderr,"%g\t%g\t%d\t%d\n",Ei/(1.0*N),Ef/(1.0*N),attemptcutjoin,acceptcutjoin);
		
	
}



