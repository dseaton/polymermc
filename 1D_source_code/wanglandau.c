#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rand.h"
#include "model.h"
#include "metrop.h"
#include "wanglandau.h"

//NEW VERSION for production run
//uses an intermediate (10th) iteration to generate the average order parameter values
void production_run(void)
{
	int i,k;
	double Hmin=0.0,Hmax=0.0,tmp=0.0;
	FILE *ifp1,*ifp2,*ifp3,*ofp;
	char s1[512],s3[512];
	
	//initialize the production run histograms
	init_production_run();

	//Reading in the final iteration number from the restart function
	ifp1=fopen("DOS_H_iter012.dat","r");
	if( ifp1 == NULL )
	{
		fprintf(stderr,"Error in production_run(), the file could not be opened!\n\n"); 
		exit(1);
	};
	fscanf(ifp1,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&lnwlf,&TotalSweeps,&IterSweeps,&tmp);
	//fprintf(stderr,"#%d\t%g\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL());
    fclose(ifp1);
	
	//Reading in the current confinguration for the restart
	ifp2=fopen("restart.XYZ","r");
	//fscanf(ifp2,"#\t%lg\t%lg\n",&tmp1,&tmp2);
	for(i=0;i<N;++i)
    {
		fscanf(ifp2,"%lg\t%lg\t%lg\n",&(x[i]),&(y[i]),&(z[i]));    
		//fprintf(stderr,"%g\t%g\t%g\n",x[i],y[i],z[i]); 
    };
	fclose(ifp2);

	//Calculate the energy of the restart configuration
	currEtot=Etot();
	currnonbondE=nonbondEnergy();
	
	//Introduce a metropolis routine to ensure the configuration is not stuck
	//Seq(2.5,0,1000,10);	

	//currEtot=Etot();
	//currnonbondE=nonbondEnergy();

	//read in the DOS where initial simulation left off
	readg();  //reads the 'numf' DOS and H
	resetWL();  //resets the histogram
	
	//number of hybrid moves - N diff, one of each other
	IterSweeps=0;
	
	//the moves below are the same as those in wlhybrid
	while(Hmin < ProductionBinSamps)
	{
		wlhybrid();
		IterSweeps+=1;
		//This loop finds the flatness and the number of unsampled bins
		
		if( (IterSweeps % 20000) == 0 )
		{
			//Initialize the max and min variables
			Hmin=1.0e300;
			Hmax=-1.0e300;
			for(i=0;i<D1BINS;++i)
			{			
				//minimum of the histogram
				if( wlH[i]<Hmin )
				{
					Hmin=wlH[i];
				};
				//maximum of the histogram
				if( wlH[i]>Hmax )
				{
					Hmax=wlH[i];
				};
			};	
	
			//Writing out the production run histogram
			ofp=fopen("prodrun_H.dat","w");

			//Label each iteration with the mod. factor, number of sweeps, and flatness
			fprintf(ofp,"#Production Run: %d\t%g\t%d\t%d\t%g\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL(),Hmin);
			for(i=0;i<D1BINS;++i)
			{	
				//Printing Out 1D Data
				fprintf(ofp,"%g\t%18.10e\n",i/invdWLD1+WLD1min,wlH[i]);
			};
			fflush(ofp);
			fclose(ofp);
			//fprintf(stderr,"%g\t%g\t%g\n",Hmin,Hmax,wlH[0]);
		
			//writeEERgyr();		
		};
		
		//fprintf(stderr,"%g\t%g\t%d\n",Hmin,Hmax,ProductionBinSamps);

	}

	//Now we divide the order parameter histograms by the number of measurements 
	for(i=0;i<D1BINS;++i)
		if(wlH[i] > 0.0)
		{
			HRg[i] = HRg[i]/(1.0*wlH[i]);
			HEEdist[i] = HEEdist[i]/(1.0*wlH[i]);
			Hcore[i] = Hcore[i]/(1.0*wlH[i]);
			for(k=0;k<GRBINS;k++)
				gr[i][k] = gr[i][k]/(1.0*wlH[i]);
			//fprintf(stderr,"%g\t%g\n",i/invdWLD1+WLD1min,HRg[i][0]);
		}
		else
		{
			HRg[i] = 0.0;
			HEEdist[i] = 0.0;
			Hcore[i] = 0.0;
			for(k=0;k<GRBINS;k++)
				gr[i][k] = gr[i][k]/(1.0*wlH[i]);
		};

	//Keeps track of the total number of sweeps - Original + Production
	TotalSweeps+=IterSweeps;

	//Writing out the production run histogram
	ofp=fopen("prodrun_H.dat","w");

	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#Production Run: %d\t%g\t%d\t%d\t%g\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL(),Hmin);
	for(i=0;i<D1BINS;++i)
		//Printing Out 1D Data
		fprintf(ofp,"%g\t%18.10e\n",i/invdWLD1+WLD1min,wlH[i]);
	  
	fflush(ofp);
	fclose(ofp);

	//Read the final DOS (since restart.DOS_H has the final numf (iteration) we can read it in from there)
	//this allows us to not worry about what the final DOS was, or that all runs have the same final DOS
	//numf=20;
	
	ifp3=fopen("restart.DOS_H","r");
	fscanf(ifp3,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&tmp,&tmp,&tmp,&tmp);
	//fprintf(stderr,"%03d\n",numf);
	readg();  //Note:  this erases the production run histogram, but it has already been written out
       
	//Print out the results
	orderqs();
	writeEERgyr();
}

//bins the order parameters into their arrays - accepts the new bin from the Wang-Landau function
void histfill_prun(int pbin)
{
	//pbin is the energy index
	HEEdist[pbin]+=EEdist();
	HRg[pbin]+=Rgyr2();
	Hcore[pbin]+=coredensity();	  
	grcalculate(pbin);
}


void init_production_run(void)
{
	int i,k;
	
	//Binning for pair distribution function g(r)
	dGR = 0.01;
	GRmax = 3.0;
	GRmin = 0.0;
	GRBINS = (int)((GRmax - GRmin)/(dGR));
	
	//Dynamic Memory Allocation - this has to be done because an input file is used.
	//WL 2D Arrays - Histogram, Density of States, Mask, and Rawmask
	//allocate storage for an array of pointers
	HRg = malloc( D1BINS * sizeof(double ) );
	HEEdist = malloc( D1BINS * sizeof(double ) );
	Hcore = malloc( D1BINS * sizeof(double ) ); 
	gr = malloc( D1BINS * sizeof(double *) );
	//Check to see if memory was allocated properly
	if ( (HRg == NULL) || (HEEdist == NULL) || (Hcore == NULL) )
    {
        fprintf(stderr,"\nFailure to allocate memory for 'Wang-Landau 2D Arrays'.  See 'initWL()'.\n");
        exit(1);
    };
	    
	//Initialize Arrays
	for(i=0;i<D1BINS;++i)
	{
		HRg[i]=0.0;
		HEEdist[i]=0.0;
		Hcore[i]=0.0;
	};

	for(i=0;i<D1BINS;i++)
	{
		gr[i] = malloc( GRBINS * sizeof(double) );
	}

	//Initialize Arrays
	for(i=0;i<D1BINS;++i)
		for(k=0;k<GRBINS;k++)
			gr[i][k] = 0.0;

}

//returns 1) wlH_min/wlH_avg and 2) the number of states below the flatness criteria
double flatWL(void)
{
	int i,k;
	double avg,min;
  
	k=0;  //number of sampled bins, used to calculate average
	min=1.0e300;  //the minimum sampled bin in the histogram wlH[][]
	avg=0.0;  //average height of the histogram wlH[][]
	numbelow_flat=0.0;  //number of bins below the flatness criteria
  
	//This loop finds the flatness and the number of unsampled bins
	for(i=0;i<D1BINS;++i)
		if(mask[i]==1)
		{			
			//minimum of the histogram
			if( (wlH[i]<min) )
			{
				min=wlH[i];
			};
			//average of the histogram
			if( (wlH[i]>0.0) )
			{
				avg+=wlH[i];
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
    if(mask[i]==1)
    {
    if( wlH[i]<Flatness*avg )
    numbelow_flat+=1.0;
    };
    //Stores the percentage of states below the flatness criteria	
    numbelow_flat = (numbelow_flat)/(1.0*D1BINS*D2BINS);
  */
	numbelow_flat = 0.0;  //Turned off for now
	//returns the current flatenss of the histogram		
	return min/avg;
  
}

//reset accumulated histogram
void resetWL()
{
	int i;
  
	for(i=0;i<D1BINS;++i)
		wlH[i]=0.0;
  
}

//reads in a mask.dat file
void readmask(void)
{
  int i,k,tmp;
  FILE *ifp;
  
  ifp=fopen("mask.dat","r");
  
  for(k=0;k<D1BINS;++k)
    {
      fscanf(ifp,"%d\t%d\n",&i,&tmp);
      mask[i]=tmp;
    };
}

//reads in a mask.dat file
void writemask(void)
{
  int i,k,tmp;
  FILE *ofp;
  
  ofp=fopen("mask.dat","w");
  
  for(i=0;i<D1BINS;++i)
      fprintf(ofp,"%g\t%d\n",i/invdWLD1+WLD1min,mask[i]);
}

//reads in the g.dat file
void readg(void)
{
	int i;
	FILE *ifp;
	char s1[512];
	double tmp;
  
	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	ifp=fopen(s1,"r");
	//Reading in the DOS and histogram from the restart function
	fscanf(ifp,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&lnwlf,&TotalSweeps,&IterSweeps,&tmp);
	//fprintf(stderr,"#%d\t%lg\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,tmp);
	for(i=0;i<D1BINS;++i)
    {
		fscanf(ifp,"%lg\t%lg\t%lg\n",&tmp,&(wllng[i]),&(wlH[i]));
		//fprintf(stderr,"%g\t%18.10e\t%g\n",tmp,wllng[i],wlH[i]);
    };
	fclose(ifp);
  
}

//writes out the DOS and the histogram
void write_DOS_H(void)
{
	int i;
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
	{
		if((wllng[i]>maxg) && (mask[i]==1))
			maxg=wllng[i];
		if((wllng[i]<ming) && (mask[i]==1))
			ming=wllng[i];
	};
  
	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#%d\t%g\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL());
	for(i=0;i<D1BINS;++i)
	{
		//Printing Out 1D Data
		//fprintf(ofp,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i]);
		fprintf(ofp,"%g\t%18.10e\t%18.10e\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i]);
	};
  
	fflush(ofp);
	fclose(ofp);
}

//Write out the normalized DOS //Currently only works for smaller chain sizes, this is because of underflow
void write_normDOS(void)
{
	int i;
	double maxg=-1.0e300;
	double area;
	FILE *ofp;
	char s1[512];
  
	sprintf(s1,"gnorm.dat");
	ofp=fopen(s1,"w");
  
	//Find the maximum of the DOS
	for(i=0;i<D1BINS;++i)
	{
		if((wllng[i]>maxg) && (mask[i]==1))
			maxg=wllng[i];
	};
  
	area = 0.0;
	
	//Taking the area under the curve (DOS)
	for(i=0;i<D1BINS-1;++i)
    {
		area = area + 0.5*(exp(wllng[i]-maxg) + exp(wllng[i+1]-maxg))*dWLD1*invN;	
	};
  
  
	for(i=0;i<D1BINS;++i)
    {
		//Printing Out 1D Data
		fprintf(ofp,"%g\t%g\n",i/invdWLD1+WLD1min,exp(wllng[i]-maxg)/area);
		//fprintf(stderr,"%g\t%g\t%g\n",wllng[i],area,exp(wllng[i]-maxg)/area);
    };
	fflush(ofp);
	fclose(ofp);
  
}

//writes out the accumlated histogram to H.dat
void writeEERgyr(void)
{
	int i,k;
	FILE *ofp1,*ofp2,*ofp3,*ofp4;
	char s1[512],s2[512],s3[512],s4[512];
  
	sprintf(s1,"Rgyr2_E.dat");
	ofp1=fopen(s1,"w");
    
	sprintf(s2,"EEdist_E.dat");
	ofp2=fopen(s2,"w");
  
	sprintf(s3,"CoreD_E.dat");
	ofp3=fopen(s3,"w");
	
	sprintf(s4,"gR_E.dat");
	ofp4=fopen(s4,"w");
    
	for(i=0;i<D1BINS;++i)
    {
		//Printing Out 1D Data
		
		fprintf(ofp1,"%g\t%g\n",i/invdWLD1+WLD1min,HRg[i]);
		fprintf(ofp2,"%g\t%g\n",i/invdWLD1+WLD1min,HEEdist[i]);
		fprintf(ofp3,"%g\t%g\n",i/invdWLD1+WLD1min,Hcore[i]);
	
		for(k=0;k<GRBINS;k++)
		{
			fprintf(ofp4,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,k*dGR,gr[i][k]);
		}
		fprintf(ofp4,"\n");
	};
    
	fflush(ofp1);
	fflush(ofp2);
	fclose(ofp3);
	fclose(ofp4);
}

//Function that writes out appropriate files used in the restart function "read_restart"
void write_restart(void)
{
	int i;
	FILE *ofp1,*ofp2,*ofp3;
	char s1[512];
  
	ofp1=fopen("restart.DOS_H","w");
  
	//Writing the DOS and histogram for the restart function
	fprintf(ofp1,"#%d\t%18.10e\t%d\t%d\t%18.10e\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL());
	for(i=0;i<D1BINS;++i)
    {
		fprintf(ofp1,"%g\t%18.10e\t%18.10e\n",i/invdWLD1+WLD1min,wllng[i],wlH[i]);
    };
	fflush(ofp1);
	fclose(ofp1);
  
	//Writing out the current confinguration for the restart
	ofp2=fopen("restart.XYZ","w");
	//fprintf(ofp2,"#\t%g\t%g\n",Etot(),fastnonbondEnergy());
	for(i=0;i<N;++i)
    {
		fprintf(ofp2,"%18.10e\t%18.10e\t%18.10e\n",x[i],y[i],z[i]);
    };
	fflush(ofp2);
	fclose(ofp2);
  
	//If in the middle of a production run, write out the order parameter histograms
	if( ProductionRun == 1 )
    {	
		ofp3=fopen("restart.ProdRun","w");
		for(i=0;i<D1BINS;++i)
		{
			fprintf(ofp3,"%18.10e\t%18.10e\t%18.10e\t%18.10e\n",i/invdWLD1+WLD1min,HRg[i],HEEdist[i],Hcore[i]);
		};	
		fflush(ofp3);
		fclose(ofp3);
    };
  
}

//Restart 'reading' function that allows simulations to be started from any given checkpoint
void read_restart(void)
{
	int i;
	double tmp=0.0,tmp1=0.0,tmp2=0.0;
	FILE *ifp1,*ifp2,*ifp3;
	ifp1=fopen("restart.DOS_H","r");

	//Reading in the DOS and histogram from the restart function
	fscanf(ifp1,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&lnwlf,&TotalSweeps,&IterSweeps,&tmp);
	//fprintf(stderr,"#%d\t%18.10e\t%d\t%d\t%18.10e\n",numf,lnwlf,TotalSweeps,IterSweeps,tmp);
	for(i=0;i<D1BINS;++i)
    {
		fscanf(ifp1,"%lg\t%lg\t%lg\n",&tmp,&(wllng[i]),&(wlH[i]));
		//fprintf(stderr,"%g\t%18.10e\t%18.10e\n",tmp,wllng[i],wlH[i]);
    };
	fclose(ifp1);
  
	//if(flatWL() > Flatness)
		//lnwlf = lnwlf/IterationFactor;
 
	//Reading in the current confinguration for the restart
	ifp2=fopen("restart.XYZ","r");
	//fscanf(ifp2,"#\t%lg\t%lg\n",&tmp1,&tmp2);
	for(i=0;i<N;++i)
    {
		fscanf(ifp2,"%lg\t%lg\t%lg\n",&(x[i]),&(y[i]),&(z[i]));    
		//fprintf(stderr,"%18.10e\t%18.10e\t%18.10e\n",x[i],y[i],z[i]); 
    };
	fclose(ifp2);

	//Calculate the energy of the restart configuration
	currEtot=Etot();
	currnonbondE=nonbondEnergy();
	//fprintf(stderr,"%g\t%d\t%f\t%f\t%f\t%f\n",lnwlf,numf,currEtot*invN,tmp1*invN,currnonbondE*invN,tmp2*invN);

	if( (currEtot < WLD1min*N) || (currEtot > WLD1max*N) )
	{
		fprintf(stderr,"Error: Total energy - currEtot=%f - is out of bounds!\n",currEtot);
		exit(1);
	};
 
	//If in the middle of a production run, read in the order parameter histograms
	if( ProductionRun == 1 )
    {	
		ifp3=fopen("restart.ProdRun","r");
		for(i=0;i<D1BINS;++i)
		{
			fscanf(ifp3,"%lg\t%lg\t%lg\t%lg\n",&tmp,&(HRg[i]),&(HEEdist[i]),&(Hcore[i]));
	    };
	
		fclose(ifp3);
    };
  
}

//routine to initialize things for the Wang-Landau simulation
void initWL(void)
{
	int i,k;

	//make sure that the polymer is initialized
	initialize();
	numbelow_flat = 0.0;  //Gives the percentage of states below the flatness criteria ( see function flatWL() )
	TotalSweeps=0;		//Total number of sweeps
	IterSweeps=0;			//Number of sweeps per iteration
	ProductionRun=0;		//Initially the production run is turned off
  
	//Primary Binning Direction for WL Simulation
	D1BINS = (int)((WLD1max - WLD1min)/(dWLD1*invN));
		  
	//Inverse Bin Width
	invdWLD1=1.0/(dWLD1*invN);
    
	//this value times Ebond-WLEbondmin, gives the bin number;
	//WLinvdEbond=(1.0*Lbond)/(WLEbondmax-WLEbondmin);
	//WLinvdEnonbond=(1.0*Lnonbond)/(WLEnonbondmax-WLEnonbondmin);

	//Dynamic Memory Allocation - this has to be done because an input file is used.
	//WL 2D Arrays - Histogram, Density of States, Mask, and Rawmask
	//allocate storage for an array of pointers
	wlH = malloc( D1BINS * sizeof(double ) );
	wllng = malloc( D1BINS * sizeof(double ) );
	mask = malloc( D1BINS * sizeof(int ) );
	//Check to see if memory was allocated properly
	if ( (wlH == NULL) || (wllng == NULL) || (mask == NULL) )
    {
        fprintf(stderr,"\nFailure to allocate memory for 'Wang-Landau 2D Arrays'.  See 'initWL()'.\n");
        exit(1);
    };

	//Run the standard MC routine until the configuration has energy within the WL simulation energy range
	T=0.5;
	while( (currEtot*invN>WLD1max) || (currEtot*invN<WLD1min) )
    {
		fprintf(stderr,"Relaxing:  %g of %g \n",currEtot*invN,WLD1max);
		mchybrid();
    };

	LOWESTE=currEtot*invN;
  
	//Initialize Arrays
	for(i=0;i<D1BINS;++i)
	{
		wllng[i]=0.0;
		wlH[i]=0.0;
		mask[i]=1;
		//fprintf(stderr,"%d\t%g\t%d\n",i,wllng[i],mask[i]);
	};
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
		wldiff();
  
    wlreptation();
	wlcutjoin();
	//wlcrankshaft();
	//wlrandpivot();
	
}


//routine to encapsulate the WangLandau algorithm, returns 1 if accepted and 0 if rejected
int WangLandau(double Ei, double Ef, double Enbi, double Enbf)
{
	int iti;//indices of initial config
	int fti;//indices of final config

	//short hand for energies
	//initial
	double D1Valuei=Ei*invN;
	//double D1Valuei=(Ei-Efasti)*invN;
	double D2Valuei=Enbi*invN;
	//final
	double D1Valuef=Ef*invN;
	//double D1Value=(Ef-Efastf)*invN;
	double D2Valuef=Enbf*invN;

	double lngi,lngf;  //values for the density of states, intial and final
	double R;

	//Primary direction index
	iti=(int) ((D1Valuei-WLD1min)*invdWLD1);
	fti=(int) ((D1Valuef-WLD1min)*invdWLD1);

	//fprintf(stderr,"%d\t%d\t%g\t%g\t%d\n",iti,fti,D1Valuei,D1Valuef,D1BINS);

	//fprintf(stderr,"%d\t%d\t%g\t%g\t%g\n",iti,fti,D1Valuei,D1Valuef,Etot());
	//fprintf(stderr,"%20.15f\t%20.15f\t%20.15f\t%20.15f\n",Ei*invN,Efasti*invN,Ef*invN,Efastf*invN);
  
	//This statement simply prints out the lowest energy configuration
	if(Ef*invN<LOWESTE)
	{
		writeinc(0,Ef*invN);
		write_mol2(0);
		LOWESTE=Ef*invN;
		//fprintf(stderr,"New Lowest E config %g\n\n",LOWESTE);
	};

	//  fprintf(stderr,"%d,%d >> %d,%d\n",iti,itj,fti,ftj);

	if((fti<0)||(fti>=D1BINS)||(mask[fti]==0))
    {
		//reject, outside of bounds
		//reject
		if(ProductionRun == 0)
		{
			wlH[iti]+=1.0;
			wllng[iti]+=lnwlf;
		}	  
		else if(ProductionRun == 1)
		{
			wlH[iti]+=1.0;
			histfill_prun(iti);
			wllng[iti]+=lnwlf;
		};
		return 0;
    }
	else
    {
		//inside of bounds
		lngi=wllng[iti];
		lngf=wllng[fti];
      
		R=exp(lngi-lngf);

		if(randd1()<R)
		{
			//accept
			currEtot=Ef;
			currnonbondE=Enbf;

			if(ProductionRun == 0)
			{	
				wlH[fti]+=1.0;
				wllng[fti]+=lnwlf;
			}
			else if(ProductionRun == 1)
			{
				wlH[fti]+=1.0;
				histfill_prun(fti);
				wllng[fti]+=lnwlf;
			};
		
			return 1;
		}
		else
		{
			//reject
			if(ProductionRun == 0)
			{
				wlH[iti]+=1.0;
				wllng[iti]+=lnwlf;
			}
			else if(ProductionRun == 1)
			{
				wlH[iti]+=1.0;
				histfill_prun(iti);
				wllng[iti]+=lnwlf;
			};
				  
			return 0;
		};
    };

}

//////////  Wang-Landau Diffusion Move  //////////
//Displaces a random monomer by a small amount 
void wldiff(void)
{
	int i=0;
	double ox=0.0,oy=0.0,oz=0.0; //old positions of the monomer
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double Enonbondi=0.0,Enonbondf=0.0; //overall non-bonded energy variables

	//randomly choose a monomer
	i=(int) (randd1()*N);
	
	//Set initial energies
	//Ei=Etot();
	Ei=currEtot;
	Ebondi=localE(i);
	Enonbondi=currnonbondE;
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
	Enonbondf=Enonbondi + (Enblocalf - Enblocali);
	Ef = Ei + (Ebondf - Ebondi) + (Enonbondf - Enonbondi);
		
    attemptdiff+=1;
  
	if(WangLandau(Ei,Ef,Enonbondi,Enonbondf)==1)
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
	double Enonbondi=0.0,Enonbondf=0.0; //overall non-bonded energy variables
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
	Enonbondi=currnonbondE;
	
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
	//Enonbondf=nonbondEnergy();
	Enonbondf=Enonbondi + (Enblocalf - Enblocali);
	Ef = Ei + (Ebondf - Ebondi) + (Enonbondf - Enonbondi);
	
    attemptsnake+=1;
	
	if(WangLandau(Ei,Ef,Enonbondi,Enonbondf)==1)
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
	double Enonbondi=0.0,Enonbondf=0.0; //overall non-bonded energy variables
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
	Enonbondi=currnonbondE;
	
	//Calculate the initial local bond energy (NOTE: this cannot go in the rotation loop)
	i=start;
	while( i != end1 )
	{
		i=i+direction;
		j=i-direction;
		distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
		Ebondi = Ebondi + bondpot(distr);	
		//fprintf(stderr,"%d\t%d\t%d\t%f\t%f\n",m,i,i-direction,distr,bondpot(distr));
	};
	
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
	i=start;
	while( i != end1 )
	{
		i=i+direction;
		j=i-direction;
		distr = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
		Ebondf = Ebondf + bondpot(distr);	
		//fprintf(stderr,"%d\t%d\t%d\t%f\t%f\n",m,i,i-direction,distr,bondpot(distr));
	};

	Enonbondf = Enonbondi + (Enblocalf - Enblocali);
	//Enonbondf=nonbondEnergy();
	Ef = Ei + (Ebondf - Ebondi) + (Enonbondf - Enonbondi);	
	//Ef=Etot();

	attemptpivot+=1;

	if(WangLandau(Ei,Ef,Enonbondi,Enonbondf)==1)
	{
		acceptpivot+=1;
		if( (acceptpivot%1000) == 0 )  //this insures that errors do not accumulate over time
		{
			currEtot = Etot();
			currnonbondE = nonbondEnergy();
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
	double Enonbondi=0.0,Enonbondf=0.0; //overall non-bonded energy variables
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
	Enonbondi=currnonbondE;

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
	for(i=start;i<end;i++)
	{
		distr = sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );
		Ebondi = Ebondi + bondpot(distr);
	};
	
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
	for(i=start;i<end;i++)
	{
		distr = sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );
		Ebondf = Ebondf + bondpot(distr);
	};
		
	//Enonbondf=nonbondEnergy();
	Enonbondf = Enonbondi + (Enblocalf - Enblocali);
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);

	attemptside+=1;

	if(WangLandau(Ei,Ef,Enonbondi,Enonbondf)==1)
	{
		acceptside+=1;
		if( (acceptside%1000) == 0 )
		{
			currEtot = Etot();
			currnonbondE = nonbondEnergy();
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
	int i=0,j=0,m=0;
	int end=0;
	double dist=0.0;
    double Ei=0.0,Ef=0.0; //intial and final energies (before and after diffusion)
	double Ebondi=0.0,Ebondf=0.0,Enblocali=0.0,Enblocalf=0.0; //local energy variables
	double Enonbondi=0.0,Enonbondf=0.0; //overall non-bonded energy variables

	//set temporary positions
	for(j=0;j<N;++j)
    {
		tx[j]=x[j];
		ty[j]=y[j];
		tz[j]=z[j];
    };

	//randomly choos a monomer (excluding ends)
	m=(int) (randd1()*(N-2) + 1);
	
	//These if statements decide which end to connect with
	if( randd1() < 0.5 )
		end = N-1;
	else
		end = 0;

	//Set initial energies
	//Ei=Etot();
	Ei=currEtot;
	Ebondi=localE(m);
	Enonbondi=currnonbondE;
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
	Enonbondf=Enonbondi + (Enblocalf - Enblocali);
	Ef = Ei + (Ebondf - Ebondi) + (Enonbondf - Enonbondi);
	
	//fprintf(stderr,"%g\t%g\t%g\t%g\n",Ei,Ef,Enonbondf,Ebondf);
				
    attemptcutjoin+=1;
  
	if(WangLandau(Ei,Ef,Enonbondi,Enonbondf)==1)
    {
		acceptcutjoin+=1;
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



