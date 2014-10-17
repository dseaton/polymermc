/*
  Contains the energy functions and basic arrays and variables needed to describe the model
 
 This is code to simulate the folding and unfolding of a polymer chain
 of length N monomers.  There are three basic types of interactions:
 
 1) a local harmonic potential for the bond lengths
 E_bond_length=LJpot+FENE, 
 where the preferred bond length is 1.0
 
 2) a local harmonic potential for the bond angles
 E_bond_angle=JA*(cos(ang)-COSANG)^2, 
 where ang is the angle between a monomer and its two nn monomers
 and ANG is the ground state angle, ANG=PI, COSANG=cos(ANG)
 
 3) a non-local interaction between non-nn (non-bonded) monomers
 if(r<Rcutnonbond)
 return -1.0*Jnonbond*(2.0*(r-2.0)^3 - 3.0*(r-2.0)^2 + 1.0);
 else
 return 0.0;
 
 the interaction is normalized such that it goes to zero r>=Rcutnonbond and is equal to Jnonbond at r=2.0.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "model.h"
#include "rand.h"
#include "wanglandau.h"

//Calculates thermodynamic quantities and writes modification-factor labeled files
void thermoqs()
{
	int i;	
	double T,U,Z,C,F,S;
	double centerE,lambda,Bw,Nfree;
	
	FILE *therm_op;
	char s1[512];
	
	//Opening File containing all thermodynamic quantities (in following order)
	//	T	U	Cv	freeE  Entropy	Rgyr2  EEdist
	sprintf(s1,"therm.dat");
	therm_op=fopen(s1,"w");
	
	if( (therm_op==NULL) )
	{
		fprintf(stderr, "\nHey, this file ( in thermoqs() ) could not be opened!\n\n");
		exit(1);
	};
	
	//Normalization for free energy
	Nfree = 0.0;
	
	//Main Temperature Loop
	for(T=TTi;T<=TTf+dTT;T=T+dTT)
	{			
		U = 0.0;	// Initializing the average energy <E>	
		Z = 0.0;	// Initializing the Partition Function 		
		C = 0.0;	// Initialize the Specific Heat
		F = 0.0;	// Initialize the free energy
		S = 0.0;	// Initialize the entropy
		Bw = 0.0;	//Boltzmann weight
		lambda = -1.0e300;	//Normalization shift (max value of DOS considering T)
		centerE = 0.0;	// Taking the center of the energy bin
		
		//Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
		for(i=0;i<D1BINS;i++)
		{
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
				
			if( (lambda < ((wllng[i]) - 1.0*(centerE)/T)) && (wllng[i] > 0.0) ) 
				lambda = ((wllng[i]) - 1.0*(centerE)/T);  
			
		};
		
		//Central Loop for calculating thermodynamic properties from the DOS
		for(i=0;i<D1BINS;i++)
		{			
			//Taking the center of the bin
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
			
			//Boltzmann Factor
			Bw =  exp( wllng[i]  - (centerE)/T - lambda );
		
			//Partition Function
			Z = Z + Bw;
				
			//Average Energy <E>
			U = U + (centerE)*Bw;		   		   
				
			//Average Energy Squared <E^2>
			C = C + (centerE)*(centerE)*Bw;								
		};
		
		//Internal Energy
		U = U/Z;	
		//Normalization of free energy
		if(T == TTi)
		{
			Nfree = -T*( lambda + log(Z) ) - U;
		};
		//Specific Heat
		C = ( (C/Z) - (U*U) ) / (T*T);	
		//Free Energy
		F = -T*( lambda + log(Z) ) - Nfree;
		//Entropy
		S = (U - F)/T;
			
		fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\n",T,U*invN,C*invN,F*invN,S*invN);	
	};	
		
	fclose(therm_op);
}


//Calculates various order-parameters associated with the model (Radius of Gyration, End to End distance, etc...)
void orderqs()
{
	int i;	
	double T,Z,Rg,EEd,CoreD;
	double centerE,lambda,Bw;
	
	FILE *order_op;
	char s1[512];
	
	//Opening File for printing all order parameter quantities 
	sprintf(s1,"orderqs.dat");
	order_op=fopen(s1,"w");
	
	if( (order_op==NULL) )
	{
		fprintf(stderr, "\nHey, this file ( in orderqs() ) could not be opened!\n\n");
		exit(1);
	};
	
	//Main Temperature Loop
	for(T=TTi;T<=TTf;T=T+dTT)
	{			
		Z = 0.0;  //Initialize the partition function
		Rg = 0.0;  //Initialize the average radius of gyration
		EEd = 0.0;    //Initialize the average end to end distance
		CoreD = 0.0;	//Initialize the core density
		Bw = 0.0;	//Boltzmann weight
		lambda = -1.0e300;	//Normalization shift (max value of DOS considering T)
		centerE = 0.0;	// Taking the center of the energy bin
		
		//Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
		for(i=0;i<D1BINS;i++)
		{
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
				
			if( (lambda < ((wllng[i]) - 1.0*(centerE)/T)) ) 
				lambda = ((wllng[i]) - 1.0*(centerE)/T);  
		};
		
		//Central Loop for calculating thermodynamic properties from the DOS
		for(i=0;i<D1BINS;i++)
		{			
			//Taking the center of the bin
			centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
			
			//Boltzmann Factor
			Bw =  exp( wllng[i]  - (centerE)/T - lambda );
				
			//Partition Function
			Z = Z + Bw;
			
			if(wlH[i]>0.0)		
			{
				//Radius of Gyration
				Rg = Rg + HRg[i]*Bw;
					
				//End to End distance
				EEd = EEd + HEEdist[i]*Bw;
					
				//Bond angle distribution
				CoreD = CoreD + Hcore[i]*Bw;
			};
		
		};
		//Radius of Gyration
		Rg = Rg/Z;
		//End to End distance
		EEd = EEd/Z;
		//Core Density
		CoreD = CoreD/Z;
	
		fprintf(order_op,"%g\t%g\t%g\t%g\n",T,Rg,EEd,CoreD);
	};	
		
	fclose(order_op);

	//fprintf(stderr,"What do you think?\n\n");
}

//ROUTINES TO MANAGE THE ACCEPTANCE RATES FOR THE MOVES

//reset the attempt and acceptance counters to 0
void resetcounters(void)
{
      acceptpivot=0;
      acceptdiff=0;
      attemptpivot=0;
      attemptdiff=0;
      acceptside=0;
      attemptside=0;
      attemptsnake=0;
      acceptsnake=0;
}

void reportcounters(FILE *fp)
{
  fprintf(fp,"%g\t%g\t%g\t%g\t%g\n",acceptdiff/(1.0*attemptdiff),acceptcutjoin/(1.0*attemptcutjoin),acceptpivot/(1.0*attemptpivot),acceptside/(1.0*attemptside),acceptsnake/(1.0*attemptsnake));
}

//INTERACTION POTENTIALS, Note bond angle potential does not have it's own function


/*
Non-bonded potential.  Returns the value of the nonbonded potential
which is -Jnonbond for r=2 and 0 for r>=Rcutnonbond
*/

double nonbondpot(double r)
{
	double pvalue=0.0;

	if(r<=Rcutnonbond)
	{
		double tmp=1.0/r;
		double tmp2=tmp*tmp;
		double tmp6=tmp2*tmp2*tmp2;
		double tmp12=tmp6*tmp6;
	
		pvalue = Jnonbond*(tmp12 - 2.0*tmp6 - shiftLJ);
	};

	return pvalue;

/*  
  //Standard LJ Potential (truncated and shifted)
  if(NBINTERACTION == 1)
  {
  	double tmp=1.0/r;
	double tmp2=tmp*tmp;
	double tmp6=tmp2*tmp2*tmp2;
	double tmp12=tmp6*tmp6;
	
	//the below potential is a truncated-shifted LJ potential
	if(r<=Rcutnonbond)
		pvalue = Jnonbond*(tmp12 - 2.0*tmp6 - shiftLJ);
	else
		pvalue = 0.0;  
  };
*/

/*
  //Quasi LJ Potential
  if(NBINTERACTION == 2)
  {
	double tmp = (r-2.0);
	double tmp2=tmp*tmp;

	if(r<Rcutnonbond)
		pvalue = -1.0*Jnonbond*(2.0*tmp2*tmp - 3.0*tmp2 + 1.0);
	else
		pvalue = 0.0;
  };
	
  //Square-Well Potential
  if(NBINTERACTION == 3)
  {
	if( (r<=Rcutnonbond) && (r>=1.0) )
		pvalue = -1.0*Jnonbond;
	else if( r<1.0 )
		pvalue = 1.0e10;
	else
		pvalue = 0.0;
  };	

  return pvalue;	

*/
}


/*
Bonded potential.  Returns the value of the bonded potential, LJ+FENE
which is 0 for r=1 and should be infinity for r>=Rcutbond
NOTE, returns 10^6 for r=Rcutbond

The LJ potential is Jbond*((l/r)^12-2(l/r)^6), where l=1.03412016183
is needed to shift the full potential such that the minimum is at r=1.0;
The shift value is included such that the potential is zero at r=1;
Rcutbond MUST BE 1.2, this is the cutoff for the FENE potential and
the l value WILL NOT BE CORRECT IF RCUTBOND!=1.2

l=(1/2+sqrt(23/44))^(1/6), when Rcutbond=1.2, in FENE, R0=Rcutbond
l=((1+sqrt(1+4*(-1/(12*(1-1/Rcutbond^2)))))/2)^(1/6)
l is found by letting r be 1.0, and then solving for the l factor in LJ+FENE
You get a quadratic equation that you can solve for the 6th root of l
*/
double bondpot(double r)
{
	if(r<1.2)
	{
		double tmp=sigmaFENE/r;
		//double tmp=1.0/r;
		double tmp2=tmp*tmp;
		double tmp6=tmp2*tmp2*tmp2;
		double tmp12=tmp6*tmp6;
	
		double fene_invR = (1.0/Rcutbond);
		double temp=r*fene_invR; // r^2/Rcutbond^2, r^2/(1.2*1.2)
		double temp2=temp*temp;
		double shift=0.09662249348037610375;
		
		return Jbond*(-0.72*log(1.0-temp2) + (tmp12 - 2.0*tmp6) + shift );
	}
	else
		return 1.0e6;
}

//ENERGY CALCULATION FUNCTIONS

//Calculates the total energy of the polymer chain
double Etot(void)
{
	int i,j;
	double dxl,dxr,dyl,dyr,dzl,dzr,ll,lr;
	double Eretval=0.0;
	double tmp;	
	double Eb=0.0;
	double Enb=0.0;
	double dist=0.0;
	
	for(i=0;i<N-1;i++)
	{
		//Calculate the bonded energy
		dist=sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );	
		Eb = Eb + bondpot(dist);
		
		for(j=i+2;j<N;j++)
		{
			if(i!=j)
			{
				//Calculate the nonbonded energy
				dist=sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
				Enb = Enb + nonbondpot(dist);
			};
		}
	};

	//Return the total energy
	return (Eb + Enb);

/*
  //calculate nn harmonic bond angle (3 body) local interaction
  if(JA > 0.0)
  {
	for(i=1;i<N-1;++i)//loop over all nn bonds
		{
			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i+1]-x[i]; 
			dyr=y[i+1]-y[i]; 
			dzr=z[i+1]-z[i]; 
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=((dxl*dxr+dyl*dyr+dzl*dzr)/(tmp)-COSANG);
			Eretval+=JA*tmp*tmp;
		}
  };
  return Eretval+fastnonbondEnergy();//+slownonbondEnergy();
*/
}

double bondEnergy(void)
{
	int i;
	double Eb=0.0;
	double dist=0.0;
	
	for(i=0;i<N-1;i++)
	{
		//Calculate the bonded energy
		dist=sqrt( (x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]) + (z[i]-z[i+1])*(z[i]-z[i+1]) );	
		Eb = Eb + bondpot(dist);
	};

	//Return the total energy
	return (Eb);

}

//returns the nonbonded energy, calculated using a fast and efficient algorithm
double fastnonbondEnergy(void)
{
  int i,j,k;
  double Esum=0.0;
  double r,r2;

  //nonbond interaction energy
  for(i=0;i<N-1;++i)
    {
      k=2;
      j=i+k;
      
      while(j<N)
	  {
		r2=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
		r=sqrt(r2);
	  
		if(r<=Rcutnonbond)
			Esum+=nonbondpot(r);
	  
		k=(int) (r-Rcutnonbond+1);
	  
	    //with this set to k=k-5, you should be OK up to T=30.0, tested for Jbond=JA=30.0, N=25
	    //smaller values than five (integer greater than or equal to 1) might run faster, but
		//don't change this unless you check first to make sure fast and slow nonbondEnergy give the same value
	    if(k>1)
			k=k-1;
	  
	    if(k<=0)
			k=1;
	  
	    j=j+k;
		//	  fprintf(stderr,"%d\t%d\t%d\n",i,j,k);
	  };
    };

  return Esum;
}

//returns the nonbonded energy, using a slow, but easily debugged algorithm
double nonbondEnergy(void)
{
	int i,j;
	double Enb=0.0;
	double dist=0.0;
	
	for(i=0;i<N-1;i++)
	{
		for(j=i+2;j<N;j++)
		{
			if(i!=j)
			{
				//Calculate the nonbonded energy
				dist=sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
				Enb = Enb + nonbondpot(dist);
				//fprintf(stderr,"%d\t%d\t%20.16f\t%20.16f\t%20.16f\n",i,j,dist,nonbondpot(dist),Enb);
			};
		}
	};
	
	//Return the total energy
	return (Enb);
}

//Calculates the total energy of the polymer chain (using 2N calculations)
double localnonbondEnergy(int m)
{
	int i;
	double local_Enb=0.0;
	double dist=0.0;

	for(i=0;i<N;i++)
	{
		if( (m!=i) && (m!=(i-1)) && (m!=(i+1)) )
		{
			//Calculate the nonbonded energy
			dist=sqrt( (x[m]-x[i])*(x[m]-x[i]) + (y[m]-y[i])*(y[m]-y[i]) + (z[m]-z[i])*(z[m]-z[i]) );
			local_Enb = local_Enb + nonbondpot(dist);
		};
	}

/*
	//Calculate the nonbonded energy - excluding monomer m and its neighbors - m+1, m-1
	//Originally this was one loop, but with this implementation, no if statement is needed - previously the if was checked for each i
	for(i=0;i<m-1;i++)
	{
		dist=sqrt( (x[m]-x[i])*(x[m]-x[i]) + (y[m]-y[i])*(y[m]-y[i]) + (z[m]-z[i])*(z[m]-z[i]) );
		local_Enb = local_Enb + nonbondpot(dist);
	}

	for(i=m+2;i<N;i++)
	{
		dist=sqrt( (x[m]-x[i])*(x[m]-x[i]) + (y[m]-y[i])*(y[m]-y[i]) + (z[m]-z[i])*(z[m]-z[i]) );
		local_Enb = local_Enb + nonbondpot(dist);
	}
*/
	//Return the total energy
	return (local_Enb);

}


//returns the bond length and bond angle energies envolving a change in position for monomer #ii
double localE(int m)
{
	double Esum=0.0;
	double dist=0.0;
	
	//calculates right neighbor bond (excludes (N-1)th monomer)
	if( m != (N-1) )
	{
		dist = sqrt( (x[m]-x[m+1])*(x[m]-x[m+1]) + (y[m]-y[m+1])*(y[m]-y[m+1]) + (z[m]-z[m+1])*(z[m]-z[m+1]));
		Esum+=bondpot(dist);
	};
	//calculates left neighbor bond (excludes 0th monomer)	
	if( m != 0 )  
	{
		dist = sqrt( (x[m]-x[m-1])*(x[m]-x[m-1]) + (y[m]-y[m-1])*(y[m]-y[m-1]) + (z[m]-z[m-1])*(z[m]-z[m-1]));
		Esum+=bondpot(dist);
	};

	return Esum;
/*
  //calculate nn harmonic bond angle (3 body) local interaction
  //potentially three bond angles that could have changed, <ii-2,ii-1,ii>,<ii-1,ii,ii+1>,and <ii,ii+1,ii+2>
  if(JA > 0.0)
  {
	//left angle
	if(ii-2>=0)
		{
			i=ii-1;

			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 
	
			//dr vector points from i to i+1, right
			dxr=x[i+1]-x[i]; 
			dyr=y[i+1]-y[i]; 
			dzr=z[i+1]-z[i]; 
      
			//length of left and right bond length vectors
			ll=sqrt(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=sqrt(dxr*dxr+dyr*dyr+dzr*dzr);

			tmp=((dxl*dxr+dyl*dyr+dzl*dzr)/(ll*lr)-COSANG);
			Eretval+=JA*tmp*tmp;
		}

		//center angle
		if((ii-1>=0)&&(ii+1<N))
		{
			i=ii;

			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i+1]-x[i]; 
			dyr=y[i+1]-y[i]; 
			dzr=z[i+1]-z[i]; 
      
			//length of left and right bond length vectors
			ll=sqrt(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=sqrt(dxr*dxr+dyr*dyr+dzr*dzr);

			tmp=((dxl*dxr+dyl*dyr+dzl*dzr)/(ll*lr)-COSANG);
			Eretval+=JA*tmp*tmp;
		}

	  //right angle
	  if(ii+2<N)
		{
			i=ii+1;
		
			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i+1]-x[i]; 
			dyr=y[i+1]-y[i]; 
			dzr=z[i+1]-z[i]; 
      
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);

			tmp=sqrt(ll*lr);
			tmp=((dxl*dxr+dyl*dyr+dzl*dzr)/(tmp)-COSANG);
			Eretval+=JA*tmp*tmp;
		}
	};
	return Eretval;
*/
}


//  OUTPUT AND GRAPHICS ROUTINES

//writes out a povray include file with the file number by frame
void writeinc(int frame, double energy)
{
  int i,j;
  char s[512];
  FILE *ofp;
  double cx=0.0,cy=0.0,cz=0.0;

  for(i=0;i<N;++i)
    {
      cx+=x[i];
      cy+=y[i];
      cz+=z[i];
    };

  cx*=invN;
  cy*=invN;
  cz*=invN;

  sprintf(s,"snap_%d.inc",frame);
  ofp=fopen(s,"w");
  
  fprintf(ofp,"//N=%d, Etot=%18.10e\n",N,currEtot*invN);
  fprintf(ofp,"#declare polychain=union{\n");

  for(i=0;i<N-1;++i)
    {
      fprintf(ofp,"sphere{<%g,%g,%g>, Rs texture{stex }}\n",x[i],y[i],z[i]);
      fprintf(ofp,"cylinder{<%g,%g,%g>,<%g,%g,%g>, Rc texture{ctex }}\n",
	      x[i],y[i],z[i],
	      x[i+1],y[i+1],z[i+1]);
    };
  fprintf(ofp,"sphere{<%g,%g,%g>, Rs texture{stex}}\n",x[N-1],y[N-1],z[N-1]);
  
  fprintf(ofp,"finish{ambient 0.25 diffuse 0.75 specular 0.3}\ntranslate -<%g,%g,%g>\n}\n",cx,cy,cz);
  //fprintf(ofp,"finish{ambient 0.25 diffuse 0.75 specular 0.3}\n}\n");
  fflush(ofp);
  fclose(ofp);
}

//similar to writeinc, but writes out a list of the x,y,z monomer positions
void writepositions(int frame)
{
  int i,j;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"positions_%d.dat",frame);
  ofp=fopen(s,"w");

  for(i=0;i<N;++i)
    {
      fprintf(ofp,"%g\t%g\t%g\n",x[i],y[i],z[i]);
    
    };
 
  fclose(ofp);
}

//similar to writeinc, but writes out a list of the x,y,z monomer positions
void write_mol2(int frame)
{
  int i,j;
  char s[512];
  FILE *ofp;
 
  sprintf(s,"snap_%d.mol2",frame);
  ofp=fopen(s,"w");

  fprintf(ofp,"@<TRIPOS>MOLECULE\n");
  fprintf(ofp,"POLYMER01\n");
  fprintf(ofp," %d %d\n",N,N-1);
  fprintf(ofp,"PROTEIN\n");
  fprintf(ofp,"NO_CHARGES\n");
  fprintf(ofp,"@<TRIPOS>ATOM\n");

  for(i=0;i<N;++i)
    {
      fprintf(ofp,"%d MON%d %g %g %g H\n",i+1,i+1,x[i],y[i],z[i]);
    
    };
  fprintf(ofp,"@<TRIPOS>BOND\n");
  
  for(i=0;i<N-1;++i)
    {
      fprintf(ofp,"%d %d %d %d\n",i+1,i+1,i+2,1);
    
    };

  fclose(ofp);
}


//INITIALIZATION and MEASUREMENT


//initialize the polymer chain with one end at the origin
//and all nn bonds at the ground state, in the x-y plane
void initialize(void)
{
	int i;
	double a=1.0;

	//for relatively big jumps
	D=0.05;  DC=1.0;  DD=0.1;

	//precalculations for chain length
	invN=1.0/(1.0*N);
	rootinvN = 1.0/sqrt(1.0*N);
  
	//term used to set the FENE+LJ potential's minimum at 1.0
	sigmaFENE = pow((0.5+sqrt(23.0/44.0)),(1.0/6.0));

	//precalculation for the LJ potential
	sigmaLJ=pow(0.5,1.0/6.0);
	
	//term used to shift the LJ potential such that it is zero at r=1.0
	shiftLJ = pow((1.0/Rcutnonbond),12.0) - 2.0*pow((1.0/Rcutnonbond),6.0);
  

	//precalculations for the nn bond angle potentials
	PI=2.0*acos(0.0);
	ANG=PI;
	COSANG=cos(ANG);
	 
	T=5.0;

	//Dynamic Memory Allocation - this has to be done because an input file is used.
	//x,y,z coordinates 
	x = malloc( N * sizeof(double) );
	y = malloc( N * sizeof(double) );
	z = malloc( N * sizeof(double) );
	if ( (x == NULL) || (y == NULL) || (z == NULL) )
    {
        fprintf(stderr,"\nFailure to allocate memory for 'x,y, or z'.  See 'initialize()'.\n");
        exit(1);
    };
  
	//tx,ty,tz - used as temporary coordinates in many of the "move" functions 
	tx = malloc( N * sizeof(double) );
	ty = malloc( N * sizeof(double) );
	tz = malloc( N * sizeof(double) );
	if ( (tx == NULL) || (ty == NULL) || (tz == NULL) )
    {
        fprintf(stderr,"\nFailure to allocate memory for 'tx,ty, or tz'.  See 'initialize()'.\n");
        exit(1);
    };

	//Initialize the tmp coordinates
	for(i=0;i++;i<N)
	{
		tx[i]=0.0;
		ty[i]=0.0;
		tz[i]=0.0;
	};
	
	//Initializing the chain
	//One end starts at the origin
	x[0]=0.0;
	y[0]=0.0;
	z[0]=0.0;
  
	//Initialize the rest of the chain
	for(i=1;i<N;++i)
	{
		x[i]=x[i-1]+1.0;
		y[i]=0.0;
		z[i]=0.0;
    };

	currEtot=Etot();
	currnonbondE=nonbondEnergy();
}


void reset_XYZ(void)
{
	int i;
	double cx=0.0,cy=0.0,cz=0.0;

	//Find the center of mass for each component
	for(i=0;i<N;i++)
    {
		cx+=x[i];
		cy+=y[i];
		cz+=z[i];
    };

	cx*=invN;
	cy*=invN;
	cz*=invN;

	//fprintf(stderr,"%g\t%g\t%g\n",cmass(),currEtot,currnonbondE);
	for(i=0;i<N;i++)
    {
		x[i]=x[i]-cx;
		y[i]=y[i]-cy;
		z[i]=z[i]-cz;
    };
	
	//Recalculate the overall energies
	currEtot=Etot();
	currnonbondE=nonbondEnergy();
	
	//fprintf(stderr,"%g\t%g\t%g\n",cmass(),currEtot,currnonbondE);
}

//calculates the center of mass (magnitude)
double cmass(void)
{
	int i;
	double cx=0.0,cy=0.0,cz=0.0;

	for(i=0;i<N;i++)
    {
		cx+=x[i];
		cy+=y[i];
		cz+=z[i];
    };

	cx*=invN;
	cy*=invN;
	cz*=invN;

	return sqrt( cx*cx + cy*cy + cz*cz );
}

//returns the radius of gyration
double Rgyr2(void)
{
  int i;
  double cx=0.0,cy=0.0,cz=0.0;
  double dx,dy,dz;
  double sum=0.0;

  for(i=0;i<N;++i)
    {
      cx+=x[i];
      cy+=y[i];
      cz+=z[i];
    };

  cx*=invN;
  cy*=invN;
  cz*=invN;

  for(i=0;i<N;++i)
    {
      dx=x[i]-cx;
      dy=y[i]-cy;
      dz=z[i]-cz;

      sum+=dx*dx+dy*dy+dz*dz;
    };

  //sum=sqrt(sum);
  sum=sum*invN;

  return sum;

}

//returns the end to end distance
double EEdist(void)
{
  int i;
  double dx,dy,dz;
  double dist=0.0;
  
  dx=x[0]-x[N-1];
  dy=y[0]-y[N-1];
  dz=z[0]-z[N-1];

  dist=dx*dx+dy*dy+dz*dz;

  dist=sqrt(dist);  

  return dist;
}


double coredensity(void)
{
	int i;
	double cx=0.0,cy=0.0,cz=0.0;
	double dx,dy,dz;
	double density=0.0;
  
	for(i=0;i<N;++i)
    {
		cx+=x[i];
		cy+=y[i];
		cz+=z[i];
    };

    cx*=invN;
    cy*=invN;
    cz*=invN;	
  
    for(i=0;i<N;i++)
	{
		dx=x[i]-cx;
		dy=y[i]-cy;
		dz=z[i]-cz;
		if( sqrt(dx*dx+dy*dy+dz*dz) < 3.0 )
			density+=1.0;	
	};
  
    return density;	
};

//Pair distribution function - calculation
void grcalculate(int pbin)
{
   	int i,r;
	int CM=-1;
	double cx=0.0,cy=0.0,cz=0.0;
	double dx,dy,dz;
	double dist,MINdist=1.0e300;
	double density=0.0;
  
	for(i=0;i<N;++i)
    	{
		cx+=x[i];
		cy+=y[i];
		cz+=z[i];
    	};

    	cx*=invN;
    	cy*=invN;
    	cz*=invN;	
  
    	for(i=0;i<N;i++)
	{
		dx=x[i]-cx;
		dy=y[i]-cy;
		dz=z[i]-cz;
		dist = sqrt(dx*dx+dy*dy+dz*dz);
		if( dist < MINdist )
		{
			MINdist = dist;	
			CM = i;
		}
	}

	for(i=0;i<N;i++)
	{
		if(i != CM)
		{
			dx=x[i]-x[CM];
			dy=y[i]-y[CM];
			dz=z[i]-z[CM];
			dist = sqrt(dx*dx+dy*dy+dz*dz);
		
			if(dist < GRmax)
			{
				r = (int)(dist/dGR);
				gr[pbin][r] = gr[pbin][r] + 1.0;
				//fprintf(stderr,"%d\t%d\t%g\t%g\n",i,r,gr[pbin][r]*dGR,dist);
			}
		}
	}

	//for(i=0;i<GRBINS;i++)
		//fprintf(stderr,"%g\t%g\n",i*dGR,gr[pbin][i]);
	
	//exit(1);

};


//ROTATION OPERATIONS


//pivots the polymer around a randomly chosen monomer by a small angle about a random axis
void randrotate(int ii, int leftright)
{
  int i,j;
  double xi1,xi2,xi;

  i=ii;
  
  //generate a random direction using the marsaglia (1973) method
  xi1=1.0-2.0*randd1();
  xi2=1.0-2.0*randd1();
  xi=xi1*xi1+xi2*xi2;
  
  while(xi>1.0)
    {
      xi1=1.0-2.0*randd1();
      xi2=1.0-2.0*randd1();
      xi=xi1*xi1+xi2*xi2;
    };

  rotationaxisx=2.0*xi1*sqrt(1.0-xi);
  rotationaxisy=2.0*xi2*sqrt(1.0-xi);
  rotationaxisz=1.0-2.0*xi;

  rotate(i,leftright,2.0*M_PI*D*(2.0*randd1()-1.0));  //rotates by a small angle

}


/*rotate the x,y,z coordinates of the chain about monomer ii on the left or right side, 
with the axis of rotation specified through the global variables rotatationaxisx, rotateaxisy, and rotateaxisz
(leftright==-1) or the right bond (leftright==1) of ii
angle is the rotation angle in radians

If you do not want to change the angle of a bond, then you should specify the axis in the direction of one of the bonds,
either to the left or right of ii

If you DO want to change the angle of the bond, but not the bond length,
then you should specify the axis of rotation perpindicular to the two bonds attached
to ii
*/
void rotate(int ii, int leftright, double angle)
{
  int i,j;
  double c,s,t,xdir,ydir,zdir;
  double transx,transy,transz,r;
  double ttx,tty,ttz;

  if((leftright!=-1)&&(leftright!=1))
  {
    fprintf(stderr,"ERROR:  leftright=%d for the pivot funcion\n",leftright);
    exit(1);
  };

  if((ii==0)||(ii==N-1))
  {
    //Nothing to do, rotation is beyond the edge of the chain
  }
  else
  {
    //translate such that x[ii]=y[ii]=z[ii]=0
    transx=x[ii];
    transy=y[ii];
    transz=z[ii];
    
    for(i=0;i<N;++i)
    {
      x[i]=x[i]-transx;
      y[i]=y[i]-transy;
      z[i]=z[i]-transz;
    };
    
    xdir=rotationaxisx;
    ydir=rotationaxisy;
    zdir=rotationaxisz;
    
    c=cos(angle);
    s=sin(angle);
    t=1.0-c;
    
    if(leftright==1)
    {
      //rotation direction is to the right
      for(i=ii+1;i<N;++i)
	{
	  ttx=x[i];
	  tty=y[i];
	  ttz=z[i];
	  
	  x[i]=(t*xdir*xdir+c)*ttx+(t*xdir*ydir+s*zdir)*tty+(t*xdir*zdir-s*ydir)*ttz;
	  y[i]=(t*xdir*ydir-s*zdir)*ttx+(t*ydir*ydir+c)*tty+(t*ydir*zdir+s*xdir)*ttz;
	  z[i]=(t*xdir*zdir+s*ydir)*ttx+(t*ydir*zdir-s*xdir)*tty+(t*zdir*zdir+c)*ttz;
	};
    }
    else
    {
      //rotation direction is to the left
      for(i=0;i<ii;++i)
      {
	ttx=x[i];
	tty=y[i];
	ttz=z[i];
	
	x[i]=(t*xdir*xdir+c)*ttx+(t*xdir*ydir+s*zdir)*tty+(t*xdir*zdir-s*ydir)*ttz;
	y[i]=(t*xdir*ydir-s*zdir)*ttx+(t*ydir*ydir+c)*tty+(t*ydir*zdir+s*xdir)*ttz;
	z[i]=(t*xdir*zdir+s*ydir)*ttx+(t*ydir*zdir-s*xdir)*tty+(t*zdir*zdir+c)*ttz;
      };
    };
     
    for(i=0;i<N;++i)
    {
      x[i]=x[i]+transx;
      y[i]=y[i]+transy;
      z[i]=z[i]+transz;
    };
  
  };
}

/*rotate the x,y,z coordinates of the chain between monomers ii and jj on an axis connecting ii and jj, 
angle is the rotation angle in radians

If you do not want to change the angle of a bond, then you should specify the axis in the direction of one of the bonds,
either to the left or right of ii

*/
void rotate2pnt(int ii, int jj, double angle)
{
  int i,j;
  double c,s,t,xdir,ydir,zdir;
  double transx,transy,transz,r;
  double ttx,tty,ttz;

  //define the rotation axis between the iith and jjth monomers, normalized rotation axis
  rotationaxisx=x[ii]-x[jj];
  rotationaxisy=y[ii]-y[jj];
  rotationaxisz=z[ii]-z[jj];

  c=sqrt(rotationaxisx*rotationaxisx+rotationaxisy*rotationaxisy+rotationaxisz*rotationaxisz);

  rotationaxisx/=c;
  rotationaxisy/=c;
  rotationaxisz/=c;


  //translate such that x[ii]=y[ii]=z[ii]=0
  transx=x[ii];
  transy=y[ii];
  transz=z[ii];
  
  for(i=0;i<N;++i)
    {
      x[i]=x[i]-transx;
      y[i]=y[i]-transy;
      z[i]=z[i]-transz;
    };
  
  xdir=rotationaxisx;
  ydir=rotationaxisy;
  zdir=rotationaxisz;
  
  c=cos(angle);
  s=sin(angle);
  t=1.0-c;
  
  //rotation direction is to the right
  for(i=ii+1;i<jj;++i)
    {
      ttx=x[i];
      tty=y[i];
      ttz=z[i];
      
      x[i]=(t*xdir*xdir+c)*ttx+(t*xdir*ydir+s*zdir)*tty+(t*xdir*zdir-s*ydir)*ttz;
      y[i]=(t*xdir*ydir-s*zdir)*ttx+(t*ydir*ydir+c)*tty+(t*ydir*zdir+s*xdir)*ttz;
      z[i]=(t*xdir*zdir+s*ydir)*ttx+(t*ydir*zdir-s*xdir)*tty+(t*zdir*zdir+c)*ttz;
    };
     
  for(i=0;i<N;++i)
    {
      x[i]=x[i]+transx;
      y[i]=y[i]+transy;
      z[i]=z[i]+transz;
    };
  
}

