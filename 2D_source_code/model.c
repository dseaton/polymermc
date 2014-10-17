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
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "model.h"
#include "rand.h"
#include "wanglandau.h"

//Calculates thermodynamic quantities and writes modification-factor labeled files
void thermoqs()
{
	int i,j;	
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
	
	for(j=0;j<D2BINS;j++)
	{
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
					
				if( (lambda < ((wllng[i][j]) - 1.0*(centerE)/T)) && (wllng[i][j] > 0.0) && mask[i][j]==1) 
					lambda = ((wllng[i][j]) - 1.0*(centerE)/T);  
			};
		
			//Central Loop for calculating thermodynamic properties from the DOS
			for(i=0;i<D1BINS;i++)
			{			
				if(mask[i][j]==1)
				{
					//Taking the center of the bin
					centerE = N*( (i/invdWLD1+WLD1min) + 0.5*(WLD1max - WLD1min)/(1.0*D1BINS) );
			
					//Boltzmann Factor
					Bw =  exp( wllng[i][j]  - (centerE)/T - lambda );
		
					//Partition Function
					Z = Z + Bw;
				
					//Average Energy <E>
					U = U + (centerE)*Bw;		   		   
				
					//Average Energy Squared <E^2>
					C = C + (centerE)*(centerE)*Bw;								
				};
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
			
			fprintf(therm_op,"%g\t%g\t%g\t%g\t%g\t%g\n",j/invdWLD2+WLD2min,T,U*invN,C*invN,F*invN,S*invN);	
		}	
		
		fprintf(therm_op,"\n");
	};		
	fclose(therm_op);
}

//Calculates various order-parameters associated with the model (Radius of Gyration, End to End distance, etc...)
void orderqs()
{
	int i,j;	
	double T,Z,Rg,Pn,CoreD,K1tmp,K2tmp,EEd;
	double centerE,lambda,Bw,Ewidth;
	
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
	
	for(j=0;j<D2BINS;j++)
	{
		//Main Temperature Loop
		for(T=TTi;T<=TTf;T=T+dTT)
		{			
			Z = 0.0;  //Initialize the partition function
			Rg = 0.0;  //Initialize the average radius of gyration
			K1tmp = 0.0;
			K2tmp = 0.0;
			Pn = 0.0;    //Initialize the average pair number
			EEd = 0.0;	//Initialize the average end to end distance
			CoreD = 0.0;	//Initialize the core density
			Bw = 0.0;	//Boltzmann weight
			lambda = -1.0e300;	//Normalization shift (max value of DOS considering T)
			centerE = 0.0;	// Taking the center of the energy bin
			Ewidth = 0.5*(WLD1max - WLD1min)/(1.0*D1BINS);
			
			//Finds the max and min of exp( wllng[][] )*exp(Etot*N/T)
			for(i=0;i<D1BINS;i++)
			{
				centerE = N*( (i/invdWLD1+WLD1min) + Ewidth );
				
				if( (lambda < ((wllng[i][j]) - 1.0*(centerE)/T)) && (wllng[i][j] > 0.0) && mask[i][j]==1) 
					lambda = ((wllng[i][j]) - 1.0*(centerE)/T);  
			};
		
			//Central Loop for calculating thermodynamic properties from the DOS
			for(i=0;i<D1BINS;i++)
			{			
				if(mask[i][j]==1 && prun_H[i][j]>0.0)
				{
					//Taking the center of the bin
					centerE = N*( (i/invdWLD1+WLD1min) + Ewidth );
			
					//Boltzmann Factor
					Bw =  exp( wllng[i][j]  - (centerE)/T - lambda );
			
					//Partition Function
					Z = Z + Bw;
			
					//if(prun_H[i][j]>0.0)		
					//{
						//Radius of Gyration
						Rg = Rg + (prun_Rgyr2[i][j]/prun_H[i][j])*Bw;
						//fprintf(stderr,"%g\t%18.10e\n",Rg,(prun_Rgyr2[i][j]/prun_H[i][j])*Bw);
						
						//1st Khalatur Parameter
						K1tmp = K1tmp + (prun_K1P[i][j]/prun_H[i][j])*Bw;
						
						//2nd Khalatur Parameter
						K2tmp = K2tmp + (prun_K2P[i][j]/prun_H[i][j])*Bw;
						
						//Pair Number
						Pn = Pn + (prun_pairnum[i][j]/prun_H[i][j])*Bw;
						
						//End to End distance
						EEd = EEd + (prun_EEd[i][j]/prun_H[i][j])*Bw;
				
						//Bond angle distribution
						//CoreD = CoreD + Hcore[i]*Bw;
					//};
					//fprintf(stderr,"%g\t%g\t%g\t%g\n",T,Z,Rg,Pn);
		
				};
				
			};
			
			//Radius of Gyration
			Rg = Rg/Z;
			//1st Khalatur Parameter
			K1tmp = K1tmp/Z;			
			//2nd Khalatur Parameter
			K2tmp = K2tmp/Z;
			//Pair number
			Pn = Pn/Z;
			//End to end distance
			EEd = EEd/Z;
			
			//fprintf(stderr,"%g\t%g\n",Rg,Pn);
			fprintf(order_op,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",j/invdWLD2+WLD2min,T,Rg,K1tmp,K2tmp,Pn,EEd);
				
		}
		fprintf(order_op,"\n");
	}
	fclose(order_op);
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
	if(r<Rcutbond)
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
	double dxl=0.0,dyl=0.0,dzl=0.0,ll=0.0;
	double dxr=0.0,dyr=0.0,dzr=0.0,lr=0.0;
	double tmp;	
	double Eb=0.0;
	double Enb=0.0;
	double Ea=0.0;
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
	//return (Eb + Enb);


	//calculate nn bond angle (3 body) local interaction
	if(JA != 0.0)
	{
		for(i=1;i<N-1;++i)//loop over all nn bonds
		{
			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i]-x[i+1]; 
			dyr=y[i]-y[i+1]; 
			dzr=z[i]-z[i+1]; 
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=( 1.0 - (dxl*dxr+dyl*dyr+dzl*dzr)/(tmp) );
			
			Ea = Ea + JA*tmp;
		}
	};
	return (Eb + Enb + Ea);

}


//returns the bond angle energy of all bond angles
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

//returns the bonded energy between neighboring monomers
double bondangleEnergy(void)
{
	int i;
	double Ea=0.0;
	double dxl=0.0,dyl=0.0,dzl=0.0,ll=0.0;
	double dxr=0.0,dyr=0.0,dzr=0.0,lr=0.0;
	double tmp;	
	
	if(JA != 0.0)
	{
		for(i=1;i<N-1;++i)//loop over all nn bonds
		{
			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i]-x[i+1]; 
			dyr=y[i]-y[i+1]; 
			dzr=z[i]-z[i+1]; 
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=( 1.0 - (dxl*dxr+dyl*dyr+dzl*dzr)/(tmp) );
			
			Ea = Ea + JA*tmp;
		}
	}

	//Return the total energy
	return (Ea);

}

//Calculates the total energy of the polymer chain (using 2N calculations)
double localnonbondEnergy(int m)
{
	int i;
	double local_Enb=0.0;
	double dist=0.0;
	
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

	//Return the total energy
	return (local_Enb);

/*
	//Simple loop over all monomers interacting via the non-bonded interaction
	//This function should be used in conjunction with a MC move, e.g., diffusion
	for(i=0;i<N;i++)
	{
		if( (m!=i) && (m!=(i-1)) && (m!=(i+1)) )
		{
			//Calculate the nonbonded energy
			dist=sqrt( (x[m]-x[i])*(x[m]-x[i]) + (y[m]-y[i])*(y[m]-y[i]) + (z[m]-z[i])*(z[m]-z[i]) );
			local_Enb = local_Enb + nonbondpot(dist);
		};
	}
*/
}


//returns the bond length and bond angle energies envolving a change in position for monomer #ii
double localE(int m)
{
	int i;
	double Esum=0.0;
	double dist=0.0;
	double dxl=0.0,dyl=0.0,dzl=0.0,ll=0.0;
	double dxr=0.0,dyr=0.0,dzr=0.0,lr=0.0;
	double tmp;
	
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

	//return Esum;

	//calculate nn harmonic bond angle (3 body) local interaction
	//potentially three bond angles that could have changed, <m-2,m-1,m>,<m-1,m,m+1>,and <m,m+1,m+2>
	if(JA != 0.0)
	{
		//left angle
		if(m-2>=0)
		{
			i=m-1;

			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i]-x[i+1]; 
			dyr=y[i]-y[i+1]; 
			dzr=z[i]-z[i+1]; 
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=( 1.0 - (dxl*dxr+dyl*dyr+dzl*dzr)/(tmp) );
			Esum = Esum + JA*tmp;
		}

		//center angle
		if((m-1>=0)&&(m+1<N))
		{
			i=m;

			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i]-x[i+1]; 
			dyr=y[i]-y[i+1]; 
			dzr=z[i]-z[i+1];
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=( 1.0 - (dxl*dxr+dyl*dyr+dzl*dzr)/(tmp) );
			Esum = Esum + JA*tmp;
		}

		//right angle
		if(m+2<N)
		{
			i=m+1;
		
			//dl vector points from i to i-1, left
			dxl=x[i-1]-x[i]; 
			dyl=y[i-1]-y[i]; 
			dzl=z[i-1]-z[i]; 

			//dr vector points from i to i+1, right
			dxr=x[i]-x[i+1]; 
			dyr=y[i]-y[i+1]; 
			dzr=z[i]-z[i+1];
				
			//length of left and right bond length vectors
			ll=(dxl*dxl+dyl*dyl+dzl*dzl);
			lr=(dxr*dxr+dyr*dyr+dzr*dzr);
	
			tmp=sqrt(ll*lr);
			tmp=( 1.0 - (dxl*dxr+dyl*dyr+dzl*dzr)/(tmp) );
			Esum = Esum + JA*tmp;
		}
	};
	
	return Esum;

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
  
  fprintf(ofp,"//N=%d, Etot=%18.10e, JA=%g\n",N,currEtot*invN,JA);
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

	//Initialize order parameters used in "gyrationtensor()"
	K1P=0.0;
	K2P=0.0;
	
	//Initialize JA
	JA=WLD2min;
	
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

	Nlist = malloc( N * sizeof(int) );
	
	//Initialize the tmp coordinates
	for(i=0;i<N;i++)
	{
		tx[i]=0.0;
		ty[i]=0.0;
		tz[i]=0.0;
		Nlist[i]=-1;
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

//returns number of interacting pairs
double icos_core_count(void)
{
	int i,j,k;
	double dx=0.0,dy=0.0,dz=0.0,dist=0.0;
	double edges=0.0,numcores=0.0;
	double NEIGHBOR_DIST=1.25;
	int neighs=0,NUMNEIGHS=12;  //12 neighbors occur when a 13 monomer icosahedral core is present
	int list[N];
	
	for(i=0;i<N;++i)
	{
		edges=0.0;
		neighs=0;
		
		for(j=0;j<N;j++)
		{
			list[j]=-1;
			
			if( i != j )
			{
				dx=x[i]-x[j];
				dy=y[i]-y[j];
				dz=z[i]-z[j];
			
				dist = dx*dx+dy*dy+dz*dz;
				if(dist < NEIGHBOR_DIST)
				{
					list[neighs]=j;
					//fprintf(stderr,"%d\t%d\n",neighs,list[neighs]);
					neighs+=1;
				}
			}
		}
		
		if(neighs >= NUMNEIGHS)
		{
			for(j=0;j<NUMNEIGHS;j++)
				for(k=j+1;k<NUMNEIGHS;k++)
				{
					dx=x[ list[j] ]-x[ list[k] ];
					dy=y[ list[j] ]-y[ list[k] ];
					dz=z[ list[j] ]-z[ list[k] ];
					
					dist = dx*dx+dy*dy+dz*dz;
					if(dist < NEIGHBOR_DIST)
					{
						edges+=1.0;
					}
					//fprintf(stderr,"%g\t%g\n",neighs,edges);
				}
			
			if(edges > 29.0 )//Chose > 29 to account for "near" icos cores (icosahedron has 30 edges)
				numcores+=1.0;
		
			fprintf(stderr,"%d\t%g\n",neighs,edges);
		}
	
		
	};
	
	return numcores;
	
}

//returns the bonded energy between neighboring monomers
double aveBAngle(void)
{
	int i;
	double BAngle=0.0,BAngle2=0.0;
	double dxl=0.0,dyl=0.0,dzl=0.0,ll=0.0;
	double dxr=0.0,dyr=0.0,dzr=0.0,lr=0.0;
	double tmp;	
/*
	for(i=0;i<N-2;i++)//loop over all nn bonds
	{
		//dl vector points from i to i-1, left
		dxl=x[i]-x[i+2]; 
		dyl=y[i]-y[i+2]; 
		dzl=z[i]-z[i+2]; 
		
		BAngle += sqrt(dxl*dxl + dyl*dyl + dzl*dzl);
	}
	
	BAngle /= (N-2);
	
	for(i=0;i<N-2;i++)//loop over all nn bonds
	{
		//dl vector points from i to i-1, left
		dxl=x[i]-x[i+2]; 
		dyl=y[i]-y[i+2]; 
		dzl=z[i]-z[i+2]; 
		
		tmp = (sqrt(dxl*dxl + dyl*dyl + dzl*dzl) - BAngle)*(sqrt(dxl*dxl + dyl*dyl + dzl*dzl) - BAngle);
		BAngle2 += tmp;
	}
	
	return BAngle2;
*/	
	for(i=1;i<N-1;++i)//loop over all nn bonds
	{
		//dl vector points from i to i-1, left
		dxl=x[i-1]-x[i]; 
		dyl=y[i-1]-y[i]; 
		dzl=z[i-1]-z[i]; 
			
		//dr vector points from i to i+1, right
		dxr=x[i]-x[i+1]; 
		dyr=y[i]-y[i+1]; 
		dzr=z[i]-z[i+1]; 
			
		//length of left and right bond length vectors
		ll=(dxl*dxl+dyl*dyl+dzl*dzl);
		lr=(dxr*dxr+dyr*dyr+dzr*dzr);
		
		tmp=sqrt(ll*lr);
		tmp=(dxl*dxr+dyl*dyr+dzl*dzr)/(tmp);
			
		BAngle = BAngle + fabs(tmp);
	}
	
	//Return the total energy
	return (BAngle/(N-2));	

}

//returns the radius of gyration
double Rgyr2(void)
{
	int i;
	double cx=0.0,cy=0.0,cz=0.0;
	double dx=0.0,dy=0.0,dz=0.0;
	double sum=0.0;
	
	for(i=0;i<N;i++)
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
		sum += dx*dx + dy*dy + dz*dz;
    };

	//sum=sqrt(sum);
	sum=sum*invN;

	
	return sum;

}

//returns the radius of gyration
double gyration_tensor(void)
{
	int i;
	double cx=0.0,cy=0.0,cz=0.0;
	double dx=0.0,dy=0.0,dz=0.0;
	double Sxx=0.0,Sxy=0.0,Sxz=0.0;
	double Syx=0.0,Syy=0.0,Syz=0.0;
	double Szx=0.0,Szy=0.0,Szz=0.0;
	double Lx2=0.0,Ly2=0.0,Lz2=0.0;
	double sum=0.0;
	
	for(i=0;i<N;i++)
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
		
		Sxx += dx*dx;
		Sxy += dx*dy;
		Sxz += dx*dz;
		
		Syx += dy*dx;
		Syy += dy*dy;
		Syz += dy*dz;
		
		Szx += dz*dx;
		Szy += dz*dy;
		Szz += dz*dz;
    };

	//Uses Gnu Scientific Library to find eigenvalues and vectors of the gyration tensor
	gsl_matrix *M;
	M = gsl_matrix_alloc(3,3);
	
	gsl_matrix_set(M,0,0,Sxx);
	gsl_matrix_set(M,0,1,Sxy);
	gsl_matrix_set(M,0,2,Sxz);
	
	gsl_matrix_set(M,1,0,Syx);
	gsl_matrix_set(M,1,1,Syy);
	gsl_matrix_set(M,1,2,Syz);
	
	gsl_matrix_set(M,2,0,Szx);
	gsl_matrix_set(M,2,1,Szy);
	gsl_matrix_set(M,2,2,Szz);
	
	//gsl_matrix_fprintf(stdout,M,"%g");
	
	gsl_vector *eval = gsl_vector_alloc(3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);
	
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	
	gsl_eigen_symmv(M,eval,evec,w);
	
	gsl_eigen_symmv_free(w);
	
	gsl_eigen_symmv_sort (eval,M,GSL_EIGEN_SORT_VAL_DESC);
	
	//for(i=0;i<3;i++)
	//	fprintf(stderr,"%d\t%g\n",i,gsl_vector_get(eval,i));
	
	Lz2 = gsl_vector_get(eval,0);
	Ly2 = gsl_vector_get(eval,1);
	Lx2 = gsl_vector_get(eval,2);
	
	gsl_matrix_free(M);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	//Sphere - K1P=K2P=1, Rod - KP1=0,K2P=1, Disk - K1P=K2P=0.5
	//The first Khalatur parameter
	K1P = ( Lx2 + Ly2 )/( Ly2 + Lz2 );
	//The second Khalatur parameter
	K2P = ( Lx2 + Lz2 )/( Ly2 + Lz2 );
	
	//Returns radius of gyration
	return ( Lx2 + Ly2 + Lz2 );
	
}

//returns number of interacting pairs
double pairnum(void)
{
	int i,j;
	double dx=0.0,dy=0.0,dz=0.0,dist=0.0;
	double sum=0.0;
	
	for(i=0;i<N-1;++i)
	{
		for(j=i+2;j<N;j++)
		{
			dx=x[i]-x[j];
			dy=y[i]-y[j];
			dz=z[i]-z[j];
		
			dist = dx*dx+dy*dy+dz*dz;
			if(dist < 1.2)
				sum+=1.0;
		}
	};
	
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
	int i,cm;
	double cx=0.0,cy=0.0,cz=0.0;
	double dx=0.0,dy=0.0,dz=0.0,dist=0.0,MINDIST=1.0e300;
	double density=0.0;
	double numcores=0.0;
	
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
		if( dist < MINDIST )
		{
			MINDIST = dist;	
			cm = i;
		};
	};
	
    for(i=0;i<N;i++)
	{
		if(i!=cm)
		{
			dx=x[i]-x[cm];
			dy=y[i]-y[cm];
			dz=z[i]-z[cm];
			if( sqrt(dx*dx+dy*dy+dz*dz) < 1.2 )
				density+=1.0;	
		};
	};
  
	//if(density >= 12.0)
	//	numcores += 1.0;
	
    return density;	
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

