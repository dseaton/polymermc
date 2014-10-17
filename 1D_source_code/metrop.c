#include <math.h>
#include <stdio.h>

#include "wanglandau.h"
#include "metrop.h"
#include "model.h"
#include "rand.h"

//ROUTINES TO PERFORM STANDARD METROPOLIS MONTE CARLO


//routine to encapsulate the Metropolis algorithm, returns 1 if accepted and 0 if rejected
int Metropolis(double Ei, double Ef, double Efastf)
{
  double R;

  if(Ef<Ei)
    {
      //accept
      currEtot=Ef;
      currnonbondE=Efastf;
      return 1;
    }
  else
    {
      R=exp(-(Ef-Ei)/T);

      if(randd1()<R)
	{
	  //accept
	  currEtot=Ef;
	  currnonbondE=Efastf;
	  return 1;
	}
      else
	{
	  //reject
	  return 0;
	};
    };
}


//excutes all of the combined MC moves, performs one step of each type of move, and one sweep of diffusion moves
void mchybrid(void)
{
  int i;
  
  mcreptation();
  mcrandpivot();
  mccrankshaft();
  
  for(i=0;i<N;++i)
    mcdiff();


}

//performs one slithering snake MC attempt
void mcreptation(void)
{
  int i,j;
  double Ei,Ef,R,r;
  double xi1,xi2,xi;
  double rx,ry,rz;
  double Eloci,Elocf,Efasti,Efastf;

  for(j=0;j<N;++j)
    {
      tx[j]=x[j];
      ty[j]=y[j];
      tz[j]=z[j];
    };

  Ei=currEtot;
  Efasti=currnonbondE;

  //generate a random direction using the marsaglia (1973) method
  //new relative position of the moved endpoint
  xi1=1.0-2.0*randd1();
  xi2=1.0-2.0*randd1();
  xi=xi1*xi1+xi2*xi2;
  
  while(xi>1.0)
    {
      xi1=1.0-2.0*randd1();
      xi2=1.0-2.0*randd1();
      xi=xi1*xi1+xi2*xi2;
    };

  rx=2.0*xi1*sqrt(1.0-xi);
  ry=2.0*xi2*sqrt(1.0-xi);
  rz=1.0-2.0*xi;

  //chose whether to move the N-1 monomer or the 0th monomer
  if(randd1()<0.5)
    {  
      Eloci=localE(0);
      //move 0th monomer to end

      //measure the length of the 0-1 bond
      r=sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0])+(z[1]-z[0])*(z[1]-z[0]));
      //r=1.0;

      //more x[1] to x[0], etc...
      for(i=0;i<N-1;++i)
	{
	  x[i]=tx[i+1];
	  y[i]=ty[i+1];
	  z[i]=tz[i+1];
	};

      //new tail has random direction, but same bond length as old head had
      x[N-1]=x[N-2]+r*rx;
      y[N-1]=y[N-2]+r*ry;
      z[N-1]=z[N-2]+r*rz;

      Elocf=localE(N-1);
    }
  else
    {  
      Eloci=localE(N-1);
      //move N-1th monomer to head

      //measure the length of the N-1 to N-2 bond
      r=sqrt((x[N-1]-x[N-2])*(x[N-1]-x[N-2])+(y[N-1]-y[N-2])*(y[N-1]-y[N-2])+(z[N-1]-z[N-2])*(z[N-1]-z[N-2]));
      //r=1.0;

      //more x[0] to x[1], etc...
      for(i=1;i<N;++i)
	{
	  x[i]=tx[i-1];
	  y[i]=ty[i-1];
	  z[i]=tz[i-1];
	};

      //new tail has random direction, but same bond length as old head had
      x[0]=x[1]+r*rx;
      y[0]=y[1]+r*ry;
      z[0]=z[1]+r*rz;
      
      Elocf=localE(0);
    }

  Efastf=fastnonbondEnergy();
  Ef=Ei-Eloci-Efasti+Elocf+Efastf;

  attemptsnake+=1;

  if(Metropolis(Ei,Ef,Efastf)==1)
    acceptsnake+=1;
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


/*
Monte Carlo pivot around a random direction
Random direction generated using the Marsaglia method

NOTE this algorithm should NOT change the bond lengths, use checkbonds() to check, but should change one bond angle at a time
*/
void mcrandpivot(void)
{
  int i,j,leftright;
  double Ei,Ef,R,r;
  double Eloci,Efasti,Efastf;

  for(j=0;j<N;++j)
    {
      tx[j]=x[j];
      ty[j]=y[j];
      tz[j]=z[j];
    };

  i=(int) (randd1()*(N-2)+1.0);  //note, never choses monomers on the ends, as this would do nothing to the internal degrees of freedom

  //choose the left or right portion of the chain to rotate
  if(randd1()<0.5)
	leftright=1;  //used in randrotate( )
  else
	leftright=-1; //used in randrotate( )

  Ei=currEtot;
  Eloci=localE(i);
  Efasti=currnonbondE;

  randrotate(i,leftright);

  Efastf=fastnonbondEnergy();
  Ef=Ei-Eloci-Efasti+localE(i)+Efastf;

  attemptpivot+=1;

  if(Metropolis(Ei,Ef,Efastf)==1)
    acceptpivot+=1;
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

/*Monte Carlo bond fluctuation  algorithm
Similar to mcpivot, but rotates monomers betwee two randomly chosen monomers by a small amount 
axis of rotation is between the two chosen monomers

NOTE this algorithm should NOT change the bond lengths, use checkbonds() to check, but should change one bond angle at a time
*/
void mccrankshaft(void)
{
  int i,j;
  double Ei,Ef,R,r;
  double Eloci,Efasti,Efastf;

  for(j=0;j<N;++j)
    {
      tx[j]=x[j];
      ty[j]=y[j];
      tz[j]=z[j];
    };

  i=(int) (randd1()*N);
  j=(int) (randd1()*N);

  Ei=currEtot;
  Eloci=localE(i)+localE(j);
  Efasti=currnonbondE;
    
  rotate2pnt(i,j,DC*2.0*PI*(2.0*randd1()-1.0));

  Efastf=fastnonbondEnergy();
  Ef=Ei-Eloci-Efasti+localE(i)+localE(j)+Efastf;

  attemptside+=1;

  if(Metropolis(Ei,Ef,Efastf)==1)
    acceptside+=1;
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


//Diffuses a single monomer
void mcdiff(void)
{
  int i,j; //randomly chosen monomer
  double ox,oy,oz; //old positions of the monomer
  double Ei,Ef; //intial and final energies (before and after diffusion)
  double R;
  double Eloci,Efasti,Efastf,Eotheri,Eotherf;
  double r;

  //randomly choose a monomer
  i=(int) (randd1()*N);

  Ei=currEtot;
  Eloci=localE(i);
  Efasti=currnonbondE;
  
  Eotheri=0.0;
  for(j=0;j<N;++j)
    if((i!=j)&&(j!=(i-1))&&(j!=(i+1)))
      {
	r=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
	Eotheri+=nonbondpot(r);
      };
  
  //diffuse the ith monomer
  ox=x[i];
  oy=y[i];
  oz=z[i];

  x[i]+=DD*(2.0*randd1()-1.0);
  y[i]+=DD*(2.0*randd1()-1.0);
  z[i]+=DD*(2.0*randd1()-1.0);
  
  Eotherf=0.0;
  for(j=0;j<N;++j)
    if((i!=j)&&(j!=(i-1))&&(j!=(i+1)))
      {
	r=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
	Eotherf+=nonbondpot(r);
      };
  
  //old way, slower
  //  Efastf=fastnonbondEnergy();
  //Ef=Ei-Eloci-Efasti+localE(i)+Efastf;
  Efastf=Efasti-Eotheri+Eotherf;
  Ef=Ei-Eloci-Eotheri+localE(i)+Eotherf;
  
  attemptdiff+=1;

  if(Metropolis(Ei,Ef,Efastf)==1)
    acceptdiff+=1;
  else
    {
      //reject
      x[i]=ox;
      y[i]=oy;
      z[i]=oz;
    };
}






//EASE OF USE ROUTINES



//using Metropolis, begins with a straight polymer and runs SAMPS, reporting every SEP
void Seq(double Tf,int DROPI,int SAMPS, int SEP)
{
  int i,j;
  FILE *ofp;
  double Etemp=0.0;
  initialize();
  
  T=Tf;
  ofp=fopen("run.dat","w");

  //system("rm *.inc");

  //relax
  for(i=0;i<DROPI;++i)
    {
		mchybrid();
    };

  //sample
  for(i=0;i<SAMPS;++i)
    {
      fprintf(ofp,"%d\t%g\t%g\t%g\t%g\n",i*SEP,currEtot*invN,(currEtot-currnonbondE)*invN,currnonbondE*invN,Rgyr2());
	  //fprintf(stderr,"%g\t%g\t%g\n",Etot()-currEtot,Etot()-slowbondEnergy()-fastnonbondEnergy(),fastnonbondEnergy());
	  //fprintf(stderr,"%d\t%g\t%g\t%g\t%g\n",i*SEP,currEtot*invN,Etot()*invN,currnonbondE*invN,fastnonbondEnergy()*invN);
      if(currEtot*invN < Etemp)
      	writeinc(0,currEtot*invN);

      fflush(ofp);
      for(j=0;j<SEP;++j)
		mchybrid();
    };

  close(ofp);
}
      

//using Metropolis, scans up and then back down in Temperature
//before initiating scan, runs DROPI mchybrid steps
//then drops DROP, before recording SAMPS at each T
void Tscan(double Ti, double Tf, double dT,int DROPI, int DROP, int SAMPS, int SEP)
{
  int i,j;
  double avgE,avgE2,avgEnonbond,avgRgyr2;
  FILE *ofp;

  initialize();

  ofp = fopen("run.dat","w");

  T=Ti;
  for(i=0;i<DROPI;++i)
    {
      mchybrid();
	};

  //UP SCAN
  for(T=Ti;T<=Tf;T=T+dT)
    {
      for(i=0;i<DROP;++i)
		{
			mchybrid();
		};

      avgE=0.0;
      avgEnonbond=0.0;
      avgRgyr2=0.0;

	  for(i=0;i<SAMPS;++i)
		{
			for(j=0;j<SEP;++j)
			  mchybrid();
				
			avgE+=currEtot;
			avgE2+=currEtot*currEtot;
			avgEnonbond+=currnonbondE;
			avgRgyr2+=Rgyr2();
			
		};

      avgE/=1.0*SAMPS;
      avgE2/=1.0*SAMPS;
	  avgEnonbond/=1.0*SAMPS;
      avgRgyr2/=1.0*SAMPS;

      fprintf(ofp,"%g\t%g\t%g\t%g\t%g\t%d\n",T,avgE*invN,avgEnonbond*invN,(1.0/(T*T))*(avgE2 - avgE*avgE)*invN,avgRgyr2,SAMPS);
	  fflush(ofp);
    };

  //DOWN SCAN
  for(T=Tf;T>=Ti;T=T-dT)
    {
      for(i=0;i<DROP;++i)
		{
			mchybrid();
	  
		};

      avgE=0.0;
	  avgE2=0.0;
      avgEnonbond=0.0;
      avgRgyr2=0.0;

      for(i=0;i<SAMPS;++i)
		{
			for(j=0;j<SEP;++j)
			  mchybrid();
			
			avgE+=currEtot;
			avgE2+=currEtot*currEtot;
			avgEnonbond+=currnonbondE;
			avgRgyr2+=Rgyr2();
	
		};

      avgE/=1.0*SAMPS;
	  avgE2/=1.0*SAMPS;
      avgEnonbond/=1.0*SAMPS;
      avgRgyr2/=1.0*SAMPS;

      fprintf(ofp,"%g\t%g\t%g\t%g\t%g\t%d\n",T,avgE*invN,avgEnonbond*invN,(1.0/(T*T))*(avgE2 - avgE*avgE)*invN,avgRgyr2,SAMPS);
      fflush(ofp);
    };

  fclose(ofp);
}






