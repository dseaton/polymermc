#include <math.h>
#include <stdio.h>

#include "wanglandau.h"
#include "metrop.h"
#include "model.h"
#include "rand.h"

//ROUTINES TO PERFORM STANDARD METROPOLIS MONTE CARLO


//routine to encapsulate the Metropolis algorithm, returns 1 if accepted and 0 if rejected
int Metropolis(double Ei, double Ef)
{
	double R;

	if(Ef<Ei)
    {
		//accept
		currEtot=Ef;
		return 1;
    }
	else
    {
		R=exp(-(Ef-Ei)/T);

		if(randd1()<R)
		{
			//accept
			currEtot=Ef;
			return 1;
		}
		else
		{
			//reject
			return 0;
		};
    };
}

//EASE OF USE ROUTINES
//using Metropolis, begins with a straight polymer and runs SAMPS, reporting every SEP
void Seq(double Tf,int DROPI,int SAMPS, int SEP)
{
	int i,j;
	FILE *ofp;
	double Etemp=0.0;
	double aveE=0.0;
	double aveE2=0.0;
	double aveEerr=0.0;
	initialize();
	
	//fprintf(stderr,"%g\t%g\n",currEtot,Etot());
	
	//JA = WLD2max;
	ofp=fopen("EMAX_METROPvalues.input","w");
	
	for(JA=WLD2min;JA<WLD2max;JA=JA+dWLD2)
	{
		T=Tf;
		aveE=0.0;
		aveE2=0.0;
		aveEerr=0.0;
		//system("rm *.inc");
	
		//relax
		for(i=0;i<DROPI;++i)
		{
			mchybrid();
		};
		
		//sample
		for(i=0;i<SAMPS;++i)
		{
			//fprintf(stderr,"%g\t%g\n",JA,currEtot*invN);
			//fprintf(stderr,"%g\t%g\t%g\n",Etot()-currEtot,Etot()-slowbondEnergy()-nonbondEnergy(),nonbondEnergy());
			//fprintf(stderr,"%d\t%g\t%g\t%g\t%g\n",i*SEP,currEtot*invN,Etot()*invN,currnonbondE*invN,nonbondEnergy()*invN);
			if(currEtot*invN < Etemp)
				writeinc(0,currEtot*invN);
		
			aveE+=currEtot*invN;
			aveE2+=currEtot*invN*currEtot*invN;
			
			for(j=0;j<SEP;++j)
				mchybrid();
		};
		
		aveE/=(1.0*SAMPS);
		aveE2/=(1.0*SAMPS);
		aveEerr = sqrt( fabs(aveE2 - aveE*aveE) );
		
		fprintf(ofp,"%g\t%g\t%g\n",JA,aveE,aveEerr);
		//fprintf(ofp,"\n");
		fflush(ofp);
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
	JA = WLD2max;
	
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


//excutes all of the combined MC moves, performs one step of each type of move, and one sweep of diffusion moves
void mchybrid(void)
{
	int i;
	
	for(i=0;i<N;++i)
		mcdiff();

	mcreptation();
	mcrandpivot();
	mccrankshaft();

}

//Diffuses a single monomer
void mcdiff(void)
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
	
	if(Metropolis(Ei,Ef)==1)
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


//performs one slithering snake MC attempt
void mcreptation(void)
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
	
	if(Metropolis(Ei,Ef)==1)
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


/*
Monte Carlo pivot around a random direction
Random direction generated using the Marsaglia method

NOTE this algorithm should NOT change the bond lengths, use checkbonds() to check, but should change one bond angle at a time
*/
void mcrandpivot(void)
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
	//Local bond and bond angle energy - assumes bond lengths don't change (unlike above commented out section)
	Ebondi = Ebondi + localE(m);
	
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
	//Local bond and bond angle energy - assumes bond lengths don't change (unlike above commented out section)
	Ebondf = Ebondf + localE(m);
	
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);	
	//Ef=Etot();
	
	attemptpivot+=1;
	
	if(Metropolis(Ei,Ef)==1)
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

/*Monte Carlo bond fluctuation  algorithm
Similar to mcpivot, but rotates monomers betwee two randomly chosen monomers by a small amount 
axis of rotation is between the two chosen monomers

NOTE this algorithm should NOT change the bond lengths, use checkbonds() to check, but should change one bond angle at a time
*/
void mccrankshaft(void)
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
	//Local bond and bond angle energy - assumes the above in unecessary
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
	//Local bond and bond angle energy - assumes the above in unecessary
	Ebondf = localE(m1) + localE(m2);
	
	//Bond-angles should be unaffected by this move
	
	Ef = Ei + (Ebondf - Ebondi) + (Enblocalf - Enblocali);
	
	attemptside+=1;
	
	if(Metropolis(Ei,Ef)==1)
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


