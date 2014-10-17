#include "rand.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double invttrand; //inverse parameter for quick calculations
struct seed_type seed = {314159265, 362436069, 521288629};

/*
  XXXXXXXXXXXXXXXXXXXXXXXXXX
  KISS random number generator and associated functions
  XXXXXXXXXXXXXXXXXXXXXXXXXX
*/  

/*function to set seeds:*/
void kisset(ii,jj,kk) unsigned int ii,jj,kk;
{ seed.i=ii; seed.j=jj; seed.k=kk; }

/*function to produce random number:*/
unsigned int kiss(void)
{ seed.j =  seed.j ^ (seed.j<<17);
  seed.k = (seed.k ^ (seed.k<<18)) & 0x7FFFFFFF;
  return (seed.i=69069*seed.i+23606797) +
  (seed.j ^= (seed.j>>15)) + (seed.k ^= (seed.k>>13));
}

/*sets one of the kiss seeds from the time*/
void shelltimeseed(unsigned int tseed)
{
   invttrand=1.0/((double) UINT_MAX);
  
  //sets the new random number sequence based on the shell time seed
  //tseed is given by $RANDOM from the unix command line  
  kisset(314159265,362436069,tseed);
  //fprintf(stderr,"#seeds: %d,%d,%d\n",314159265,362436069,tseed);
}

/*sets one of the kiss seeds from the time*/
void gettimeseed(void)
{
  int i;

  time_t rseed;
  unsigned int tseed, high,low;

  rseed=time(&rseed);
  tseed=(unsigned int) rseed;
  high=tseed % 1000;
  low=(tseed /1000) % 1000000;
  low=((low+high)*611953)%1000000;
  tseed=high*1000000+low;

  invttrand=1.0/((double) UINT_MAX);

  //fprintf(stderr,"#seeds: %d,%d,%d\n",314159265,362436069,tseed);
  RSEED = tseed; //Purely for printing reasons in main.c
  
  kisset(314159265,362436069,tseed);

}

/*uniformly distributed random number (0,1) in double*/
double randd1(void)
{ 
  double trand;
  trand=(double) kiss();
  trand=trand*invttrand;

  while((trand ==0.0) || (trand ==1.0))
    {
       trand=(double) kiss();
       trand=trand*invttrand;
    };

  return trand;
}

/*Gaussian distributed random number*/
void gaussian(double m, double sigma, double *x1, double *x2)
{
  double pi=2*acos(0.0);
  double r1=randd1(), r2=randd1();
  double s,mtheta;
  double t1,t2;

  s=sigma*sqrt(-2.0*log(r1));
  mtheta=2.0*pi*r2;
  t1=m+s*cos(mtheta);
  t2=m+s*sin(mtheta);
  /*  while((fabs(t1)>2.0)||(fabs(t2)>2.0))
    {
    */  r1=randd1();
      r2=randd1();
      s=sigma*sqrt(-2.0*log(r1));
      mtheta=2.0*pi*r2;
      t1=m+s*cos(mtheta);
      t2=m+s*sin(mtheta);
      /*    };*/

  *x1=t1;
  *x2=t2;

  /*  fprintf(stderr,"%g\n%g\n",t1,t2);
   */
}



