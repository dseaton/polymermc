#ifndef METROP_H
#define METROP_H

//Input File Metropolis Variables
int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
double MTi;		//Initial Temperature used in simple Seq() loop
double MTf;		//Final Tempertaure used in simple Seq() loop
double MdT;		//Tempertaure Increment used in simple Seq() loop
double MSAMPS;	//Number of Samples taking in the Seq() function
double MSEP;		//Number of Point seperating actual taken data
double MDROP;	//Number of MC Steps dropped before each temperature run

//Metropolis routines
int Metropolis(double Ei, double Ef);
void mchybrid(void);

void mcreptation(void);
void mcrandpivot(void);
void mccrankshaft(void);
void mcdiff(void);


void Seq(double Tf,int DROPI,int SAMPS, int SEP);
void Tscan(double Ti, double Tf, double dT,int DROPI, int DROP, int SAMPS, int SEP);

#endif




