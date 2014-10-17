#ifndef WANGLANDAU_H
#define WANGLANDAU_H
//Input Parameters
double dWLD1;  //Bin Widths for Primary sampling direction
double WLD1max,WLD1min;  //sampling boundaries for PRIMARY (Energy) direction
int numf;  //Used to label DOS for each modification factor
int TotalSweeps,IterSweeps;  //Number of sweeps (total and per iteration)
int nummoves;

double Flatness;
int NumSweepsFlat;
double ModFactorInit;
double IterationFactor;
double ModFactorFinal;
int ProductionRun;
int ProductionBinSamps;

double numbelow_flat;  //gives the number of bins below the flatness criteria

//Parameters which are found from input parameters
int D1BINS;  //Number of Bins in the DOS and Histogram
double invdWLD1;  //Stores inverse bin widths

//For 2D WL simulations
double *wlH;  //the Wang-Landau accumulated histogram
double *wllng;  //the Wang-Landau natural log of the density of states
double *HRg;  //Array used to keep track of radius gyration (eventually holds the average)
double *HEEdist;  //Array used to keep track of end to end distance
double *Hcore; 
double **gr;  
double lnwlf;  //natural log of the Wang-Landau update factor, f
int *mask;

//Associated with the production run
int GRBINS;
double dGR;
double GRmax,GRmin;

double LOWESTE;

void production_run(void);
void histfill_prun(int pbin);
void init_production_run(void);
void initWL(void); //initializes for WL
void sweepWL(double sweeeps);//runs sweeps WL attempts
void resetWL(void);//resets H histogram array
double flatWL(void); //returns H_min/H_avg

void write_DOS_H(void);
void read_DOS_H(void);
void write_normDOS(void);
void write_mov_g(int frame);
void readg(void);
void writeEERgyr(void);
void readmask(void);
void writemask(void);
void write_restart(void);
void read_restart(void);

//Wang-Landau routines
int WangLandau(double Ei, double Ef, double Enbi, double Enbf);

void wlhybrid(void);

void wlrandpivot(void);
void wlcrankshaft(void);
void wldiff(void);
void wlreptation(void);
void wlcutjoin(void);

#endif
