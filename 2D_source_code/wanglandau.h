#ifndef WANGLANDAU_H
#define WANGLANDAU_H
//Input Parameters
int EXPLORE;    //Switch that allows the method to add new bins if = to 1, if 0, samples the input energy range
double dWLD1;  //Bin Widths for Primary sampling direction
double WLD1max,WLD1min;  //sampling boundaries for PRIMARY (Energy) direction
double dWLD2;  //Bin Widths for Primary sampling direction
double WLD2max,WLD2min;  //sampling boundaries for PRIMARY (Energy) direction
int numf;  //Used to label DOS for each modification factor
int TotalSweeps,IterSweeps;  //Number of sweeps (total and per iteration)
int MASKCUTLOWE;
double PERCENT_DOS_CUT;  //Cuts and keeps this portion of the DOS, Emax being the point of reference, so only low E bins are cut

double Flatness;
int WLBINVISITS;
int NUM_UNSAMPLED_BINS;
int NumSweepsFlat;
double ModFactorInit;
double IterationFactor;
double ModFactorFinal;
int ProductionRun;
int ProductionBinSamps;

//Parameters which are found from input parameters
int D1BINS;  //Number of Bins in the first direction of DOS and Histogram
int D2BINS;  //Number of Bins in the second direction of DOS and Histogram
double invdWLD1;  //Stores inverse bin widths
double invdWLD2;  //Stores inverse bin widths

//For 2D WL simulations
double **wlH;  //the Wang-Landau accumulated histogram
double **wllng;  //the Wang-Landau natural log of the density of states
int **mask;
double *LOWESTE;
int *LOWESTESTATE;
double *HIGHESTE;
int *HIGHESTESTATE;
double **diffvalues;
double **diffrate;
double **diffattempts;
double **JArate;
double **JAattempts;
double **crankvalues;
double **crankrate;
double **crankattempts;
double **cutjoinrate;
double **cutjoinattempts;

void initialize_prodrun(void);
void prunbin( int Estate, int JAstate );
void MUCAstep(void);
double **prun_pairnum;
double **prun_EEd;
double **prun_Rgyr2;
double **prun_K1P;
double **prun_K2P;
double **prun_H;

double lnwlf;  //natural log of the Wang-Landau update factor, f

//double LOWESTE;

void initWL(void); //initializes for WL
void sweepWL(double sweeeps);//runs sweeps WL attempts
void resetWL(void);//resets H histogram array
double flatWL(void); //returns H_min/H_avg
int sampledBins(void);
void resize_mask(double percentage, int currMF);

void write_diff(void);
void read_diff(void);
void write_EMIN(void);
void write_DOS_H(void);
void readg(void);

//Wang-Landau routines
void production_run();
int WangLandau(double Ei, double Ef, double DBvalue, double JAi, double JAf);

void wlhybrid(void);

void wlJAchange(void);
void wldiffE(void);
void wlsinglecrank(void);
double EPSILON;

void wlrandpivot(void);
void wlcrankshaft(void);
void wldiff(void);
void wlreptation(void);
void wlcutjoin(void);

#endif
