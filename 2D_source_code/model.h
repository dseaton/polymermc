#ifndef MODEL_H
#define MODEL_H

//Input Parameters
int N;  //number of monomers
int NBINTERACTION;  //Selects the non-bonded interaction (1 = LJ potential, 2 = Quasi LJ potential)
double JA;  //Semiflexible coefficient (strength), if greater than zero, it turns on (0.5 applied in code, so this is K only)

//General Parameters
#define MAXNEIGHS 500 //maximum number of neighbors per monomer in the nonbonded neighbor table, code terminates if this is ever violated
double D,DC,DD;//D is for pivot angle (1.0 is 2pi), DC is for crank angle, and DD is for diffussion between +/-DD in each dir
//above moves set in the "intialize()" function

#define L 1.0 //ground state nn bond length
double PI,ANG,COSANG;//PI,ANG,COSANG are defined with machine precision in main

#define Rcutnonbond 3.0  //cutoff distance for nonbonded potential, DO NOT CHANGE
#define Rcutbond 1.2 //DO NOT CHANGE, bondpot assumes Rcutbond=1.2, for normalization
double sigmaLJ;  //set in the "initialize()" function
double shiftLJ;  //set in the "initialize()" function
double sigmaFENE;  //set in the "initialize()" function


#define Jbond 2.0 //bond length stiffness, 2 gives a similar potential to the previous harmonic potential which had JL=30
#define Jnonbond 1.0 //strength of nonbonded potential

double currEtot,currnonbondE,T; //current total energy, total bonbonded energy, and temperature
double TTi,TTf,dTT;

double *x,*y,*z; //x,y,z coordinates of the N monomers
double *tx,*ty,*tz;//temporary storage of the x,y,z arrays, STORES OLD CONFIG
int *Nlist;//stores neighbor list in cut-and-join moves

double rotationaxisx,rotationaxisy,rotationaxisz;  //the x,y, and z values for the axis of rotation

double invN,rootinvN; //inverse parameters, 1/N, sqrt(1/N)

int acceptpivot,acceptdiff,attemptpivot,attemptdiff;
int acceptside,attemptside;
int acceptsnake,attemptsnake;
int acceptcutjoin,attemptcutjoin;

double nonbondpot(double r);
double bondpot(double r);

double Etot(void);
double nonbondEnergy(void);
double localE(int ii);
double localnonbondEnergy(int m);
double bondEnergy(void);
double bondangleEnergy(void);

void writeinc(int frame, double energy);
void writepositions(int frame);
void write_mol2(int frame);
void thermoqs(void);
void orderqs(void);

void initialize(void);
void reset_XYZ(void);
double cmass(void);
double icos_core_count(void);
double aveBAngle(void);
double Rgyr2(void);
double gyration_tensor(void);
double K1P,K2P;
double pairnum(void);
double EEdist(void);
double coredensity(void);

void randrotate(int ii, int leftright);
void rotate(int ii, int leftright, double angle);
void rotate2pnt(int ii, int jj, double angle);

#endif
