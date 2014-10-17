#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "InputOutput.h"
#include "model.h"
#include "wanglandau.h"
#include "metrop.h"

/*
NOTE - when adding a new input parameter, it is necessary to update all of 
the functions in this file
*/

//This function reads in the input parameters from the input file given at the command line
////////////////////////////////////////////////////////////////////////////
void ReadInput(const char* filename)
{

  char line[200];
  char pname[51];

  FILE* f = fopen(filename, "r");
  if (f == NULL) ErrorMsg(1, filename);

  while (fgets(line, sizeof(line), f) != NULL) {

    if ((line[0] != '#') && (line[0] != '\n')) {

      if (sscanf(line, "%50s", pname) == 1) {

//Number of Monomers in a Single Chain
	if (!strcmp(pname, "N")) {
	  if (sscanf(line, "%*50s %d", &N) != 1) ErrorMsg(2, pname);
	  if (N < 1) ErrorMsg(3, pname);
	  continue;
	}
	
//Nonbonded Potential Selection
	if (!strcmp(pname, "NBINTERACTION")) {
	  if (sscanf(line, "%*50s %d", &NBINTERACTION) != 1) ErrorMsg(2, pname);
	  if ( (NBINTERACTION < 1) || (NBINTERACTION > 4) ) ErrorMsg(3, pname);
	  continue;
	}

//Bin Width for Primary Sampling Direction	
	if (!strcmp(pname, "dWLD1")) {
	  if (sscanf(line, "%*50s %lg", &dWLD1) != 1) ErrorMsg(2, pname);
	  if (dWLD1 > 10.0 || (dWLD1 < 0.0001)) ErrorMsg(3, pname);  //Generic Boundary Limits
	  continue;
	}

//Maximum Boundary for Primary Sampling Direction	
	if (!strcmp(pname, "WLD1min")) {
	  if (sscanf(line, "%*50s %lg", &WLD1min) != 1) ErrorMsg(2, pname);
	  if (WLD1min > 20.0) ErrorMsg(3, pname);  
	  continue;
	}
	
//Minimum Boundary for Primary Sampling Direction	
	if (!strcmp(pname, "WLD1max")) {
	  if (sscanf(line, "%*50s %lg", &WLD1max) != 1) ErrorMsg(2, pname);
	  if (WLD1max < -20.0) ErrorMsg(3, pname);  
	  continue;
	}

//Flatness Criteria
	if (!strcmp(pname, "Flatness")) {
	  if (sscanf(line, "%*50s %lg", &Flatness) != 1) ErrorMsg(2, pname);
	  if ( (Flatness <= 0.0) || (Flatness >= 1.0) ) ErrorMsg(3, pname);
	  continue;
	}

//Initial Modification Factor
	if (!strcmp(pname, "ModFactorInit")) {
	  if (sscanf(line, "%*50s %lg", &ModFactorInit) != 1) ErrorMsg(2, pname);
	  if (ModFactorInit <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Modification Factor Iterator
	if (!strcmp(pname, "IterationFactor")) {
	  if (sscanf(line, "%*50s %lg", &IterationFactor) != 1) ErrorMsg(2, pname);
	  if (IterationFactor <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Final Modification Factor	
	if (!strcmp(pname, "ModFactorFinal")) {
	  if (sscanf(line, "%*50s %lg", &ModFactorFinal) != 1) ErrorMsg(2, pname);
	  if (ModFactorFinal <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Turns on the production run (0 = off, 1 = on)	
	if (!strcmp(pname, "ProductionBinSamps")) {
	  if (sscanf(line, "%*50s %d", &ProductionBinSamps) != 1) ErrorMsg(2, pname);
	  if ( ProductionBinSamps < 1 ) ErrorMsg(3, pname);
	  continue;
	}
		
//Initial Temperture used in Thermoqs()
	if (!strcmp(pname, "TTi")) {
	  if (sscanf(line, "%*50s %lg", &TTi) != 1) ErrorMsg(2, pname);
	  if ( TTi < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Final Temperture used in Thermoqs()
	if (!strcmp(pname, "TTf")) {
	  if (sscanf(line, "%*50s %lg", &TTf) != 1) ErrorMsg(2, pname);
	  if ( TTf < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Temperture Increment used in Thermoqs()
	if (!strcmp(pname, "dTT")) {
	  if (sscanf(line, "%*50s %lg", &dTT) != 1) ErrorMsg(2, pname);
	  if ( dTT < 0 ) ErrorMsg(3, pname);
	  continue;
	}

//Turns on Metropolis Sampling
	if (!strcmp(pname, "MetropolisSampling")) {
	  if (sscanf(line, "%*50s %d", &MetropolisSampling) != 1) ErrorMsg(2, pname);
	  if ( (MetropolisSampling < 0) || (MetropolisSampling > 1) ) ErrorMsg(3, pname);
	  continue;
	}
	
//Initial Temperture used in simple Seq() loop
	if (!strcmp(pname, "MTi")) {
	  if (sscanf(line, "%*50s %lg", &MTi) != 1) ErrorMsg(2, pname);
	  if ( MTi < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Final Temperture used in simple Seq() loop
	if (!strcmp(pname, "MTf")) {
	  if (sscanf(line, "%*50s %lg", &MTf) != 1) ErrorMsg(2, pname);
	  if ( MTf < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Temperture Increment used in simple Seq() loop
	if (!strcmp(pname, "MdT")) {
	  if (sscanf(line, "%*50s %lg", &MdT) != 1) ErrorMsg(2, pname);
	  if ( MdT < 0 ) ErrorMsg(3, pname);
	  continue;
	}

//Number of MC Samples taken in the Seq() function
	if (!strcmp(pname, "MSAMPS")) {
	  if (sscanf(line, "%*50s %lg", &MSAMPS) != 1) ErrorMsg(2, pname);
	  if ( MSAMPS < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Number of Samples separating data in Seq()
	if (!strcmp(pname, "MSEP")) {
	  if (sscanf(line, "%*50s %lg", &MSEP) != 1) ErrorMsg(2, pname);
	  if ( MSEP < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Number of Dropped Samples before each T run
	if (!strcmp(pname, "MDROP")) {
	  if (sscanf(line, "%*50s %lg", &MDROP) != 1) ErrorMsg(2, pname);
	  if ( MDROP < 0 ) ErrorMsg(3, pname);
	  continue;
	}

	ErrorMsg(4, pname);

      }

      else ErrorMsg(5, "");

    }

  }

  fclose(f);

}

//Writes the input parameters to see if they were read in correctly
////////////////////////////////////////////////////////////////////////////
void WriteInput()
{

  printf("\n### Program Parameters #################################\n\n");

  //if(MetropolisSampling == 0)	
  //{	
	printf("#  Number of Monomers  %d\n", N);
	/*
	printf("#  Interaction Constant  %g\n", Jkb);
	printf("#  Number of Energy Bins  %d\n", EBINS);
	printf("#  Maximum Energy Bound  %g\n", WLEmax);
	printf("#  Minimum Energy Bound  %g\n", WLEmin);
	printf("#  Frontier Sampling  %d\n", FrontierSampling);
	printf("#  Maximum Possible Boost Factor  %g\n", MaxBoost);
	printf("#  Number of Sweeps in between Boosts  %d\n", BoostSweeps);
	printf("#  Top Portion of the DOS kept after Boosting  %g\n", Mountain);
	printf("#  Steps Between Smoothing  %d\n", SmoothSkip);
	printf("#  Flatness Criteria  %g\n", Flatness);
	printf("#  Initial Modification Factor  %g\n", ModFactorInit);
	printf("#  Modification Factor Divider  %g\n", IterationFactor);
	printf("#  Final Modification Factor  %g\n", ModFactorFinal);
	*/
  //};
  /*
  if(MetropolisSampling == 1)
  {
	printf("#  Metropolis Sampling  %d\n", MetropolisSampling);
	printf("#  Initial Temperature  %g\n", MTi);
	printf("#  Final Temperature  %g\n", MTf);
	printf("#  Temperature Increment  %g\n", MdT);
	printf("#  Number of Samples in Metropolis  %g\n", SeqSAMPS);
	printf("#  Number of Samples Separating Taken Data  %g\n", SeqSEP);
	printf("#  Number of Point Initially Dropped  %g\n", SeqDROP);
  };
  */
  printf("\n########################################################\n\n");

}


//Error Messages
////////////////////////////////////////////////////////////////////////////
void ErrorMsg(int message, const char* arg)
{
  switch (message) {
  case  0 : {
    printf("Syntax: error in or command line arguments [input-file]\n");
    break;
  }
  case  1 : {
    printf("Error: Input-file '%s' not found.\n", arg);
    break;
  }
  case  2 : {
    printf("Error: Parameter '%s' has invalid format.\n", arg);
    break;
  }
  case  3 : {
    printf("Error: Parameter '%s' is out of range.\n", arg);
    break;
  }
  case  4 : {
    printf("Error: Parameter '%s' is unknown.\n", arg);
    break;
  }
  
  case  5 : {
    printf("Error: Unreadable parameter.\n");
    break;
  }
  
  default : printf("Error.\n");
  }
  printf("Program aborted.\n");
  exit(1);
}

