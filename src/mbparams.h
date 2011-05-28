
/** \file mbparams.h */
#ifndef MBPARAMS_H
#define MBPARAMS_H 1

#include <stdlib.h>
#include <string>
#include <list>

using namespace::std;

/// Data structure to hold parameters for running Mr Bayes
/** This data structure simply holds all parameters one might need to
    run Mr Bayes */
class MrBayesParameters {
public:
    list        <string> inputNexFiles; // List of files to process
    int         MBP_nst; // GTR rate variance model.
    string      MBP_rates; // gamma shaped rate variation.
    int         MBP_nruns; // Number of independent runs. 
    int         MBP_ngen; // perform 10,000 runs.
    int         MBP_sampleFreq; // Sample every 10 chains.
    string      MBP_pathToMB;// = "~/Desktop/mrbayes-3.1.2/mb"; // location of mb
    int         MBP_saveOutput;// = 0; // 1 save output, 0 not 
    // Parameters for MPI implementation of mr bayes.
    int         MPI_np;// = 4; // Number of processors to use.
    string      parametersBlockFileName;
private:
};
#endif
