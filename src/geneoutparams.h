
/** \file geneoutparams.h */

#ifndef GENEOUT_H
#define GENEOUT_H 1

#include <string>
#include <list>
#include <set>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "debugoutput.h"
#include "nexus.h"
#include "kernelmethod.h"
#include "bootstrapparams.h"
#include "mbparams.h"
#include "sampleparams.h"
#include "multindparams.h"
#include "tempprefix.h"


using namespace::std;

/// This class will read in the <NAME> NEXUS block which will
/// specify all the parameters to run our software to identify
/// gene outliers.
class GeneOutParameters 
{
public: 
    GeneOutParameters();
    GeneOutParameters(string fileName);
    GeneOutParameters(std::istream &in);
    friend std::ostream& operator << (std::ostream &out, const GeneOutParameters &someGeneOutParameters);
    friend std::istream& operator >> (std::istream &in, GeneOutParameters &someGeneOutParameters);
    void processArg(int argc, char **argv);

    list <string> inputFileNames; // For backwards compatibility
    list <string> inputFileNamesGroupOne;
    list <string> inputFileNamesGroupTwo;

    // Parameters for mrbayes with their default values.
    MrBayesParameters MrBayesParam;
    BootstrapParameters BootstrapParam;
    MultIndParameters MultIndParam;
    SampleParameters SampleParam;

    time_t      randGenSeed;
    int         doMB; // By default, do not do mr bayes
    int         doBootstrap;      // Default is to do bootstrap
    int         doMultInd;
    int         numInd; // Default is 1 individual per species
    int         noTreeCalc; // Do not skip mr bayes calculation.
    int         concatGroups;   // Default is 0, do not concat groups
    int         numGroupOne; // First file is group one and thats all.
    int         statTest; // 0 mean no stat test, 1 means stat test
    int         geneStatTest; // Do the statistical test based on the hyp test on gene trees
    int         numInitCalc; // 
    int         numStatTests;
    int         skipInitCalc; // Dont skip initial calculation for stat test
    double      mySeparation;
    int         doSim; 
    string      simCommand;
                                    // 1 is bootstrap-align, 2 is tree, 3 is no-bootstrap
    int         multSVMSepMethod;   // When doing multiple svm sep test, which method: bootstrap-align, tree, no-bootstrap
                            // 1 is singlesvmsep, 2 is multiplesvmsep
    int         testType;   // What test to run: single svm calc or multiple svm calc.
    string      tempPrefix;
    int         BSStepOne;
    int         indBootstrap; // 1 means use randomly selected alignments to bootstrap/bootstrap from.
                              // 0 means to use the concatenated alignments. Default is 0.
    int         permuteOrig;  // 1 means permute original alignments and bootstrap/bootstrap to the correct size
                              // 0 means not.
    int         allowAnyPermuation; //1 means allow anything. 0 means dissallow new group one to be all from group one.

protected:
    void init ();
};

// Old usage message
//    string usageMessage = "main <file1> <file2> ... <filen> -(options)\n\n \
//          OPTIONS:\n \
//          -doMB                     <1 run mr bayes (default), 0 not>\n \
//          -MBP_nst                  <GTR rate variance model>\n \
//          -MBP_rates                <gamma shape>\n \
//          -MBP_nruns                <number of runs>\n \
//          -MBP_ngen                 <number of generations>\n \
//          -MBP_sampleFreq           <sample frequency>\n \
//          -MBP_pathToMB             <path to mb file>\n \
//          -MBP_saveOutput           <1 saveoutput, 0 not (default)>\n \
//          -MPI_np                   <Number of processes to use for My Bayes>\n \
//          -burninPercent            <burn-in percentage>\n \
//          -burninNumber             <fixed burn-in number>\n \
//          -randGenSeed              <random generator seed>\n \
//          -doBootstrap              <1 do bootstrap, 0 not>\n \
//          -bootstrapColSize         <number of columns (k) for bootstrap >\n \
//          -bootstrapCount           <number of bootstraps to perform>\n \
//          -doMultInd                <1 do multInd, 0 not>\n \
//          -numInd                   <number of ind per species, assumed grouped sequentially>\n \
//          -numTreesReconstruct      <number of trees to reconstruct>\n \
//          -alignToTreeCommand       <Command which reads in an alignment in arg1, outputs tree to stdout>\n \
//          -doSVM                    <Do SVM calculation>\n \
//          -SVM_sampleSize           <SVM sample size for each set\n \
//          -SVM_resampleSampleSize   <SVM sample size for each set\n \
//          -noTreeCalc               <1 skip Tree calculations, 0 not>\n \
//          -numTreesPerFile          <Number of trees per file is noTreeCalc == 1>\n \
//          -numGroupOne              <Number of files in group one>\n \
//          -concatGroups             <1 concat alignments for group one and two, 0 not>\n \
//          -projectSVD               <1 project the points via SVD, 0 not>\n \
//          -projectSVDcutOff         <cutOff percentage for singular values>\n \
//          -distanceOne              <0 use newick distance, 1 use distance 1>\n \
//          -scaleToOne               <0 no scaline, 1 scaleToOne\n \
//          -doDiffMeans              <calculate difference of Means>\n \
//          -statTest                 <0 no stat test, 1 yes>\n \
//          -geneStatTest             <1 do gene stat test, 0 no>\n \
//          -numStatTests             <Number of stat tests to do>\n \
//          -numInitCalc              <Number of init calculations for separation value>\n \
//          -skipInitCalc             <1 skip init separation calc for stat test, 0 not>\n \
//          -mySeparation             <Separation value if skiping init calculation>\n \
//          -doSim                    <In statTest, simulate data instead of jackknifing>\n \
//          -simCommand               <Command to simulate new data if using statTest>\n \
//          -tempPrefix               <prefix where to write/read temporary files>\n";
 
#endif
