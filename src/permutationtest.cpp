// $Rev: 545 $ $Date: 2010-03-07 18:52:38 -0500 (Sun, 07 Mar 2010) $

/// \file permutationtest.cpp 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include "svm.h"
#include <map>
#include <list>
#include <set>
#include <math.h>
#include "alignment.h"
#include "debugoutput.h"
#include "tempprefix.h"
#include "kernelmethod.h"


// Used to control debug output. 
int DEBUG_OUTPUT = 0;

// Used to control prefix location of temporary files. Default is nothing,
// i.e. the current directory.
char emptyPrefix[] = "";
string tempPrefix = emptyPrefix;

using namespace::std;

int main (int argc, char **argv)
{
    int mainStartTime = time(0);
    string usageMessage = "permutationtest <file1> <file2> ... <filen> -(options)\n\n \
          OPTIONS:\n \
          -doMB                     <1 run mr bayes (default), 0 not>\n \
          -MBP_nst                  <GTR rate variance model>\n \
          -MBP_rates                <gamma shape>\n \
          -MBP_nruns                <number of runs>\n \
          -MBP_ngen                 <number of generations>\n \
          -MBP_sampleFreq           <sample frequency>\n \
          -MBP_pathToMB             <path to mb file>\n \
          -MBP_saveOutput           <1 saveoutput, 0 not (default)>\n \
          -MPI_np                   <Number of processes to use for My Bayes>\n \
          -burninPercent            <burn-in percentage>\n \
          -burninNumber             <fixed burn-in number>\n \
          -randGenSeed              <random generator seed>\n \
          -doJackknife              <1 do jackknife, 0 not>\n \
          -jackknifeColSize         <number of columns (k) for jackknife >\n \
          -jackknifeCount           <number of jackknifes to perform>\n \
          -distanceOne              <0 use newick distance, 1 use distance 1>\n \
          -alignToTreeCommand       <Command which reads in an alignment in arg1, outputs tree to stdout>\n \
          -noTreeCalc               <1 skip Tree calculations, 0 not>\n \
          -numTreesPerFile          <Number of trees per file is noTreeCalc == 1>\n \
          -numGroupOne              <Number of files in group one>\n \
          -concatGroups             <1 concat alignments for group one and two, 0 not>\n \
          -projectSVD               <1 project the points via SVD, 0 not>\n \
          -projectSVDcutOff         <cutOff percentage for singular values>\n \
          -distanceOne              <0 use newick distance, 1 use distance 1>\n \
          -scaleToOne               <0 no scaline, 1 scaleToOne\n \
          -numStatTests             <Number of stat tests to do>\n \
          -doSim                    <In statTest, simulate data instead of jackknifing>\n \
          -simCommand               <Command to simulate new data if using statTest>\n \
          -bootstrapAlignments      <1 means add extra alignments using bootstrapping, 0 not>\n \
          -bootstrapCount           <Number of bootstrap alignemnts to create per input alignments>\n \
          -bootstrapConcat          <1 means when bootstraping, make new alignments from concat of groups>\n \
          -bootstrapColSize         <Number of columns to choose>\n";

            
          
    // Parse command line input.
    // Anything not preceded by '-' we will consider as an input file for mrbayes.
    // We assume each input file contains the same taxa, but each file corresponds
    // to a different gene.
    list <string> inputFileNames;

    // Parameters for mrbayes with their default values.
    MrBayesParameters myMrBayesParam;

        myMrBayesParam.MBP_nst = 6; // GTR rate variance model.
        myMrBayesParam.MBP_rates = "invgamma"; // gamma shaped rate variation.
        myMrBayesParam.MBP_nruns = 2; // Number of independent runs. 
        myMrBayesParam.MBP_ngen = 10000; // perform 10,000 runs.
        myMrBayesParam.MBP_sampleFreq = 10; // Sample every 10 chains.

        myMrBayesParam.MBP_pathToMB = "mb";
        myMrBayesParam.MBP_saveOutput = 0; // 1 save output, 0 not 

        // Parameters for MPI implementation of mr bayes.
        myMrBayesParam.MPI_np = 4; // Number of processors to use.

    JackknifeParameters myJackknifeParam;
        myJackknifeParam.jackknifeColSize = 10; // Default
        myJackknifeParam.jackknifeCount = 1000; // Default
        myJackknifeParam.alignToTreeCommand = "./run_dnadist_neighbor_mult"; // Default
    
    SampleParameters mySampleParam;
        mySampleParam.doSVM = 1;
        mySampleParam.doDiffMeans = 1;
        mySampleParam.SVM_sampleSize = 10; // Default
        mySampleParam.SVM_resampleSize = 10; // Default
        mySampleParam.projectSVD = 0; // Do no project with svd by default
        mySampleParam.projectSVDcutOff = 1; // Default value for svd singular value projection
        mySampleParam.distanceOne = 0; // 0 means use distance specified in newick, 1 means use 1.
        mySampleParam.scaleToOne = 0; // 0 means no scaling, 1 means scale treeLength to 1
        mySampleParam.modelType = LINEAR; // DEFAULT
        mySampleParam.numTreesPerFile = -1; // Default is unspecified!
        // Parameters for the burnin type and amount.
        mySampleParam.burninFormat = -1;    // 0 will mean percentage, 1 will mean fixed number.
        mySampleParam.burninPercent = 0.25; // Assume 25% burnin.
        mySampleParam.burninNumber = 0;    // Default is 0 since we assume default of 25% burnin


    // These are all variables used by the local program
    time_t      randGenSeed             = 0;
    int         doMB                    = 1; // By default, do mr bayes
    int         doJackknife             = 0;      // Default is not to do jackknife
    int         noTreeCalc              = 0; // Do not skip mr bayes calculation.
    int         concatGroups            = 0;   // Default is 0, do not concat groups
    int         numGroupOne             = 1; // First file is group one and thats all.

    int         numStatTests            = 10;
    int         doSim                   = 0; 
    string      simCommand;
    int         bootstrapAlignments     = 0;
    int         bootstrapCount          = 10; // How many bootstrap alignments to add per alignments read
    int         bootstrapColSize        = 0;  // col size of bootstrapping
    int         bootstrapConcat         = 0;  // means concat all alignments for a group and when adding
                                              // a new boostrap alignment, add from concat alignments
    list <list <string> > TreeListListFileNames; // For getTreeMrBayes
    list <string> treeFileNamesGroupOne; // This should be the exact and literal filenames 
    list <string> treeFileNamesGroupTwo; // This should be the exact and literal filenames 
    list <string> treeFileNames; // This should be the exact and literal filenames 

    set <unsigned> groupOneTrees;
    set <unsigned> groupTwoTrees;

    cout << " ***** Reading in file names and command-line arguments *****" << endl;
    for (int i=1;i<argc;i++)
    {
        if (argv[i][0] != '-'){
            // This argument is an input file for mrbayes
            string newMrbayesInput(argv[i]);
            inputFileNames.push_back(newMrbayesInput);
        }
        else {
            // This argument is an option. Parse according to case.
            string tempString = argv[i];
            string tempString2;
            tempString = tempString.substr(1,tempString.size());
            if (DEBUG_OUTPUT >= 1){
                cout << "parameter " << tempString << endl;
            }

            if (tempString == "doMB")
            {
                doMB = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doMB = " << doMB << endl;
                }
            }
            else if (tempString == "MBP_nst")
            {
                myMrBayesParam.MBP_nst = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_nst = " << myMrBayesParam.MBP_nst << endl;
                }
            }
            else if (tempString == "MBP_rates")
            {
                myMrBayesParam.MBP_rates = argv[++i];
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_rates = " << myMrBayesParam.MBP_rates << endl;
                }
            }
            else if (tempString == "MBP_nruns")
            {
                myMrBayesParam.MBP_nruns = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_nruns = " << myMrBayesParam.MBP_nruns << endl;
                }
            }
            else if (tempString == "MBP_ngen")
            {
                myMrBayesParam.MBP_ngen = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_ngen = " << myMrBayesParam.MBP_ngen << endl;
                }
            }
            else if (tempString == "MBP_sampleFreq")
            {
                myMrBayesParam.MBP_sampleFreq = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_sampleFreq = " << myMrBayesParam.MBP_sampleFreq << endl;
                }
            }
            else if (tempString == "MBP_pathToMB")
            {
                myMrBayesParam.MBP_pathToMB = argv[++i];
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_pathToMB = " << myMrBayesParam.MBP_pathToMB << endl;
                }
            }
            else if (tempString == "MBP_saveOutput")
            {
                myMrBayesParam.MBP_saveOutput = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_saveOutput = " << myMrBayesParam.MBP_saveOutput << endl;
                }
            }
            else if (tempString == "MPI_np")
            {
                myMrBayesParam.MPI_np = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MPI_np = " << myMrBayesParam.MPI_np << endl;
                }
            }
            else if (tempString == "burninPercent")
            {
                mySampleParam.burninFormat = 0;
                sscanf(argv[++i],"%lf",&mySampleParam.burninPercent);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: burninFormat = " << mySampleParam.burninFormat << endl;
                    cout << "       Setting: burninPercent = " << mySampleParam.burninPercent << endl;
                }
            }
            else if (tempString == "burninNumber")
            {
                mySampleParam.burninFormat = 1;
                mySampleParam.burninNumber = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: burninFormat = " << mySampleParam.burninFormat << endl;
                    cout << "       Setting: burninNumber = " << mySampleParam.burninNumber << endl;
                }
            }
            else if (tempString == "randGenSeed")
            {
                randGenSeed = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: randGenSeed = " << randGenSeed << endl;
                }
            }
            else if (tempString == "doJackknife")
            {
                doJackknife = atoi(argv[++i]); 
                if (doJackknife == 1){
                    doMB = 0;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doJackknife = " << doJackknife << endl;
                }
                mySampleParam.burninFormat = 1;
                mySampleParam.burninNumber = 0;
            }
            else if (tempString == "jackknifeColSize")
            {
                myJackknifeParam.jackknifeColSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: jackknifeColSize = " << myJackknifeParam.jackknifeColSize << endl;
                }
            }
            else if (tempString == "jackknifeCount")
            {
                myJackknifeParam.jackknifeCount = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: jackknifeCount = " << myJackknifeParam.jackknifeCount << endl;
                }
            }
            else if (tempString == "alignToTreeCommand")
            {
                myJackknifeParam.alignToTreeCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: alignToTreeCommand = " << myJackknifeParam.alignToTreeCommand << endl;
                }
            }
            else if (tempString == "SVM_sampleSize")
            {
                mySampleParam.SVM_sampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: SVM_sampleSize = " << mySampleParam.SVM_sampleSize << endl;
                }
            }
            else if (tempString == "SVM_resampleSize")
            {
                mySampleParam.SVM_resampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: SVM_resampleSize = " << mySampleParam.SVM_resampleSize << endl;
                }
            }
            else if (tempString == "noTreeCalc")
            {
                noTreeCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: noTreeCalc = " << noTreeCalc << endl;
                    //cout << "       Setting: groupOne
                }
            }
            else if (tempString == "numTreesPerFile")
            {
                mySampleParam.numTreesPerFile = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numTreesPerFile = " << mySampleParam.numTreesPerFile << endl;
                }
            }
            else if (tempString == "numGroupOne")
            {
                numGroupOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numGroupOne = " << numGroupOne << endl;
                }
            }
            else if (tempString == "numStatTests")
            {
                numStatTests = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numStatTests = " << numStatTests << endl;
                }
            }
            else if (tempString == "concatGroups")
            {
                concatGroups = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: concatGroups = " << concatGroups << endl;
                }
            }
            else if (tempString == "projectSVD")
            {
                mySampleParam.projectSVD = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: projectSVD = " << mySampleParam.projectSVD << endl;
                }
            }
            else if (tempString == "projectSVDcutOff")
            {
                mySampleParam.projectSVDcutOff = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: projectSVDcutOff = " << mySampleParam.projectSVDcutOff << endl;
                }
            }
            else if (tempString == "distanceOne")
            {
                mySampleParam.distanceOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: distanceOne = " << mySampleParam.distanceOne << endl;
                }
            }
            else if (tempString == "scaleToOne")
            {
                mySampleParam.scaleToOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: scaleToOne = " << mySampleParam.scaleToOne << endl;
                }
            }
            else if (tempString == "doSim")
            {
                doSim = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doSim = " << doSim << endl;
                }
            }
            else if (tempString == "doSVM")
            {
                mySampleParam.doSVM = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doSVM = " << mySampleParam.doSVM << endl;
                }
            }
            else if (tempString == "doDiffMeans")
            {
                mySampleParam.doDiffMeans = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doDiffMeans = " << mySampleParam.doDiffMeans << endl;
                }
            }
            else if (tempString == "simCommand")
            {
                simCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: simCommand = " << simCommand << endl;
                }
            }
            else if (tempString == "bootstrapAlignments")
            {
                bootstrapAlignments = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapAlignments = " << bootstrapAlignments << endl;
                }
            }
            else if (tempString == "bootstrapCount")
            {
                bootstrapCount = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapCount = " << bootstrapCount << endl;
                }
            }
            else if (tempString == "bootstrapConcat")
            {
                bootstrapConcat = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapConcat = " << bootstrapConcat << endl;
                }
            }
            else if (tempString == "bootstrapColSize")
            {
                bootstrapColSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapColSize = " << bootstrapColSize << endl;
                }
            }
            else if (tempString == "v")
            {
                DEBUG_OUTPUT = 1;
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: DEBUG_OUTPUT = 1" << endl;
                }
            }
            else if (tempString == "vv")
            {
                DEBUG_OUTPUT = 2;
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: DEBUG_OUTPUT = 2" << endl;
                }
            }
            else if (tempString == "vvv")
            {
                DEBUG_OUTPUT = 3;
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: DEBUG_OUTPUT = 3" << endl;
                }
            }
            else if (tempString == "q")
            {
                DEBUG_OUTPUT = -1;
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: DEBUG_OUTPUT = -1" << endl;
                }
            }
            else {
                cout << "Parameter: " << tempString << " not understood." << endl;
                cout << usageMessage;
                exit(1);
            }
        }
    }
    if (inputFileNames.size() == 0){
        cout << usageMessage;
        cout << "No input files!" << endl;
        exit(0);
    }
    cout << endl;
    // Output all parameters.
    cout << " ***** Parameters ***** " << endl;
    cout << " Mr Bayes Parameters: " << endl;
    cout << "   doMB: " << doMB << endl;
    cout << "   MBP_nst: " << myMrBayesParam.MBP_nst << endl;
    cout << "   MBP_rates: " << myMrBayesParam.MBP_rates << endl;
    cout << "   MBP_nruns: " << myMrBayesParam.MBP_nruns << endl;
    cout << "   MBP_ngen: " << myMrBayesParam.MBP_ngen << endl;
    cout << "   MBP_sampleFreq: " << myMrBayesParam.MBP_sampleFreq << endl;
    cout << "   MBP_pathToMB: " << myMrBayesParam.MBP_pathToMB << endl;
    cout << "   MBP_saveOutput: " << myMrBayesParam.MBP_saveOutput << endl;
    cout << "   MPI_np: " << myMrBayesParam.MPI_np << endl;
    cout << " Jackknife Parameters: " << endl;
    cout << "   doJackknife: " << doJackknife << endl;
    cout << "   jackknifeColSize: " << myJackknifeParam.jackknifeColSize << endl;
    cout << "   jackknifeCount: " << myJackknifeParam.jackknifeCount << endl;
    cout << "   alignToTreeCommand: " << myJackknifeParam.alignToTreeCommand << endl;
    cout << "   burninFormat: " << mySampleParam.burninFormat << endl;
    cout << "   burninPercent: " << mySampleParam.burninPercent << endl;
    cout << "   burninNumber: " << mySampleParam.burninNumber << endl;
    cout << "   projectSVD " << mySampleParam.projectSVD << endl;
    cout << "   projectSVDcutOff " << mySampleParam.projectSVDcutOff << endl;
    cout << "   scaleToOne: " << mySampleParam.scaleToOne << endl;
    cout << "   distanceOne: " << mySampleParam.distanceOne << endl;
    cout << " Misc:" << endl;
    cout << "   doDiffMeans: " << mySampleParam.doDiffMeans << endl;
    cout << "   randGenSeed: " << randGenSeed << endl;
    cout << "   noTreeCalc: " << noTreeCalc << endl;
    cout << "   numTreesPerFile: " << mySampleParam.numTreesPerFile << endl;
    cout << "   numGroupOne: " << numGroupOne << endl;
    cout << "   concatGroups: " << concatGroups << endl;
    cout << "   bootstrapAlignments: " << bootstrapAlignments << endl;
    cout << "   bootstrapColSize: " << bootstrapColSize << endl;
    cout << "   bootstrapCount: " << bootstrapCount << endl;
    cout << "   bootstrapConcat: " << bootstrapConcat << endl;
    cout << "   numStatTests: " << numStatTests << endl;
    cout << "   doSim: " << doSim << endl;
    cout << "   simCommand: " << simCommand << endl;


    // Do some error checking
    if (doMB + doJackknife % 2 == 0){
        cout << "Must set either doMB or doJackknife, not both" << endl;
        exit (0);
    }
    if ((int)inputFileNames.size () < numGroupOne){
        cout << "Not enough files for numGroupOne" << endl;
        exit (0);
    }

    if (noTreeCalc == 1 && mySampleParam.numTreesPerFile == -1){
        cout << "mySampleParam.numTreesPerFile == -1 AND noTreeCalc == 1" << endl;
        exit(0);
    }

    // Error check for SVM vs number of trees that are going to be available.

    if (randGenSeed == 0){
        srand(time(0));
    }
    else {
        srand(time(&randGenSeed));
    }


    if (noTreeCalc != 1 && doMB == 1){
        mySampleParam.numTreesPerFile = (myMrBayesParam.MBP_ngen/myMrBayesParam.MBP_sampleFreq);
    }
    if (noTreeCalc != 1 && doJackknife == 1){
        mySampleParam.numTreesPerFile = myJackknifeParam.jackknifeCount;
    }

    int burninPerTreeFile; // This is the burnin number for each file.

    if (mySampleParam.burninFormat == -1){ // No burnin
        burninPerTreeFile = 0;
    }
    else if (mySampleParam.burninFormat == 0){ // Percentage burnin
        burninPerTreeFile = (int)floor((double)mySampleParam.numTreesPerFile*(mySampleParam.burninPercent));
    }
    else if (mySampleParam.burninFormat == 1){ // Fixed number burnin
        burninPerTreeFile = mySampleParam.burninNumber;
    }

    if (DEBUG_OUTPUT >= 0){
        cout << "Input files:" << endl;
        for (list<string>::const_iterator lit=inputFileNames.begin();lit!=inputFileNames.end();lit++){
             cout << *lit << endl; 
        }
        cout << endl;
    }

    // If we are working on alignments, we read them all in here.
    list <int> inputAlignmentsLengths;
    list <Alignment *> inputAlignments;
    list <Alignment *> AlignmentsOne, AlignmentsTwo;
    if (noTreeCalc != 1){
        if (concatGroups != 1) {

            for (list<string>::const_iterator lsit=inputFileNames.begin();lsit!=inputFileNames.end();lsit++){
                Alignment *newAlignments = new Alignment(*lsit);
                inputAlignmentsLengths.push_back(newAlignments->get_nchar()); 
                inputAlignments.push_back(newAlignments);
                // Here we will add extra bootstrap alignments if called for
                if (bootstrapAlignments == 1 && bootstrapConcat != 1){
                    for (int i=0;i<bootstrapCount;i++){
                        Alignment *tmpAlignments = new Alignment;
                        *tmpAlignments = newAlignments->getJackknife(bootstrapColSize);
                        inputAlignments.push_back(tmpAlignments);
                    }
                }
                //cout << *newAlignments;
            }
            // If we are adding bootstrapAlignments, then we need to increase numGroupOne 
            if (bootstrapAlignments == 1 && bootstrapConcat != 1){
                numGroupOne *= (bootstrapCount + 1);
            }

            Alignment *groupOneAlignmentsConcat;
            Alignment *groupTwoAlignmentsConcat;
            int firstGroupOne = 1;
            int firstGroupTwo = 1;

            // Again, default is to compare the first input vs the rest.
            // Lets us numGroupOne
            list <Alignment *>::iterator lait=inputAlignments.begin();
            int groupAddCount = 0;
            while (groupAddCount < numGroupOne){
                AlignmentsOne.push_back(*lait);
                if (firstGroupOne == 1){
                    firstGroupOne = 0;
                    groupOneAlignmentsConcat = new Alignment;
                    *groupOneAlignmentsConcat = *(*lait);
                }
                else {
                    *groupOneAlignmentsConcat += *(*lait);
                }
                lait++;
                groupAddCount++;
            }
            list <int>::iterator lit = inputAlignmentsLengths.begin();
            if (bootstrapAlignments == 1 && bootstrapConcat == 1){
                for (int j=0;j<numGroupOne;j++){
                    for (int i=0;i<bootstrapCount;i++){
                        Alignment *tmpAlignments = new Alignment;
                        *tmpAlignments = groupOneAlignmentsConcat->getJackknife(bootstrapColSize,(int)floor((double)*lit/bootstrapColSize));
                        AlignmentsOne.push_back(tmpAlignments);
                    }
                    lit++;
                }
                numGroupOne *= (bootstrapCount + 1);
            }

            // Put rest in AlignmentsTwo
            while(lait!=inputAlignments.end()){
                AlignmentsTwo.push_back(*lait);
                if (firstGroupTwo == 1){
                    firstGroupTwo = 0;
                    groupTwoAlignmentsConcat = new Alignment;
                    *groupTwoAlignmentsConcat = *(*lait);
                }
                else {
                    *groupTwoAlignmentsConcat += *(*lait);
                }
                lait++;
            }
            if (bootstrapAlignments == 1 && bootstrapConcat == 1){
                while (lit!=inputAlignmentsLengths.end()) {
                    for (int i=0;i<bootstrapCount;i++){
                        Alignment *tmpAlignments = new Alignment;
                        *tmpAlignments = groupTwoAlignmentsConcat->getJackknife(bootstrapColSize,(int)floor((double)*lit/bootstrapColSize));
                        AlignmentsTwo.push_back(tmpAlignments);
                    }
                    lit++;
                }
            }
            delete groupOneAlignmentsConcat;
            delete groupTwoAlignmentsConcat;
            if (bootstrapAlignments == 1 && bootstrapConcat == 1){
                inputAlignments.clear();
                for (list <Alignment *>::iterator lait=AlignmentsOne.begin();lait!=AlignmentsOne.end();lait++){
                    inputAlignments.push_back(*lait);
                }
                for (list <Alignment *>::iterator lait=AlignmentsTwo.begin();lait!=AlignmentsTwo.end();lait++){
                    inputAlignments.push_back(*lait);
                }
            }
        }
        else if (concatGroups == 1) {
            int fileCount = 0;
            Alignment *groupOneAlignmentsConcat;
            Alignment *groupTwoAlignmentsConcat;
            int firstGroupOne = 1;
            int firstGroupTwo = 1;
            for (list<string>::const_iterator lsit=inputFileNames.begin();lsit!=inputFileNames.end();lsit++){
                if (fileCount < numGroupOne) {
                    if (firstGroupOne == 1){
                        firstGroupOne = 0;
                        groupOneAlignmentsConcat = new Alignment(*lsit);
                    }
                    Alignment *newAlignments = new Alignment(*lsit);
                    *groupOneAlignmentsConcat += *newAlignments;
                    delete newAlignments;
                }
                else {
                    if (firstGroupTwo == 1){
                        firstGroupTwo = 0;
                        groupTwoAlignmentsConcat = new Alignment(*lsit);
                    }
                    Alignment *newAlignments = new Alignment(*lsit);
                    *groupTwoAlignmentsConcat += *newAlignments;
                    delete newAlignments;
                }
                fileCount++;
            }
            inputAlignmentsLengths.push_back(groupOneAlignmentsConcat->get_nchar()); 
            inputAlignmentsLengths.push_back(groupTwoAlignmentsConcat->get_nchar()); 
            inputAlignments.push_back(groupOneAlignmentsConcat);
            inputAlignments.push_back(groupTwoAlignmentsConcat);
            AlignmentsOne.push_back(groupOneAlignmentsConcat);
            AlignmentsTwo.push_back(groupTwoAlignmentsConcat);
        }
        if (DEBUG_OUTPUT >= 1){
            for (list <Alignment *>::iterator lait=inputAlignments.begin();lait!=inputAlignments.end();lait++){
                cout << *(*lait);
            }
        }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "inputAlignments.size() = " << inputAlignments.size() << endl;
        cout << "AlignmentsOne.size() = " << AlignmentsOne.size() << endl;
        cout << "AlignmentsTwo.size() = " << AlignmentsTwo.size() << endl;
    }

    int maxTreeIndex=0;
    // Compute the trees, and form up sets for groupOne and groupTwo
    if (noTreeCalc == 1){
        treeFileNames = inputFileNames;
        fillUnsignedSet(groupOneTrees,0,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne - 1);
        fillUnsignedSet(groupTwoTrees,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne,(mySampleParam.numTreesPerFile-burninPerTreeFile)*inputFileNames.size()-1);
        maxTreeIndex = (mySampleParam.numTreesPerFile-burninPerTreeFile)*inputFileNames.size()-1;
    }
    else if (noTreeCalc != 1){
        if (doMB == 1){
            MrBayesResults MB_results;
            getTreesMrBayes(myMrBayesParam,inputAlignments,TreeListListFileNames, MB_results);
            treeFileNames = listListStringToListString(TreeListListFileNames);

            fillUnsignedSet(groupOneTrees,0,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne*myMrBayesParam.MBP_nruns - 1);
            fillUnsignedSet(groupTwoTrees,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne*myMrBayesParam.MBP_nruns,(mySampleParam.numTreesPerFile-burninPerTreeFile)*inputAlignments.size()*myMrBayesParam.MBP_nruns-1);
            maxTreeIndex=(mySampleParam.numTreesPerFile-burninPerTreeFile)*inputAlignments.size()*myMrBayesParam.MBP_nruns-1;
        }
        else if (doJackknife == 1){
            for (int i=0;i<inputAlignments.size();i++){
                string baseName = "temptreesfile_ID";
                char tmpStr[80];
                sprintf(tmpStr,"%d%d%d_file_%d.nex",rand(),rand(),rand(),i);
                baseName.append(tmpStr);
                treeFileNames.push_back(baseName);
            }
            getTreesJackknife(myJackknifeParam,inputAlignments,treeFileNames);
            fillUnsignedSet(groupOneTrees,0,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne - 1);
            fillUnsignedSet(groupTwoTrees,(mySampleParam.numTreesPerFile-burninPerTreeFile)*numGroupOne,(mySampleParam.numTreesPerFile-burninPerTreeFile)*inputAlignments.size()-1);
            maxTreeIndex=(mySampleParam.numTreesPerFile-burninPerTreeFile)*inputAlignments.size()-1;
        }
    }

    if (DEBUG_OUTPUT >= 0){
        cout << "groupOneTrees.size() = " << groupOneTrees.size() << endl;
        cout << "groupTwoTrees.size() = " << groupTwoTrees.size() << endl;
        cout << "treeFileNames:" << endl;
        for (list <string>::iterator sit=treeFileNames.begin();sit!=treeFileNames.end();sit++){
            cout << "   " << *sit << endl;
        }
    }
    if (DEBUG_OUTPUT >= 2){
        for (set <unsigned>::const_iterator sit=groupOneTrees.begin();sit!=groupOneTrees.end();sit++){
            cout << *sit << endl;
        }
        cout << endl;
        for (set <unsigned>::const_iterator sit=groupTwoTrees.begin();sit!=groupTwoTrees.end();sit++){
            cout << *sit << endl;
        }
    }
    
    NewickWorker                *groupOneNewickWorker;
    NewickWorker                *groupTwoNewickWorker;
    NewickWorker_mean           *groupOneNewickWorker_mean;
    NewickWorker_mean           *groupTwoNewickWorker_mean;
    NewickWorker_mean_topology  *groupOneNewickWorker_mean_topology;
    NewickWorker_mean_topology  *groupTwoNewickWorker_mean_topology;
    if (mySampleParam.distanceOne == 1){
        if (DEBUG_OUTPUT >= 0){
            cout << "Selecting NewickWorker_mean_topology:" << endl;
        }
        groupOneNewickWorker_mean_topology = new NewickWorker_mean_topology;
        groupTwoNewickWorker_mean_topology = new NewickWorker_mean_topology; 
        groupOneNewickWorker = groupOneNewickWorker_mean_topology;
        groupTwoNewickWorker = groupTwoNewickWorker_mean_topology;
    }
    else{
        if (DEBUG_OUTPUT >= 0){
            cout << "Selecting NewickWorker_mean:" << endl;
        }
        groupOneNewickWorker_mean = new NewickWorker_mean;
        groupTwoNewickWorker_mean = new NewickWorker_mean; 
        groupOneNewickWorker = groupOneNewickWorker_mean;
        groupTwoNewickWorker = groupTwoNewickWorker_mean;
    }

    
    readNexusTreesUniform(treeFileNames, mySampleParam, groupOneTrees, *groupOneNewickWorker);
    readNexusTreesUniform(treeFileNames, mySampleParam, groupTwoTrees, *groupTwoNewickWorker);

    long double initDiffMeans;
    if (mySampleParam.distanceOne == 1){
        initDiffMeans = (groupOneNewickWorker_mean_topology->mean() - groupTwoNewickWorker_mean_topology->mean()).twoNorm();
        if (DEBUG_OUTPUT >= 0){
            cout << "l2 Mean group one: " << groupOneNewickWorker_mean_topology->l2norm() << endl;
            cout << "l2 Mean group two: " << groupTwoNewickWorker_mean_topology->l2norm() << endl;
            cout << "Difference of means: " << initDiffMeans << endl;
        }
    }
    else {
        initDiffMeans = (groupOneNewickWorker_mean->mean() - groupTwoNewickWorker_mean->mean()).twoNorm();
        if (DEBUG_OUTPUT >= 0){
            cout << "l2 Mean group one: " << groupOneNewickWorker_mean->l2norm() << endl;
            cout << "l2 Mean group two: " << groupTwoNewickWorker_mean->l2norm() << endl;
            cout << "Difference of means: " << initDiffMeans << endl;
        }
    }
    int origGroupOneTreesSize = groupOneTrees.size();
    int origGroupTwoTreesSize = groupTwoTrees.size();

    double pvalue = 0;
    for (int i=0;i<numStatTests;i++){
        groupOneTrees.clear();
        groupTwoTrees.clear();
        groupOneNewickWorker->clear();
        groupTwoNewickWorker->clear();

        while (groupOneTrees.size() != origGroupOneTreesSize){
            groupOneTrees.insert(rand() % (maxTreeIndex + 1));
        }
        if (DEBUG_OUTPUT >= 1){
            cout << "groupOneTrees.size() = " << groupOneTrees.size() << endl;
        }
        if (DEBUG_OUTPUT >= 2){
            for (set <unsigned>::const_iterator sit=groupOneTrees.begin();sit!=groupOneTrees.end();sit++){
                cout << *sit << endl;
            }
        }
        fillUnsignedSetExclude(groupTwoTrees,0,maxTreeIndex,groupOneTrees);

        if (DEBUG_OUTPUT >= 1){
            cout << "groupTwoTrees.size() = " << groupTwoTrees.size() << endl;
        }
        if (DEBUG_OUTPUT >= 2){
            for (set <unsigned>::const_iterator sit=groupTwoTrees.begin();sit!=groupTwoTrees.end();sit++){
                cout << *sit << endl;
            }
        }
        readNexusTreesUniform(treeFileNames, mySampleParam, groupOneTrees, *groupOneNewickWorker);
        readNexusTreesUniform(treeFileNames, mySampleParam, groupTwoTrees, *groupTwoNewickWorker);
        long double thisDiffMeans;
        if (mySampleParam.distanceOne == 1){
            thisDiffMeans = (groupOneNewickWorker_mean_topology->mean() - groupTwoNewickWorker_mean_topology->mean()).twoNorm();
            if (DEBUG_OUTPUT >= 0){
                //cout << "l2 Mean group one: " << groupOneNewickWorker_mean_topology->l2norm() << endl;
                //cout << "l2 Mean group two: " << groupTwoNewickWorker_mean_topology->l2norm() << endl;
                cout << "   This difference of means(" << i+1 << "): " << thisDiffMeans << endl;
            }
        }
        else {
            thisDiffMeans = (groupOneNewickWorker_mean->mean() - groupTwoNewickWorker_mean->mean()).twoNorm();
            if (DEBUG_OUTPUT >= 0){
                //cout << "l2 Mean group one: " << groupOneNewickWorker_mean->l2norm() << endl;
                //cout << "l2 Mean group two: " << groupTwoNewickWorker_mean->l2norm() << endl;
                cout << "   This difference of means(" << i+1 << "): " << thisDiffMeans << endl;
            }
        }
        if (thisDiffMeans > initDiffMeans){
            pvalue++;
            if (DEBUG_OUTPUT >= 0){
                cout << "       pvalue++"  << endl;
            }
        }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "pvalue/numStatTests " << pvalue/(long double)numStatTests << endl;
    }

    cout << endl;
    int mainEndTime = time(0);
    cout << "Total run time (wall clock): " << mainEndTime - mainStartTime << endl;
    cout << endl;
    for (list <Alignment *>::iterator lait=inputAlignments.begin();lait!=inputAlignments.end();lait++){
        delete *lait;
    }
}
