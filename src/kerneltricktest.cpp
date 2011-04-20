// $Rev: 402 $ $Date: 2009-11-19 23:10:48 -0500 (Thu, 19 Nov 2009) $

/// \file kerneltricktest.cpp 

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
#include "alignments.h"
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
    string usageMessage = "main <file1> <file2> ... <filen> -(options)\n\n \
          OPTIONS:\n \
          -samplePercent            <percent of trees to randomly sample>\n \
          -sampleNumber             <number of trees to randomly sample\n \
          -sampleSeq                <Sample sequential trees 0 to ...>\n \
          -doLHS                      <LHS only>\n \
          -doRHS                      <RHS only>\n \
          -burninPercent            <burn-in percentage>\n \
          -burninNumber             <fixed burn-in number>\n \
          -randGenSeed              <random generator seed>\n \
          -numTreesPerFile          <Number of trees per file is noTreeCalc == 1>\n \
          -distanceOne              <0 use newick distance, 1 use distance 1>\n \
          -numGroupOne              <Number of files in group one>\n \
          -scaleToOne               <0 no scaline, 1 scaleToOne\n";

            
          
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
        myMrBayesParam.MBP_ngen = 1000; // perform 10,000 runs.
        myMrBayesParam.MBP_sampleFreq = 1; // Sample every 10 chains.

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
        mySampleParam.burninFormat = -1;    // 0 will mean percentage, 1 will mean fixed number, -1 means no burnin.
        mySampleParam.burninPercent = 0; // Assume 25% burnin.
        mySampleParam.burninNumber = 0;    // Default is 0 since we assume default of 25% burnin


    // These are all variables used by the local program
    time_t      randGenSeed  = 0;
    int         doMB = 1; // By default, do mr bayes
    int         doJackknife = 0;      // Default is not to do jackknife
    int         noTreeCalc = 0; // Do not skip mr bayes calculation.
    int         concatGroups = 0;   // Default is 0, do not concat groups
    int         numGroupOne = 1; // First file is group one and thats all.
    int         statTest = 0; // 0 mean no stat test, 1 means stat test
    int         geneStatTest = 1; // Do the statistical test based on the hyp test on gene trees

    int         numInitCalc = 1; // 
    int         numStatTests = 10;
    int         skipInitCalc = 0; // Dont skip initial calculation for stat test
    double      mySeparation = -1;
    int         doSim = 0; 
    string      simCommand;
    string      outputFileName = "";
    int         sampleType = 0;  //This means sample all. If samplePercent,1. sampleNumber 2.
    int         sampleNumber = 0;
    int         samplePercent = 0;
    int         sampleSeq = 0;
    int         doLHS = 1;
    int         doRHS = 1;

    list <list <string> > TreeListListFileNames; // For getTreeMrBayes
    list <string> treeFileNamesGroupOne; // This should be the exact and literal filenames 
    list <string> treeFileNamesGroupTwo; // This should be the exact and literal filenames 

    list <svm_node *> *projectionMatrix; 
    svm_node    *empiricalMean;

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
            else if (tempString == "statTest")
            {
                statTest = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: statTest = " << statTest << endl;
                }
            }
            else if (tempString == "geneStatTest")
            {
                geneStatTest = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: geneStatTest = " << geneStatTest << endl;
                }
            }
            else if (tempString == "numStatTests")
            {
                numStatTests = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numStatTests = " << numStatTests << endl;
                }
            }
            else if (tempString == "numInitCalc")
            {
                numInitCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numInitCalc = " << numInitCalc << endl;
                }
            }
            else if (tempString == "skipInitCalc")
            {
                skipInitCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: skipInitCalc = " << skipInitCalc << endl;
                }
            }
            else if (tempString == "mySeparation")
            {
                mySeparation = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: mySeparation = " << mySeparation << endl;
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
            else if (tempString == "samplePercent")
            {
                samplePercent = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: samplePercent = " << samplePercent << endl;
                }
                sampleType = 1;
            }
            else if (tempString == "sampleNumber")
            {
                sampleNumber = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: sampleNumber = " << sampleNumber << endl;
                }
                sampleType = 2;
            }
            else if (tempString == "sampleSeq")
            {
                sampleSeq = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: sampleSeq = " << sampleSeq << endl;
                }
            }
            else if (tempString == "doLHS")
            {
                doLHS = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doLHS = " << doLHS << endl;
                }
                doLHS = 1;
                doRHS = 0;
            }
            else if (tempString == "doRHS")
            {
                doRHS = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doRHS = " << doRHS << endl;
                }
                doRHS = 1;
                doLHS = 0;
            }
            else if (tempString == "simCommand")
            {
                simCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: simCommand = " << simCommand << endl;
                }
            }
            else if (tempString == "outputFileName")
            {
                outputFileName = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: outputFileName = " << outputFileName << endl;
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
    cout << "   outputFileName: " << outputFileName << endl;
    cout << "   sampleType: " << sampleType << endl;
    cout << "   samplePercent: " << samplePercent << endl;
    cout << "   sampleNumber: " << sampleNumber << endl;
    cout << "   sampleSeq: " << sampleSeq << endl;
    cout << "   burninFormat: " << mySampleParam.burninFormat << endl;
    cout << "   burninPercent: " << mySampleParam.burninPercent << endl;
    cout << "   burninNumber: " << mySampleParam.burninNumber << endl;
    cout << "   doLHS: " << doLHS << endl;
    cout << "   doRHS: " << doRHS << endl;
    cout << "   projectSVD " << mySampleParam.projectSVD << endl;
    cout << "   projectSVDcutOff " << mySampleParam.projectSVDcutOff << endl;
    cout << "   scaleToOne: " << mySampleParam.scaleToOne << endl;
    cout << "   distanceOne: " << mySampleParam.distanceOne << endl;
    cout << "   randGenSeed: " << randGenSeed << endl;
    cout << "   numTreesPerFile: " << mySampleParam.numTreesPerFile << endl;


    // Do some error checking
    if (doMB + doJackknife % 2 == 0){
        cout << "Must set either doMB or doJackknife, not both" << endl;
        exit (0);
    }
    if ((int)inputFileNames.size () < numGroupOne){
        cout << "Not enough files for numGroupOne" << endl;
        exit (0);
    }
    if (skipInitCalc == 1 && mySeparation == -1){
        cout << "skipInitCalc selected but no mySeparation set." << endl;
        exit(0);
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

    if (DEBUG_OUTPUT >= 0){
        cout << "Input files:" << endl;
        for (list<string>::const_iterator lit=inputFileNames.begin();lit!=inputFileNames.end();lit++){
             cout << *lit << endl; 
        }
        cout << endl;
    }

    list <int>      numTreesPerFileGroupOne;
    list <int>      burninPerTreeFileGroupOne;
    list <int>      numTreesPerFileGroupTwo;
    list <int>      burninPerTreeFileGroupTwo;
    int             numTreesPostBurninGroupOne = 0;
    int             numTreesPostBurninGroupTwo = 0;
    list <string>   treesListGroupOne;
    list <string>   treesListGroupTwo;
    list <string>   treesList;
    int fileCount = 0;
    for (list<string>::const_iterator lit=inputFileNames.begin();lit!=inputFileNames.end();lit++){
        int tempInt = numTreesNexus(*lit);
        int tempBurnin;
        if (mySampleParam.burninFormat == -1){ // No burnin
            tempBurnin = 0;
        }
        else if (mySampleParam.burninFormat == 0){ // Percentage burnin
            tempBurnin = (int)floor((double)tempInt*(mySampleParam.burninPercent));
        }
        else if (mySampleParam.burninFormat == 1){ // Fixed number burnin
            tempBurnin = mySampleParam.burninNumber;
        }
        if (fileCount < numGroupOne){
            numTreesPerFileGroupOne.push_back(tempInt);
            numTreesPostBurninGroupOne += tempInt - tempBurnin;
            burninPerTreeFileGroupOne.push_back(tempBurnin);
            treesListGroupOne.push_back(*lit);
        }
        else {
            numTreesPerFileGroupTwo.push_back(tempInt);
            numTreesPostBurninGroupTwo += tempInt - tempBurnin;
            burninPerTreeFileGroupTwo.push_back(tempBurnin);
            treesListGroupTwo.push_back(*lit);
        }
        treesList.push_back(*lit);
        cout << "File " << *lit << " has " << tempInt << " trees." << endl;
        cout << "   Burnin = " << tempBurnin << endl;
        if (fileCount < numGroupOne){
            cout << "               Group one." << endl;
        } else {
            cout << "               Group two." << endl;
        }
        mySampleParam.numTreesPerFile = tempInt;
        fileCount++;
    }

    ostream *output;
    fstream outputFile;
    if (outputFileName != ""){
        outputFile.open(outputFileName.c_str(),fstream::out);
        output = &outputFile;
    }
    else {
        output = &cout;
    }

    NewickWorker *myNewickWorkerMeanGroupOne;
    NewickWorker *myNewickWorkerMeanGroupTwo;
    NewickWorker *myNewickWorker_l2normGroupOne_vs_GroupTwo;
    NewickWorker *myNewickWorker_l2normGroupTwo;
    NewickWorker *myNewickWorker_l2normGroupOne;

    NewickWorker_mean                   myNewickWorker_meanGroupOne;
    NewickWorker_mean                   myNewickWorker_meanGroupTwo;
    NewickWorker_mean_l2norm            myNewickWorker_mean_l2normGroupOne_vs_GroupTwo;
    NewickWorker_mean_l2norm            myNewickWorker_mean_l2normGroupOne;
    NewickWorker_mean_l2norm            myNewickWorker_mean_l2normGroupTwo;

    NewickWorker_mean_topology          myNewickWorker_mean_topologyGroupOne;
    NewickWorker_mean_topology          myNewickWorker_mean_topologyGroupTwo;
    NewickWorker_mean_topology_l2norm   myNewickWorker_mean_topology_l2normGroupOne_vs_GroupTwo;
    NewickWorker_mean_topology_l2norm   myNewickWorker_mean_topology_l2normGroupOne;
    NewickWorker_mean_topology_l2norm   myNewickWorker_mean_topology_l2normGroupTwo;

    if (mySampleParam.distanceOne == 1)
    {
        myNewickWorkerMeanGroupOne                     = &myNewickWorker_mean_topologyGroupOne;
        myNewickWorkerMeanGroupTwo                     = &myNewickWorker_mean_topologyGroupTwo;
        myNewickWorker_l2normGroupOne_vs_GroupTwo      = &myNewickWorker_mean_topology_l2normGroupOne_vs_GroupTwo;
        myNewickWorker_l2normGroupOne                  = &myNewickWorker_mean_topology_l2normGroupOne;
        myNewickWorker_l2normGroupTwo                  = &myNewickWorker_mean_topology_l2normGroupTwo;
    }
    else {
        myNewickWorkerMeanGroupOne                     = &myNewickWorker_meanGroupOne;
        myNewickWorkerMeanGroupTwo                     = &myNewickWorker_meanGroupTwo;
        myNewickWorker_l2normGroupOne_vs_GroupTwo      = &myNewickWorker_mean_l2normGroupOne_vs_GroupTwo;
        myNewickWorker_l2normGroupTwo                  = &myNewickWorker_mean_l2normGroupOne;
        myNewickWorker_l2normGroupOne                  = &myNewickWorker_mean_l2normGroupTwo;
    }

    int     sampleNumberGroupOne;
    int     sampleNumberGroupTwo;
    // Set sample number based off numTreesPostBurninGroup{One,Two}
    if (sampleType == 0){ //All
        sampleNumberGroupOne = numTreesPostBurninGroupOne;
        sampleNumberGroupTwo = numTreesPostBurninGroupTwo;
    }
    else if (sampleType == 1) { //Percent
        sampleNumberGroupOne = samplePercent*numTreesPostBurninGroupOne;
        sampleNumberGroupTwo = samplePercent*numTreesPostBurninGroupTwo;
    }
    else if (sampleType == 2) { //Number
        sampleNumberGroupOne = sampleNumber;
        sampleNumberGroupTwo = sampleNumber;
    }

    set <unsigned> treesToReadGroupOne;
    set <unsigned> treesToReadGroupTwo;
    if (sampleNumberGroupOne != numTreesPostBurninGroupOne){
        if (sampleSeq == 1) {
            fillUnsignedSet(treesToReadGroupOne,0,sampleNumberGroupOne - 1);
            fillUnsignedSet(treesToReadGroupTwo,numTreesPostBurninGroupOne,numTreesPostBurninGroupOne + sampleNumberGroupOne - 1);
        }
        else {
            while (treesToReadGroupOne.size() < sampleNumberGroupOne) {
                treesToReadGroupOne.insert(rand() % numTreesPostBurninGroupOne);
            }
            while (treesToReadGroupTwo.size() < sampleNumberGroupTwo) {
                treesToReadGroupTwo.insert(numTreesPostBurninGroupOne + (rand() % numTreesPostBurninGroupTwo));
            }
        }
    }
    else {
        fillUnsignedSet(treesToReadGroupOne,0,numTreesPostBurninGroupOne-1);
        fillUnsignedSet(treesToReadGroupTwo,numTreesPostBurninGroupOne,numTreesPostBurninGroupOne + numTreesPostBurninGroupTwo-1);
    }
    cout << "treesToReadGroupOne.size()  = " << treesToReadGroupOne.size() << endl;
    cout << "treesToReadGroupTwo.size()  = " << treesToReadGroupTwo.size() << endl;
    if (DEBUG_OUTPUT >= 2){
        for (set <unsigned>::iterator uit=treesToReadGroupOne.begin();uit!=treesToReadGroupOne.end();uit++){
            cout << setw(3) << *uit << " ";
        }
        cout << endl;
        for (set <unsigned>::iterator uit=treesToReadGroupTwo.begin();uit!=treesToReadGroupTwo.end();uit++){
            cout << setw(3) << *uit << " ";
        }
        cout << endl;
    }

    // LHS \| \mu_x - \mu_y \|
    if (doLHS == 1) {
        cout << "Computing LHS...";
        readNexusTreesUniform(treesList, mySampleParam, treesToReadGroupOne, *myNewickWorkerMeanGroupOne);
        readNexusTreesUniform(treesList, mySampleParam, treesToReadGroupTwo, *myNewickWorkerMeanGroupTwo);
        cout << "done." << endl;
    }

    if (doRHS == 1) {
        cout << "Computing RHS...";
        readNexusTreesUniformAllPairs(treesList, mySampleParam, treesToReadGroupOne, treesToReadGroupTwo, *myNewickWorker_l2normGroupOne_vs_GroupTwo);
        readNexusTreesUniformAllPairs(treesList, mySampleParam, treesToReadGroupOne, treesToReadGroupOne, *myNewickWorker_l2normGroupOne);
        readNexusTreesUniformAllPairs(treesList, mySampleParam, treesToReadGroupTwo, treesToReadGroupTwo, *myNewickWorker_l2normGroupTwo);
        cout << "done." << endl;
    }


    if (mySampleParam.distanceOne == 1){
        if (doLHS == 1) {
            cout << "LHS = || mu_one - mu_two ||^2 = ";// <<  " || " << myNewickWorker_mean_topologyGroupOne.mean() << " - " << myNewickWorker_mean_topologyGroupTwo.mean() << " || = ";
            cout << (myNewickWorker_mean_topologyGroupOne.mean() - myNewickWorker_mean_topologyGroupTwo.mean()).twoNormSquared() << endl;
        }
        if (doRHS == 1) {
            cout << "RHS = " << myNewickWorker_mean_topology_l2normGroupOne_vs_GroupTwo.mean() << " - (1/2)" << myNewickWorker_mean_topology_l2normGroupOne.mean() << " - (1/2)" << myNewickWorker_mean_topology_l2normGroupTwo.mean() << " = ";
            cout << myNewickWorker_mean_topology_l2normGroupOne_vs_GroupTwo.mean() - (1.0/2.0)*myNewickWorker_mean_topology_l2normGroupOne.mean() - (1.0/2.0)*myNewickWorker_mean_topology_l2normGroupTwo.mean() << endl;
        }
    }
    else {
        if (doLHS == 1) {
            cout << "LHS = || mu_one - mu_two ||^2 = ";// << "|| " << myNewickWorker_meanGroupOne.mean() << " - " << myNewickWorker_meanGroupTwo.mean() << " || = ";
            cout << (myNewickWorker_meanGroupOne.mean() - myNewickWorker_meanGroupTwo.mean()).twoNormSquared() << endl;
        }
        if (doRHS == 1) {
            cout << "RHS = " << myNewickWorker_mean_l2normGroupOne_vs_GroupTwo.mean() << " - (1/2)" << myNewickWorker_mean_l2normGroupOne.mean() << " - (1/2)" << myNewickWorker_mean_l2normGroupTwo.mean() << " = ";
            cout << myNewickWorker_mean_l2normGroupOne_vs_GroupTwo.mean() - (1.0/2.0)*myNewickWorker_mean_l2normGroupOne.mean() - (1.0/2.0)*myNewickWorker_mean_l2normGroupTwo.mean() << endl;
        }
    }


    cout << endl;
    int mainEndTime = time(0);
    cout << "Total run time (wall clock): " << mainEndTime - mainStartTime << endl;
    cout << endl;
}

