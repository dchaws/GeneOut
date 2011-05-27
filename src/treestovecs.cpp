// $Rev: 808 $ $Date: 2011-04-10 13:32:07 -0400 (Sun, 10 Apr 2011) $

/// \file treestovecs.cpp 

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
    string usageMessage = "main <file1> <file2> ... <filen> -(options)\n\n \
          OPTIONS:\n \
          -outputFileName           <file to output trees to>\n \
          -samplePercent            <percent of trees to randomly sample>\n \
          -sampleNumber             <number of trees to randomly sample\n \
          -burninPercent            <burn-in percentage>\n \
          -burninNumber             <fixed burn-in number>\n \
          -randGenSeed              <random generator seed>\n \
          -numTreesPerFile          <Number of trees per file is noTreeCalc == 1>\n \
          -projectSVD               <1 project the points via SVD, 0 not>\n \
          -projectSVDcutOff         <cutOff percentage for singular values>\n \
          -projectSVDnumVecs        <number of singular vectors to take>\n \
          -distanceOne              <0 use newick distance, 1 use distance 1>\n \
          -scaleToOne               <0 no scaline, 1 scaleToOne\n \
          -printSplits              <0 no, 1 print splits>\n";

            
          
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

    BootstrapParameters myBootstrapParam;
        myBootstrapParam.bootstrapColSize = 1; // Default
        myBootstrapParam.bootstrapCount = 1000; // Default
        myBootstrapParam.alignToTreeCommand = "./run_dnadist_neighbor_mult"; // Default
    
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
    int         doBootstrap = 0;      // Default is not to do bootstrap
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
    int         projectSVDsampleSize = 100;
    int         projectSVDnumVecs = 0;
    int         printSplits = 0;

    list <list <string> > TreeListListFileNames; // For getTreeMrBayes
    list <string> treeFileNamesGroupOne; // This should be the exact and literal filenames 
    list <string> treeFileNamesGroupTwo; // This should be the exact and literal filenames 

    list <svm_node *> *projectionMatrix; 
    svm_node    *empiricalMean;

    if (DEBUG_OUTPUT >= 1){
        cout << " ***** Reading in file names and command-line arguments *****" << endl;
    }
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
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: doMB = " << doMB << endl;
                }
            }
            else if (tempString == "MBP_nst")
            {
                myMrBayesParam.MBP_nst = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_nst = " << myMrBayesParam.MBP_nst << endl;
                }
            }
            else if (tempString == "MBP_rates")
            {
                myMrBayesParam.MBP_rates = argv[++i];
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_rates = " << myMrBayesParam.MBP_rates << endl;
                }
            }
            else if (tempString == "MBP_nruns")
            {
                myMrBayesParam.MBP_nruns = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_nruns = " << myMrBayesParam.MBP_nruns << endl;
                }
            }
            else if (tempString == "MBP_ngen")
            {
                myMrBayesParam.MBP_ngen = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_ngen = " << myMrBayesParam.MBP_ngen << endl;
                }
            }
            else if (tempString == "MBP_sampleFreq")
            {
                myMrBayesParam.MBP_sampleFreq = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_sampleFreq = " << myMrBayesParam.MBP_sampleFreq << endl;
                }
            }
            else if (tempString == "MBP_pathToMB")
            {
                myMrBayesParam.MBP_pathToMB = argv[++i];
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_pathToMB = " << myMrBayesParam.MBP_pathToMB << endl;
                }
            }
            else if (tempString == "MBP_saveOutput")
            {
                myMrBayesParam.MBP_saveOutput = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MBP_saveOutput = " << myMrBayesParam.MBP_saveOutput << endl;
                }
            }
            else if (tempString == "MPI_np")
            {
                myMrBayesParam.MPI_np = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: MPI_np = " << myMrBayesParam.MPI_np << endl;
                }
            }
            else if (tempString == "burninPercent")
            {
                mySampleParam.burninFormat = 0;
                sscanf(argv[++i],"%lf",&mySampleParam.burninPercent);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: burninFormat = " << mySampleParam.burninFormat << endl;
                    cout << "       Setting: burninPercent = " << mySampleParam.burninPercent << endl;
                }
            }
            else if (tempString == "burninNumber")
            {
                mySampleParam.burninFormat = 1;
                mySampleParam.burninNumber = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: burninFormat = " << mySampleParam.burninFormat << endl;
                    cout << "       Setting: burninNumber = " << mySampleParam.burninNumber << endl;
                }
            }
            else if (tempString == "randGenSeed")
            {
                randGenSeed = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: randGenSeed = " << randGenSeed << endl;
                }
            }
            else if (tempString == "doBootstrap")
            {
                doBootstrap = atoi(argv[++i]); 
                if (doBootstrap == 1){
                    doMB = 0;
                }
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: doBootstrap = " << doBootstrap << endl;
                }
                mySampleParam.burninFormat = 1;
                mySampleParam.burninNumber = 0;
            }
            else if (tempString == "bootstrapColSize")
            {
                myBootstrapParam.bootstrapColSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: bootstrapColSize = " << myBootstrapParam.bootstrapColSize << endl;
                }
            }
            else if (tempString == "bootstrapCount")
            {
                myBootstrapParam.bootstrapCount = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: bootstrapCount = " << myBootstrapParam.bootstrapCount << endl;
                }
            }
            else if (tempString == "alignToTreeCommand")
            {
                myBootstrapParam.alignToTreeCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: alignToTreeCommand = " << myBootstrapParam.alignToTreeCommand << endl;
                }
            }
            else if (tempString == "SVM_sampleSize")
            {
                mySampleParam.SVM_sampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: SVM_sampleSize = " << mySampleParam.SVM_sampleSize << endl;
                }
            }
            else if (tempString == "SVM_resampleSize")
            {
                mySampleParam.SVM_resampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: SVM_resampleSize = " << mySampleParam.SVM_resampleSize << endl;
                }
            }
            else if (tempString == "noTreeCalc")
            {
                noTreeCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: noTreeCalc = " << noTreeCalc << endl;
                    //cout << "       Setting: groupOne
                }
            }
            else if (tempString == "numTreesPerFile")
            {
                mySampleParam.numTreesPerFile = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: numTreesPerFile = " << mySampleParam.numTreesPerFile << endl;
                }
            }
            else if (tempString == "numGroupOne")
            {
                numGroupOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: numGroupOne = " << numGroupOne << endl;
                }
            }
            else if (tempString == "projectSVDnumVecs")
            {
                projectSVDnumVecs = atoi(argv[++i]); 
                mySampleParam.projectSVD = 1;
                mySampleParam.projectSVDcutOff = -1;
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: projectSVDnumVecs = " << projectSVDnumVecs << endl;
                }
            }
            else if (tempString == "concatGroups")
            {
                concatGroups = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: concatGroups = " << concatGroups << endl;
                }
            }
            else if (tempString == "projectSVD")
            {
                mySampleParam.projectSVD = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: projectSVD = " << mySampleParam.projectSVD << endl;
                }
            }
            else if (tempString == "projectSVDcutOff")
            {
                mySampleParam.projectSVDcutOff = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: projectSVDcutOff = " << mySampleParam.projectSVDcutOff << endl;
                }
            }
            else if (tempString == "distanceOne")
            {
                mySampleParam.distanceOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: distanceOne = " << mySampleParam.distanceOne << endl;
                }
            }
            else if (tempString == "scaleToOne")
            {
                mySampleParam.scaleToOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: scaleToOne = " << mySampleParam.scaleToOne << endl;
                }
            }
            else if (tempString == "statTest")
            {
                statTest = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: statTest = " << statTest << endl;
                }
            }
            else if (tempString == "geneStatTest")
            {
                geneStatTest = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: geneStatTest = " << geneStatTest << endl;
                }
            }
            else if (tempString == "numStatTests")
            {
                numStatTests = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: numStatTests = " << numStatTests << endl;
                }
            }
            else if (tempString == "numInitCalc")
            {
                numInitCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: numInitCalc = " << numInitCalc << endl;
                }
            }
            else if (tempString == "projectSVDsampleSize")
            {
                projectSVDsampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: projectSVDsampleSize = " << projectSVDsampleSize << endl;
                }
            }
            else if (tempString == "skipInitCalc")
            {
                skipInitCalc = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: skipInitCalc = " << skipInitCalc << endl;
                }
            }
            else if (tempString == "mySeparation")
            {
                mySeparation = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: mySeparation = " << mySeparation << endl;
                }
            }
            else if (tempString == "doSim")
            {
                doSim = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: doSim = " << doSim << endl;
                }
            }
            else if (tempString == "doSVM")
            {
                mySampleParam.doSVM = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: doSVM = " << mySampleParam.doSVM << endl;
                }
            }
            else if (tempString == "doDiffMeans")
            {
                mySampleParam.doDiffMeans = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: doDiffMeans = " << mySampleParam.doDiffMeans << endl;
                }
            }
            else if (tempString == "samplePercent")
            {
                samplePercent = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: samplePercent = " << samplePercent << endl;
                }
                sampleType = 1;
            }
            else if (tempString == "sampleNumber")
            {
                sampleNumber = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: sampleNumber = " << sampleNumber << endl;
                }
                sampleType = 2;
            }
            else if (tempString == "simCommand")
            {
                simCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: simCommand = " << simCommand << endl;
                }
            }
            else if (tempString == "outputFileName")
            {
                outputFileName = argv[++i]; 
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: outputFileName = " << outputFileName << endl;
                }
            }
            else if (tempString == "v")
            {
                DEBUG_OUTPUT = 1;
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: DEBUG_OUTPUT = 1" << endl;
                }
            }
            else if (tempString == "vv")
            {
                DEBUG_OUTPUT = 2;
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: DEBUG_OUTPUT = 2" << endl;
                }
            }
            else if (tempString == "vvv")
            {
                DEBUG_OUTPUT = 3;
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: DEBUG_OUTPUT = 3" << endl;
                }
            }
            else if (tempString == "q")
            {
                DEBUG_OUTPUT = -1;
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: DEBUG_OUTPUT = -1" << endl;
                }
            }
            else if (tempString == "printSplits")
            {
                printSplits = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 1){
                    cout << "       Setting: printSplits = " << printSplits << endl;
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

    if (DEBUG_OUTPUT >= 1){
        cout << endl;
        // Output all parameters.
        cout << " ***** Parameters ***** " << endl;
        cout << "   outputFileName: " << outputFileName << endl;
        cout << "   sampleType: " << sampleType << endl;
        cout << "   samplePercent: " << samplePercent << endl;
        cout << "   sampleNumber: " << sampleNumber << endl;
        cout << "   burninFormat: " << mySampleParam.burninFormat << endl;
        cout << "   burninPercent: " << mySampleParam.burninPercent << endl;
        cout << "   burninNumber: " << mySampleParam.burninNumber << endl;
        cout << "   projectSVD " << mySampleParam.projectSVD << endl;
        cout << "   projectSVDcutOff " << mySampleParam.projectSVDcutOff << endl;
        cout << "   projectSVDsampleSize " << projectSVDsampleSize << endl;
        cout << "   scaleToOne: " << mySampleParam.scaleToOne << endl;
        cout << "   distanceOne: " << mySampleParam.distanceOne << endl;
        cout << "   randGenSeed: " << randGenSeed << endl;
        cout << "   numTreesPerFile: " << mySampleParam.numTreesPerFile << endl;
        cout << "   printSplits: " << printSplits << endl;
    }

    // Do some error checking
    if (doMB + doBootstrap % 2 == 0){
        cout << "Must set either doMB or doBootstrap, not both" << endl;
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



    if (DEBUG_OUTPUT >= 1){
        cout << "Input files:" << endl;
        for (list<string>::const_iterator lit=inputFileNames.begin();lit!=inputFileNames.end();lit++){
             cout << *lit << endl; 
        }
        cout << endl;
    }

    list <int> numTreesPerFile;
    list <int> burninPerTreeFile;
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
        numTreesPerFile.push_back(tempInt);
        burninPerTreeFile.push_back(tempBurnin);
        if (DEBUG_OUTPUT >= 1){
            cout << "File " << *lit << " has " << tempInt << " trees." << endl;
            cout << "   Burnin = " << tempBurnin << endl;
        }
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

    NewickWorker *myNewickWorker;

    NewickWorker_print myNewickWorker_print;
    myNewickWorker_print.init (*output,0);

    NewickWorker_print_topology myNewickWorker_print_topology;
    myNewickWorker_print_topology.init (*output,0);

    NewickWorker_getSVMList myNewickWorker_getSVMList;
    NewickWorker_getSVMList_topology myNewickWorker_getSVMList_topology;
    
    NewickWorker_print_splits myNewickWorker_print_splits;
    myNewickWorker_print_splits.init (*output,0);

    if (printSplits == 1)
    {
            myNewickWorker = &myNewickWorker_print_splits;
    }
    else{
        if (mySampleParam.projectSVD == 1) {
            if (mySampleParam.distanceOne == 1)
            {
                myNewickWorker = &myNewickWorker_getSVMList_topology;
            }
            else {
                myNewickWorker = &myNewickWorker_getSVMList;
            }
        }
        else {
            if (mySampleParam.distanceOne == 1)
            {
                myNewickWorker = &myNewickWorker_print_topology;
            }
            else {
                myNewickWorker = &myNewickWorker_print;
            }
        }
    }

    list <int>::const_iterator numTrees_it = numTreesPerFile.begin();
    list <int>::const_iterator burnin_it = burninPerTreeFile.begin();
    for (list<string>::const_iterator lit=inputFileNames.begin();lit!=inputFileNames.end();lit++){
        list <string> tempStringList;
        list <svm_node * > groupOne;
        list <svm_node * > groupTwo;
        list <svm_node *> *projectionMatrix; 
        svm_node    *empiricalMean;
        set <unsigned> treesToRead;
        tempStringList.clear();
        tempStringList.push_back(*lit);
        treesToRead.clear();
        groupOne.clear();
        groupTwo.clear();

        if (mySampleParam.projectSVD == 1) {
            // First sample to do pca
            while (treesToRead.size() < projectSVDsampleSize) {
                treesToRead.insert(rand() % (*numTrees_it - *burnin_it ));
            }
            readNexusTreesUniform(tempStringList, mySampleParam, treesToRead, *myNewickWorker);
            treesToRead.clear();

            if (mySampleParam.distanceOne == 1)
            {
                groupOne = myNewickWorker_getSVMList_topology.getNodeList();
                myNewickWorker_getSVMList_topology.clear();
            }
            else {
                groupOne = myNewickWorker_getSVMList.getNodeList();
                myNewickWorker_getSVMList.clear();
            }
            groupTwo.clear();
            getPCA_data(groupOne, groupTwo, mySampleParam.projectSVDcutOff, &projectionMatrix, &empiricalMean,projectSVDnumVecs);
            for (list <svm_node * >::iterator svmit=projectionMatrix->begin();svmit!=projectionMatrix->end();svmit++){
                cout << "PROJECTION MATRIX: ";
                int curIndex=0;
                while ((*svmit)[curIndex].index != -1) {
                    cout << (*svmit)[curIndex].value << " "; 
                    curIndex++;
                }
                cout << endl;
            }
            cout << "EMPIRICAL MEAN: ";
            int curIndexEmp=0;
            while (empiricalMean[curIndexEmp].index != -1) {
                cout << empiricalMean[curIndexEmp].value << " "; 
                curIndexEmp++;
            }
            cout << endl;
        }
        if (sampleType == 0){
            fillUnsignedSet(treesToRead,0,*numTrees_it - *burnin_it - 1);
        }
        else if (sampleType == 1){ // Percent
            int numTreesToRead = (*numTrees_it - *burnin_it)*samplePercent;
            while (treesToRead.size() != numTreesToRead){
                treesToRead.insert(rand() % (*numTrees_it - *burnin_it));
            }
        }
        else if (sampleType == 2){ // Number
            while (treesToRead.size() != sampleNumber){
                treesToRead.insert(rand() % (*numTrees_it - *burnin_it));
            }
        }
        readNexusTreesUniform(tempStringList, mySampleParam, treesToRead, *myNewickWorker);
        // Apply PCA and print out svm list.
        if (mySampleParam.projectSVD == 1) {
            if (mySampleParam.distanceOne == 1)
            {
                groupOne = myNewickWorker_getSVMList_topology.getNodeList();
                groupTwo.clear();
                doPCA(groupOne, groupTwo, projectionMatrix, empiricalMean);
            }
            else {
                groupOne = myNewickWorker_getSVMList.getNodeList();
                groupTwo.clear();
                doPCA(groupOne, groupTwo, projectionMatrix, empiricalMean);
            }
            //delete projectionMatrix;
            for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
                int curIndex=0;
                while ((*svmit)[curIndex].index != -1) {
                    cout << (*svmit)[curIndex].value << " "; 
                    curIndex++;
                }
                cout << endl;
            }
        }
        
        numTrees_it++;
        burnin_it++;
    }

}

