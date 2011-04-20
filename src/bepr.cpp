// $Rev: 402 $ $Date: 2009-11-19 23:10:48 -0500 (Thu, 19 Nov 2009) $

/// \file bepr.cpp 

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
          -burninPercent            <burn-in percentage>\n \
          -burninNumber             <fixed burn-in number>\n \
          -randGenSeed              <random generator seed>\n \
          -numTreesPerFile          <Number of trees per file is noTreeCalc == 1>\n \
          -numGroupOne              <Number of files in group one>\n \
          -outputFileName           <Name of output file, else stdout>\n \
          -distanceOne              <0 use newick distance, 1 use distance 1>\n";

            
          
    // Parse command line input.
    // Anything not preceded by '-' we will consider as an input file.
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
        mySampleParam.burninFormat = 0;    // 0 will mean percentage, 1 will mean fixed number.
        mySampleParam.burninPercent = 0.25; // Assume 25% burnin.
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
    fstream     *outputFile;
    ostream     *outStream;
    string      outputFileName = "";
    

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
    cout << "   burninFormat: " << mySampleParam.burninFormat << endl;
    cout << "   burninPercent: " << mySampleParam.burninPercent << endl;
    cout << "   burninNumber: " << mySampleParam.burninNumber << endl;
    cout << "   projectSVD " << mySampleParam.projectSVD << endl;
    cout << "   projectSVDcutOff " << mySampleParam.projectSVDcutOff << endl;
    cout << "   scaleToOne: " << mySampleParam.scaleToOne << endl;
    cout << "   distanceOne: " << mySampleParam.distanceOne << endl;
    cout << "   numTreesPerFile: " << mySampleParam.numTreesPerFile << endl;
    cout << "   numGroupOne: " << numGroupOne << endl;
    cout << "   outputFileName: " << outputFileName << endl;

    if (outputFileName == "") {
        outStream = &cout;
    }
    else {
        outputFile = new fstream;
        outputFile->open(outputFileName.c_str(),fstream::out);
        outStream = outputFile;
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

    // If N input files are given, then this program treats the first N-1 files 
    // as nexus files which each contains numTreesPerFile. The last file
    // is a nexus file which holds one tree.

    list <string> treeFileNames;
    list <string> origTreeName;
    list <string>::iterator sit = inputFileNames.begin();
    for (int i=0;i< inputFileNames.size() - 1;i++){
        treeFileNames.push_back(*sit);
        sit++;
    }
    origTreeName.push_back(*sit);

    NewickWorker_mean_topology topologicalMean;
    NewickWorker_mean_topology origTree;
    
    set <unsigned> treesToRead;
    set <unsigned> readOrigTree;
    // All trees (minus burnin) are to be read
    fillUnsignedSet(treesToRead,0,(mySampleParam.numTreesPerFile-burninPerTreeFile)*(inputFileNames.size() - 1) - 1);
    if (DEBUG_OUTPUT >= 0) {
        cout << "Number of trees to read total: " << treesToRead.size() << endl;
    }
    readOrigTree.insert(0);

    readNexusTreesUniform(treeFileNames, mySampleParam, treesToRead, topologicalMean);
    if (DEBUG_OUTPUT >= 1) {
        cout << "Topological mean: " << endl << topologicalMean.mean();
    }

    SampleParameters tsp = mySampleParam;
    tsp.burninFormat = 1; 
    tsp.burninNumber = 0;

    readNexusTreesUniform(origTreeName, tsp, readOrigTree, origTree);
    if (DEBUG_OUTPUT >= 1) {
        cout << "Original tree: " << endl << origTree.mean();
    }

    NewickWorker_print_l2_diff_topology diffWorker;
    diffWorker.init(origTree.mean(), topologicalMean.mean(), *outStream);
    *outStream << "Orig Mean" << endl;
    readNexusTreesUniform(treeFileNames, mySampleParam, treesToRead, diffWorker);

    int mainEndTime = time(0);
    cout << "Total run time (wall clock): " << mainEndTime - mainStartTime << endl;
    cout << endl;

}


