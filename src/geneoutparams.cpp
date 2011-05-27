    // $Rev: 767 $ $Date: 2010-10-10 15:02:33 -0400 (Sun, 10 Oct 2010) $

/** \file geneoutparams.cpp */
#include "geneoutparams.h"

GeneOutParameters::GeneOutParameters()
{
    init();
}

GeneOutParameters::GeneOutParameters(string fileName)
{
    init();
    fstream inputFile;
    inputFile.open(fileName.c_str(),fstream::in);
    if (inputFile.bad()) {
        cout << "Could not open file " << fileName << endl;
        exit(0);
    }
}

GeneOutParameters::GeneOutParameters(std::istream &in)
{
    init();
    in >> *this;
}

std::ostream& operator << (std::ostream &out, const GeneOutParameters &sGOP)
{
    //cout << "GeneOut Class(" << &sGOP << ")" << endl;
    cout << "Group one files: ";
    for (list <string>::const_iterator lsit=sGOP.inputFileNamesGroupOne.begin();lsit!=sGOP.inputFileNamesGroupOne.end();lsit++){
        cout << *lsit << " ";
    }
    cout << endl;
    cout << "Group two files: ";
    for (list <string>::const_iterator lsit=sGOP.inputFileNamesGroupTwo.begin();lsit!=sGOP.inputFileNamesGroupTwo.end();lsit++){
        cout << *lsit << " ";
    }
    cout << endl;

    if (DEBUG_OUTPUT >= 1)
    {
        cout << "MrBayesParameters: " << endl;
        cout << "   inputNexFiles: " << endl;
        for (list <string>::const_iterator lsit=sGOP.MrBayesParam.inputNexFiles.begin();lsit!=sGOP.MrBayesParam.inputNexFiles.end();lsit++){
            cout << *lsit << " ";
        }
        cout << endl;
        cout << "   MBP_nst: " << sGOP.MrBayesParam.MBP_nst << endl;
        cout << "   MBP_rates: " << sGOP.MrBayesParam.MBP_rates << endl;
        cout << "   MBP_nruns: " << sGOP.MrBayesParam.MBP_nruns << endl;
        cout << "   MBP_ngen: " << sGOP.MrBayesParam.MBP_ngen << endl;
        cout << "   MBP_sampleFreq: " << sGOP.MrBayesParam.MBP_sampleFreq << endl;
        cout << "   MBP_pathToMB: " << sGOP.MrBayesParam.MBP_pathToMB << endl;
        cout << "   MBP_saveOutput: " << sGOP.MrBayesParam.MBP_saveOutput << endl;
        cout << "   MPI_np: " << sGOP.MrBayesParam.MPI_np << endl;
        cout << "   parametersBlockFileName: " << sGOP.MrBayesParam.parametersBlockFileName << endl;
        cout << "BootstrapParameters: " << endl;
        cout << "   inputNexFiles: " << endl;
        for (list <string>::const_iterator lsit=sGOP.BootstrapParam.inputNexFiles.begin();lsit!=sGOP.BootstrapParam.inputNexFiles.end();lsit++){
            cout << *lsit << " ";
        }
        cout << endl;
        cout << "   bootstrapColSize: " << sGOP.BootstrapParam.bootstrapColSize << endl;
        cout << "   bootstrapCount: " << sGOP.BootstrapParam.bootstrapCount << endl;
        cout << "   alignToTreeCommand: " << sGOP.BootstrapParam.alignToTreeCommand << endl;
        cout << "SampleParameters: " << endl;
        cout << "   doSVM: " << sGOP.SampleParam.doSVM << endl;
        cout << "   doDiffMeans: " << sGOP.SampleParam.doDiffMeans << endl;
        cout << "   SVM_sampleUniform: " << sGOP.SampleParam.SVM_sampleUniform << endl;
        cout << "   SVM_sampleSize: " << sGOP.SampleParam.SVM_sampleSize << endl;
        cout << "   SVM_resampleSize: " << sGOP.SampleParam.SVM_resampleSize << endl;
        cout << "   projectSVD: " << sGOP.SampleParam.projectSVD << endl;
        cout << "   projectSVDcutOff: " << sGOP.SampleParam.projectSVDcutOff << endl;
        cout << "   distanceOne: " << sGOP.SampleParam.distanceOne << endl;
        cout << "   scaleToOne: " << sGOP.SampleParam.scaleToOne << endl;
        cout << "   modelType: " << sGOP.SampleParam.modelType << endl;
        cout << "   numTreesPerFile: " << sGOP.SampleParam.numTreesPerFile << endl;
        cout << "   burninFormat: " << sGOP.SampleParam.burninFormat << endl;
        cout << "   burninPercent: " << sGOP.SampleParam.burninPercent << endl;
        cout << "   burninNumber: " << sGOP.SampleParam.burninNumber << endl;

        cout << "randGenSeed : " << sGOP.randGenSeed  << endl;
        cout << "doMB        : " << sGOP.doMB         << endl;
        cout << "doBootstrap : " << sGOP.doBootstrap  << endl;
        cout << "doMultInd   : " << sGOP.doMultInd    << endl;
        cout << "numInd      : " << sGOP.numInd       << endl;
        cout << "noTreeCalc  : " << sGOP.noTreeCalc   << endl;
        cout << "concatGroups: " << sGOP.concatGroups << endl;
        cout << "numGroupOne : " << sGOP.numGroupOne  << endl;
        cout << "statTest    : " << sGOP.statTest     << endl;
        cout << "geneStatTest: " << sGOP.geneStatTest << endl;
        cout << "numInitCalc : " << sGOP.numInitCalc  << endl;
        cout << "numStatTests: " << sGOP.numStatTests << endl;
        cout << "skipInitCalc: " << sGOP.skipInitCalc << endl;
        cout << "mySeparation  : " << sGOP.mySeparation   << endl;
        cout << "doSim       : " << sGOP.doSim        << endl;
        cout << "simCommand  : " << sGOP.simCommand   << endl;
        cout << "multSVMSepMethod  : " << sGOP.multSVMSepMethod   << endl;
        cout << "testType    : " << sGOP.testType << endl;
        cout << "tempPrefix  : " << sGOP.tempPrefix << endl;
        cout << "permuteOrig : " << sGOP.permuteOrig << endl;
        cout << "allowAnyPermuation : " << sGOP.allowAnyPermuation << endl;
        cout << "sGOP.BSStepOne : " << sGOP.BSStepOne << endl;
    }

}

std::istream& operator >> (std::istream &in, GeneOutParameters &sGOP)
{
    // If we are going to read in a new alignment, delete previous alignments
    // and clear everything out.
    if (in.bad()){
        cout << "GeneOutParameters friend operator >>, input bad bit set" << endl;
        exit (0);
    }

    string input;
    // First look for '#NEXUS'
    input = nextNexusToken(in);
    if (input != "#NEXUS"){
        if (DEBUG_OUTPUT >= 0){
            cout << "Warning: GeneOutParameters friend operator >>, '#NEXUS' not first line" << endl;
        }
        // Lets not be so strict about it.
        // Seek to begining
        in.seekg(0,ios::beg);
        //exit(0);
    }
    else {
        if (DEBUG_OUTPUT >= 1) {
            cout << "GeneOutParameters friend operator >>, '#NEXUS' found" << endl;
        }
    }


    while (in.good()){
        input = nextNexusToken(in);
        transform(input.begin(), input.end(),input.begin(), ::toupper);
        if (input == "BEGIN"){
            input = nextNexusTokenUpper(in);
            if (input == "SINGLESVMSEP" || input == "MULTSVMSEP"){
                if (input == "SINGLESVMSEP"){
                    sGOP.testType = 1; //single svm sep 
                }
                if (input == "MULTSVMSEP"){
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "Setting \"MULTSVMSEP\" method." << endl;
                    }
                    sGOP.testType = 2; //multiple svm sep 
                    sGOP.statTest = 1;
                }
                readUntilSemiColonNexus(in);
                input = nextNexusTokenUpper(in);
                while (input != "END" && input != "ENDBLOCK"){
                    if (input == "DISTMETHOD"){
                        input = nextNexusTokenUpper(in);
                        if (input == "MB"){
                            sGOP.doMB = 1;
                            sGOP.doBootstrap = 0;
                            sGOP.noTreeCalc = 0;
                        }
                        else if (input == "JK" || input == "BS"){
                            sGOP.doMB = 0;
                            sGOP.doBootstrap = 1;
                            sGOP.noTreeCalc = 0;
                        }
                        else if (input == "TREES"){
                            sGOP.doMB = 0;
                            sGOP.doBootstrap = 0;
                            sGOP.noTreeCalc = 1;
                        }
                        readUntilSemiColonNexus(in);
                    }
                    if (input == "MBP"){
                        input = nextNexusTokenUpper(in);
                        while (input != ";"){
                            if (input == "MBP_NST"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.MrBayesParam.MBP_nst = atoi(input.c_str());
                            }
                            else if (input == "MBP_RATES"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                //input = nextNexusNumber(in);
                                input = nextQuoteBlock(in);
                                sGOP.MrBayesParam.MBP_rates = input;
                            }
                            else if (input == "MBP_NRUNS"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.MrBayesParam.MBP_nruns = atoi(input.c_str());
                            }
                            else if (input == "MBP_NGEN"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.MrBayesParam.MBP_ngen = atoi(input.c_str());
                            }
                            else if (input == "MBP_SAMPLEFREQ"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.MrBayesParam.MBP_sampleFreq = atoi(input.c_str());
                            }
                            else if (input == "MBP_PATHTOMB"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                // Read everything in between quotes.
                                input = nextQuoteBlock(in);
                                sGOP.MrBayesParam.MBP_pathToMB = atoi(input.c_str());
                            }
                            else if (input == "MBP_SAVEOUTPUT"){
                                if (DEBUG_OUTPUT >= 1) {
                                    cout << "   MBP_SAVEOUTPUT detected." << endl;
                                }
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusTokenUpper(in);
                                if (input == "YES") {
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "   MBP_SAVEOUTPUT \"YES\" detected." << endl;
                                    }
                                    sGOP.MrBayesParam.MBP_saveOutput = 1;
                                }
                                else if (input == "NO") {
                                    sGOP.MrBayesParam.MBP_saveOutput = 0;
                                }
                            }
                            else if (input == "MPI_NP"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.MrBayesParam.MPI_np = atoi(input.c_str());
                            }
                            else if (input == "MBP_PARAMFILE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                // Read everything in between quotes.
                                input = nextQuoteBlock(in);
                                sGOP.MrBayesParam.parametersBlockFileName = input;
                            }
                            input = nextNexusTokenUpper(in);
                        }
                    }
                    if (input == "JKP" || input == "BSP" ){ // JKP is for legacy support.
                        input = nextNexusTokenUpper(in);
                        while (input != ";"){
                            if (input == "COLSIZE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.BootstrapParam.bootstrapColSize = atoi(input.c_str());
                            }
                            else if (input == "COUNT"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.BootstrapParam.bootstrapCount = atoi(input.c_str());
                            }
                            else if (input == "ALIGNTOTREECOMMAND"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextQuoteBlock(in);
                                sGOP.BootstrapParam.alignToTreeCommand = input;
                            }
                            input = nextNexusTokenUpper(in);
                        }
                    }
                    if (input == "SVMP"){
                        input = nextNexusTokenUpper(in);
                        while (input != ";"){
                            if (input == "SAMPLESIZE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.SampleParam.SVM_sampleSize = atoi(input.c_str());
                            }
                            else if (input == "RESAMPLESIZE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.SampleParam.SVM_resampleSize = atoi(input.c_str());
                            }
                            else if (input == "BURNINPERCENT"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.SampleParam.burninPercent = atof(input.c_str());
                                sGOP.SampleParam.burninFormat = 0;
                            }
                            else if (input == "BURNINNUMBER"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.SampleParam.burninNumber = atoi(input.c_str());
                                sGOP.SampleParam.burninFormat = 1;
                            }
                            else if (input == "PCACUTOFF"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusNumber(in);
                                sGOP.SampleParam.projectSVDcutOff = atof(input.c_str());
                                sGOP.SampleParam.projectSVD = 1;
                            }
                            else if (input == "PATHDIFFERENCE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusTokenUpper(in);
                                if (input == "YES") {
                                    sGOP.SampleParam.distanceOne = 1;
                                    sGOP.SampleParam.scaleToOne = 0;
                                }
                            }
                            else if (input == "SCALETOONE"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusTokenUpper(in);
                                if (input == "YES") {
                                    sGOP.SampleParam.distanceOne = 0;
                                    sGOP.SampleParam.scaleToOne = 1;
                                }
                            }
                            else if (input == "PRINTMEANTOP"){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusTokenUpper(in);
                                if (input == "YES") {
                                    sGOP.SampleParam.doDiffMeans = 1;
                                }
                            }
                            else if (input == "SAMPLEUNIFORM    "){
                                input = nextNexusToken(in);
                                errorCheckNexusToken(input,"=");
                                input = nextNexusTokenUpper(in);
                                if (input == "YES") {
                                    sGOP.SampleParam.SVM_sampleUniform = 1;
                                }
                                else if (input == "NO") {
                                    sGOP.SampleParam.SVM_sampleUniform = 0;
                                }

                            }
                            input = nextNexusTokenUpper(in);
                        }
                    }
                    if (input == "RANDSEED"){
                        input = nextNexusNumber(in);
                        sGOP.randGenSeed = atoi(input.c_str());
                        readUntilSemiColonNexus(in);
                    }
                    if (input == "CONCATGROUPS"){
                        input = nextNexusTokenUpper(in);
                        if (input == "YES") {
                            sGOP.concatGroups = 1;
                        }
                        readUntilSemiColonNexus(in);
                    }
                    if (input == "JKSTEPONE" || input == "BSSTEPONE"){
                        input = nextNexusTokenUpper(in);
                        if (input == "YES") {
                            sGOP.BSStepOne = 1;
                        }
                        else if (input == "NO") {
                            sGOP.BSStepOne = 0;
                        }
                        readUntilSemiColonNexus(in);
                    }
                    if (input == "TEMPPREFIX"){
                        input = nextQuoteBlock(in);
                        sGOP.tempPrefix = input;
                        tempPrefix = sGOP.tempPrefix;
                        readUntilSemiColonNexus(in);
                    }
                    if (input == "GROUPONEFILES"){
                        input = nextQuoteBlock(in);
                        while (input != ";"){
                            sGOP.inputFileNames.push_back(input); // This assume the block GROUPONEFILES is seen first
                            sGOP.inputFileNamesGroupOne.push_back(input);
                            input = nextQuoteBlock(in);
                        }
                        sGOP.numGroupOne = sGOP.inputFileNamesGroupOne.size();
                    }
                    if (input == "GROUPTWOFILES"){
                        input = nextQuoteBlock(in);
                        while (input != ";"){
                            sGOP.inputFileNames.push_back(input);
                            sGOP.inputFileNamesGroupTwo.push_back(input);
                            input = nextQuoteBlock(in);
                        }
                    }
                    // Process commands in MULTSVMSEP block
                    if (sGOP.testType == 2){
                        if (input == "METHOD"){
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Detecting \"METHOD\"." << endl;
                            }
                            input = nextNexusTokenUpper(in);
                            if (input == "JK" || input == "BS") {
                                sGOP.multSVMSepMethod = 1;
                                sGOP.noTreeCalc = 0;
                                sGOP.statTest = 1;
                                sGOP.geneStatTest = 0;
                            }
                            else if (input == "TREE") {
                                if (DEBUG_OUTPUT >= 1) {
                                    cout << "Setting \"TREE\" method" << endl;
                                }
                                sGOP.multSVMSepMethod = 2;
                                sGOP.noTreeCalc = 1;
                                sGOP.statTest = 1;
                                sGOP.geneStatTest = 0;
                            }
                            else if (input == "NOJK" || input == "NOBS") {
                                sGOP.multSVMSepMethod = 3;
                                sGOP.noTreeCalc = 0;
                                sGOP.statTest = 1;
                                sGOP.geneStatTest = 1;
                            }
                        }
                        if (input == "READALIGNMENTS"){
                            input = nextQuoteBlock(in);
                            sGOP.simCommand = input;
                        }
                        if (input == "TESTPARAMETERS"){
                            input = nextNexusTokenUpper(in);
                            while (input != ";"){
                                if (input == "NUMINITCALC"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.numInitCalc = atoi(input.c_str());
                                    sGOP.skipInitCalc = 0;;
                                }
                                else if (input == "NUMSTATTESTS"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.numStatTests = atoi(input.c_str());
                                }
                                else if (input == "SETINITCALC"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.mySeparation = atof(input.c_str());
                                    sGOP.skipInitCalc = 0;;
                                }
                                else if (input == "BOOTSTRAPIND"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.indBootstrap = atoi(input.c_str());
                                }
                                else if (input == "PERMUTE"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.permuteOrig = atoi(input.c_str());
                                }
                                else if (input == "ALLOWANYPERMUATION"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusNumber(in);
                                    sGOP.allowAnyPermuation = atoi(input.c_str());
                                }
                                input = nextNexusTokenUpper(in);
                            }

                        }
                    }
                    input = nextNexusTokenUpper(in);
                }
                readUntilSemiColonNexus(in);
            }
            if (input == "MULTSVMSEP"){
                readUntilSemiColonNexus(in);
            }
        }
    }

    // Should do some error checking here.

    return in;


}

void GeneOutParameters::init ()
{

    MrBayesParam.MBP_nst = 6; // GTR rate variance model.
    MrBayesParam.MBP_rates = "invgamma"; // gamma shaped rate variation.
    MrBayesParam.MBP_nruns = 2; // Number of independent runs. 
    MrBayesParam.MBP_ngen = 1000; // perform 10,000 runs.
    MrBayesParam.MBP_sampleFreq = 1; // Sample every 10 chains.
    MrBayesParam.MBP_pathToMB = "mb";
    MrBayesParam.MBP_saveOutput = 0; // 1 save output, 0 not 

    // Parameters for MPI implementation of mr bayes.
    MrBayesParam.MPI_np = 1; // Number of processors to use. 1 means no mpi

    BootstrapParam.bootstrapColSize = 1; // Default
    BootstrapParam.bootstrapCount = 1000; // Default
    BootstrapParam.alignToTreeCommand = "./run_dnadist_neighbor_mult"; // Default

    MultIndParam.numTreesReconstruct = 1000; 
    MultIndParam.alignToTreeCommand = "./run_dnadist_neighbor_mult"; // Default

    SampleParam.doSVM = 1;
    SampleParam.doDiffMeans = 0;
    SampleParam.SVM_sampleUniform = 1; // Default. Should sample uniformly from each alignment.
    SampleParam.SVM_sampleSize = -1; // Default. Should calculate based on dimension.
    SampleParam.SVM_resampleSize = -1; // Default. Should calculate based on dimension.
    SampleParam.projectSVD = 0; // Do no project with svd by default
    SampleParam.projectSVDcutOff = 1; // Default value for svd singular value projection
    SampleParam.distanceOne = 0; // 0 means use distance specified in newick, 1 means use 1.
    SampleParam.scaleToOne = 0; // 0 means no scaling, 1 means scale treeLength to 1
    SampleParam.modelType = LINEAR; // DEFAULT
    SampleParam.numTreesPerFile = -1; // Default is unspecified!
    // Parameters for the burnin type and amount.
    SampleParam.burninFormat = -1;    // 0 will mean percentage, 1 will mean fixed number.
    SampleParam.burninPercent = 0.25; // Assume 25% burnin.
    SampleParam.burninNumber = 0;    // Default is 0 since we assume default of 25% burnin


    randGenSeed  = 0;
    doMB = 0; // By default, do not do mr bayes
    doBootstrap = 1;      // Default is to do bootstrap
    doMultInd = 0;
    numInd = 1; // Default is 1 individual per species
    noTreeCalc = 0; // Do not skip mr bayes calculation.
    concatGroups = 0;   // Default is 0, do not concat groups
    numGroupOne = 1; // First file is group one and thats all.
    statTest = 0; // 0 mean no stat test, 1 means stat test
    geneStatTest = 0; // Do the statistical test based on the hyp test on gene trees

    numInitCalc = 1; // 
    numStatTests = 10;
    skipInitCalc = 0; // Dont skip initial calculation for stat test
    mySeparation = -1;
    doSim = 0; 
    simCommand = "";
    multSVMSepMethod = -1;
    testType = -1;
    tempPrefix = "";
    BSStepOne = 0;    // Default is to use the original alignments
    indBootstrap = 1; // 0 use concatenated alignments. 1 means use randomly selected alignments to bootstrap data.
    permuteOrig = 0;
    allowAnyPermuation = 0; //1 means allow anything. 0 means dissallow new group one to be all from group one.
}

void GeneOutParameters::processArg(int argc, char **argv)
{
    if (DEBUG_OUTPUT >= 1) {
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
                if (doMB == 1){
                    doBootstrap = 0;
                    doMultInd = 0;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doMB = " << doMB << endl;
                }
            }
            else if (tempString == "MBP_nst")
            {
                MrBayesParam.MBP_nst = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_nst = " << MrBayesParam.MBP_nst << endl;
                }
            }
            else if (tempString == "MBP_rates")
            {
                MrBayesParam.MBP_rates = argv[++i];
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_rates = " << MrBayesParam.MBP_rates << endl;
                }
            }
            else if (tempString == "MBP_nruns")
            {
                MrBayesParam.MBP_nruns = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_nruns = " << MrBayesParam.MBP_nruns << endl;
                }
            }
            else if (tempString == "MBP_ngen")
            {
                MrBayesParam.MBP_ngen = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_ngen = " << MrBayesParam.MBP_ngen << endl;
                }
            }
            else if (tempString == "MBP_sampleFreq")
            {
                MrBayesParam.MBP_sampleFreq = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_sampleFreq = " << MrBayesParam.MBP_sampleFreq << endl;
                }
            }
            else if (tempString == "MBP_pathToMB")
            {
                MrBayesParam.MBP_pathToMB = argv[++i];
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_pathToMB = " << MrBayesParam.MBP_pathToMB << endl;
                }
            }
            else if (tempString == "MBP_saveOutput")
            {
                MrBayesParam.MBP_saveOutput = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MBP_saveOutput = " << MrBayesParam.MBP_saveOutput << endl;
                }
            }
            else if (tempString == "MPI_np")
            {
                MrBayesParam.MPI_np = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: MPI_np = " << MrBayesParam.MPI_np << endl;
                }
            }
            else if (tempString == "burninPercent")
            {
                SampleParam.burninFormat = 0;
                sscanf(argv[++i],"%lf",&SampleParam.burninPercent);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: burninFormat = " << SampleParam.burninFormat << endl;
                    cout << "       Setting: burninPercent = " << SampleParam.burninPercent << endl;
                }
            }
            else if (tempString == "burninNumber")
            {
                SampleParam.burninFormat = 1;
                SampleParam.burninNumber = atoi(argv[++i]);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: burninFormat = " << SampleParam.burninFormat << endl;
                    cout << "       Setting: burninNumber = " << SampleParam.burninNumber << endl;
                }
            }
            else if (tempString == "randGenSeed")
            {
                randGenSeed = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: randGenSeed = " << randGenSeed << endl;
                }
            }
            else if (tempString == "doBootstrap")
            {
                doBootstrap = atoi(argv[++i]); 
                if (doBootstrap == 1){
                    doMB = 0;
                    doMultInd = 0;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doBootstrap = " << doBootstrap << endl;
                }
                SampleParam.burninFormat = 1;
                SampleParam.burninNumber = 0;
            }
            else if (tempString == "bootstrapColSize")
            {
                BootstrapParam.bootstrapColSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapColSize = " << BootstrapParam.bootstrapColSize << endl;
                }
            }
            else if (tempString == "bootstrapCount")
            {
                BootstrapParam.bootstrapCount = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: bootstrapCount = " << BootstrapParam.bootstrapCount << endl;
                }
            }
            else if (tempString == "doMultInd")
            {
                doMultInd = atoi(argv[++i]); 
                if (doMultInd == 1){
                    doBootstrap = 0;
                    doMB = 0;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doMultInd = " << doMultInd << endl;
                }
            }
            else if (tempString == "numInd")
            {
                numInd = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numInd = " << numInd << endl;
                }
            }
            else if (tempString == "numTreesReconstruct")
            {
                MultIndParam.numTreesReconstruct = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numTreesReconstruct = " << MultIndParam.numTreesReconstruct << endl;
                }
            }
            else if (tempString == "alignToTreeCommand")
            {
                BootstrapParam.alignToTreeCommand = argv[++i]; 
                MultIndParam.alignToTreeCommand = BootstrapParam.alignToTreeCommand;
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: alignToTreeCommand = " << BootstrapParam.alignToTreeCommand << endl;
                }
            }
            else if (tempString == "tempPrefix")
            {
                tempPrefix = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    if ( (tempPrefix.substr(tempPrefix.size()-1,1))[0] != '/'){
                        if (DEBUG_OUTPUT >= 0){
                            cout << "Appending \"/\" to tempPrefix" << endl;
                        }
                        tempPrefix.append("/");
                    }
                    cout << "       Setting: tempPrefix = " << tempPrefix << endl;
                }
            }
            else if (tempString == "SVM_sampleSize")
            {
                SampleParam.SVM_sampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: SVM_sampleSize = " << SampleParam.SVM_sampleSize << endl;
                }
            }
            else if (tempString == "SVM_resampleSize")
            {
                SampleParam.SVM_resampleSize = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: SVM_resampleSize = " << SampleParam.SVM_resampleSize << endl;
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
                SampleParam.numTreesPerFile = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: numTreesPerFile = " << SampleParam.numTreesPerFile << endl;
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
                SampleParam.projectSVD = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: projectSVD = " << SampleParam.projectSVD << endl;
                }
            }
            else if (tempString == "projectSVDcutOff")
            {
                SampleParam.projectSVDcutOff = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: projectSVDcutOff = " << SampleParam.projectSVDcutOff << endl;
                }
            }
            else if (tempString == "distanceOne")
            {
                SampleParam.distanceOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: distanceOne = " << SampleParam.distanceOne << endl;
                }
            }
            else if (tempString == "scaleToOne")
            {
                SampleParam.scaleToOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: scaleToOne = " << SampleParam.scaleToOne << endl;
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
            else if (tempString == "Separation")
            {
                mySeparation = atof(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: Separation = " << mySeparation << endl;
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
                SampleParam.doSVM = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doSVM = " << SampleParam.doSVM << endl;
                }
            }
            else if (tempString == "doDiffMeans")
            {
                SampleParam.doDiffMeans = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: doDiffMeans = " << SampleParam.doDiffMeans << endl;
                }
            }
            else if (tempString == "BSStepOne")
            {
                BSStepOne = atoi(argv[++i]); 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: BSStepOne = " << BSStepOne << endl;
                }
            }
            else if (tempString == "simCommand")
            {
                simCommand = argv[++i]; 
                if (DEBUG_OUTPUT >= 0){
                    cout << "       Setting: simCommand = " << simCommand << endl;
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
                exit(1);
            }
        }
    }
    // Now use numGroupOne to place files in tempFileNames
    int tmpFileCount = 0;
    for (list <string>::const_iterator lsit=inputFileNames.begin();lsit!=inputFileNames.end();lsit++){
        if (tmpFileCount < numGroupOne){
            inputFileNamesGroupOne.push_back(*lsit);
        }
        else{
            inputFileNamesGroupTwo.push_back(*lsit);
        }
        tmpFileCount++;
    }
}
