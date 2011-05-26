// $Rev: 776 $ $Date: 2011-02-10 12:00:19 -0500 (Thu, 10 Feb 2011) $

/// \file geneout.cpp 

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
#include "mww.h"
#include "geneoutparams.h"


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
    string usageMessage = "geneout -nex <config file>\n";

    GeneOutParameters   myGeneOutParam;

    list <list <string> > TreeListListFileNames; // For getTreeMrBayes
    list <string> treeFileNamesGroupOne; // This should be the exact and literal filenames 
    list <string> treeFileNamesGroupTwo; // This should be the exact and literal filenames 

    // Look for -nex command and if it is present, read in all options using nexus file.
    int nexArgFound = 0;
    string nexBlockFileName;
    for (int i=1;i<argc;i++)
    {
        if (argv[i][0] == '-'){
            string tempString = argv[i];
            string tempString2;
            tempString = tempString.substr(1,tempString.size());
            if (tempString == "nex")
            {
                nexArgFound = 1;
                nexBlockFileName = argv[++i];
            }
            // We should also parse DEBUG_OUTPUT related argument here
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
        }
    }

    if (nexArgFound == 1) {
        fstream nexBlockFile;
        nexBlockFile.open(nexBlockFileName.c_str(),fstream::in);
        if (nexBlockFile.bad()){
            cout << "Could not open " << nexBlockFileName << endl;
            exit (0);
        }
        nexBlockFile >> myGeneOutParam;
        nexBlockFile.close();
    }

    if (nexArgFound == 0){ // Then read in parameters from the command line
    cout << " ***** Reading in command-line arguments *****" << endl;
    // Lets let GeneOutParameters handle it.
    myGeneOutParam.processArg(argc,argv);
    }
    if (myGeneOutParam.inputFileNames.size() == 0){
        cout << "No input files!" << endl;
        cout << usageMessage;
        exit(0);
    }
    cout << endl;
    // Output all parameters.
    cout << " ***** Parameters ***** " << endl;
    cout << myGeneOutParam;


    // Do some error checking

    // Maybe this case we dont exit, just increase the sample and resampleSize.
    if (myGeneOutParam.SampleParam.SVM_sampleUniform == 1 && (myGeneOutParam.SampleParam.SVM_sampleSize != -1 || myGeneOutParam.SampleParam.SVM_resampleSize != 1))
    {
        if ( ((myGeneOutParam.SampleParam.SVM_sampleSize % myGeneOutParam.inputFileNamesGroupOne.size()) != 0 ) || ((myGeneOutParam.SampleParam.SVM_sampleSize % myGeneOutParam.inputFileNamesGroupTwo.size()) != 0 ) )
        {
            if ( ((myGeneOutParam.SampleParam.SVM_resampleSize % myGeneOutParam.inputFileNamesGroupOne.size()) != 0 ) || ((myGeneOutParam.SampleParam.SVM_resampleSize % myGeneOutParam.inputFileNamesGroupTwo.size()) != 0 ) )
            {
                // cout << "sampleSize or resampleSize does not divide inputFileNamesGroupOne or Two" << endl;
                // exit(0);
                // Lets instead up the sample and resample size.
                int fileSizeLCM = lcm(myGeneOutParam.inputFileNamesGroupOne.size(), myGeneOutParam.inputFileNamesGroupTwo.size());
                int newCoeffOne = ceil(double(myGeneOutParam.SampleParam.SVM_sampleSize) / double(fileSizeLCM));
                int newCoeffTwo = ceil(double(myGeneOutParam.SampleParam.SVM_resampleSize) / double(fileSizeLCM));

                if (DEBUG_OUTPUT >= 0){
                    cout << "sampleSize or resampleSize does not divide inputFileNamesGroupOne or Two" << endl;
                    cout << "Correcting..." << endl;
                    cout << "fileSizeLCM = " << fileSizeLCM << endl;
                    cout << "newCoeffOne = " << newCoeffOne << endl;
                    cout << "newCoeffTwo = " << newCoeffTwo << endl;
                    cout << "new myGeneOutParam.SampleParam.SVM_sampleSize = " << newCoeffOne*fileSizeLCM << endl;
                    cout << "myGeneOutParam.SampleParam.SVM_resampleSize = " << newCoeffTwo*fileSizeLCM << endl;
                }

                myGeneOutParam.SampleParam.SVM_sampleSize = newCoeffOne*fileSizeLCM;
                myGeneOutParam.SampleParam.SVM_resampleSize = newCoeffTwo*fileSizeLCM;
            }
        }
    }


    if (myGeneOutParam.doMB + myGeneOutParam.doJackknife + myGeneOutParam.doMultInd % 3 != 1){
        cout << "Must set either myGeneOutParam.doMB or myGeneOutParam.doJackknife or myGeneOutParam.doMultInd, not all or any two." << endl;
        exit (0);
    }
    if ((int)myGeneOutParam.inputFileNames.size () < myGeneOutParam.numGroupOne){
        cout << "Not enough files for myGeneOutParam.numGroupOne" << endl;
        exit (0);
    }
    if (myGeneOutParam.skipInitCalc == 1 && myGeneOutParam.mySeparation == -1){
        cout << "myGeneOutParam.skipInitCalc selected but no myGeneOutParam.mySeparation set." << endl;
        exit(0);
    }

    if (myGeneOutParam.noTreeCalc == 1 && myGeneOutParam.SampleParam.numTreesPerFile == -1){
        if (DEBUG_OUTPUT >= 0) {
            cout << "myGeneOutParam.SampleParam.numTreesPerFile == -1 AND myGeneOutParam.noTreeCalc == 1" << endl;
            cout << "Computing number of trees per file." << endl;
        }
        int firstTreeRun = 1;
        for (list<string>::const_iterator lit=myGeneOutParam.inputFileNames.begin();lit!=myGeneOutParam.inputFileNames.end();lit++){
            if (firstTreeRun == 1){
                myGeneOutParam.SampleParam.numTreesPerFile = numTreesNexus(*lit);
                firstTreeRun = 0;
            }
            else {
                if (numTreesNexus(*lit) != myGeneOutParam.SampleParam.numTreesPerFile) {
                    cout << "Not all tree files have same number of trees." << endl;
                    exit (0);
                }
            }
        }
    }


    // Error check for SVM vs number of trees that are going to be available.

    if (myGeneOutParam.randGenSeed == 0){
        srand(time(0));
    }
    else {
        srand(time(&myGeneOutParam.randGenSeed));
    }

    if (DEBUG_OUTPUT >= 0){
        cout << "Input files:" << endl;
        for (list<string>::const_iterator lit=myGeneOutParam.inputFileNames.begin();lit!=myGeneOutParam.inputFileNames.end();lit++){
             cout << *lit << endl; 
        }
        cout << endl;
    }

    // If we are working on alignments, we read them all in here.
    int numTaxa = 0;
    list <int> inputAlignmentsLengths;
    list <int> inputAlignmentsLengthsOne;
    list <int> inputAlignmentsLengthsTwo;
    list <Alignments *> inputAlignments;
    list <Alignments *> origInputAlignments;
    list <Alignments *> AlignmentsOne, AlignmentsTwo;
    list <Alignments *> origAlignmentsOne, origAlignmentsTwo;
    Alignments         *AlignmentsOneConcat, *AlignmentsTwoConcat;
    if (myGeneOutParam.noTreeCalc != 1){
        if (myGeneOutParam.concatGroups != 1) {
            for (list<string>::const_iterator lsit=myGeneOutParam.inputFileNames.begin();lsit!=myGeneOutParam.inputFileNames.end();lsit++){
                Alignments *newAlignments = new Alignments(*lsit);
                numTaxa = newAlignments->get_ntax ();
                inputAlignmentsLengths.push_back(newAlignments->get_nchar()); 
                inputAlignments.push_back(newAlignments);
                //cout << *newAlignments;
            }
            // Again, default is to compare the first input vs the rest.
            // Lets us myGeneOutParam.numGroupOne
            list <Alignments *>::iterator lait=inputAlignments.begin();
            int groupAddCount = 0;
            AlignmentsOneConcat = new Alignments;
            int firstConcatRunOne = 1;
            while (groupAddCount < myGeneOutParam.numGroupOne){
                if (firstConcatRunOne == 1){
                    firstConcatRunOne=0;
                    *AlignmentsOneConcat = *(*lait);
                }
                else {
                    *AlignmentsOneConcat += *(*lait);
                }
                inputAlignmentsLengthsOne.push_back((*lait)->get_nchar());
                AlignmentsOne.push_back(*lait);
                // Make copy of data and store in origAlignmentsOne.
                // Don't want to just put in *lait since they would both point
                // to the same data which is cleared later.
                Alignments *newAlignments = new Alignments;
                (*newAlignments) = *(*lait);
                origAlignmentsOne.push_back(newAlignments);
                origInputAlignments.push_back(newAlignments);

                lait++;
                groupAddCount++;
            }

            // Put rest in AlignmentsTwo
            AlignmentsTwoConcat = new Alignments;
            int firstConcatRunTwo = 1;
            while(lait!=inputAlignments.end()){
                if (firstConcatRunTwo == 1){
                    firstConcatRunTwo=0;
                    *AlignmentsTwoConcat = *(*lait);
                }
                else {
                    *AlignmentsTwoConcat += *(*lait);
                }
                inputAlignmentsLengthsTwo.push_back((*lait)->get_nchar());
                AlignmentsTwo.push_back(*lait);
                // Make copy of data and store in origAlignmentsTwo.
                // Don't want to just put in *lait since they would both point
                // to the same data which is cleared later.
                Alignments *newAlignments = new Alignments;
                *newAlignments = *(*lait);
                origAlignmentsTwo.push_back(newAlignments);
                origInputAlignments.push_back(newAlignments);
                lait++;
            }

            
            ////Concat all the alignments together, except the first one.
            //list <Alignments *>::iterator lait=inputAlignments.begin();
            //int newColSizeCount = (int)floor((*lait)->get_nchar()/myGeneOutParam.JackknifeParam.jackknifeColSize);
            ////cout << "newColSizeCount " << newColSizeCount << endl;
            ////cout << *(*lait);
            //// Why just the first! Use myGeneOutParam.numGroupOne
            //lait++;
            ////cout << *(*lait);
            //Alignments *allConcatMinusFirst = new Alignments;
            //*allConcatMinusFirst = *(*lait); // Copy second alignment
            //lait++;
            //while(lait!=inputAlignments.end()){
            //    //cout << *(*lait);
            //    *allConcatMinusFirst += *(*lait);
            //    lait++;
            //}
            ////cout << "All alignments minus first: " << endl;
            ////cout << *allConcatMinusFirst;

        }
        else if (myGeneOutParam.concatGroups == 1) {
            int fileCount = 0;
            Alignments *groupOneAlignmentsConcat;
            Alignments *groupTwoAlignmentsConcat;
            int firstGroupOne = 1;
            int firstGroupTwo = 1;
            for (list<string>::const_iterator lsit=myGeneOutParam.inputFileNames.begin();lsit!=myGeneOutParam.inputFileNames.end();lsit++){
                if (fileCount < myGeneOutParam.numGroupOne) {
                    if (firstGroupOne == 1){
                        firstGroupOne = 0;
                        groupOneAlignmentsConcat = new Alignments(*lsit);
                    }
                    Alignments *newAlignments = new Alignments(*lsit);
                    *groupOneAlignmentsConcat += *newAlignments;
                    delete newAlignments;
                }
                else {
                    if (firstGroupTwo == 1){
                        firstGroupTwo = 0;
                        groupTwoAlignmentsConcat = new Alignments(*lsit);
                    }
                    Alignments *newAlignments = new Alignments(*lsit);
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
            Alignments *newAlignments = new Alignments;
            *newAlignments = (*groupOneAlignmentsConcat);
            origAlignmentsOne.push_back(newAlignments);
            AlignmentsTwo.push_back(groupTwoAlignmentsConcat);
            newAlignments = new Alignments;
            *newAlignments = (*groupTwoAlignmentsConcat);
            origAlignmentsTwo.push_back(newAlignments);
        }
        if (DEBUG_OUTPUT >= 1){
            for (list <Alignments *>::iterator lait=inputAlignments.begin();lait!=inputAlignments.end();lait++){
                cout << *(*lait);
            }
        }
    }


    // Set species groupings according to myGeneOutParam.numInd
    set <unsigned> currGrouping;
    for (int i=0;i<numTaxa;i++){
        currGrouping.insert(i);
        if((i+1) % (myGeneOutParam.numInd) == 0) {
            myGeneOutParam.MultIndParam.speciesGroupings.push_back(currGrouping);
            currGrouping.clear();
        }
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "Number of taxa: " << numTaxa << endl;
        for (list <set <unsigned> >::iterator lsit=myGeneOutParam.MultIndParam.speciesGroupings.begin();lsit!=myGeneOutParam.MultIndParam.speciesGroupings.end();lsit++){
            cout << "   +" << endl;
            for (set <unsigned>::iterator sit=(*lsit).begin();sit!=(*lsit).end();sit++){
                cout << setw(3) << *sit << " ";
            }
            cout << endl;
        }
    }

    int groupOneCount=0;
    int groupTwoCount=0;

    // This is the method without permutation.
    if (myGeneOutParam.statTest != 1){
        if (myGeneOutParam.noTreeCalc != 1 && myGeneOutParam.doMB == 1){
            if (DEBUG_OUTPUT >= 0){
                cout << endl << " ***** Performing Mr Bayes calculations ***** " << endl << endl;
            }
            myGeneOutParam.MrBayesParam.inputNexFiles = myGeneOutParam.inputFileNames;
            // Unecessary. numTreesPerFile set in calcSVM....
            myGeneOutParam.SampleParam.numTreesPerFile = (myGeneOutParam.MrBayesParam.MBP_ngen/myGeneOutParam.MrBayesParam.MBP_sampleFreq);
            SVM_separationResults mySeparationResults;
            MrBayesResults MB_results;
            calcSVMseparationMrBayes(myGeneOutParam.MrBayesParam,MB_results, myGeneOutParam.SampleParam,mySeparationResults,myGeneOutParam.numGroupOne);

        }
        else if (myGeneOutParam.noTreeCalc != 1 && myGeneOutParam.doJackknife == 1){
            if (DEBUG_OUTPUT >= 0){
                cout << endl << " ***** Performing Jackknife ***** " << endl << endl;
            }
            SVM_separationResults mySeparationResults;
            calcSVMseparationJackknife(myGeneOutParam.inputFileNames,myGeneOutParam.numGroupOne,myGeneOutParam.JackknifeParam, myGeneOutParam.SampleParam,mySeparationResults);
        }
        else if (myGeneOutParam.noTreeCalc == 1) {
            if (DEBUG_OUTPUT >= 0){
                cout << endl << " ***** Skipping Tree Calculations ***** " << endl << endl;
                cout << "Forming files for group one." << endl;
            }
            list <string>::const_iterator lsit=myGeneOutParam.inputFileNames.begin();
            for (int i=0;i<myGeneOutParam.numGroupOne;i++){
                if (DEBUG_OUTPUT >= 0){
                    cout << "   " << *lsit << endl;
                }
                treeFileNamesGroupOne.push_back(*lsit);
                lsit++;
            }
            if (DEBUG_OUTPUT >= 0){
                cout << "Forming files for group two." << endl;
            }
            while(lsit != myGeneOutParam.inputFileNames.end()){
                if (DEBUG_OUTPUT >= 0){
                    cout << "   " << *lsit << endl;
                }
                treeFileNamesGroupTwo.push_back(*lsit);
                lsit++;
            }

            //myGeneOutParam.SampleParam.numTreesPerFile = ;
            SVM_separationResults mySeparationResults;
            calcSVMseparation(treeFileNamesGroupOne,treeFileNamesGroupTwo,myGeneOutParam.SampleParam,mySeparationResults);
            
        }
    } // END NO STATISTICAL TEST
    //  _____ _        _     _____         _   
    // /  ___| |      | |   |_   _|       | |  
    // \ `--.| |_ __ _| |_    | | ___  ___| |_ 
    //  `--. \ __/ _` | __|   | |/ _ \/ __| __|
    // /\__/ / || (_| | |_    | |  __/\__ \ |_ 
    // \____/ \__\__,_|\__|   \_/\___||___/\__|
    else if (myGeneOutParam.statTest == 1) // STAT TEST
    {
        if (myGeneOutParam.noTreeCalc != 1 && myGeneOutParam.geneStatTest != 1){
            if (DEBUG_OUTPUT >= 0){
                cout << "**** Performing alignments statistical test ****" << endl;
                if (myGeneOutParam.doJackknife == 1){
                    cout << "    Using Jackknifing." << endl;
                }
                if (myGeneOutParam.doMB == 1){
                    cout << "    Using Mr Bayes." << endl;
                }
            }

            int groupOneCount;
            int groupTwoCount;
            list <double> sepValuesOne; // Holds the separation value for step one.
            list <double> sepValuesTwo; // Holds the separation value for step two.
            double initSeparation = 0;
            if (myGeneOutParam.skipInitCalc != 1){
                int firstRun = 1;
                SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
                if (DEBUG_OUTPUT >= 0){
                    cout << "   Calculating primary separation value between group one and group two." << endl;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "Performing " << myGeneOutParam.numInitCalc << " initial separation calculations and averaging." << endl;
                }
                // Get a cutoff/baseline, pc.
                // Why am I repeatedly calling Mr Bayes with the same input. If we bootstrap the alignments in step 2,
                // we should do that here also!
                for (int i=0;i<myGeneOutParam.numInitCalc;i++){
                    if (myGeneOutParam.JKStepOne == 1) {
                        cout << "   Jackknifing new alignments from original." << endl;
                        // We should jackknife to create new alignments. Not use the original. Basically, we should copy step 2!
                        // Cant just clear. There is dynamically allocated data here.
                        AlignmentsOne.clear();
                        AlignmentsTwo.clear();
                        getSomeJackknifeAlignments(*AlignmentsOneConcat,AlignmentsOne,inputAlignmentsLengthsOne,myGeneOutParam.JackknifeParam.jackknifeColSize);
                        getSomeJackknifeAlignments(*AlignmentsTwoConcat,AlignmentsTwo,inputAlignmentsLengthsTwo,myGeneOutParam.JackknifeParam.jackknifeColSize);
                    }
                    else {
                        cout << "   Using original alignments." << endl;
                    }
                    // Generate two new alignments from allConcatMinusFirst and try svm.
                    // Dangerous, because we cant see any output, but necessary.
                    int initPrimSepStartTime = time(0);
                    int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                    SVM_separationResults mySeparationResults;
                    MrBayesResults MB_results;
                    MB_results.clear();
                    if (DEBUG_OUTPUT <= 0) {
                        DEBUG_OUTPUT = -1;
                    }
                    if (myGeneOutParam.doJackknife == 1){
                        calcSVMseparationJackknife(AlignmentsOne,AlignmentsTwo, myGeneOutParam.JackknifeParam, myGeneOutParam.SampleParam,mySeparationResults);
                    }
                    if (myGeneOutParam.doMB == 1){
                        // Unecessary. numTreesPerFile set in calcSVM...
                        myGeneOutParam.SampleParam.numTreesPerFile = (myGeneOutParam.MrBayesParam.MBP_ngen/myGeneOutParam.MrBayesParam.MBP_sampleFreq);
                        calcSVMseparationMrBayes(AlignmentsOne,AlignmentsTwo, myGeneOutParam.MrBayesParam,MB_results, myGeneOutParam.SampleParam,mySeparationResults);
                    }
                    if (myGeneOutParam.doMultInd == 1){
                        calcSVMseparationMultInd(AlignmentsOne,AlignmentsTwo, myGeneOutParam.MultIndParam, myGeneOutParam.SampleParam,mySeparationResults);
                    }
                    DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                    if (firstRun == 1){
                        averagedResults = mySeparationResults;
                        firstRun = 0;
                    }
                    else {
                        averagedResults += mySeparationResults;
                    }
                    int initPrimSepEndTime = time(0);
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   " << initPrimSepEndTime - initPrimSepStartTime << " seconds." << endl;
                        //cout << "   " << mySeparationResults.secondsToCompute << " seconds." << endl;
                    }

                    //double thisInitSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                    double thisInitSeparation = mySeparationResults.resampleSeparationPercentage;
                    if (DEBUG_OUTPUT >= 0){
                        cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                        cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                        //cout << "       Separation percentage: " << setw(6) << 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize) << endl;
                        cout << "       This Separation percentage(" << i + 1 << "): " << setw(6) << thisInitSeparation << endl;
                        if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                            mySeparationResults.print("This Initial");
                            if (myGeneOutParam.SampleParam.projectSVD == 1){
                                mySeparationResults.printPCA("This Initial");
                            }
                        }
                        if (myGeneOutParam.doMB == 1){
                            cout << "       Split frequencies: ";
                            double  tempDouble = 0;
                            for(list <double>::const_iterator dlit=MB_results.splitFreqs.begin();dlit!=MB_results.splitFreqs.end();dlit++){
                                cout << *dlit <<  " ";
                                tempDouble += *dlit;
                            }
                            cout << endl;
                            cout << "       Average split frequency: " << tempDouble/MB_results.splitFreqs.size() << endl; 
                        }
                    }
                    initSeparation+= thisInitSeparation;
                    sepValuesOne.push_back(thisInitSeparation);
                }
                averagedResults /= myGeneOutParam.numInitCalc;
                initSeparation /= myGeneOutParam.numInitCalc;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Averaged Initial Separation Percentage: " << setw(6) << initSeparation << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        averagedResults.print("Averaged Initial");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            averagedResults.printPCA("Averaged Initial");
                        }
                    }
                }
            }
            else if (myGeneOutParam.skipInitCalc){
                initSeparation = myGeneOutParam.mySeparation;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Setting initSeparation = " << myGeneOutParam.mySeparation << endl;
                }
            }
            
            
            if (DEBUG_OUTPUT >= 0){
                cout << endl;
                cout << "Now performing " << myGeneOutParam.numStatTests << " statistical tests." << endl;
            }

            Alignments *tempAlignments;
            // for 1 to myGeneOutParam.numStatTests create a new set of alignments from groupTwo using jacknifing
            //      then for each new set of alignments, run mr bayes or jackknifing
            //      and svm to find the separation 
            //      if (separation > pc) pval++
            // return pval/myGeneOutParam.numStatTests
            double pval = 0;
            double averageSeparation=0;
            SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
            int firstRun = 1;
            for (int i=0;i<myGeneOutParam.numStatTests;i++){
                SVM_separationResults mySeparationResults;
                MrBayesResults MB_results;
                int innerStatTestStartTime = time(0);
                // Clear previous Alignments.
                for (list <Alignments *>::iterator lait=AlignmentsOne.begin();lait!=AlignmentsOne.end();lait++){
                    delete (*lait);
                }
                AlignmentsOne.clear();
                for (list <Alignments *>::iterator lait=AlignmentsTwo.begin();lait!=AlignmentsTwo.end();lait++){
                    delete (*lait);
                }
                AlignmentsTwo.clear();

                // Create an alignment to compare the rest to
                // We need to be able to either jackknife OR if we are simulating
                // data we should have option to generate from parameters instead
                // of jackknifing.
                if (myGeneOutParam.doSim == 1){
                    userGenerateNewAlignmentsUniform(myGeneOutParam.simCommand, AlignmentsOne, AlignmentsTwo,inputAlignmentsLengths.size(),myGeneOutParam.numGroupOne); 
                }
                else {
                    if (myGeneOutParam.indJackknife == 0 && myGeneOutParam.permuteOrig == 0) {
                        if (DEBUG_OUTPUT >= 0){
                            cout << "Bootstrapping alignments from concat of second group of alignments." << endl;
                        }
                        getSomeJackknifeAlignments(*AlignmentsTwoConcat,AlignmentsOne,AlignmentsTwo,inputAlignmentsLengths,myGeneOutParam.JackknifeParam.jackknifeColSize,myGeneOutParam.numGroupOne);
                    }
                    if (myGeneOutParam.indJackknife == 1 && myGeneOutParam.permuteOrig == 0) {
                        if (DEBUG_OUTPUT >= 0){
                            cout << "Bootstrapping alignments from second group of alignments." << endl;
                        }
                        getSomeJackknifeAlignments(origAlignmentsTwo,AlignmentsOne,AlignmentsTwo,inputAlignmentsLengths,myGeneOutParam.JackknifeParam.jackknifeColSize,myGeneOutParam.numGroupOne);
                    }
                    if (myGeneOutParam.permuteOrig == 1) {
                        if (DEBUG_OUTPUT >= 0){
                            cout << "Permuting all alignments and bootstrapping to appropriate size." << endl;
                        }
                        getSomeJackknifeAlignmentsPermute(origInputAlignments,AlignmentsOne,AlignmentsTwo,inputAlignmentsLengths,myGeneOutParam.JackknifeParam.jackknifeColSize,myGeneOutParam.numGroupOne,myGeneOutParam.allowAnyPermuation);
                    }

                }

                int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                if (DEBUG_OUTPUT <= 0) {
                    DEBUG_OUTPUT = -1;
                }
                if (myGeneOutParam.doJackknife == 1){
                    calcSVMseparationJackknife(AlignmentsOne,AlignmentsTwo, myGeneOutParam.JackknifeParam, myGeneOutParam.SampleParam,mySeparationResults);
                }
                if (myGeneOutParam.doMB == 1){
                    // Unnecessary.
                    myGeneOutParam.SampleParam.numTreesPerFile = (myGeneOutParam.MrBayesParam.MBP_ngen/myGeneOutParam.MrBayesParam.MBP_sampleFreq);
                    calcSVMseparationMrBayes(AlignmentsOne,AlignmentsTwo, myGeneOutParam.MrBayesParam,MB_results, myGeneOutParam.SampleParam,mySeparationResults);
                }
                if (myGeneOutParam.doMultInd == 1){
                    calcSVMseparationMultInd(AlignmentsOne,AlignmentsTwo, myGeneOutParam.MultIndParam, myGeneOutParam.SampleParam,mySeparationResults);
                }
                DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                if (firstRun == 1){
                    averagedResults = mySeparationResults;
                    firstRun = 0;
                }
                else {
                    averagedResults += mySeparationResults;
                }
                //double thisSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                double thisSeparation = mySeparationResults.resampleSeparationPercentage;
                
                int innerStatTestEndTime = time(0);
                if (DEBUG_OUTPUT >= 0){
                    cout << endl;
                    cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                    cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                    cout << "       Separation percentage: " << setw(6) << thisSeparation << endl;
                    cout << "   " << innerStatTestEndTime - innerStatTestStartTime << " seconds. Stat test " << i+1 << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        mySeparationResults.print("This Stat Test");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            mySeparationResults.printPCA("This Stat Test");
                        }
                    }
                    if (myGeneOutParam.doMB == 1){
                        cout << "       Split frequencies: ";
                        double  tempDouble = 0;
                        for(list <double>::const_iterator dlit=MB_results.splitFreqs.begin();dlit!=MB_results.splitFreqs.end();dlit++){
                            cout << *dlit <<  " ";
                            tempDouble += *dlit;
                        }
                        cout << endl;
                        cout << "       Average split frequency: " << tempDouble/MB_results.splitFreqs.size() << endl; 
                    }
                }
                averageSeparation += thisSeparation;
                sepValuesTwo.push_back(thisSeparation);
                if (initSeparation <= thisSeparation){
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value higher than or equal to initSeparation." << endl;
                        cout << "       pval++" << endl;
                    }
                    pval++;
                }
                else {
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value lower than initSeparation." << endl;
                    }
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "   ### Current(" << i+1 << ") pval = " << setw(3) << pval/(i+1) << " ###" << endl;
                }

            }
            averageSeparation/= myGeneOutParam.numStatTests;
            averagedResults /= myGeneOutParam.numStatTests;
            if (DEBUG_OUTPUT >= 0){
                cout << "pval/myGeneOutParam.numStatTests = " << pval/myGeneOutParam.numStatTests << endl;
                cout << "p-value: " << pval/myGeneOutParam.numStatTests << endl;
                cout << "Averaged Separation Percentage: " << averageSeparation << endl;
                if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                    averagedResults.print("Averaged Stat Test");
                    if (myGeneOutParam.SampleParam.projectSVD == 1){
                        averagedResults.printPCA("Averaged Stat Test");
                    }
                }
                cout << "Init separation values: ";
                for (list <double>::const_iterator lit=sepValuesOne.begin();lit!=sepValuesOne.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;
                cout << "Step two separation values: ";
                for (list <double>::const_iterator lit=sepValuesTwo.begin();lit!=sepValuesTwo.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;

//                cout << "MWW pvalue: " << MannWhitneyWilcoxon(sepValuesOne,sepValuesTwo) << endl;
//                cout << "Min pvalue: " << pvalueMin(sepValuesOne, sepValuesTwo) << endl;
            }
        }
//        _______                  ____        _       
//       |__   __|                / __ \      | |      
//          | |_ __ ___  ___  ___| |  | |_ __ | |_   _ 
//          | | '__/ _ \/ _ \/ __| |  | | '_ \| | | | |
//          | | | |  __/  __/\__ \ |__| | | | | | |_| |
//          |_|_|  \___|\___||___/\____/|_| |_|_|\__, |
//                                                __/ |
//                                               |___/ 
        if (myGeneOutParam.noTreeCalc == 1 && myGeneOutParam.geneStatTest != 1){
            if (DEBUG_OUTPUT >= 0){
                cout << "**** Performing trees only statistical test ****" << endl;
            }
            // Here we are only working with a set of tree files.
            list <string>::const_iterator lsit=myGeneOutParam.inputFileNames.begin();
            if (DEBUG_OUTPUT >= 0){
                cout << "Forming files for group one." << endl;
            }
            for (int i=0;i<myGeneOutParam.numGroupOne;i++){
                if (DEBUG_OUTPUT >= 0){
                    cout << "   " << *lsit << endl;
                }
                treeFileNamesGroupOne.push_back(*lsit);
                lsit++;
            }
            if (DEBUG_OUTPUT >= 0){
                cout << "Forming files for group two." << endl;
            }
            while(lsit != myGeneOutParam.inputFileNames.end()){
                if (DEBUG_OUTPUT >= 0){
                    cout << "   " << *lsit << endl;
                }
                treeFileNamesGroupTwo.push_back(*lsit);
                lsit++;
            }

            int groupOneCount;
            int groupTwoCount;
            list <double> sepValuesOne; // Holds the separation value for step one.
            list <double> sepValuesTwo; // Holds the separation value for step two.

            double initSeparation = 0;
            if (myGeneOutParam.skipInitCalc != 1){
                int firstRun = 1;
                SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
                if (DEBUG_OUTPUT >= 0){
                    cout << "   Calculating primary separation value between group one and group two." << endl;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "Performing " << myGeneOutParam.numInitCalc << " initial separation calculations and averaging." << endl;
                }
                // Get a cutoff/baseline, pc.
                // I would also like the ability to get an average of some number of runs
                for (int i=0;i<myGeneOutParam.numInitCalc;i++){

                    SVM_separationResults mySeparationResults;
                    int initPrimSepStartTime = time(0);
                    // Dangerous, because we cant see any output, but necessary.
                    int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                    if (DEBUG_OUTPUT <= 0) {
                        DEBUG_OUTPUT = -1;
                    }
                    calcSVMseparation(treeFileNamesGroupOne,treeFileNamesGroupTwo,myGeneOutParam.SampleParam,mySeparationResults);
                    DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                    if (firstRun == 1){
                        averagedResults = mySeparationResults;
                        firstRun = 0;
                    }
                    else {
                        averagedResults += mySeparationResults;
                    }
                    int initPrimSepEndTime = time(0);
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   " << initPrimSepEndTime - initPrimSepStartTime << " seconds." << endl;
                    }

                    //double thisInitSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                    double thisInitSeparation = mySeparationResults.resampleSeparationPercentage;
                    if (DEBUG_OUTPUT >= 0){
                        cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                        cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                        //cout << "       Separation percentage: " << setw(6) << 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize) << endl;
                        cout << "       This Separation percentage(" << i + 1 << "): " << setw(6) << thisInitSeparation << endl;
                        if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                            mySeparationResults.print("This Initial");
                            if (myGeneOutParam.SampleParam.projectSVD == 1){
                                mySeparationResults.printPCA("This Initial");
                            }
                        }
                    }
                    initSeparation+= thisInitSeparation;
                    sepValuesOne.push_back(thisInitSeparation);
                }
                averagedResults /= myGeneOutParam.numInitCalc;
                initSeparation /= myGeneOutParam.numInitCalc;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Averaged Initial Separation percentage: " << setw(6) << initSeparation << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        averagedResults.print("Averaged Initial");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            averagedResults.printPCA("Averaged Initial");
                        }
                    }
                }
            }
            else if (myGeneOutParam.skipInitCalc){
                initSeparation = myGeneOutParam.mySeparation;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Setting initSeparation = " << myGeneOutParam.mySeparation << endl;
                }
            }

            if (DEBUG_OUTPUT >= 0){
                cout << endl;
                cout << "Now performing " << myGeneOutParam.numStatTests << " statistical tests." << endl;
            }

            // if (separation > pc) pval++
            // return pval/myGeneOutParam.numStatTests
            double pval = 0;
            double averageSeparation=0;
            SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
            int firstRun = 1;
            for (int i=0;i<myGeneOutParam.numStatTests;i++){
                SVM_separationResults mySeparationResults;
                int innerStatTestStartTime = time(0);
                // Clear previous Alignments.

                int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                if (DEBUG_OUTPUT <= 0) {
                    DEBUG_OUTPUT = -1;
                }
                // Compute separation in treeFileNamesGroupTwo
                calcSVMseparation(treeFileNamesGroupTwo,treeFileNamesGroupTwo,myGeneOutParam.SampleParam,mySeparationResults);
                DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                if (firstRun == 1){
                    averagedResults = mySeparationResults;
                    firstRun = 0;
                }
                else {
                    averagedResults += mySeparationResults;
                }
                //double thisSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                double thisSeparation = mySeparationResults.resampleSeparationPercentage;
                
                int innerStatTestEndTime = time(0);
                if (DEBUG_OUTPUT >= 0){
                    cout << endl;
                    cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                    cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                    cout << "       Separation percentage: " << setw(6) << thisSeparation << endl;
                    cout << "   " << innerStatTestEndTime - innerStatTestStartTime << " seconds. Stat test " << i+1 << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        mySeparationResults.print("This Stat Test");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            mySeparationResults.printPCA("This Stat Test");
                        }
                    }
                }
                averageSeparation += thisSeparation;
                sepValuesTwo.push_back(thisSeparation);
                if (initSeparation <= thisSeparation){
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value higher than or equal to initSeparation." << endl;
                        cout << "       pval++" << endl;
                    }
                    pval++;
                }
                else {
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value lower than initSeparation." << endl;
                    }
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "   ### Current(" << i+1 << ") pval = " << setw(3) << pval/(i+1) << " ###" << endl;
                }

            }
            averageSeparation/= myGeneOutParam.numStatTests;
            averagedResults /= myGeneOutParam.numStatTests;
            if (DEBUG_OUTPUT >= 0){
                cout << "pval/myGeneOutParam.numStatTests = " << pval/myGeneOutParam.numStatTests << endl;
                cout << "p-value: " << pval/myGeneOutParam.numStatTests << endl;
                cout << "Averaged Separation Percentage: " << averageSeparation << endl;
                if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                    averagedResults.print("Averaged Stat Test");
                    if (myGeneOutParam.SampleParam.projectSVD == 1){
                        averagedResults.printPCA("Averaged Stat Test");
                    }
                }
                cout << "Init separation values: ";
                for (list <double>::const_iterator lit=sepValuesOne.begin();lit!=sepValuesOne.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;
                cout << "Step two separation values: ";
                for (list <double>::const_iterator lit=sepValuesTwo.begin();lit!=sepValuesTwo.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;
//                cout << "MWW pvalue: " << MannWhitneyWilcoxon(sepValuesOne,sepValuesTwo) << endl;
//                cout << "Min pvalue: " << pvalueMin(sepValuesOne, sepValuesTwo) << endl;
            }
        }
//                                _        _   
//                               | |      | |  
//      __ _  ___ _ __   ___  ___| |_ __ _| |_ 
//     / _` |/ _ \ '_ \ / _ \/ __| __/ _` | __|
//    | (_| |  __/ | | |  __/\__ \ || (_| | |_ 
//     \__, |\___|_| |_|\___||___/\__\__,_|\__|
//      __/ |                                  
//     |___/ 
        else if (myGeneOutParam.noTreeCalc != 1 && myGeneOutParam.geneStatTest == 1){ // Gene stat test BEGIN
            MrBayesResults MB_results;
            if (DEBUG_OUTPUT >= 0){
                cout << "**** Performing genes statistical test ****" << endl;
                if (myGeneOutParam.doJackknife == 1){
                    cout << "    Using Jackknifing." << endl;
                }
                if (myGeneOutParam.doMB == 1){
                    cout << "    Using Mr Bayes." << endl;
                }
            }
            // This is the statistical test where we only compute the distribution of
            // trees once.
            list <list <string> > treeFileNames;

            if (myGeneOutParam.doJackknife == 1){
                //myGeneOutParam.JackknifeParam.inputNexFiles = myGeneOutParam.inputFileNames;
                myGeneOutParam.SampleParam.numTreesPerFile = myGeneOutParam.JackknifeParam.jackknifeCount;
                //getTreesJackknife(myGeneOutParam.JackknifeParam, treeFileNames);
                //getTreesJackknife(myGeneOutParam.JackknifeParam, treeFileNames);
                getTreesJackknife(myGeneOutParam.JackknifeParam, inputAlignments, treeFileNames);
            }
            if (myGeneOutParam.doMB == 1){
                //list <Alignments *> myAlignments;
                //for (list <string>::const_iterator lsit=myGeneOutParam.inputFileNames.begin();lsit!=myGeneOutParam.inputFileNames.end();lsit++){
                //    Alignments *newAlignments = new Alignments(*lsit);
                //    myAlignments.push_back(newAlignments);
                //}
                // Let getTreesJackknife open up the inputfiles
                myGeneOutParam.SampleParam.numTreesPerFile = (myGeneOutParam.MrBayesParam.MBP_ngen/myGeneOutParam.MrBayesParam.MBP_sampleFreq);
                getTreesMrBayes(myGeneOutParam.MrBayesParam, inputAlignments, treeFileNames, MB_results);
            }

            list <list <string> >::iterator llit = treeFileNames.begin();
            cout << "treeFileNames:" << endl;
            for (int i=0;i<myGeneOutParam.numGroupOne;i++){
                for (list <string>::iterator lsit = (*llit).begin();lsit!= (*llit).end();lsit++){
                    cout << "*lsit: " << *lsit << endl;
                    treeFileNamesGroupOne.push_back(*lsit);
                }
                llit++;
            }
            while(llit != treeFileNames.end()){
                for (list <string>::iterator lsit = (*llit).begin();lsit!= (*llit).end();lsit++){
                    cout << "*lsit: " << *lsit << endl;
                    treeFileNamesGroupTwo.push_back(*lsit);
                }
                llit++;
            }

            int groupOneCount;
            int groupTwoCount;
            list <double> sepValuesOne; // Holds the separation value for step one.
            list <double> sepValuesTwo; // Holds the separation value for step two.

            double initSeparation = 0;
            if (myGeneOutParam.skipInitCalc != 1){
                int firstRun = 1;
                SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
                if (DEBUG_OUTPUT >= 0){
                    cout << "   Calculating primary separation value between group one and group two." << endl;
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "Performing " << myGeneOutParam.numInitCalc << " initial separation calculations and averaging." << endl;
                }
                // Get a cutoff/baseline, pc.
                // I would also like the ability to get an average of some number of runs
                for (int i=0;i<myGeneOutParam.numInitCalc;i++){

                    SVM_separationResults mySeparationResults;
                    int initPrimSepStartTime = time(0);
                    // Dangerous, because we cant see any output, but necessary.
                    int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                    if (DEBUG_OUTPUT <= 0) {
                        DEBUG_OUTPUT = -1;
                    }
                    calcSVMseparation(treeFileNamesGroupOne,treeFileNamesGroupTwo,myGeneOutParam.SampleParam,mySeparationResults);
                    DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                    if (firstRun == 1){
                        averagedResults = mySeparationResults;
                        firstRun = 0;
                    }
                    else {
                        averagedResults += mySeparationResults;
                    }
                    int initPrimSepEndTime = time(0);
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   " << initPrimSepEndTime - initPrimSepStartTime << " seconds." << endl;
                    }

                    //double thisInitSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                    double thisInitSeparation = mySeparationResults.resampleSeparationPercentage;
                    if (DEBUG_OUTPUT >= 0){
                        cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                        cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                        //cout << "       Separation percentage: " << setw(6) << 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize) << endl;
                        cout << "       This Separation percentage(" << i + 1 << "): " << setw(6) << thisInitSeparation << endl;
                        if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                            mySeparationResults.print("This Initial");
                            if (myGeneOutParam.SampleParam.projectSVD == 1){
                                mySeparationResults.printPCA("This Initial");
                            }
                        }
                    }
                    initSeparation+= thisInitSeparation;
                    sepValuesOne.push_back(thisInitSeparation);
                }
                initSeparation /= myGeneOutParam.numInitCalc;
                averagedResults /= myGeneOutParam.numInitCalc;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Averaged Initial Separation percentage: " << setw(6) << initSeparation << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        averagedResults.print("Averaged Initial");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            averagedResults.printPCA("Averaged Initial");
                        }
                    }
                }
            }
            else if (myGeneOutParam.skipInitCalc){
                initSeparation = myGeneOutParam.mySeparation;
                if (DEBUG_OUTPUT >= 0){
                    cout << "Setting initSeparation = " << myGeneOutParam.mySeparation << endl;
                }
            }

            if (DEBUG_OUTPUT >= 0){
                cout << endl;
                cout << "Now performing " << myGeneOutParam.numStatTests << " statistical tests." << endl;
            }

            // if (separation > pc) pval++
            // return pval/myGeneOutParam.numStatTests
            double pval = 0;
            double averageSeparation=0;
            int firstRun = 1;
            SVM_separationResults averagedResults; //Construct makes all non-matrix items 0.
            for (int i=0;i<myGeneOutParam.numStatTests;i++){
                SVM_separationResults mySeparationResults;
                int innerStatTestStartTime = time(0);
                // Clear previous Alignments.

                int old_DEBUG_OUTPUT = DEBUG_OUTPUT;
                if (DEBUG_OUTPUT <= 0) {
                    DEBUG_OUTPUT = -1;
                }
                // Compute separation in treeFileNamesGroupTwo
                calcSVMseparation(treeFileNamesGroupTwo,treeFileNamesGroupTwo,myGeneOutParam.SampleParam,mySeparationResults);
                DEBUG_OUTPUT = old_DEBUG_OUTPUT;
                if (firstRun == 1){
                    averagedResults = mySeparationResults;
                    firstRun = 0;
                }
                else {
                    averagedResults += mySeparationResults;
                }
                //double thisSeparation = 100*(double)(groupOneCount + groupTwoCount)/(2*myGeneOutParam.SampleParam.SVM_resampleSize);
                double thisSeparation = mySeparationResults.resampleSeparationPercentage;
                
                int innerStatTestEndTime = time(0);
                if (DEBUG_OUTPUT >= 0){
                    cout << endl;
                    cout << "       groupOneCount " << mySeparationResults.resampleGroupOneCount << endl;
                    cout << "       groupTwoCount " << mySeparationResults.resampleGroupTwoCount << endl;
                    cout << "       Separation percentage: " << setw(6) << thisSeparation << endl;
                    cout << "   " << innerStatTestEndTime - innerStatTestStartTime << " seconds. Stat test " << i+1 << endl;
                    if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                        mySeparationResults.print("This Stat Test");
                        if (myGeneOutParam.SampleParam.projectSVD == 1){
                            mySeparationResults.printPCA("This Stat Test");
                        }
                    }
                }
                averageSeparation += thisSeparation;
                sepValuesTwo.push_back(thisSeparation);
                if (initSeparation <= thisSeparation){
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value higher than or equal to initSeparation." << endl;
                        cout << "       pval++" << endl;
                    }
                    pval++;
                }
                else {
                    if (DEBUG_OUTPUT >= 0){
                        cout << "   Separation value lower than initSeparation." << endl;
                    }
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "   ### Current(" << i+1 << ") pval = " << setw(3) << pval/(i+1) << " ###" << endl;
                }

            }
            averageSeparation/= myGeneOutParam.numStatTests;
            averagedResults /= myGeneOutParam.numStatTests;
            if (DEBUG_OUTPUT >= 0){
                cout << "pval/myGeneOutParam.numStatTests = " << pval/myGeneOutParam.numStatTests << endl;
                cout << "p-value: " << pval/myGeneOutParam.numStatTests << endl;
                cout << "Averaged Separation Percentage: " << averageSeparation << endl;
                if (myGeneOutParam.SampleParam.doDiffMeans == 1){
                    averagedResults.print("Averaged Stat Test");
                    if (myGeneOutParam.SampleParam.projectSVD == 1){
                        averagedResults.printPCA("Averaged Stat Test");
                    }
                }
                cout << "Init separation values: ";
                for (list <double>::const_iterator lit=sepValuesOne.begin();lit!=sepValuesOne.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;
                cout << "Step two separation values: ";
                for (list <double>::const_iterator lit=sepValuesTwo.begin();lit!=sepValuesTwo.end();lit++){
                    cout << *lit << setw(5) << " ";
                }
                cout << endl;
//                cout << "MWW pvalue: " << MannWhitneyWilcoxon(sepValuesOne,sepValuesTwo) << endl;
//                cout << "Min pvalue: " << pvalueMin(sepValuesOne, sepValuesTwo) << endl;
            }
        }
    } // END STAT TEST
    cout << endl;
    int mainEndTime = time(0);
    cout << "Total run time (wall clock): " << mainEndTime - mainStartTime << endl;
    cout << endl;
}
