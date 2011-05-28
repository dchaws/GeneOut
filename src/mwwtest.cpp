
/// \file mwwtest.cpp 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <map>
#include <list>
#include <set>
#include <math.h>
#include "mww.h"
#include "debugoutput.h"


// Used to control debug output. 
int DEBUG_OUTPUT = 3;

using namespace::std;

/// Usage: This program first reads in two positive integers, which specify
/// the number of individuals in sample one and number of individuals in
/// sample two respectively. What should follow are sample one and sample two.
int main (int argc, char **argv)
{
    int numSampleOne, numSampleTwo;
    cin >> numSampleOne;
    cin >> numSampleTwo;

    list <double> sampleOne;
    double  tmpDouble;
    for (int i=0;i<numSampleOne;i++){
        cin >> tmpDouble;
        sampleOne.push_back(tmpDouble);
    }

    list <double> sampleTwo;
    for (int i=0;i<numSampleTwo;i++){
        cin >> tmpDouble;
        sampleTwo.push_back(tmpDouble);
    }

    cout << "P-value: " << MannWhitneyWilcoxon(sampleOne,sampleTwo) << endl;
}

