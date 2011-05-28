
/// \file geneoutparameterstest.cpp 

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
#include "debugoutput.h"
#include "geneoutparams.h"

// Used to control debug output. 
int DEBUG_OUTPUT = 3;

string tempPrefix = "";

using namespace::std;

int main (int argc, char **argv)
{
    GeneOutParameters myGeneOutParameters(cin);
    cout << myGeneOutParameters << endl;
}


