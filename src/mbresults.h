#ifndef MBRESULTS_H
#define MBRESULTS_H 1

#include <stdlib.h>
#include <string>
#include <list>

using namespace::std;

/// Data structure to hold results for running Mr Bayes
class MrBayesResults {
public:
    MrBayesParameters       parametersCopy;
    list <double>           splitFreqs; 

    void clear() {splitFreqs.clear();};

private:
};


#endif

