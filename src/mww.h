
/** \file mww.h */
#ifndef MWW_H
#define MWW_H 1

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <list>
#include <set>
#include <math.h>
#include "debugoutput.h"

using namespace::std;

typedef struct doubleIndexStruct {
    double  value;
    int     index;
    double  rank;
} doubleIndex;

struct ltDoubleIndex {
    bool operator()(const doubleIndex doubleIndexOne, const doubleIndex doubleIndexTwo) const {
        if ((double)doubleIndexOne.value < (double)doubleIndexTwo.value){
            return 1;
        }
        else {
            return 0;
        }
    }
};

/// Returns the sum of integers 1 to n using formula n(n+1)/2. Uses integer math.
int sumUpTo(int n);

/// Returns the cumulative normal distribution using "The most used algorithm is algorithm 26.2.17
/// from Abromowitz and Stegun, Handbook of Mathematical Functions. It has a maximum absolute error of 7.5e^-8." 
double normalCumDist(const double x);

/// Impliments Mann-Whitney-Wilcoxon non-parametric statistical test. Returns the p-value
double MannWhitneyWilcoxon(const list <double> sampleOne, const list <double> sampleTwo);
#endif
