// $Rev: 433 $ $Date: 2009-11-24 01:44:21 -0500 (Tue, 24 Nov 2009) $

/** \file svmresults.h */
#ifndef SVMRESULTS_H
#define SVMRESULTS_H 1

#include <stdlib.h>
#include <iomanip>
#include <string>
#include <list>
#include <math.h>
#include "matrix.h"

using namespace::std;

///
/// Data structure to hold the results of SVM separation.
/// Also contains information on difference of means, variance.
///
class SVM_separationResults{
public:
    SVM_separationResults();
    void        print(); //< This will print all contained non-PCA data
    void        print(const string prefix);
    void        print(const char *);
    void        printPCA(); //< This will print all contained PCA data
    void        printPCA(const string prefix);
    void        printPCA(const char *);
    int         groupOneCount;
    int         groupTwoCount;

    int         resampleGroupOneCount;
    int         resampleGroupTwoCount;

    double      separationPercentage;
    double      resampleSeparationPercentage;

    int         numSingularValuesTaken;

    Matrix      diffMeans;
    double      diffMeans_l2_norm;

    Matrix      varianceGroupOne;
    double      traceVarianceGroupOne;
    long double      detVarianceGroupOne;
    long double      logDetVarianceGroupOne;

    Matrix      varianceGroupTwo;
    double      traceVarianceGroupTwo;
    long double      detVarianceGroupTwo;
    long double      logDetVarianceGroupTwo;

    Matrix      diffMeansPCA;
    double      diffMeans_l2_normPCA;

    Matrix      varianceGroupOnePCA;
    double      traceVarianceGroupOnePCA;
    long double      detVarianceGroupOnePCA;
    long double      logDetVarianceGroupOnePCA;

    Matrix      varianceGroupTwoPCA;
    double      traceVarianceGroupTwoPCA;
    long double      detVarianceGroupTwoPCA;
    long double      logDetVarianceGroupTwoPCA;

    int         secondsToCompute;

    SVM_separationResults & operator+=(const SVM_separationResults &rhs);
    SVM_separationResults & operator*=(const double &rhs); 
    SVM_separationResults & operator/=(const double &rhs);
    SVM_separationResults & operator= (const SVM_separationResults& rhs);   // Assignment operator

private:
};

#endif
