// $Rev: 757 $ $Date: 2010-09-22 15:39:42 -0400 (Wed, 22 Sep 2010) $

/** \file sampleparams.h */
#ifndef SAMPLEPARAMS_H
#define SAMPLEPARAMS_H 1

#include <stdlib.h>

using namespace::std;

/// Data structure to hold parameters for SVM sampling of trees
class SampleParameters {
public:
    // Paramters to control if SVM and/or difference of means is done
    int         doSVM;
    int         doDiffMeans;
    int         SVM_sampleUniform; // 0 means sample uniformly frmo mixture model. 1 means sample sample uniformly from each file.
                                    // If this is 1, we need to check SVM_sampleSize and SVN_resampleSize.
    int         SVM_sampleSize;// -1 // Default
    int         SVM_resampleSize;// -1 // Default
    int         projectSVD;// = 0; // Do no project with svd by default
    double      projectSVDcutOff;// = 1; // Default value for svd singular value projection
    int         distanceOne;// = 0; // 0 means use distance specified in newick, 1 means use 1.
    int         scaleToOne;// = 0; // 0 means no scaling, 1 means scale treeLength to 1
    int         modelType; // default should be LINEAR 
    int         numTreesPerFile; // The number of trees to expect in each file
    // Parameters for the burnin type and amount. This really has nothing to do with mr bayes
    int         burninFormat;// = 0;    // 0 will mean percentage, 1 will mean fixed number.
    double      burninPercent;// = 0.25; // Assume 25% burnin.
    int         burninNumber;// = 0;    // Default is 0 since we assume default of 25% burnin
private:
};

#endif
