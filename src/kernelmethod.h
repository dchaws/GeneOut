// $Rev: 776 $ $Date: 2011-02-10 12:00:19 -0500 (Thu, 10 Feb 2011) $

/** \file kernelmethod.h */
#ifndef KERNELMETHOD_H
#define KERNELMETHOD_H 1

#include <string>
#include <stdlib.h>
#include "svm.h"
#include <map>
#include <list>
#include <set>
#include <math.h>
#include "alignments.h"
#include "debugoutput.h"
#include "matrix.h"
#include "newickworker.h"
#include "tempprefix.h"
#include "mbparams.h"
#include "jkparams.h"
#include "multindparams.h"
#include "sampleparams.h"
#include "svmresults.h"
#include "nexustreefile.h"
#include "treefuncs.h"
#include "mbresults.h"

extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double A[], int *lda, double s[], double U[], int *ldu, double VT[], int *ldtv, double Work[], int *lwork, int *info);

//// Max number of nodes, internal and leaves. If the number of nodes exceeds
//// this amount, the behavior is unknown and will most likely cause a seg fault.
#define MAX_NUM_NODES 300

// Used to track the recursion level for NewickToDistanceMatrix
extern int recurseLevel;

using namespace::std;

/// This takes in a nodeCount x nodeCount square distance matrix and a translation
/// array leafNumberToNodeNumber and leafCount. It returns a leafCount*(leafCount-1)/2 
/// size vector of the distances indexed as 
/// d[leaf 1][leaf 2], d[leaf 1][leaf 3], ... ,  d[leaf 1][leaf leafCount], d[leaf 2][leaf 3], ...
/// i.e. lexicographic order
double *MatrixToVector(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount);

/// This takes in a nodeCount x nodeCount square distance matrix and a translation
/// array leafNumberToNodeNumber and leafCount. It returns a leafCount*(leafCount-1)/2 
/// size vector of the distances indexed as 
/// d[leaf 1][leaf 2], d[leaf 1][leaf 3], ... ,  d[leaf 1][leaf leafCount], d[leaf 2][leaf 3], ...
/// i.e. lexicographic order. This returns the distance vector as a sparse vector in
/// svm_node format.
svm_node *MatrixTo_svm_node_Vector(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount);

/// This takes in a nodeCount x nodeCount square distance matrix and a translation
/// array leafNumberToNodeNumber and leafCount. It returns a leafCount*(leafCount-1)/2 
/// size vector of the distances indexed as 
/// d[leaf 1][leaf 2], d[leaf 1][leaf 3], ... ,  d[leaf 1][leaf leafCount], d[leaf 2][leaf 3], ...
/// i.e. lexicographic order. This returns the distance vector as a Matrix.
Matrix distanceMatrixToMatrix(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount);

/// This recursive function takes in the string treeString which is in Newick format.
/// It parses the string accordingly and recurses when it encounters nested '(' ')'.
/// beg points to the location of the treeString for which this call of NewickToDistanceMatrix
/// is to parse.
/// nodeCount should be 0 when first called and it will track the total number of nodes seen.
/// leafcount should be 0 and it will track the total number of leaves.
/// transition and D are assumed to contain all -1's. 
/// D will be the distance matrix once the treeString is fully parsed. D[i][j] will be the distance
/// from node i to node j. transition[i] is the parent of node i.
/// leafNumbertoNodeNumber is a matrix where leafNumberToNodeNumber[i] is the index in D of leaf i.
/// distanceOne = 0 means use the distances specified in the newick file, 1 means
/// assume distance 1.
/// If treeLengthFactor is non-zero, then assume that we want to scale each distance seen
/// by (distance)*treeLengthFactor;
void NewickToDistanceMatrix (string &treeString, int beg, double D[MAX_NUM_NODES][MAX_NUM_NODES],\
        int &nodeCount, int &leafCount, int transition[MAX_NUM_NODES],\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int distanceOne, double &treeLengthFactor,\
        list <set <string> > &splits);

/** permutes the elements of an svm_node */
void permute_svm_node(svm_node *tempNode);

/// This does the hard work. It takes in a list of file names which it expects
/// to be in mr bayes output format. It then chooses SVM_samplesize many trees
/// uniformly from the files, considering the burn in for each file.
/// This will return a list of svm_nodes compatible with svn_train input.
 list <svm_node *> *readNexusTreesUniform(list <string> treeFileNames, int SVM_sampleSize, int numTrees, int burninFormat, double burninPercent, int burninNumber, int distanceOne, int scaleToOne);

///
/// This will read through treeFileNames all tree indexed by treesToRead,
/// considering of course any burn in specified in tsp. For every tree,
/// this function will call someNewickWorker. Trees are indexed starting
/// with 0.
///
void readNexusTreesUniform(list <string> treeFileNames, SampleParameters &tsp,\
        set <unsigned> treesToRead, NewickWorker &someNewickWorker);

///
/// This will read through treeFileNames all tree tuples (T1, T2) where
/// T1 in treesToReadOne and T2 in treesToReadTwo,
/// considering of course any burn in specified in tsp. For every tuple,
/// this function will call somNewickWorker. Trees are indexed starting
/// with 0.
///
void readNexusTreesUniformAllPairs(list <string> treeFileNames, SampleParameters &tsp,\
        set <unsigned> treesToReadOne, set <unsigned> treesToReadTwo, NewickWorker &someNewickWorker);

/// This now can handle the case when we are sampling every single tree. In this case we use
/// NewickWorkers's and set the trees to sample
void Sample_SVM(list <string> treeFileNamesGroupOne, list <string> treeFileNamesGroupTwo, SampleParameters &tsp, list <svm_node *> &groupOne, list <svm_node *> &groupTwo);

/// This function will take in a list of svm_nodes and will project them using
/// Principal component analysis via SVD calculation. This function calls
/// the dgesvd function from LAPACK to calculate the SVD.
/// It will return a list of svm_node pointers of the projected points.
/// The vectors are placed as the rows of a matrix A.
/// The SVD decomposition is calculated: A = USV^T
/// This will return U_l S_l where U_l is the first l columns of U and
/// S_l is the first l rows of S. l is determined by the number of
/// singular values greater than cutOffValue
/// This function returns V_l for use as a projection matrix
/// To project the original matrix A, perform the operation A*V_l = U_l*S_l
/// Every element of the return list is a row of V
list <svm_node *> *project_SVM_nodesMatrix (list <svm_node *> &svmVecs, double cutOffValue, int numSingValues);

/// This function will take in a list of svm_nodes and will project them using
/// Principal component analysis via SVD calculation. This function calls
/// the dgesvd function from LAPACK to calculate the SVD.
/// It will return a list of svm_node pointers of the projected points.
/// The vectors are placed as the rows of a matrix A.
/// The SVD decomposition is calculated: A = USV^T
/// This will return U_l S_l where U_l is the first l columns of U and
/// S_l is the first l rows of S. l is determined by the number of
/// singular values greater than cutOffValue
list <svm_node *> *project_SVM_nodes (list <svm_node *> &svmVecs, double cutOffValue, int numSingValues);

/// This will actually call and run Mr Bayes using the parameters in tmbp.
/// This function will fill treeFileNames with a list of a list of filenames for
/// each file, since for each file there can be multiple run files
void getTreesMrBayes(MrBayesParameters &tmbp, list <list <string> > &treeFileNames, MrBayesResults &MB_results);

void getTreesMrBayes(MrBayesParameters &tmbp,list <Alignments *> &inputAlignments, list <list <string> > &treeFileNames, MrBayesResults &MB_results);

/// This takes in tjp (ignoring nexInputFiles) and uses inputAlignments.
/// It will create the trees using jackknife for each Alignments and place
/// them in outputFileNames
void getTreesJackknife(JackknifeParameters &tjp,list <Alignments *> &inputAlignments, list <string> outputFileNames);

/// This will only sample tsp.sampleSize many trees uniformly from the inputAlignments
/// This will create outputFileNames, which will be one file.
void getTreesJackknife(JackknifeParameters &tjp, SampleParameters &tsp, list <Alignments *> &inputAlignments, string &outputFileName);

void getTreesJackknife(JackknifeParameters &tjp, list <list <string> > &treeFileNames);

/// Fill treeFileNames with random names
void getTreesJackknife(JackknifeParameters &tjp, list <Alignments *> &inputAlignments, list <list <string> > &treeFileNames);

/// This takes in tjp (ignoring nexInputFiles) and uses inputAlignments.
/// It will create trees by randomly selecting representative individuals and place
/// the resulting trees in outputFileNames
void getTreesMultInd(MultIndParameters &tmip,list <Alignments *> &inputAlignments, list <string> outputFileNames);

/// This function actually does the PCA projection using the empirical mean and the projection
/// matrix.
void doPCA(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, list <svm_node *> *projectionMatrix, svm_node *empiricalMean);

/// This functions gets the empirical mean and the projection matrix.
/// It does not do the actual calculation
void getPCA_data(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, double projectSVDcutOff, list <svm_node *> **projectionMatrix, svm_node **empiricalMean);

/// This is the actual function. If projectSVDcutOff == -1 then we take numSingValues
void getPCA_data(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, double projectSVDcutOff, list <svm_node *> **projectionMatrix, svm_node **empiricalMean, int numSingValues);

svm_problem form_svm_problem(list <svm_node * > &groupOne, list <svm_node * > &groupTwo);


svm_parameter getLinear_svm_parameter();

void delete_svm_problem(svm_problem &prob);

void get_svm_predictions(svm_problem &myProblem, svm_model &myModel, int &groupOnePredictCorrectCount, int &groupTwoPredictCorrectCount);

void calcSVMseparationJackknife(list <string> inputFileNames, int numGroupOne, JackknifeParameters &tjp, SampleParameters &tsvmp, SVM_separationResults &results);

void calcSVMseparationJackknife(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, JackknifeParameters &tjp, SampleParameters &tsvmp, SVM_separationResults &results);

/// This will use AlignementOne and AlignemntsTwo for input
void calcSVMseparationMrBayes(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, MrBayesParameters &tmbp, MrBayesResults &MB_results, SampleParameters &tsvmp, SVM_separationResults &results);

/// This will use the filenames specified in MrBayesParamters as input
void calcSVMseparationMrBayes(MrBayesParameters &tmbp, MrBayesResults &MB_results, SampleParameters &tsvmp, SVM_separationResults &results, int numGroupOne);

/// This will use AlignementOne and AlignemntsTwo for input
void calcSVMseparationMultInd(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, MultIndParameters &tmip, SampleParameters &tsvmp, SVM_separationResults &results);

/// This will use the filenames specified in MrBayesParamters as input
void calcSVMseparationMultInd(list <string> inputFileNames, MultIndParameters &tmip, SampleParameters &tsvmp, SVM_separationResults &results, int numGroupOne);

/// This does the work of actually sampling and resampling and giving the groupOneCount and groupTwo count according to resampleSize
void calcSVMseparation(list <string> &treeFileNamesGroupOne, list <string> &treeFileNamesGroupTwo,list <string> &treeFileNamesGroupOneResample, list <string> &treeFileNamesGroupTwoResample, SampleParameters &tsp, SVM_separationResults &results);

/// This calls the above function, but using the same treeFileNames for the sample and resample
void calcSVMseparation(list <string> &treeFileNamesGroupOne, list <string> &treeFileNamesGroupTwo, SampleParameters &tsp, SVM_separationResults &results);

/// This will will use origAlignments to fill AlignmentsOne with groupOneSize jackknifes, and 
/// AlignmentsLengths.size() - groupOneSize jackknifes in AlignmentsTwo. It will create jackknife 
/// alignments of length specified in AlignmentsLengths.
/// The new alignments will be added to to AlignmentsOne and AlignmentsTwo, according
/// to groupOneSize, i.e. groupOneSize will go into AlignmentsOne.
void getSomeJackknifeAlignments(Alignments &origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize);

/// This will use origAlignments to fill AlignmentsOne with jackknifes using
/// AligmentsLengths as a guide for the number of characters and the number of
/// alignments total;
void getSomeJackknifeAlignments(Alignments &origAlignments, list <Alignments *> &AlignmentsOne, list <int> &AlignmentsLengths, int colSize);

/// This will use origAlignments to fill AlignmentsOne with jackknifes using
/// AligmentsLengths as a guide for the number of characters and the number of
/// alignments total. When generating each new alignment, only one randomly
/// selected alignment of origAlignments is used.
void getSomeJackknifeAlignments(list <Alignments *> &origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize);


/// 9/1/10 New function to take in all original alignments, permute them (sample without replacement) then bootstrap to the appropriate size
/// This function assumes that the first groupOneSize alignments in origAlignments are the original first group.
/// allowAnyPermutation 1, means allow any permuation, 0 means do not allow AlignmentsOne to be alignments all from the original group one.
void getSomeJackknifeAlignmentsPermute(list <Alignments *> origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize, int allowAnyPermutation);

/// This will call the command specified by input to generate new alignments
/// The word 'uniform' refers to the fact that the program will be used for
/// all new alignments created.
/// We will fill AlignmentsOne with the first alignment then the rest will
/// go into AlignmentsTwo
/// We expect to use genCommand as 'genCommand <filename>'
/// The new alignments will be added to to AlignmentsOne and AlignmentsTwo, according
/// to groupOneSize, i.e. groupOneSize will go into AlignmentsOne.
void userGenerateNewAlignmentsUniform(string genCommand, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, int numNewAlignments, int gropuOneSize);

/// This calculates the mean of the vectors in groupOne
svm_node *Mean(list <svm_node *> &groupOne);

/// This calculates the difference of means between groupOne and groupTwo
/// This actually calculates Mean(groupOne) - Mean(groupTwo);
svm_node *DifferenceOfMeans(list <svm_node *> &groupOne, list <svm_node *> &groupTwo);
Matrix DifferenceOfMeansMatrix(list <svm_node *> &groupOne, list <svm_node *> &groupTwo);

double DifferenceOfMeans_l2(list <svm_node *> &groupOne, list <svm_node *> &groupTwo);

/// This returns the variance of vectors in groupOne
svm_node *Variance(list <svm_node *> &groupOne);

Matrix SampleVariance(list <svm_node *> &groupOne);

/// Returns a column vector in the form of Matrix
Matrix svm_node_to_Matrix(svm_node *);

/// Returns a svm_node copied from Matrix
svm_node *Matrix_to_svm_node(Matrix &);

/// Adds all values low to high, inclusive to someSet.
void fillUnsignedSet(set <unsigned> &someSet, unsigned low, unsigned high);

/// Adds all values low to high, inclusive to someSet, excludes values from excludeSet
void fillUnsignedSetExclude(set <unsigned> &someSet, unsigned low, unsigned high, set <unsigned> &excludeSet);

/// Takes in a list of list of strings and returns a list of strings
list <string> listListStringToListString(list <list <string> > &TreeListListFileNames);

/// This will take in two lists and return #{x > min(groupOne) | x \in groupTwo }/groupTwo.size()
double pvalueMin(list <double> groupOne, list <double> groupTwo);

/// GCD function copied from web.
int GCD(int x,int y);

/// LCM function copied from web.
int lcm(int a,int b);

#endif
