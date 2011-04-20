// $Rev: 29 $ $Date: 2009-08-28 14:27:49 -0400 (Fri, 28 Aug 2009) $
#ifndef JACKKNIFE_H
#define JACKKNIFE_H 1

#include <string>
#include <list>

using namespace::std;


char **getAlignments(string nexInputFileName,int &ntax, int &nchar, list <string> &taxaNames);

// This will take in char[m][n] where char[i] is the ith columns (i.e. transposed)
// and delete any columns that do not contain 'A','T','G', or 'C'. The number of rows
// of the returned data will be the same, but it will have newM columns.
char **delColumnsWithMissingData(char **alignments,int m, int n, int &newM);

void writeAlignmentsToNex(ostream &out, char **alignments, int ntax, int nchar);

void writeAlignmentsTo_phylip(ostream &out, char **alignments, int ntax, int nchar);

void writeAlignmentsTo_phylip(ostream &out, char **alignments, int ntax, int nchar, list <string> taxaNames);

void writeAlignmentsToNex(ostream &out, char **alignments, int ntax, int nchar, list <string> taxaNames);

// This takes in char[nchar][ntax] table of alignments and will return.
// a randomly chooses a k-contiguous block of columns and returns it.
char **getRandAlignmentsContCol(char **alignments, int ntax, int nchar, int k);

// This glues the two alignments alignmentsOne[ncharOne][ntaxOne] and alignments[ncharTwo][ntaxTwo]
char **concatAlignments(char **alignmentsOne, char **alignmentsTwo, int ntaxOne, int ncharOne, int ntaxTwo, int ncharTwo);

// This takes in the filename for a .nex file. It will strip out any columns
// that do not contain ATCG. If there are N columns left this will perform
// numJackknifes, random jackknife calculations to create a new alignment of the approprate size 
// Usually numJackknifes = N/k
char **getJackknife(string nexFileName,int k,int numJackknifes, int &return_ntax, int &return_nchar);
char **getJackknife(char **alignments,int ntax, int nchar,int k,int numJackknifes, int &return_ntax, int &return_nchar);

// This deletes alignments[nchar][ntax]
void deleteAlignments(char **alignments,int ntax, int nchar);

// This glues the multiple alignments 
char **concatAlignments(list <char **> alignmentsList, list <int> ntaxList, list <int> ncharList, int &new_ntax, int &new_nchar);

// end jackknife
#endif
