// $Rev: 774 $ $Date: 2011-02-09 17:53:15 -0500 (Wed, 09 Feb 2011) $

/** \file alignment.h */
#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H 1

#include <string>
#include <list>
#include <set>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>

#define OF_STANDARD 1
#define OF_NEXUS 2
#define OF_PHYLIP 3

// Input format.
#define IF_NEXUS 1
#define IF_PHYLIP 2

using namespace::std;

///
/// Class to hold a set of aligned sequence of characters.
///

class Alignment 
{
public:
    Alignment ();
    //Alignment (Alignment &someAlignment);
    Alignment (int tmp_ntax, int tmp_nchar, string tmp_AlignmentName, list <string> tmp_taxaNames, char **tmp_alignment, int tmp_alignmentAllocated );
    ~Alignment ();
    Alignment (string fileName); // opens the filename string
    Alignment (string fileName, int setInputFormat); // opens filename with said input format.
    Alignment (std::istream &in);
    void readFile (string fileName); 

    friend std::ostream& operator << (std::ostream &out, const Alignment &someAlignment);
    friend std::istream& operator >> (std::istream &in, Alignment &someAlignment);

    Alignment & operator = (const Alignment &someAlignment);
    bool operator == (const Alignment &someAlignment);
    // Concat, x + y concats alignment x + y and returns the new alignment
    Alignment operator + (const Alignment &someAlignmens);
    // Concats someAlignment to this
    Alignment & operator += (const Alignment &someAlignment);
    char & operator() (unsigned row, unsigned col); 
    char operator() (unsigned row, unsigned col) const; 

    const int get_ntax();
    const int get_nchar();
    const list <string> get_taxaNames();
    const string get_AlignmentName();
    void setOutputFormat (int format);
    void setInputFormat (int format);

    // Computes the sequence divergence as min pairwise hamming distance divided by number of base pairs.
    double sequenceDivergenceMin();
    // Computes the sequence divergence as avg pairwise hamming distance divided by number of base pairs.
    double sequenceDivergenceAvg();
    // Computes the pairwise sequence divergence hamming distance divided by number of base pairs.
    void printSequenceDivergencePairs();


    // This return k contigious columns of alignment (with wraping) start at
    // position pos
    Alignment getContiguousColumns(int pos, int k);

    // This will return an alignment with ntax taxa and k*count columns.
    // It will form it by randomly selecting k contiguous columns and appending them
    // to the result
    Alignment getBootstrap(int k, int count);
    // This will do the same as above but count = floor(nchar/k);
    Alignment getBootstrap(int k);

    // This will go through all elements of alignment and delete any columns
    // that do not contain 'A', 'T', 'C' or 'G'
    Alignment delColumnsWithMissingData();

    // Returns the subset of characters given by someTaxa
    Alignment getTaxaSubset(set <unsigned> someTaxa);

private:
    void init();
    // Delete the alignment and set alignment to 0
    void deleteAlignment();

    // Allocates alignment data structure to be char[nchar][ntax]
    void allocateAlignment();
    // Allocates alignment data structure to be char[tmp_nchar][tmp_ntax]
    void allocateAlignment(int tmp_ntax, int tmp_nchar);
    string AlignmentName;
    list <string> taxaNames;
    int ntax; // Number of taxa
    int nchar; // Number of characters or base pairs
    int outputFormat; // See above for formats
    int inputFormat;  // See above for formats

    // alignment when allocated is char[nchar][ntax]
    // We store in this format to facilitate jackknifing, where we take subsets 
    // of columns;
    char **alignment;

    int alignmentAllocated; // 0 means not allocated, 1 means allocated.
                             // Should only be set by allocate and delete functions.
};

#endif
