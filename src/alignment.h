// $Rev: 774 $ $Date: 2011-02-09 17:53:15 -0500 (Wed, 09 Feb 2011) $

/** \file alignments.h */
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

class Alignments 
{
public:
    Alignments ();
    //Alignments (Alignments &someAlignments);
    Alignments (int tmp_ntax, int tmp_nchar, string tmp_AlignmentsName, list <string> tmp_taxaNames, char **tmp_alignments, int tmp_alignmentsAllocated );
    ~Alignments ();
    Alignments (string fileName); // opens the filename string
    Alignments (string fileName, int setInputFormat); // opens filename with said input format.
    Alignments (std::istream &in);
    void readFile (string fileName); 

    friend std::ostream& operator << (std::ostream &out, const Alignments &someAlignments);
    friend std::istream& operator >> (std::istream &in, Alignments &someAlignments);

    Alignments & operator = (const Alignments &someAlignments);
    bool operator == (const Alignments &someAlignments);
    // Concat, x + y concats alignments x + y and returns the new alignments
    Alignments operator + (const Alignments &someAlignmens);
    // Concats someAlignments to this
    Alignments & operator += (const Alignments &someAlignments);
    char & operator() (unsigned row, unsigned col); 
    char operator() (unsigned row, unsigned col) const; 

    const int get_ntax();
    const int get_nchar();
    const list <string> get_taxaNames();
    const string get_AlignmentsName();
    void setOutputFormat (int format);
    void setInputFormat (int format);

    // Computes the sequence divergence as min pairwise hamming distance divided by number of base pairs.
    double sequenceDivergenceMin();
    // Computes the sequence divergence as avg pairwise hamming distance divided by number of base pairs.
    double sequenceDivergenceAvg();
    // Computes the pairwise sequence divergence hamming distance divided by number of base pairs.
    void printSequenceDivergencePairs();


    // This return k contigious columns of alignments (with wraping) start at
    // position pos
    Alignments getContiguousColumns(int pos, int k);

    // This will return an alignment with ntax taxa and k*count columns.
    // It will form it by randomly selecting k contiguous columns and appending them
    // to the result
    Alignments getJackknife(int k, int count);
    // This will do the same as above but count = floor(nchar/k);
    Alignments getJackknife(int k);

    // This will go through all elements of alignments and delete any columns
    // that do not contain 'A', 'T', 'C' or 'G'
    Alignments delColumnsWithMissingData();

    // Returns the subset of characters given by someTaxa
    Alignments getTaxaSubset(set <unsigned> someTaxa);

private:
    void init();
    // Delete the alignments and set alignments to 0
    void deleteAlignments();

    // Allocates alignments data structure to be char[nchar][ntax]
    void allocateAlignments();
    // Allocates alignments data structure to be char[tmp_nchar][tmp_ntax]
    void allocateAlignments(int tmp_ntax, int tmp_nchar);
    string AlignmentsName;
    list <string> taxaNames;
    int ntax; // Number of taxa
    int nchar; // Number of characters or base pairs
    int outputFormat; // See above for formats
    int inputFormat;  // See above for formats

    // alignments when allocated is char[nchar][ntax]
    // We store in this format to facilitate jackknifing, where we take subsets 
    // of columns;
    char **alignments;

    int alignmentsAllocated; // 0 means not allocated, 1 means allocated.
                             // Should only be set by allocate and delete functions.
};

#endif
