// $Rev: 808 $ $Date: 2011-04-10 13:32:07 -0400 (Sun, 10 Apr 2011) $

// \file newickworker.h

#ifndef NEWICKWORKER_H
#define NEWICKWORKER_H 1

#include <string>
#include <stdlib.h>
#include <map>
#include <list>
#include <set>
#include <math.h>
#include "alignments.h"
#include "debugoutput.h"
#include "matrix.h"

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif 


#ifdef HAVE_LIBSVM
#include "svm.h"
#endif

/// Takes in a newick format string and does some action. Virtual class.

/** This is the base class, not meant to be used by itself. Other classes
    are to be derived from this one to perform 'work' such as, taking sum,
    and reporting the mean. There can be any number of uses. */
class NewickWorker {
public:
    // Empty constructor
    NewickWorker() {};
    // Empty deconstructor
    ~NewickWorker() {};
    // This is the function that does the real work. Declared virtual
    /** This is one of two main functions for this code. It will take in
        a string in newick format and do something to it. What it does
        is up to the derived class */
    virtual void work(string) {}; // This takes in a Matrix and does some action.
    // This is the function that does the real work. Declared virtual
    /** This is one of two main functions for this code. It will take in
        a string in newick format and do something to it. What it does
        is up to the derived class. For some NewickWorkers it is more
        natural to work on tuples. */
    virtual void work(string,string) {}; // This takes in a Matrix and does some action.
    virtual void clear() {}; /// Clear out the data structure
protected:
};

#ifdef HAVE_LIBSVM
class NewickWorker_getSVMList : public NewickWorker {
public:
    NewickWorker_getSVMList() {};
    ~NewickWorker_getSVMList() {};
    virtual void work(string);
    virtual void work(string,string) {};
    void clear(); ///< Clear out the data structure
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
    list <svm_node *> getNodeList ();
protected:
    list <svm_node *> nodeList;
};

class NewickWorker_getSVMList_topology : public NewickWorker_getSVMList {
public:
    NewickWorker_getSVMList_topology() {};
    ~NewickWorker_getSVMList_topology() {};
    void work(string);
    void work(string,string) {};
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
protected:
};

class NewickWorker_getSVMList_scaleToOne : public NewickWorker_getSVMList {
public:
    NewickWorker_getSVMList_scaleToOne() {};
    ~NewickWorker_getSVMList_scaleToOne() {};
    void work(string);
    void work(string,string) {};
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
protected:
};
#endif

class NewickWorker_print : public NewickWorker {
public:
    NewickWorker_print() {};
    ~NewickWorker_print() {};
    void work(string);
    void work(string,string) {};
    void clear() {}; ///< Clear out the data structure
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
protected:
    ostream *outputFile;
    int     printFormat; ///< 0 means whitespace deliminated, one tree per line.
};

class NewickWorker_print_topology : public NewickWorker {
public:
    NewickWorker_print_topology() {};
    ~NewickWorker_print_topology() {};
    void work(string);
    void work(string,string) {};
    void clear() {}; ///< Clear out the data structure
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
protected:
    ostream *outputFile;
    int     printFormat; ///< 0 means whitespace deliminated, one tree per line.
};

class NewickWorker_print_splits : public NewickWorker {
public:
    NewickWorker_print_splits() {};
    ~NewickWorker_print_splits() {};
    void work(string);
    void work(string,string) {};
    void clear() {}; ///< Clear out the data structure
    void init(ostream &out,int setPrintFormat); ///< At each work step, will output tree.
protected:
    ostream *outputFile;
    int     printFormat; ///< 0 means whitespace deliminated, one tree per line.
};

class NewickWorker_print_l2_diff_topology : public NewickWorker {
public:
    NewickWorker_print_l2_diff_topology() {};
    ~NewickWorker_print_l2_diff_topology() {};
    void work(string);
    void work(string,string) {};
    void clear() {}; ///< Clear out the data structure
    void init(Matrix A, Matrix B, ostream &out); ///< At each work step, will output l2 difference of item with A and B to out.
protected:
    Matrix M1, M2;
    ostream *outputFile;
};

/// This will simply sum all distance vectors it works on.
/// It can return the l2 norm
class NewickWorker_mean : public NewickWorker {
public:
    NewickWorker_mean() {numStringsProcessed = 0;};
    long double l2norm();
    Matrix mean();
    void work(string); // Add Matrix to totalSum
    void work(string A,string B);  // Add A - B to totalSum
    void clear ();
protected:
    int numStringsProcessed;
    Matrix totalSum; // This is the sum of all vectors processed.
};

/// This will store the sum of the l2 norms of distance vectors it works on.
/// It can return the mean of the sum.
class NewickWorker_mean_l2norm : public NewickWorker {
public:
    NewickWorker_mean_l2norm() {numStringsProcessed = 0;};
    long double mean ();
    void work(string A); // Add l2 norm of the Matrix to totalSum
    void work(string A,string B);  // Add l2 norm of A - B to totalSum
    void clear ();
protected:
    int numStringsProcessed;
    long double totalSum; // This is the sum of all vectors processed.
};

/// This will simply sum all topological distance vectors it works on.
/// It can return the l2 norm
class NewickWorker_mean_topology : public NewickWorker {
public:
    NewickWorker_mean_topology() {numStringsProcessed = 0;};
    long double l2norm();
    Matrix mean();
    void work(string); // Add Matrix to totalSum
    void work(string A,string B);  // Add A - B to totalSum
    void clear ();
protected:
    int numStringsProcessed;
    Matrix totalSum; // This is the sum of all vectors processed.
};

/// This will store the sum of the l2 norms of topological distance vectors it works on.
/// It can return the mean of the sum.
class NewickWorker_mean_topology_l2norm : public NewickWorker {
public:
    NewickWorker_mean_topology_l2norm() {numStringsProcessed = 0;};
    long double mean ();
    void work(string A); // Add l2 norm of the Matrix to totalSum
    void work(string A,string B);  // Add l2 norm of A - B to totalSum
    void clear ();
protected:
    int numStringsProcessed;
    long double totalSum; // This is the sum of all vectors processed.
};

#endif
