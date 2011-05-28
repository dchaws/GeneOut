// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//


/** \file matrix.h */
//. Author David Haws

#ifndef MATRIX_H
#define MATRIX_H 

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <set>

/* cannot have gmpxx without gmp */
#ifndef HAVE_LIBGMP
#ifdef HAVE_LIBGMPXX
#undef HAVE_LIBGMPXX
#endif
#endif

#ifdef HAVE_LIBGMPXX
#include <gmpxx.h>
#endif

using namespace::std;

class Matrix 
{
public:
    Matrix();   // Does not allocate Entries. Sets rows = cols = 0
    Matrix(unsigned,unsigned); // Creates mxn matrix
    Matrix(unsigned); // Creates square matrix
    Matrix(unsigned,unsigned,int, int); // Creates mxn matrix filled with random integers ranged x to y
    ~Matrix();        
    Matrix(const Matrix& m);               // Copy constructor
    Matrix& operator= (const Matrix& m);   // Assignment operator

    bool operator== (const Matrix &m);     // Equal operator

    // Indexed 0 - (rows - 1), 0 - (cols - 1)
    double& operator() (unsigned row, unsigned col);
    double  operator() (unsigned row, unsigned col) const;
    Matrix  operator* (const Matrix &B) const;
    Matrix  operator* (const double &rhs) const;
    Matrix  operator+ (const Matrix &B) const;
    Matrix  operator- (const Matrix &B) const;
    friend std::ostream& operator<< (std::ostream& o, const Matrix &someMatrix);
    friend std::istream& operator>> (std::istream& in, Matrix &someMatrix);

    // This returns the transpose of this matrix
    Matrix  transpose();

    // This will output the matrix in matlab format with label 
    void writeMatlab(std::ostream &o, string label);

    //Matrix *subColumns(set<unsigned> &colSet);
    // Returns a submatrix using colSet as columns to include
    Matrix subColumns(const set<unsigned> &colSet);

    //Matrix *subColumns(set<unsigned> &colSet);
    // Returns a submatrix using NOT colSet as columns to include
    Matrix subColumnsDiff(const set<unsigned> &colSet);

    //given a mxn matrix, returns a mx1 matrix which is the sum of each row
    Matrix rowSum();
    
    // Returns the cofactor matrix determined by i and j.
    // rows and cols > 1
    Matrix cofactor(unsigned i, unsigned j);
    
    // Returns the result of Gaussian Elimination 
    //Matrix *GE ();
    Matrix GE ();

#ifdef HAVE_LIBGMPXX 
    // Returns the rank of this matrix. Uses arbitrary precision and expects
    // all entries to be integer. Also will set rank of the matrix
    int GE_gmp ();
    int GE_gmp_verbose ();

#endif
    int rank ();    // Returns the rank of the matrix 

    double trace(); // Computes the trace

    long double det(); // Computes the determinant
    long double detCofactor(); // Computes the determinant using cofactor expansion
    long double detLU(); // Computes the determinant using SVD decomposition
    long double log_abs_det(); // Computes log(abs(det()));

    // This uses dgesvd function of lapack to compute *this = U*S*V.transpose()
    void SVD(Matrix &U, Matrix &S, Matrix &V);
    
    // This uses dgetrf function of lapack to compute *this = P*L*U.
    // returns the sign determined by the number of row transpositions in P
    double LU(Matrix &P, Matrix &L, Matrix &U);

    // This uses the SVD decomposition function above to compute the inverse.
    Matrix inverse();

    // Uses LAPACK and calls dgesvd_. Counts number of singular values
    // that are not < 10^-16 * SIGMA_MAX
    int rank_LAPACK (); 
    const int getCols ();
    const int getRows ();
    //double det();
    // Swaps columns i and j
    void swapCols(unsigned i, unsigned j);
    // Swaps rows i and j
    void swapRows(unsigned i, unsigned j);
    static int printPadLength;

    // If this matrix is a column vector, return the 2-norm. Else exit with error
    double twoNorm ();
    double twoNormSquared ();
    unsigned rows;
    unsigned cols;
private:
    double **Entries;
    void deleteEntries ();
    int matrixRank;
    void initialize(unsigned,unsigned);
};

struct ltcolvec
{
    bool operator()(const Matrix &m1, const Matrix &m2) const
    {
        int allEqual = 1;
        if ((m1.rows) < (m2.rows))
        {
            return 1;
        }
        for(unsigned i=0;i<m1.rows;i++)
        {
            if ((m1)(i,0) != (m2)(i,0))
            {
                allEqual = 0;
            }
            if ((m1)(i,0) < (m2)(i,0))
            {
                return 1;
            }
            else if ((m1)(i,0) > (m2)(i,0))
            {
                return 0;
            }
        }
        if (allEqual == 1)
        {
            return 0;
        }
        return 1;
    }
};

#endif

