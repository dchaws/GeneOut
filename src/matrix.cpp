// Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University
// of California at Davis, IBM, Inc. All rights reserved.
//
// Distributed under the Eclipse Public License v 1.0. See ../LICENSE


// matrix.h. Author David Haws
// Taken from MOCHA package, which is distributed under the Eclipse 1.0 license.
// See https://projects.coin-or.org/Mocha

#include <iostream>
#include <iomanip>
#include <math.h>
#include "matrix.h"
#include <cstdlib>
#include <ctime>
#include <set>

#define LAPACK_RANK_ZERO 0.000000000000001

using namespace::std;

extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double A[], int *lda, double s[], double U[], int *ldu, double VT[], int *ldtv, double Work[], int *lwork, int *info);

extern "C" void dgetrf_(int *M, int *N, double A[], int *LDA, int IPIV[], int *INFO);

int Matrix::printPadLength = 0;

void Matrix::initialize(unsigned r,unsigned c)
{
    rows=r;
    cols=c;

    Entries = new double *[rows];
    if (Entries == NULL)
    {
        cerr << "Could not allocate Entries in Matrix::Matrix." << endl;
        exit (0);
    }
    for(unsigned i=0;i<rows;i++)
    {
        Entries[i] = new double[cols];
        if (Entries[i] == NULL)
        {
            cerr << "Could not allocate Entries in Matrix::Matrix." << endl;
            exit (0);
        }
        for(unsigned j=0;j<cols;j++)
        {
            Entries[i][j] = 0;
        }
    }
    matrixRank = -1;
}

Matrix::Matrix(unsigned r,unsigned c)
{
    initialize(r,c);
}

Matrix::Matrix(unsigned d)
{
    initialize(d,d);
}

Matrix::Matrix(unsigned r,unsigned c,int lower,int upper)
{
    initialize(r,c);

    srand((unsigned)time(0));

    for (unsigned i=0;i<r;i++)
    {
        for (unsigned j=0;j<c;j++)
        {
            Entries[i][j] = (rand() % (upper - lower)) + lower;
        }
    }
}

Matrix::Matrix()
{
    Entries = 0;
    rows=0;
    cols=0;
    matrixRank = -1;
}

Matrix::~Matrix()
{
    deleteEntries();
}

void Matrix::deleteEntries ()
{
    for(unsigned i=0;i<rows;i++)
    {
        delete[] Entries[i];
    }

    delete[] Entries;
}

Matrix::Matrix(const Matrix& m)               // Copy constructor
{
    rows = m.rows;
    cols = m.cols;
    matrixRank = m.matrixRank;

    Entries = new double *[rows];
    for(unsigned i=0;i<rows;i++)
    {
        Entries[i] = new double[cols];
        for(unsigned j=0;j<cols;j++)
        {
            Entries[i][j] = 0;
        }
    }

    for (unsigned i=0;i<m.rows;i++){
        for(unsigned j=0;j<m.cols;j++){
            Entries[i][j] = m.Entries[i][j];
        }
    }
}

Matrix& Matrix::operator= (const Matrix& m)   // Assignment operator
{
    if (rows != 0 || cols != 0)
    {
        deleteEntries();
    }
    initialize(m.rows,m.cols);
    //rows = m.rows;
    //cols = m.cols;
    //matrixRank = m.matrixRank;

    //Entries = new double *[rows];
    //for(unsigned i=0;i<rows;i++)
    //{
    //    Entries[i] = new double[cols];
    //    for(unsigned j=0;j<cols;j++)
    //    {
    //        Entries[i][j] = 0;
    //    }
    //}

    for (unsigned i=0;i<m.rows;i++){
        for(unsigned j=0;j<m.cols;j++){
            Entries[i][j] = m.Entries[i][j];
        }
    }
    matrixRank = m.matrixRank; 
     
    return *this;
}

bool Matrix::operator== (const Matrix &m)
{
    if (rows != m.rows || cols != m.cols)
    {
        return false;
    }
    for (unsigned i=0;i<rows;i++)
    {
        for (unsigned j=0;j<cols;j++)
        {
            if (Entries[i][j] != m.Entries[i][j])
            {
                return false;
            }

        }
    }
    return true;
}

double& Matrix::operator() (unsigned row, unsigned col)
{
    if (row >= rows || col >= cols)
    {
        cerr << "Matrix operator ():Matrix subscript out of bounds" << endl;
        cerr << "(" << row << ", " << col << ")" <<endl;
        cerr << this << endl;
        cerr << *this;
        exit(0);
    }
    return Entries[row][col];
}

double Matrix::operator() (unsigned row, unsigned col) const
{
    if (row >= rows || col >= cols)
    {
        cerr << "Matrix operator ():Matrix subscript out of bounds" << endl;
        cerr << this << endl;
        exit(0);
    }
    return Entries[row][col];
} 

std::ostream& operator<< (std::ostream& o, const Matrix& someMatrix)
{
    ios::fmtflags oldFlags;
    int oldStreamSize;

    oldStreamSize = (int)o.precision(4);
    oldFlags = o.setf (ios::fixed);

    for (unsigned i=0;i<someMatrix.rows;i++){
        for (int pl=0;pl<Matrix::printPadLength;pl++)
        {
            o << " ";
        }
        for(unsigned j=0;j<someMatrix.cols;j++){
            o.width(3);
            o << setw(8) << someMatrix.Entries[i][j] << " ";
        }
        o << "\n";
    }

    o.setf (oldFlags);
    o.precision(oldStreamSize);

    return o;
}


std::istream& operator>> (std::istream& in, Matrix &M)
{
    // First delete Entries if allocated. Determinted by rows
    if (M.rows != 0 || M.cols != 0)
    {
        M.deleteEntries();
    }
    // First entry is the number of rows
    in >> M.rows;
    in >> M.cols;
    M.initialize(M.rows, M.cols);
    for (unsigned j=0;j<M.rows;j++)
    {
        for (unsigned i=0;i<M.cols;i++)
        {
            in >> M.Entries[j][i];
        }
    }

    return in;
}

//double Matrix::det ()
//{
//    if (rows != cols)
//    {
//        cerr << "det called with rows != cols\n";
//        exit(0);
//    }
//}


//Matrix *Matrix::GE ()
//{
//    Matrix  *newM = new Matrix;
//
//    *newM = *this;
//
//    unsigned i=0,j=0;
//    unsigned maxi;
//    int rankCount = 0;
//    double Entry_i_j;
//    double Entry_u_j;
//    ios::fmtflags oldFlags;
//    int oldStreamSize;
//        
//    Matrix::printPadLength=4;
//    oldStreamSize = cout.precision(2);
//    //cout.width(4);
//    //cout.fill ('x');
//    oldFlags = cout.setf (ios::fixed);
//
//
//    while (i < rows && j < cols)
//    {
//        //cout << "i=" << i << " j=" << j << endl;
//        maxi = i;
//        for (unsigned k=i+1;k < rows;k++)
//        {
//            if (fabs(newM->Entries[k][j]) > fabs(newM->Entries[maxi][j]))
//            {
//                maxi=k;
//            }
//        }
//        //cout << "maxi=" << maxi << endl;
//        if (newM->Entries[maxi][j] != 0)
//        {
//            ////cout << "Entries[maxi][j] != 0\n";
//            //cout << *newM;
//            if(maxi != i) {
//                //cout << "Swap" << endl;
//                newM->swapRows(maxi,i);
//                //cout << *newM;
//            }
//            Entry_i_j = newM->Entries[i][j];
//            for (unsigned k=0;k<cols;k++)
//            {
//                newM->Entries[i][k] = newM->Entries[i][k] / Entry_i_j;
//            }
//            //cout << "Divide row i=" << i << " Entries[i][j]=" << Entry_i_j << endl;
//            //cout << *newM;
//            for (unsigned u=i+1;u<rows;u++)
//            {
//                Entry_u_j = newM->Entries[u][j];
//                if (Entry_u_j != 0)
//                {
//                    for (unsigned l=0;l<cols;l++)
//                    {
//                        newM->Entries[u][l] = newM->Entries[u][l] - Entry_u_j*newM->Entries[i][l];
//                    }
//                    //cout << "Processing row " << u << ". Dividing by Entries[u][j]=" << Entry_u_j;// << endl;
//                    //cout << ". u=" << u << " j=" << j << endl;
//                    //cout << *newM;
//                }
//            }
//            i++;
//            rankCount++;
//        }
//        //else {
//        //    //cout << "Entries[maxi][j] == 0\n";
//        //}
//        j++;
//    }
//    Matrix::printPadLength=0;
//    cout.setf (oldFlags);
//    cout.precision(oldStreamSize);
//    //cout << "Rank " << rankCount << endl;
//    newM->matrixRank=rankCount;
//    matrixRank=rankCount;
//    return newM;
//}

Matrix Matrix::GE ()
{
    Matrix  newM;

    newM = *this;

    unsigned i=0,j=0;
    unsigned maxi;
    int rankCount = 0;
    double Entry_i_j;
    double Entry_u_j;
    //ios::fmtflags oldFlags;
    //int oldStreamSize;
        
    //Matrix::printPadLength=4;
    //oldStreamSize = cout.precision(2);
    //cout.width(4);
    ////cout.fill ('x');
    //oldFlags = cout.setf (ios::fixed);


    while (i < rows && j < cols)
    {
        //cout << "i=" << i << " j=" << j << endl;
        maxi = i;
        for (unsigned k=i+1;k < rows;k++)
        {
            if (fabs(newM.Entries[k][j]) > fabs(newM.Entries[maxi][j]))
            {
                maxi=k;
            }
        }
        //cout << "maxi=" << maxi << endl;
        if (newM.Entries[maxi][j] != 0)
        //if (fabs(newM.Entries[maxi][j]) > 0.0000000001) //Handle -0.00 hanging around
        //if (fabs(newM.Entries[maxi][j]) != 0)
        {
            //cout << "Entries[maxi][j] != 0\n";
            //cout << newM;
            if(maxi != i) {
                //cout << "Swap" << endl;
                newM.swapRows(maxi,i);
                //cout << newM;
            }
            Entry_i_j = newM.Entries[i][j];
            for (unsigned k=0;k<cols;k++)
            {
                newM.Entries[i][k] = newM.Entries[i][k] / Entry_i_j;
            }
            //cout << "Divide row i=" << i << " Entries[i][j]=" << Entry_i_j << " fabs=" << fabs(Entry_i_j) << endl;
            //cout << newM;
            for (unsigned u=i+1;u<rows;u++)
            {
                Entry_u_j = newM.Entries[u][j];
                if (Entry_u_j != 0)
                {
                    for (unsigned l=0;l<cols;l++)
                    {
                        newM.Entries[u][l] = newM.Entries[u][l] - Entry_u_j*newM.Entries[i][l];
                    }
                    //cout << "Processing row " << u << ". Dividing by Entries[u][j]=" << Entry_u_j;// << endl;
                    //cout << ". u=" << u << " j=" << j << endl;
                    //cout << newM;
                }
            }
            i++;
            rankCount++;
        }
        //else {
        //    //cout << "Entries[maxi][j] == 0\n";
        //}
        j++;
    }
    //Matrix::printPadLength=0;
    //cout.setf (oldFlags);
    //cout.precision(oldStreamSize);
    //cout << "Rank " << rankCount << endl;
    newM.matrixRank=rankCount;
    matrixRank=rankCount;
    return newM;
}

#ifdef HAVE_LIBGMPXX
int Matrix::GE_gmp_verbose ()
{
    unsigned i=0,j=0;
    unsigned maxi;
    int rankCount = 0;
    mpq_class Entry_i_j;
    mpq_class Entry_u_j;
    mpq_class tempmpq_class;

    // Allocate double array of mpq_class to hold our entries
    mpq_class   **tempEntries = new mpq_class *[rows];
    for(unsigned i=0;i<rows;i++)
    {
        tempEntries[i] = new mpq_class[cols];
        for(unsigned j=0;j<cols;j++)
        {
            tempEntries[i][j] = Entries[i][j];
        }
    }


    //ios::fmtflags oldFlags;
    //int oldStreamSize;
        
    //Matrix::printPadLength=4;
    //oldStreamSize = cout.precision(2);
    //cout.width(4);
    ////cout.fill ('x');
    //oldFlags = cout.setf (ios::fixed);


    while (i < rows && j < cols)
    {
        cout << "i=" << i << " j=" << j << endl;
        maxi = i;
        for (unsigned k=i+1;k < rows;k++)
        {
            if (abs(tempEntries[k][j]) > abs(tempEntries[maxi][j]))
            {
                maxi=k;
            }
        }
        cout << "maxi=" << maxi << endl;
        if (tempEntries[maxi][j] != 0)
        //if (fabs(newM.Entries[maxi][j]) > 0.0000000001) //Handle -0.00 hanging around
        //if (fabs(newM.Entries[maxi][j]) != 0)
        {
            cout << "Entries[maxi][j] = " << tempEntries[maxi][j] << endl;;

            //cout << newM;
            if(maxi != i) {
                // Swap rows maxi and i
                for (unsigned l=0;l<cols;l++)
                {
                    tempmpq_class = tempEntries[maxi][l];
                    tempEntries[maxi][l] = tempEntries[i][l];
                    tempEntries[i][l] = tempmpq_class;
                }
                cout << "Swap row " << maxi << "with " << i << endl;
                //newM.swapRows(maxi,i);
                //cout << newM;
            }
            Entry_i_j = tempEntries[i][j];
            for (unsigned k=0;k<cols;k++)
            {
                tempEntries[i][k] = tempEntries[i][k] / Entry_i_j;
            }
            cout << "Divided row i=" << i << " Entries[i][j]=" << Entry_i_j << " fabs=" << abs(Entry_i_j) << endl;
            cout << "*******************************************************" << endl;
            for(unsigned ib=0;ib<rows;ib++)
            {
                for(unsigned jb=0;jb<cols;jb++)
                {
                    cout << tempEntries[ib][jb];
                    if (jb < cols-1)
                    {
                        cout << " ";
                    }
                }
                cout << endl;
            }
            cout << "*******************************************************" << endl;
            //cout << newM;
            for (unsigned u=i+1;u<rows;u++)
            {
                Entry_u_j = tempEntries[u][j];
                if (Entry_u_j != 0)
                {
                    for (unsigned l=0;l<cols;l++)
                    {
                        tempEntries[u][l] = tempEntries[u][l] - Entry_u_j*tempEntries[i][l];
                    }
                    cout << "Processing row " << u << ". Dividing by Entries[u][j]=" << Entry_u_j;// << endl;
                    cout << ". u=" << u << " j=" << j << endl;
                    //cout << newM;
                }
            }
            i++;
            rankCount++;
            cout << "*******************************************************" << endl;
            for(unsigned ib=0;ib<rows;ib++)
            {
                for(unsigned jb=0;jb<cols;jb++)
                {
                    cout << tempEntries[ib][jb];
                    if (jb < cols-1)
                    {
                        cout << " ";
                    }
                }
                cout << endl;
            }
            cout << "*******************************************************" << endl;
        }
        else {
            cout << "Entries[maxi][j] == 0\n";
        }
        j++;
    }
    cout << "GE_gmp matrix." << endl;
    for(unsigned i=0;i<rows;i++)
    {
        for(unsigned j=0;j<cols;j++)
        {
            cout << tempEntries[i][j];
            if (j < cols-1)
            {
                cout << " ";
            }
        }
        cout << endl;
    }
    //Matrix::printPadLength=0;
    //cout.setf (oldFlags);
    //cout.precision(oldStreamSize);
    cout << "Rank " << rankCount << endl;
    for(unsigned i=0;i<rows;i++)
    {
        delete[] tempEntries[i];
    }
    delete[] tempEntries;

    matrixRank=rankCount;
    return rankCount;
}

int Matrix::GE_gmp ()
{
    unsigned i=0,j=0;
    unsigned maxi;
    int rankCount = 0;
    mpq_class Entry_i_j;
    mpq_class Entry_u_j;
    mpq_class tempmpq_class;

    // Allocate double array of mpq_class to hold our entries
    mpq_class   **tempEntries = new mpq_class *[rows];
    for(unsigned i=0;i<rows;i++)
    {
        tempEntries[i] = new mpq_class[cols];
        for(unsigned j=0;j<cols;j++)
        {
            tempEntries[i][j] = Entries[i][j];
        }
    }


    //ios::fmtflags oldFlags;
    //int oldStreamSize;
        
    //Matrix::printPadLength=4;
    //oldStreamSize = cout.precision(2);
    //cout.width(4);
    ////cout.fill ('x');
    //oldFlags = cout.setf (ios::fixed);


    while (i < rows && j < cols)
    {
        //cout << "i=" << i << " j=" << j << endl;
        maxi = i;
        for (unsigned k=i+1;k < rows;k++)
        {
            if (abs(tempEntries[k][j]) > abs(tempEntries[maxi][j]))
            {
                maxi=k;
            }
        }
        //cout << "maxi=" << maxi << endl;
        if (tempEntries[maxi][j] != 0)
        //if (fabs(newM.Entries[maxi][j]) > 0.0000000001) //Handle -0.00 hanging around
        //if (fabs(newM.Entries[maxi][j]) != 0)
        {
            //cout << "Entries[maxi][j] != 0\n";
            //cout << newM;
            if(maxi != i) {
                // Swap rows maxi and i
                for (unsigned l=0;l<cols;l++)
                {
                    tempmpq_class = tempEntries[maxi][l];
                    tempEntries[maxi][l] = tempEntries[i][l];
                    tempEntries[i][l] = tempmpq_class;
                }
                //cout << "Swap" << endl;
                //newM.swapRows(maxi,i);
                //cout << newM;
            }
            Entry_i_j = tempEntries[i][j];
            for (unsigned k=0;k<cols;k++)
            {
                tempEntries[i][k] = tempEntries[i][k] / Entry_i_j;
            }
            //cout << "Divide row i=" << i << " Entries[i][j]=" << Entry_i_j << " fabs=" << abs(Entry_i_j) << endl;
            //cout << newM;
            for (unsigned u=i+1;u<rows;u++)
            {
                Entry_u_j = tempEntries[u][j];
                if (Entry_u_j != 0)
                {
                    for (unsigned l=0;l<cols;l++)
                    {
                        tempEntries[u][l] = tempEntries[u][l] - Entry_u_j*tempEntries[i][l];
                    }
                    //cout << "Processing row " << u << ". Dividing by Entries[u][j]=" << Entry_u_j;// << endl;
                    //cout << ". u=" << u << " j=" << j << endl;
                    //cout << newM;
                }
            }
            i++;
            rankCount++;
        }
        //else {
        //    //cout << "Entries[maxi][j] == 0\n";
        //}
        j++;
    }
    //cout << "GE_gmp matrix." << endl;
    //for(unsigned i=0;i<rows;i++)
    //{
    //    for(unsigned j=0;j<cols;j++)
    //    {
    //        cout << tempEntries[i][j];
    //        if (j < cols-1)
    //        {
    //            cout << " ";
    //        }
    //    }
    //    cout << endl;
    //}
    //Matrix::printPadLength=0;
    //cout.setf (oldFlags);
    //cout.precision(oldStreamSize);
    //cout << "Rank " << rankCount << endl;
    for(unsigned i=0;i<rows;i++)
    {
        delete[] tempEntries[i];
    }
    delete[] tempEntries;

    matrixRank=rankCount;
    return rankCount;
}
#endif

// Swaps columns i and j
void Matrix::swapCols(unsigned i, unsigned j)
{
    double temp;

    if (i >= rows || j >= cols)
    {
        cerr << "Matrix subscript out of bounds";
        exit(0);
    }

    for(unsigned k=0;k<rows;k++)
    {
        temp = Entries[k][i];
        Entries[k][i] = Entries[k][j];
        Entries[k][j] = temp;

    }
}

// Swaps rows i and j
void Matrix::swapRows(unsigned i, unsigned j)
{
    double temp;

    if (i >= rows || j >= cols)
    {
        cerr << "Matrix subscript out of bounds";
        exit(0);
    }

    for(unsigned k=0;k<cols;k++)
    {
        temp = Entries[i][k];
        Entries[i][k] = Entries[j][k];
        Entries[j][k] = temp;
    }
}

int Matrix::rank_LAPACK()
{
    if (cols == 0 || rows == 0)
    {
        return 0;
    }
    double *AT = new double[rows*cols];
    double *U = new double[rows*cols];
    double *VT = new double[rows*cols];
    double *S = new double[min(rows,cols)];
    double *WORK = new double[2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4)];

    char JOBU, JOBVT;
    int INFO, LDA, LDU, LDVT, LWORK, M, N;

    M = rows;
    N = cols;
    JOBU = 'N';
    JOBVT = 'N';
    LDA = M;
    LDU = 1;
    LDVT = 1;
    LWORK = 2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4);

    // FORTRAN expects column major order.
    
    for(unsigned i=0;i<cols;i++)
    {
        for(unsigned j=0;j<rows;j++)
        {
            AT[j+rows*i] = Entries[j][i];
        }
    }

    //extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double A[], int *lda, double s[], double U[], int *ldu, double VT[], int *ldtv, double Work[], int *lwork, int *info);
    //extern void dgesvd_(char *, char *, long *, long *, double *, long *, double *, double *, long *, double *, long *, double *, long *, long *);

    dgesvd_(&JOBU, &JOBVT, &M, &N, AT, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);

    if (INFO != 0)
    {
        cerr << "Matrix::rank_LAPACK(): INFO = " << INFO << endl;
    }


    double SIGMA_MAX = -1;

    for (unsigned i=0;i<min(rows,cols);i++)
    {
        if (S[i] > SIGMA_MAX)
        {
            SIGMA_MAX = S[i];
        }
    }

    // Now count the elements of S that are > than 10^(-16) * SIGMA_MAX
    int non_zero_count = 0;

    for (unsigned i=0;i<min(rows,cols);i++)
    {
        if (S[i] > SIGMA_MAX*LAPACK_RANK_ZERO )
        {
            non_zero_count++;
        }
    }


    delete AT ;
    delete U;
    delete VT;
    delete S;
    delete WORK;
    return non_zero_count;
}

int Matrix::rank ()
{
    // If rank is not calculated, use GE
    if (matrixRank == -1)
    {
        Matrix tempM;
        //tempM = GE ();
        //cout << "Rank calculated. GE." << endl;
        //cout << tempM;
        //int tempMatrixRank = tempM.rank ();
        int LAPACK_rank = rank_LAPACK();
        //int gmp_rank = GE_gmp();
        //if (gmp_rank != LAPACK_rank)
        //{
        //      cerr << *this;
        //      cerr << "gmp rank: " << gmp_rank << endl;
        //      cerr << "LAPACK rank: " << LAPACK_rank << endl;
        //      cerr << "Setting to GMP rank " << endl;
        //      LAPACK_rank = gmp_rank;
        //      exit (0);
        //}
        //if (gmp_rank != tempMatrixRank)
        //{
        //    cerr << "rank not same between floating and gmp." << endl;
        //    cerr << "Matrix:" << *this << endl;
        //    cerr << "floating GE" << endl;
        //    cerr << tempM;
        //    cerr << "tempM rank: " << tempMatrixRank << endl;
        //    cerr << "gmp rank: " << gmp_rank << endl;
        //    cerr << "LAPACK rank: " << LAPACK_rank << endl;
        //    exit (0);
        //}
        //matrixRank = tempMatrixRank;
        matrixRank = LAPACK_rank;
        //matrixRank = GE_gmp ();
    }

    return matrixRank; 
}

Matrix Matrix::subColumns(const set<unsigned> &colSet)
{
    if (colSet.size() == 0)
    {
        cerr << "Matrix::subColumns called with colSet.size()==0" << endl;
    }
    int colCount=0;

    Matrix newM(rows,(unsigned)colSet.size());

    //cout << "newM.cols=" << newM.cols << endl; 
    //newM.matrixRank = -1;

    for(set <unsigned>::const_iterator p = colSet.begin(); p != colSet.end(); p++)
    {
        if ((*p < 0) || (*p > cols)) {
            cerr << "subColumns: Index given in subset does not fit dimensions.\n";
            exit(0);
        }
        //cout << "   Copying column " << *p << endl;
        for (unsigned i=0;i<newM.rows;i++)
        {
            //cout << "       (i,*p)=(" << i << "," << *p << ")\n";
            newM.Entries[i][colCount] = Entries[i][*p];
            //cout << "       Entries[i][*p]=" << Entries[i][*p] << " " << endl;
            //cout << "       newM.Entries[i][colCount]=" << newM.Entries[i][colCount] << " " << endl;
        }
        colCount++;
    }
    //cout << newM;

    return newM;
}

Matrix Matrix::subColumnsDiff(const set<unsigned> &colSet)
{
    Matrix newM(rows,(unsigned)(cols - colSet.size()));

    //cout << "newM.cols=" << newM.cols << endl; 
    newM.matrixRank = -1;

    int newMCount = 0;
    for (unsigned i=0;i<cols;i++)
    {
        // If i is NOT in colSet then copy it
        if (colSet.find(i) == colSet.end())
        {
            for (unsigned j=0;j<rows;j++)
            {
                newM(j,newMCount) = Entries[j][i];
            }
            newMCount++;
        }
    }

    return newM;
}

Matrix Matrix::operator* (const double &rhs) const
{
    Matrix C = *this;
    for (unsigned i=0;i<C.rows;i++)
    {
        for (unsigned j=0;j<C.cols;j++)
        {
                C.Entries[i][j] *= rhs;
        }
    }

    return C;
}

Matrix  Matrix::operator* (const Matrix &B) const
{
    //Dimensions must match, else give error
    if (this->cols != B.rows) 
    {
        cerr << "*operator: Dimensions of matrices do not match for multiplication" << endl;
        exit (0);
    }
    Matrix C(this->rows,B.cols);
    for (unsigned i=0;i<this->rows;i++)
    {
        for (unsigned j=0;j<B.cols;j++)
        {
            C.Entries[i][j] = 0;
            for (unsigned k=0;k<this->cols;k++)
            {
                C.Entries[i][j] += this->Entries[i][k]*B.Entries[k][j];
            }
        }
    }

    return C;
}

Matrix  Matrix::operator+ (const Matrix &B) const
{
    if (this->cols != B.cols || this->rows != B.rows) 
    {
        cerr << "+operator: Dimensions of matrices do not match" << endl;
        exit (0);
    }
    Matrix C(this->rows,this->cols);
    for (unsigned i=0;i<this->rows;i++)
    {
        for (unsigned j=0;j<this->cols;j++)
        {
            C.Entries[i][j] = this->Entries[i][j] + B.Entries[i][j];
        }
    }

    return C;
}

Matrix  Matrix::operator- (const Matrix &B) const
{
    if (this->cols != B.cols || this->rows != B.rows) 
    {
        cerr << "+operator: Dimensions of matrices do not match" << endl;
        exit (0);
    }
    Matrix C(this->rows,this->cols);
    for (unsigned i=0;i<this->rows;i++)
    {
        for (unsigned j=0;j<this->cols;j++)
        {
            C.Entries[i][j] = this->Entries[i][j] - B.Entries[i][j];
        }
    }

    //cout << *this << " - " << endl << B << " = " << C;

    return C;

}

const int Matrix::getCols ()
{
    return cols;
}

const int Matrix::getRows ()
{
    return rows;
}


//given a mxn matrix, returns a mx1 matrix which is the sum of each row
Matrix Matrix::rowSum()
{
    Matrix  newM(rows,1);

    Matrix  tempM(cols,1);

    for (unsigned i=0;i<cols;i++)
    {
        tempM(i,0) = 1;
    }

    newM = (*this)*tempM; 
    return newM;
}

// If this matrix is a column vector, return the 2-norm squared. Else exit with error
double Matrix::twoNormSquared()
{
    if (cols != 1)
    {
        cerr << "Matirx::twoNorm called on Matrix that is not a column vector." << endl;
        exit (0);
    }
    double tempD = 0;
    for (unsigned i=0;i<rows;i++)
    {
        tempD += (Entries[i][0])*(Entries[i][0]);
    }
    
    //cout << *this;
    //cout << "twoNormSquared = " << tempD << endl;
    return tempD;
}

// If this matrix is a column vector, return the 2-norm. Else exit with error
double Matrix::twoNorm()
{
    return sqrt(twoNormSquared());
}

Matrix Matrix::transpose()
{
    Matrix  tempM(cols,rows);

    for(unsigned i=0;i<rows;i++)
    {
        for(unsigned j=0;j<cols;j++)
        {
            tempM(j,i) = Entries[i][j];
        }
    }

    return tempM;
}

void Matrix::writeMatlab(std::ostream &o, string label)
{
    o << label << " = [" << endl;
    for(unsigned i=0;i<rows;i++)
    {
        for(unsigned j=0;j<cols;j++)
        {
            o << Entries[i][j] << " ";
        }
        o << ";" << endl;
    }
    o << "];" << endl;
}

// This uses dgsvd function of lapack to compute *this = U*S*V.transpose()
void Matrix::SVD(Matrix &reU, Matrix &reS, Matrix &reV)
{
    if (cols == 0 || rows == 0)
    {
        return;
    }
    double *AT = new double[rows*cols];
    reU = Matrix(rows,rows);
    double *U = new double[rows*rows];
    double *VT = new double[cols*cols];
    reV = Matrix(cols,cols);
    double *S = new double[min(rows,cols)];
    reS = Matrix(rows,cols);
    double *WORK = new double[2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4)];
    for(unsigned i=0;i<cols;i++)
    {
        for(unsigned j=0;j<cols;j++)
        {
            U[i + j*cols] = 0;
        }
    }

    char JOBU, JOBVT;
    int INFO, LDA, LDU, LDVT, LWORK, M, N;

    M = rows;
    N = cols;
    JOBU = 'A';
    JOBVT = 'A';
    LDA = M;
    LDU = M;
    LDVT = N;
    LWORK = 2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4);

    // FORTRAN expects column major order.
    for(unsigned i=0;i<cols;i++)
    {
        for(unsigned j=0;j<rows;j++)
        {
            AT[j+rows*i] = Entries[j][i];
        }
    }

    dgesvd_(&JOBU, &JOBVT, &M, &N, AT, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);

    if (INFO != 0)
    {
        cerr << "Matrix::SVD(): INFO = " << INFO << endl;
    }

    for(unsigned i=0;i<rows;i++)
    {
        for(unsigned j=0;j<rows;j++)
        {
            reU(i,j) = U[i + rows*j];
        }
    }
    for(unsigned i=0;i<min(rows,cols);i++)
    {
        reS(i,i) = S[i];
    }
    for(unsigned i=0;i<cols;i++)
    {
        for(unsigned j=0;j<cols;j++)
        {
            reV(i,j) = VT[j + cols*i];
        }
    }

    delete AT;
    delete U;
    delete VT;
    delete S;
    delete WORK;
}

double Matrix::LU(Matrix &P, Matrix &L, Matrix &U)
{
    int M = rows;
    int N = cols;
    double *A = new double[rows*cols];
    int LDA = M;
    int *IPIV = new int[min(M,N)];
    int INFO;
    double sign=1; //Sign determined by number of row transposes

    // FORTRAN expects column major order.
    for(unsigned i=0;i<cols;i++)
    {
        for(unsigned j=0;j<rows;j++)
        {
            A[j+rows*i] = Entries[j][i];
        }
    }
    dgetrf_(&M,&N,A,&LDA,IPIV,&INFO);

    P = Matrix(rows,rows);
    for (unsigned i=0;i<rows;i++){
        for (unsigned j=0;j<rows;j++){
            P(i,j)=0;
        }
        P(i,i) = 1;
    }
    for (int i=min(rows,cols)-1;i>=0;i--){
        //cout << "IPIV[" << i << "] = " << IPIV[i]-1 << endl;
        P.swapRows(i,IPIV[i]-1);
        if (i != IPIV[i]-1){
            sign *= -1;
        }
    }
    if (M > N){
        L = Matrix(M,N);
        U = Matrix(N,N);
    }
    else if (M < N){
        L = Matrix(M,M);
        U = Matrix(M,N);
    }
    else if (M == N){
        L = Matrix(M,M);
        U = Matrix(M,M);
    }

    for (unsigned i=0;i<(unsigned)min(M,N);i++){
        for (unsigned j=0;j<(unsigned)min(M,N);j++){
            if (i > j){
                L(i,j) = A[i+M*j];
            } 
            else if (i < j){
                U(i,j) = A[i+M*j];
            }
            else if (i == j){
                U(i,j) = A[i+M*j];
                L(i,j) = 1.0;
            }
        }
    }

    delete A;
    delete IPIV;
    return sign;
}

// This uses the SVD decomposition function above to compute the inverse.
Matrix Matrix::inverse()
{
    Matrix U,S,V;
    this->SVD(U,S,V);

    for (unsigned i=0;i<min(S.rows,S.cols);i++){
        S(i,i) = 1 / S(i,i);
    }

    return V*(S.transpose())*(U.transpose());
}

double Matrix::trace()
{
    double returnValue=0;
    for (unsigned i=0;i<min(rows,cols);i++){
        returnValue += Entries[i][i];
    }
    return returnValue;
}

long double Matrix::det()
{
    if (rows != cols){
        cout << "Matrix::detCofactor(): rows != cols" << endl;
        exit (0);
    }
    return detLU();
    //return detCofactor();
}

long double Matrix::detCofactor()
{
    if (rows != cols){
        cout << "Matrix::detCofactor(): rows != cols" << endl;
        exit (0);
    }
    if (rows == 1 && cols == 1){
        return Entries[0][0];
    }
    long double returnValue=0;;
    // Do cofactor expansion on first row
    long double sign = 1;
    for (unsigned i=0;i<cols;i++){
        //cout << "detCofactor i=" << i << endl;
        returnValue += sign*((this->cofactor(0,i)).det());
        sign *= -1;
    }
    return returnValue;
}

Matrix Matrix::cofactor(unsigned i, unsigned j)
{
    if (rows < 2 || cols < 2){
        cout << "Matrix::cofactor: rows < 2 || cols < 2" << endl;
        exit (0);
    }

    Matrix returnMatrix(rows-1,cols-1);

    for (unsigned k=0;k<rows;k++){
        for (unsigned l=0;l<cols;l++){
            if (k != i && l != j){
                unsigned tmp_k=k;
                unsigned tmp_l=l;
                if (k>i){
                    tmp_k--;
                }
                if (l>j){
                    tmp_l--;
                }
                returnMatrix(tmp_k,tmp_l)=Entries[k][l];
            }
        }
    }

    return returnMatrix;
}

long double Matrix::detLU()
{
    long double returnValue=1;

    if (rows != cols) {
        cout << "Matrix::detCofactor(): rows != cols" << endl;
        exit (0);
    }
    if (rows == 1 && cols == 1){
        return Entries[0][0];
    }
    Matrix P,L,U;
    double sign;
    sign = LU(P,L,U);
    for (unsigned i=0;i<rows;i++){
        returnValue *= U(i,i);
    }
    return returnValue*sign;
}

long double Matrix::log_abs_det()
{
    long double returnValue=0;

    if (rows != cols) {
        cout << "Matrix::detCofactor(): rows != cols" << endl;
        exit (0);
    }
    if (rows == 1 && cols == 1){
        return Entries[0][0];
    }
    Matrix P,L,U;
    double sign;
    sign = LU(P,L,U);
    for (unsigned i=0;i<rows;i++){
        returnValue += (double)log((double)fabs((double)U(i,i)));
    }
    return returnValue;
}
