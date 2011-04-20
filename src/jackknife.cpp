// $Rev: 287 $ $Date: 2009-11-01 17:41:14 -0500 (Sun, 01 Nov 2009) $
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <map>
#include <list>
#include <set>
#include <math.h>
#include "debugoutput.h"
#include "jackknife.h"

using namespace::std;

char **getAlignments(string nexInputFileName,int &ntax, int &nchar, list <string> &taxaNames)
{
    fstream nexInput;
     // Read in until matrix seen
 
    ntax = 0;
    nchar = 0;
    char **alignments; // This is a double index matrix where alignments[i] is the 
                       // ith column of the alignments. Hence, this is transposed
    nexInput.open(nexInputFileName.c_str(),fstream::in);
 
    if (!nexInput.good()){
        cout << "Can not open " << nexInputFileName << endl;
    }
 
    string nexString;
    while (nexInput.good()){
        getline(nexInput,nexString);

        if (nexString.find("dimensions") != string::npos && nexInput.good()){
            // This current line contains "ntax"
            size_t ntax_start = nexString.find("ntax");
            unsigned ntax_digitPos = ntax_start;
            while (!isdigit(nexString[ntax_digitPos])){
                ntax_digitPos++;
            }
            string tempStr;
            tempStr.clear();
            while (isdigit(nexString[ntax_digitPos])){
                tempStr += nexString[ntax_digitPos]; 
                ntax_digitPos++;
            }
            if (DEBUG_OUTPUT >= 1){
                cout << "NUMBER OF TAXA: " << tempStr << endl;
            }
            ntax = atoi(tempStr.c_str());

            size_t nchar_start = nexString.find("nchar");
            unsigned nchar_digitPos = nchar_start;
            while (!isdigit(nexString[nchar_digitPos])){
                nchar_digitPos++;
            }
            tempStr.clear();
            while (isdigit(nexString[nchar_digitPos])){
                tempStr += nexString[nchar_digitPos]; 
                nchar_digitPos++;
            }
            if (DEBUG_OUTPUT >= 1){
                cout << "NUMBER OF CHARACTERS: " << tempStr << endl;
            }
            nchar = atoi(tempStr.c_str());
        }
        if (nexString.find("matrix") != string::npos && nexInput.good() && ntax > 0 && nchar > 0){

            alignments = new char*[nchar];
            for (int i=0;i<nchar;i++){
                alignments[i] = new char[ntax];
            }
            unsigned curAlignment = 0;
            // Now each line until ';' seen is a taxa
            getline(nexInput,nexString);
            while (nexInput.good() && nexString.find(";") == string::npos){
                if (DEBUG_OUTPUT >= 1){
                    cout << "STRING: " << nexString << endl;
                }
                size_t firstNonSpace = nexString.find_first_not_of(" ");
                size_t endOfName = nexString.find_first_of(" ",firstNonSpace);

                string taxaName=nexString.substr(firstNonSpace, endOfName-firstNonSpace);

                if (DEBUG_OUTPUT >= 1){
                    cout << "NAME: " << taxaName << endl;
                }
                taxaNames.push_back(taxaName);

                // Now everything after, not a space, is the alignment
                string alignment;
                alignment.clear();
                unsigned nexStringPos = endOfName;
                unsigned curAlignmentPos = 0;
                while (nexString[nexStringPos] != 0){
                    if (nexString[nexStringPos] != ' '){
                        alignment += nexString[nexStringPos];
                        alignments[curAlignmentPos][curAlignment] = nexString[nexStringPos]; 
                        curAlignmentPos++;
                    }
                    nexStringPos++;
                }
                if (DEBUG_OUTPUT >= 1){
                    cout << "ALIGNMENT: " << alignment << endl;
                }
                curAlignment++;
                getline(nexInput,nexString);
            }
        }
    }

    if (DEBUG_OUTPUT >= 1){
        cout << "ntax " << ntax << endl;
        cout << "nchar " << nchar << endl;
        cout << endl << "**** alignments ****" << endl;
        for(int i=0;i<ntax;i++){
            //cout << "i=" << i << "  ";
            for(int j=0;j<nchar;j++){
                cout << alignments[j][i];
            }
            cout << endl;
        }
    }

    return alignments;
}

// This will take in char[m][n] where char[i] is the ith columns (i.e. transposed)
// and delete any columns that do not contain 'A','T','G', or 'C'. The number of rows
// of the returned data will be the same, but it will have newM columns.
char **delColumnsWithMissingData(char **alignments,int m, int n, int &newM)
{
    char **newAlignments;
    newM = 0;
    int delColumn[m];

    for(int i=0;i<m;i++){
        delColumn[i]=0;
    }

    // Go through every column
    for(int i=0;i<m;i++){
        //cout << "Checking column " << i << endl;
        for (int j=0;j<n;j++){
            if (alignments[i][j] != 'A' && alignments[i][j] != 'T' && alignments[i][j] != 'C' && alignments[i][j] != 'G'){
                //cout << "       Bad character alignments[" << i << "][" << j << "] = " << alignments[i][j] << endl;
                delColumn[i]=1; // Delete this column
            }
        }
        if (delColumn[i] == 0){
            //cout << "   Column GOOD" << endl;
            newM++;
        }
        else {
            //cout << "   Column BAD" << endl;
        }
    }

    newAlignments = new char*[newM];
    for(int i=0;i<newM;i++){
        newAlignments[i] = new char[n];
    }
    int tempCount = 0;

    for(int i=0;i<m;i++){
        if (delColumn[i] == 0){
            for (int j=0;j<n;j++){
                newAlignments[tempCount][j] = alignments[i][j];
            }

            tempCount++;
        }
    }
    

    return newAlignments;
}

void writeAlignmentsToNex(ostream &out, char **alignments, int ntax, int nchar)
{
    out << "begin data;" << endl;
    out << "dimensions ntax=" << ntax << " nchar=" << nchar << ";" << endl;
    out << "format datatype=dna missing=? gap=-;" << endl;
    out << "matrix" << endl;
    for (int i=0;i<ntax;i++){
        out << "   " << i+1 << "     ";
        for (int j=0;j<nchar;j++){
            out << alignments[j][i];
        }
        out << endl;
    }
    out << ";" << endl;
    out << "end;" << endl;

}

void writeAlignmentsTo_phylip(ostream &out, char **alignments, int ntax, int nchar)
{
    out << ntax << " " << nchar << endl;

    for (int i=0;i<ntax;i++){
        out << i+1 << setw(10);
        for (int j=0;j<nchar;j++){
            out << alignments[j][i];
        }
        out << endl;
    }
}

void writeAlignmentsTo_phylip(ostream &out, char **alignments, int ntax, int nchar, list <string> taxaNames)
{
    out << ntax << " " << nchar << endl;

    list <string>::iterator lsit = taxaNames.begin();
    for (int i=0;i<ntax;i++){
        out << *lsit << setw(10);
        lsit++;
        for (int j=0;j<nchar;j++){
            out << alignments[j][i];
        }
        out << endl;
    }
}

void writeAlignmentsToNex(ostream &out, char **alignments, int ntax, int nchar, list <string> taxaNames)
{
    out << "begin data;" << endl;
    out << "dimensions ntax=" << ntax << " nchar=" << nchar << ";" << endl;
    out << "format datatype=dna missing=? gap=-;" << endl;
    out << "matrix" << endl;
    list <string>::iterator lsit = taxaNames.begin();
    for (int i=0;i<ntax;i++){
        out << "   " << *lsit << "     ";
        lsit++;
        for (int j=0;j<nchar;j++){
            out << alignments[j][i];
        }
        out << endl;
    }
    out << ";" << endl;
    out << "end;" << endl;

}

// This takes in char[nchar][ntax] table of alignments and will return.
// a randomly chooses a k-contiguous block of columns and returns it.
char **getRandAlignmentsContCol(char **alignments, int ntax, int nchar, int k)
{
    // Choose a value between 1 and nchar - k
    int startPos = (rand () % (nchar - k + 1));
    char **newAlignments;

    newAlignments = new char*[k];
    for(int i=0;i<k;i++){
        newAlignments[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            newAlignments[i][j] = alignments[startPos + i][j];
        }
    }
    return newAlignments;

}

// This glues the two alignments alignmentsOne[ncharOne][ntaxOne] and alignments[ncharTwo][ntaxTwo]
char **concatAlignments(char **alignmentsOne, char **alignmentsTwo, int ntaxOne, int ncharOne, int ntaxTwo, int ncharTwo)
{
    if (ntaxOne != ntaxTwo){
        cout << "ntaxOne != ntaxTwo" << endl;
        exit(0);
    }

    char **newAlignments;
    newAlignments = new char*[ncharOne + ncharTwo];
    for(int i=0;i<ncharOne + ncharTwo;i++){
        newAlignments[i] = new char[ntaxOne];
    }
    for(int i=0;i<ncharOne;i++){
        for (int j=0;j<ntaxOne;j++){
            newAlignments[i][j] = alignmentsOne[i][j];
        }
    }
    for(int i=0;i<ncharTwo;i++){
        for (int j=0;j<ntaxTwo;j++){
            newAlignments[ncharOne + i][j] = alignmentsTwo[i][j];
        }
    }
   

    return newAlignments;
}

// This takes in the filename for a .nex file. It will strip out any columns
// that do not contain ATCG. If there are N columns left this will perform
// numJackknifes, random jackknife calculations to create a new alignment of the approprate size 
// Usually numJackknifes = N/k
char **getJackknife(string nexFileName,int k,int numJackknifes, int &return_ntax, int &return_nchar)
{
    
    int ntax, nchar;
    list <string> taxaNames;
    char **alignments = getAlignments(nexFileName,ntax,nchar,taxaNames);

    int pruned_nchar;
    char **prunedAlignments = delColumnsWithMissingData(alignments,nchar,ntax,pruned_nchar);
    
    // Delete alignments now.
    for (int i=0;i<nchar;i++){
        delete [] alignments[i];
    }
    delete alignments;


    if (numJackknifes != -1 && numJackknifes > (int)floor((double)pruned_nchar/k)*k){
        cout << "getJackknife:: numJackknives > (int)floor(pruned_nchar/k)*k" << endl;
        exit (0);
    }
    int new_nchar;
    if (numJackknifes == -1){
        new_nchar = (int)floor((double)pruned_nchar/k)*k;
    }
    else {
        new_nchar = numJackknifes*k;
    }
    // reallocate to the correct size
    alignments = new char*[new_nchar];
    for(int i=0;i<new_nchar;i++){
        alignments[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            alignments[i][j] = 0;
        }
    }

    if (DEBUG_OUTPUT >= 1){
        cout << "k=" << k << "   pruned_nchar=" << pruned_nchar << "    new_nchar=" << new_nchar << endl;
    }

    char **randAlignments;
    // Do pruned_nchar/k jackknives and place in alignments
    for (int i=0;i<new_nchar/k;i++){
        randAlignments = getRandAlignmentsContCol(prunedAlignments,ntax,pruned_nchar,k);
        if (DEBUG_OUTPUT >= 1){
            cout << "Random alignment " << i << endl;
            writeAlignmentsToNex(cout,randAlignments,ntax,k,taxaNames);
        }
        for (int j=0;j<k;j++){
            for (int l=0;l<ntax;l++){
                alignments[i*k + j][l] = randAlignments[j][l];
            }
        }

        for (int j=0;j<k;j++){
            delete [] randAlignments[j];
        }
        delete [] randAlignments;

    }

    return_ntax = ntax;
    return_nchar = new_nchar;
    if (DEBUG_OUTPUT >= 1){
        cout << "All alignments." << endl;
        writeAlignmentsToNex(cout,alignments,ntax,new_nchar,taxaNames);
    }
   
    return alignments;

}

char **getJackknife(char **alignments,int ntax, int nchar,int k,int numJackknifes, int &return_ntax, int &return_nchar)
{
    return (char **)0;
}

// This deletes alignments[nchar][ntax]
void deleteAlignments(char **alignments,int ntax, int nchar)
{
    for (int i=0;i<nchar;i++){
        delete [] alignments[i];
    }
    delete [] alignments;
}

// This glues the multiple alignments 
char **concatAlignments(list <char **> alignmentsList, list <int> ntaxList, list <int> ncharList, int &new_ntax, int &new_nchar)
{
    if ( alignmentsList.size() != ntaxList.size() || alignmentsList.size() != ncharList.size() || ntaxList.size() != ncharList.size()){
        cout << "concatAlignments: list sizes not the same" << endl;
        exit(0);
    }

    char **newAlignments;
    new_ntax = *(ntaxList.begin());
    new_nchar = 0;
    for (list <int>::iterator iit = ncharList.begin(); iit != ncharList.end();iit++)
    {
        new_nchar += *iit;
    }
    newAlignments = new char*[new_nchar];
    for(int i=0;i<new_nchar;i++){
        newAlignments[i] = new char[new_ntax];
    }
    int alignmentsCount = 0; // Track which alignments we are on
    list <char **>::iterator lcit = alignmentsList.begin();
    list <int>::iterator ntaxit = ntaxList.begin();
    list <int>::iterator ncharit = ncharList.begin();
    for(int k=0;k<(int)alignmentsList.size();k++){
        for(int i=0;i<*ncharit;i++){
            for (int j=0;j<*ntaxit;j++){
                newAlignments[alignmentsCount*(*ncharit) + i][j] = (*lcit)[i][j];
            }
        }

        ntaxit++;
        ncharit++;
        lcit++;
        alignmentsCount++;
    }

    return newAlignments;
}



