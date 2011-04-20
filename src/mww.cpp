// $Rev: 452 $ $Date: 2009-12-01 16:34:16 -0500 (Tue, 01 Dec 2009) $

/** \file mww.cpp */
#include "mww.h"

double normalCumDist(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

int sumUpTo(int n)
{
    if ((n % 2) == 0){ // EVEN
        return (n/2)*(n+1);
    }
    else { // ODD
        return ((n+1)/2)*n;
    }
}

double MannWhitneyWilcoxon(const list <double> sampleOne, const list <double> sampleTwo)
{
    int N = sampleOne.size() + sampleTwo.size();
    int n1 = sampleOne.size();
    int n2 = sampleTwo.size();

    multiset <doubleIndex, ltDoubleIndex> rankedIndividuals;
    if (DEBUG_OUTPUT >= 1){
        cout << "Sample one (" << sampleOne.size() << ") values: ";
        for (list <double>::const_iterator lit=sampleOne.begin();lit!=sampleOne.end();lit++){
            cout << *lit << setw(5) << " ";
        }
        cout << endl;
        cout << "Sample two (" << sampleTwo.size() << ") values: ";
        for (list <double>::const_iterator lit=sampleTwo.begin();lit!=sampleTwo.end();lit++){
            cout << *lit << setw(5) << " ";
        }
        cout << endl;
    }

    int oldRankedIndividualsSize = rankedIndividuals.size();
    for (list <double>::const_iterator lit=sampleOne.begin();lit!=sampleOne.end();lit++){
        doubleIndex newDoubleIndex;
        newDoubleIndex.value = *lit;
        newDoubleIndex.index = 1;
        rankedIndividuals.insert(newDoubleIndex);
        if (DEBUG_OUTPUT >= 1){
            cout << *lit << setw(5) << endl;
            cout << "rankedIndividuals.size() = " << rankedIndividuals.size() << endl;
        }
    }

    oldRankedIndividualsSize = rankedIndividuals.size();
    for (list <double>::const_iterator lit=sampleTwo.begin();lit!=sampleTwo.end();lit++){
        doubleIndex newDoubleIndex;
        newDoubleIndex.value = *lit;
        newDoubleIndex.index = 2;
        rankedIndividuals.insert(newDoubleIndex);
        if (DEBUG_OUTPUT >= 1){
            cout << *lit << setw(5) << endl;
            cout << "rankedIndividuals.size() = " << rankedIndividuals.size() << endl;
        }
    }
    // Asign the ranks of this multi-set. Use the rule that ties are given the average of the ranks they should be
    // under a random tie-breaking linear ordering.
    int numIndividualsSeen = 1;
    list <double> individualRanks;
    for (multiset <doubleIndex, ltDoubleIndex>::iterator sdit=rankedIndividuals.begin();sdit!=rankedIndividuals.end();sdit++){
        int thisMultiSetCount = rankedIndividuals.count((*sdit));
        if (DEBUG_OUTPUT >= 1){
            cout << "numIndividualsSeen = " << numIndividualsSeen << endl;
            cout << "(*sdit).value = " << (*sdit).value << "   " << "(*sdit).index = " << (*sdit).index << endl;
            cout << "rankedIndividuals.count((*sdit)) = " << thisMultiSetCount << endl;
        }
        if (thisMultiSetCount > 1){
            double avgRank = ((double)sumUpTo(numIndividualsSeen + thisMultiSetCount-1) - (double)sumUpTo(numIndividualsSeen-1))/(double)thisMultiSetCount;
            if (DEBUG_OUTPUT >= 1){
                cout << "   rankedIndividuals.count((*sdit)) = " << thisMultiSetCount << endl;
                cout << "   numIndividualsSeen = " << numIndividualsSeen << endl;
                cout << "   avgRank = " << avgRank << endl;
            }
            for (int i=0;i<thisMultiSetCount;i++){
                individualRanks.push_back(avgRank);
                sdit++;
                numIndividualsSeen++;
            }
            sdit--;
        }
        else {
            individualRanks.push_back(numIndividualsSeen);
            numIndividualsSeen++;
        }
    }

    if (DEBUG_OUTPUT >= 1){
        cout << "Ranked individuals: ";
        list <double>::const_iterator ldit=individualRanks.begin();
        int tmpCount = 1;
        for (multiset <doubleIndex, ltDoubleIndex>::const_iterator sdit=rankedIndividuals.begin();sdit!=rankedIndividuals.end();sdit++){
            cout << tmpCount << "[" << *ldit << "," << (*sdit).index << "," << (*sdit).value << "] ";
            ldit++;
            tmpCount++;
        }
        cout << endl;
    }

    double rankSumSampleOne = 0;
    list <double>::const_iterator ldit=individualRanks.begin();
    for (multiset <doubleIndex, ltDoubleIndex>::const_iterator sdit=rankedIndividuals.begin();sdit!=rankedIndividuals.end();sdit++){
        if ((*sdit).index == 1) {
            rankSumSampleOne += *ldit;
        }
        ldit++;
    }
    double rankSumSampleTwo = 0;
    ldit=individualRanks.begin();
    for (multiset <doubleIndex, ltDoubleIndex>::const_iterator sdit=rankedIndividuals.begin();sdit!=rankedIndividuals.end();sdit++){
        if ((*sdit).index == 2) {
            rankSumSampleTwo += *ldit;
        }
        ldit++;
    }

    int rightTerm;

    if ((sampleOne.size() % 2) == 0){
        rightTerm = (sampleOne.size()/2)*(sampleOne.size() + 1);
    }
    else {
        rightTerm = (sampleOne.size())*((sampleOne.size() + 1)/2);
    }
        
    double U = rankSumSampleOne - rightTerm;

    // Now get the z-stat.

    double meanU = (double)(sampleOne.size())*(double)(sampleTwo.size())/2.0;

    // From wikipedia and other references.
    double varU = sqrt(((double)(sampleOne.size())*(double)(sampleTwo.size())*((double)(sampleOne.size())+(double)(sampleTwo.size()) + 1))/12.0);

    double z = (U - meanU)/varU;

    // Calculation using handout Arne gave me.
    double ArneVarU=0;
    for (ldit=individualRanks.begin();ldit!=individualRanks.end();ldit++){
        ArneVarU += ((*ldit) - ((double)N + 1.0)/2.0)*((*ldit) - ((double)N + 1.0)/2.0);
    }
    ArneVarU = ArneVarU/(((double)N)*((double)N)*((double)(N-1)));
    ArneVarU = sqrt(ArneVarU);

    double pN = (1/(double)n1)*(rankSumSampleTwo/(double)n2 - ((double)n2 + 1)/2);
    double ArneZ = sqrt((double)n1*(double)n2/((double)N*ArneVarU))*(pN - 0.5);

    if (DEBUG_OUTPUT >= 0){
        cout << "n1                 = " << n1 << endl;
        cout << "n2                 = " << n2 << endl;
        cout << "N                  = " << N << endl;
        cout << "R_1                = " << rankSumSampleOne << endl;
        cout << "R_2                = " << rankSumSampleTwo << endl;
        cout << "rightTerm          = " << rightTerm << endl;
        cout << "U                  = " << U << endl;
        cout << "meanU              = " << meanU << endl;
        cout << "varU               = " << varU << endl;
        cout << "z                  = " << z << endl;
        cout << "normalCumDist(-fabs(z))   = " << normalCumDist(-1.0*fabs(z)) << endl;
        cout << "---------------------------" << endl;
        cout << "ArneVarU           = " << ArneVarU << endl;
        cout << "pN                 = " << pN << endl;
        cout << "ArneZ              = " << ArneZ << endl;
        cout << "normalCumDist(-fabs(ArneZ))   = " << normalCumDist(-1.0*fabs(ArneZ)) << endl;
        cout << "normalCumDist(ArneZ)   = " << normalCumDist(ArneZ) << endl;
    }
        
    //return normalCumDist(-1.0*fabs(z));
    //return normalCumDist(-1.0*fabs(ArneZ));
    return normalCumDist(ArneZ);
}
