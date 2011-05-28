
/** \file svmresults.cpp */
#include "svmresults.h"
SVM_separationResults::SVM_separationResults()
{
    groupOneCount = 0;
    groupTwoCount = 0;
    resampleGroupOneCount = 0;
    resampleGroupTwoCount = 0;
    separationPercentage = 0;
    resampleSeparationPercentage = 0;
    numSingularValuesTaken = 0;
    diffMeans_l2_norm = 0;
    traceVarianceGroupOne = 0;
    detVarianceGroupOne = 0;
    logDetVarianceGroupOne = 0;
    traceVarianceGroupTwo = 0;
    detVarianceGroupTwo = 0;
    logDetVarianceGroupTwo = 0;
    diffMeans_l2_normPCA = 0;
    traceVarianceGroupOnePCA = 0;
    detVarianceGroupOnePCA = 0;
    logDetVarianceGroupOnePCA = 0;
    traceVarianceGroupTwoPCA = 0;
    detVarianceGroupTwoPCA = 0;
    logDetVarianceGroupTwoPCA = 0;
    secondsToCompute = 0;
}

SVM_separationResults & SVM_separationResults::operator+=(const SVM_separationResults &rhs) 
{
    groupOneCount                =  groupOneCount                + rhs.groupOneCount               ;     
    groupTwoCount                =  groupTwoCount                + rhs.groupTwoCount               ;     
    resampleGroupOneCount        =  resampleGroupOneCount        + rhs.resampleGroupOneCount       ;             
    resampleGroupTwoCount        =  resampleGroupTwoCount        + rhs.resampleGroupTwoCount       ;             
    separationPercentage         =  separationPercentage         + rhs.separationPercentage        ;            
    resampleSeparationPercentage =  resampleSeparationPercentage + rhs.resampleSeparationPercentage;                    
    numSingularValuesTaken       =  numSingularValuesTaken       + rhs.numSingularValuesTaken      ;              
    //diffMeans                    =  diffMeans                    + rhs.diffMeans                   ; 
    diffMeans_l2_norm            =  diffMeans_l2_norm            + rhs.diffMeans_l2_norm           ;         
    //varianceGroupOne             =  varianceGroupOne             + rhs.varianceGroupOne            ;        
    traceVarianceGroupOne        =  traceVarianceGroupOne        + rhs.traceVarianceGroupOne       ;             
    detVarianceGroupOne          =  detVarianceGroupOne          + rhs.detVarianceGroupOne         ;           
    logDetVarianceGroupOne       =  logDetVarianceGroupOne       + rhs.logDetVarianceGroupOne      ;              
    //varianceGroupTwo             =  varianceGroupTwo             + rhs.varianceGroupTwo            ;        
    traceVarianceGroupTwo        =  traceVarianceGroupTwo        + rhs.traceVarianceGroupTwo       ;             
    detVarianceGroupTwo          =  detVarianceGroupTwo          + rhs.detVarianceGroupTwo         ;           
    logDetVarianceGroupTwo       =  logDetVarianceGroupTwo       + rhs.logDetVarianceGroupTwo      ;              
    //diffMeansPCA                 =  diffMeansPCA                 + rhs.diffMeansPCA                ;    
    diffMeans_l2_normPCA         =  diffMeans_l2_normPCA         + rhs.diffMeans_l2_normPCA        ;            
    //varianceGroupOnePCA          =  varianceGroupOnePCA          + rhs.varianceGroupOnePCA         ;           
    traceVarianceGroupOnePCA     =  traceVarianceGroupOnePCA     + rhs.traceVarianceGroupOnePCA    ;                
    detVarianceGroupOnePCA       =  detVarianceGroupOnePCA       + rhs.detVarianceGroupOnePCA      ;              
    logDetVarianceGroupOnePCA    =  logDetVarianceGroupOnePCA    + rhs.logDetVarianceGroupOnePCA   ;                 
    //varianceGroupTwoPCA          =  varianceGroupTwoPCA          + rhs.varianceGroupTwoPCA         ;           
    traceVarianceGroupTwoPCA     =  traceVarianceGroupTwoPCA     + rhs.traceVarianceGroupTwoPCA    ;                
    detVarianceGroupTwoPCA       =  detVarianceGroupTwoPCA       + rhs.detVarianceGroupTwoPCA      ;              
    logDetVarianceGroupTwoPCA    =  logDetVarianceGroupTwoPCA    + rhs.logDetVarianceGroupTwoPCA   ;                 
    secondsToCompute             =  secondsToCompute             + rhs.secondsToCompute            ;        
    return *this;
}
SVM_separationResults & SVM_separationResults::operator*=(const double &rhs) 
{
    groupOneCount                =  (int)floor((double)groupOneCount                * rhs);     
    groupTwoCount                =  (int)floor((double)groupTwoCount                * rhs);     
    resampleGroupOneCount        =  (int)floor((double)resampleGroupOneCount        * rhs);             
    resampleGroupTwoCount        =  (int)floor((double)resampleGroupTwoCount        * rhs);             
    separationPercentage         =  separationPercentage         * rhs;            
    resampleSeparationPercentage =  resampleSeparationPercentage * rhs;                    
    numSingularValuesTaken       =  (int)floor((double)numSingularValuesTaken       * rhs);
    diffMeans                    =  diffMeans                    * rhs; 
    diffMeans_l2_norm            =  diffMeans_l2_norm            * rhs;         
    varianceGroupOne             =  varianceGroupOne             * rhs;        
    traceVarianceGroupOne        =  traceVarianceGroupOne        * rhs;             
    detVarianceGroupOne          =  detVarianceGroupOne          * rhs;           
    logDetVarianceGroupOne       =  logDetVarianceGroupOne       * rhs;              
    varianceGroupTwo             =  varianceGroupTwo             * rhs;        
    traceVarianceGroupTwo        =  traceVarianceGroupTwo        * rhs;             
    detVarianceGroupTwo          =  detVarianceGroupTwo          * rhs;           
    logDetVarianceGroupTwo       =  logDetVarianceGroupTwo       * rhs;              
    diffMeansPCA                 =  diffMeansPCA                 * rhs;    
    diffMeans_l2_normPCA         =  diffMeans_l2_normPCA         * rhs;            
    varianceGroupOnePCA          =  varianceGroupOnePCA          * rhs;           
    traceVarianceGroupOnePCA     =  traceVarianceGroupOnePCA     * rhs;                
    detVarianceGroupOnePCA       =  detVarianceGroupOnePCA       * rhs;              
    logDetVarianceGroupOnePCA    =  logDetVarianceGroupOnePCA    * rhs;                 
    varianceGroupTwoPCA          =  varianceGroupTwoPCA          * rhs;           
    traceVarianceGroupTwoPCA     =  traceVarianceGroupTwoPCA     * rhs;                
    detVarianceGroupTwoPCA       =  detVarianceGroupTwoPCA       * rhs;              
    logDetVarianceGroupTwoPCA    =  logDetVarianceGroupTwoPCA    * rhs;                 
    secondsToCompute             =  (int)floor((double)secondsToCompute             * rhs);        
    return *this;
}
SVM_separationResults & SVM_separationResults::operator/=(const double &rhs) 
{
    *this *= (double)(1.0/rhs);
    return *this;
      
}
SVM_separationResults & SVM_separationResults::operator= (const SVM_separationResults& rhs)
{
    groupOneCount                =   rhs.groupOneCount;
    groupTwoCount                =   rhs.groupTwoCount;
    resampleGroupOneCount        =   rhs.resampleGroupOneCount;
    resampleGroupTwoCount        =   rhs.resampleGroupTwoCount;
    separationPercentage         =   rhs.separationPercentage;
    resampleSeparationPercentage =   rhs.resampleSeparationPercentage;
    numSingularValuesTaken       =   rhs.numSingularValuesTaken;
    diffMeans                    =   rhs.diffMeans;
    diffMeans_l2_norm            =   rhs.diffMeans_l2_norm;
    varianceGroupOne             =   rhs.varianceGroupOne;
    traceVarianceGroupOne        =   rhs.traceVarianceGroupOne;
    detVarianceGroupOne          =   rhs.detVarianceGroupOne;
    logDetVarianceGroupOne       =   rhs.logDetVarianceGroupOne;
    varianceGroupTwo             =   rhs.varianceGroupTwo;
    traceVarianceGroupTwo        =   rhs.traceVarianceGroupTwo;
    detVarianceGroupTwo          =   rhs.detVarianceGroupTwo;
    logDetVarianceGroupTwo       =   rhs.logDetVarianceGroupTwo;
    diffMeansPCA                 =   rhs.diffMeansPCA;
    diffMeans_l2_normPCA         =   rhs.diffMeans_l2_normPCA;
    varianceGroupOnePCA          =   rhs.varianceGroupOnePCA;
    traceVarianceGroupOnePCA     =   rhs.traceVarianceGroupOnePCA;
    detVarianceGroupOnePCA       =   rhs.detVarianceGroupOnePCA;
    logDetVarianceGroupOnePCA    =   rhs.logDetVarianceGroupOnePCA;
    varianceGroupTwoPCA          =   rhs.varianceGroupTwoPCA;
    traceVarianceGroupTwoPCA     =   rhs.traceVarianceGroupTwoPCA;
    detVarianceGroupTwoPCA       =   rhs.detVarianceGroupTwoPCA;
    logDetVarianceGroupTwoPCA    =   rhs.logDetVarianceGroupTwoPCA;
    secondsToCompute             =   rhs.secondsToCompute;
    return *this;
}

void SVM_separationResults::print(){
    string tempString = "";
    print(tempString);
}
void SVM_separationResults::print(const string prefix){
    print(prefix.c_str());
}
void SVM_separationResults::print(const char *prefix){
    cout << prefix << " diffMeans_l2_norm: "        << setw(6) << diffMeans_l2_norm << endl;
    cout << prefix << " traceVarianceGroupOne: "    << setw(6) << traceVarianceGroupOne << endl;
    cout << prefix << " detVarianceGroupOne: "      << setw(6) << detVarianceGroupOne << endl;
    cout << prefix << " logDetVarianceGroupOne: "   << setw(6) << logDetVarianceGroupOne << endl;
    cout << prefix << " traceVarianceGroupTwo: "    << setw(6) << traceVarianceGroupTwo << endl;
    cout << prefix << " detVarianceGroupTwo: "      << setw(6) << detVarianceGroupTwo << endl;
    cout << prefix << " logDetVarianceGroupTwo: "   << setw(6) << logDetVarianceGroupTwo << endl;
}

void SVM_separationResults::printPCA(){
    string tempString = "";
    printPCA(tempString);
}

void SVM_separationResults::printPCA(const string prefix){
    printPCA(prefix.c_str());
}
void SVM_separationResults::printPCA(const char *prefix){
    cout << prefix << " numSingularValuesTaken: "       << setw(6) << numSingularValuesTaken << endl;
    cout << prefix << " diffMeans_l2_normPCA: "         << setw(6) << diffMeans_l2_normPCA << endl;
    cout << prefix << " traceVarianceGroupOnePCA: "     << setw(6) << traceVarianceGroupOnePCA << endl;
    cout << prefix << " detVarianceGroupOnePCA: "       << setw(6) << detVarianceGroupOnePCA << endl;
    cout << prefix << " logDetVarianceGroupOnePCA: "    << setw(6) << logDetVarianceGroupOnePCA << endl;
    cout << prefix << " traceVarianceGroupTwoPCA: "     << setw(6) << traceVarianceGroupTwoPCA << endl;
    cout << prefix << " detVarianceGroupTwoPCA: "       << setw(6) << detVarianceGroupTwoPCA << endl;
    cout << prefix << " logDetVarianceGroupTwoPCA: "    << setw(6) << logDetVarianceGroupTwoPCA << endl;
}

