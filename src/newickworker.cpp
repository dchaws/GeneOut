#include "newickworker.h"
#include "kernelmethod.h"

void NewickWorker_print_topology::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);
    for (unsigned i=0;i<tempM.getRows();i++){
        *outputFile << tempM(i,0) << setw(2) << " ";
    }
   *outputFile << endl;
}

void NewickWorker_print_topology::init(ostream &out,int setPrintFormat)
{
    outputFile = &out;
    printFormat = setPrintFormat;
}

void NewickWorker_print::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,0, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);
    for (unsigned i=0;i<tempM.getRows();i++){
        *outputFile << tempM(i,0) << setw(12) << " ";
    }
    *outputFile << endl;
}

void NewickWorker_print::init(ostream &out,int setPrintFormat)
{
    outputFile = &out;
    printFormat = setPrintFormat;
}

void NewickWorker_print_splits::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,0, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);
    list <set <string> >::const_iterator lsit;
    for (lsit=splits.begin();lsit!=splits.end();lsit++)
    {
        set <string>::const_iterator sit = (*lsit).begin();
        set <string>::const_iterator temp_sit;
        *outputFile << "{";
        for (;sit!=(*lsit).end();sit++)
        {
            *outputFile << *sit;
            temp_sit = sit;
            temp_sit++;
            if (temp_sit != (*lsit).end())
            {
                *outputFile << ",";
            }
        }
        *outputFile << "} ";
    }
    *outputFile << endl;
}

void NewickWorker_print_splits::init(ostream &out,int setPrintFormat)
{
    outputFile = &out;
    printFormat = setPrintFormat;
}


void NewickWorker_print_l2_diff_topology::init(Matrix A, Matrix B, ostream &out)
{
    M1 = A;
    M2 = B;
    outputFile = &out;
}

void NewickWorker_print_l2_diff_topology::work (string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);

    *outputFile << (tempM - M1).twoNorm() << " " << (tempM - M2).twoNorm() << endl;
}

void NewickWorker_mean::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,0, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);

    if (numStringsProcessed == 0){
        totalSum = tempM;
    }
    else {
        totalSum = totalSum + tempM;
    }
    numStringsProcessed++;
}

void NewickWorker_mean::work(string A, string B)
{
    Matrix tempM1, tempM2;
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,0, treeLengthFactor,splits); 
        tempM1 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(B,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,0, treeLengthFactor,splits); 
        tempM2 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }

    if (numStringsProcessed == 0){
        totalSum = tempM1 - tempM2;
    }
    else {
        totalSum = totalSum + (tempM1 - tempM2);
    }
    numStringsProcessed++;
}

void NewickWorker_mean::clear()
{
    numStringsProcessed=0;
}

long double NewickWorker_mean::l2norm()
{
    Matrix tempM = totalSum;
    tempM = tempM * (double)(1.0/(double)numStringsProcessed);
    return tempM.twoNorm();
}

Matrix NewickWorker_mean::mean()
{
    Matrix tempM = totalSum;
    tempM = tempM * (double)(1.0/(double)numStringsProcessed);
    return tempM;
}

void NewickWorker_mean_l2norm::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,0, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);

    if (numStringsProcessed == 0){
        totalSum = tempM.twoNorm();
    }
    else {
        totalSum = totalSum + tempM.twoNormSquared();
    }
    numStringsProcessed++;
}

void NewickWorker_mean_l2norm::work(string A, string B)
{
    Matrix tempM1, tempM2;
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,0, treeLengthFactor,splits); 
        tempM1 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(B,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,0, treeLengthFactor,splits); 
        tempM2 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    if (numStringsProcessed == 0){
        totalSum = (tempM1 - tempM2).twoNormSquared();
    }
    else {
        totalSum = totalSum + (tempM1 - tempM2).twoNormSquared();
    }
    numStringsProcessed++;
}

void NewickWorker_mean_l2norm::clear()
{
    numStringsProcessed=0;
}

long double NewickWorker_mean_l2norm::mean()
{
    double returnValue = totalSum;
    returnValue *= (double)(1.0/(double)numStringsProcessed);
    return returnValue;
}


void NewickWorker_mean_topology::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);

    if (numStringsProcessed == 0){
        totalSum = tempM;
    }
    else {
        totalSum = totalSum + tempM;
    }
    numStringsProcessed++;
}

void NewickWorker_mean_topology::work(string A, string B)
{
    Matrix tempM1, tempM2;
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,1, treeLengthFactor,splits); 
        tempM1 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(B,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,1, treeLengthFactor,splits); 
        tempM2 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }

    if (numStringsProcessed == 0){
        totalSum = tempM1 - tempM2;
    }
    else {
        totalSum = totalSum + (tempM1 - tempM2);
    }
    numStringsProcessed++;
}

void NewickWorker_mean_topology::clear()
{
    numStringsProcessed=0;
}

long double NewickWorker_mean_topology::l2norm()
{
    Matrix tempM = totalSum;
    tempM = tempM * (double)(1.0/(double)numStringsProcessed);
    return tempM.twoNorm();
}

Matrix NewickWorker_mean_topology::mean()
{
    Matrix tempM = totalSum;
    tempM = tempM * (double)(1.0/(double)numStringsProcessed);
    return tempM;
}

void NewickWorker_mean_topology_l2norm::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    Matrix tempM = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
            leafNumberToNodeNumber, leafCount);

    if (numStringsProcessed == 0){
        totalSum = tempM.twoNorm();
    }
    else {
        totalSum = totalSum + tempM.twoNormSquared();
    }
    numStringsProcessed++;
}

void NewickWorker_mean_topology_l2norm::work(string A, string B)
{
    Matrix tempM1, tempM2;
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,1, treeLengthFactor,splits); 
        tempM1 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    {
        double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
        int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        int leafCount = 0; 
        double treeLength = 0;
        double treeLengthFactor = 0;
        int transition[MAX_NUM_NODES];
        int leafNumberToNodeNumber[MAX_NUM_NODES];
        list <set <string> > splits;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }
        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        NewickToDistanceMatrix(B,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,1, treeLengthFactor,splits); 
        tempM2 = distanceMatrixToMatrix(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
    }
    if (numStringsProcessed == 0){
        totalSum = (tempM1 - tempM2).twoNormSquared();
    }
    else {
        totalSum = totalSum + (tempM1 - tempM2).twoNormSquared();
    }
    numStringsProcessed++;
}

void NewickWorker_mean_topology_l2norm::clear()
{
    numStringsProcessed=0;
}

long double NewickWorker_mean_topology_l2norm::mean()
{
    double returnValue = totalSum;
    returnValue *= (double)(1.0/(double)numStringsProcessed);
    return returnValue;
}



#ifdef HAVE_LIBSVM
void NewickWorker_getSVMList::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,0, treeLengthFactor,splits); 
    svm_node    *tempSVMNode = MatrixTo_svm_node_Vector(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);

    nodeList.push_back(tempSVMNode);

}

list <svm_node *> NewickWorker_getSVMList::getNodeList ()
{
    return nodeList;
}

void NewickWorker_getSVMList::clear()
{
    nodeList.clear();
}

void NewickWorker_getSVMList_topology::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    list <set <string> > splits;
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    svm_node    *tempSVMNode = MatrixTo_svm_node_Vector(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);

    nodeList.push_back(tempSVMNode);

}

void NewickWorker_getSVMList_scaleToOne::work(string A)
{
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];
    for (int i=0;i<MAX_NUM_NODES;i++){
        for (int j=0;j<MAX_NUM_NODES;j++){
            distanceMatrix[i][j] = -1; // -1 means infinity or no connection
        }
    }
    // This is a matrix where transition[i] is the parent of node i
    for (int i=0;i<MAX_NUM_NODES;i++)
    {
        transition[i] = -1; //This signifies no parent
        leafNumberToNodeNumber[i] = -1; //This signifies no parent
    }
    //Go through the distanceMatrix
    // Before we call NewickToDistanceMatrix we need to compute the treeLength
    // and the treeLengthFactor
    int curIndex = 0; // Current index in TreeInput
    int distanceCount = 0; // Number of distances encountered.
    list <set <string> > splits;
    double curTreeLength = 0;
    // Go through TreeInput and when we encounter ':', read until no more 
    // digits are read. This is a distance.

    //cout << "TreeInput: " << TreeInput << endl;
    while ((unsigned)curIndex <= A.size()){
        // Find next ':'
        while (A[curIndex] != ':' && (unsigned)curIndex <= A.size()){
            curIndex++;
        }
        //cout << A[curIndex]; 
        curIndex++;
        distanceCount++;
        string tmpDistance = "";

        // while we have a digit
        while (((A[curIndex] >= '0' && A[curIndex] <= '9') || A[curIndex] == '.') && (unsigned)curIndex <= A.size()){ 
            //cout << A[curIndex]; 
            tmpDistance.append(A.substr(curIndex,1));
            curIndex++;
        }
        //cout << "tmpDistance: " << tmpDistance << endl;
        curTreeLength+= atof(tmpDistance.c_str());
        //cout << "curIndex " << curIndex << endl;
    }
    distanceCount--;

    treeLengthFactor = 1/curTreeLength; 
    if (DEBUG_OUTPUT >= 3){
        cout << "curTreeLength = " << curTreeLength << endl;
        cout << "treeLengthFactor = " << treeLengthFactor << endl;
        cout << "distanceCount = " << distanceCount << endl;
    }
            
    NewickToDistanceMatrix(A,0,distanceMatrix,nodeCount,leafCount,transition,\
            leafNumberToNodeNumber,1, treeLengthFactor,splits); 
    svm_node    *tempSVMNode = MatrixTo_svm_node_Vector(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);

    nodeList.push_back(tempSVMNode);
}
#endif
