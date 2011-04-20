// $Rev: 776 $ $Date: 2011-02-10 12:00:19 -0500 (Thu, 10 Feb 2011) $

/** \file kernelmethod.cpp */

#include "newickworker.h"
#include "kernelmethod.h"

using namespace::std;

// Used to track the recursion level for NewickToDistanceMatrix
int recurseLevel = 0;

double *MatrixToVector(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount)
{
    int length;
    if ((leafCount % 2) == 0){
        length = (leafCount/2)*(leafCount-1);
    }
    else {
        length = ((leafCount-1)/2)*leafCount;
    }
    double *returnMatrix = new double[length];
    int count = 0;

    // Assume leafNumberToNodeNumber has non -1 entries
    // in positions 1 ... leafCount
    for (int i=1;i<=leafCount;i++)
    {
        for (int j=i+1;j<=leafCount;j++)
        {
            returnMatrix[count] = D[min(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])]\
                                  [max(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])];
            count++;
        }
    }
    return returnMatrix;
}

svm_node *MatrixTo_svm_node_Vector(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount)
{
    int length;
    if ((leafCount % 2) == 0){
        length = (leafCount/2)*(leafCount-1);
    }
    else {
        length = ((leafCount-1)/2)*leafCount;
    }
    int nonZeroCount = 0;

    // Assume leafNumberToNodeNumber has non -1 entries
    // in positions 1 ... leafCount
    for (int i=1;i<=leafCount;i++)
    {
        for (int j=i+1;j<=leafCount;j++)
        {
            if (D[min(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])]\
                                 [max(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])] != -1){

                nonZeroCount++;
            }
        }
    }
    svm_node *returnVector = new svm_node[nonZeroCount+1];

    int svmNodeCount = 0;
    int count = 0;;
    for (int i=1;i<=leafCount;i++)
    {
        for (int j=i+1;j<=leafCount;j++)
        {
            if (D[min(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])]\
                               [max(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])] != -1){
                returnVector[svmNodeCount].index = count;
                if (DEBUG_OUTPUT >=3){
                    cout << "Setting node " << svmNodeCount << " to " << D[min(leafNumberToNodeNumber[i],\
                        leafNumberToNodeNumber[j])][max(leafNumberToNodeNumber[i],\
                            leafNumberToNodeNumber[j])] << endl;
                }
                returnVector[svmNodeCount].value = D[min(leafNumberToNodeNumber[i],\
                        leafNumberToNodeNumber[j])][max(leafNumberToNodeNumber[i],\
                            leafNumberToNodeNumber[j])];
                svmNodeCount++;
            }
            count++;
        }
    }
    returnVector[nonZeroCount].index = -1;

    return returnVector;
}

Matrix distanceMatrixToMatrix(double D[MAX_NUM_NODES][MAX_NUM_NODES], int nodeCount,\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int leafCount)
{
    int length;
    if ((leafCount % 2) == 0){
        length = (leafCount/2)*(leafCount-1);
    }
    else {
        length = ((leafCount-1)/2)*leafCount;
    }
    int nonZeroCount = 0;
    Matrix returnMatrix(length,1);

    // Assume leafNumberToNodeNumber has non -1 entries
    // in positions 1 ... leafCount
    //for (int i=1;i<=leafCount;i++)
    //{
    //    for (int j=i+1;j<=leafCount;j++)
    //    {
    //        if (D[min(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])]\
    //                             [max(leafNumberToNodeNumber[i],leafNumberToNodeNumber[j])] != -1){
    //            nonZeroCount++;
    //        }
    //    }
    //}
    //svm_node *returnVector = new svm_node[nonZeroCount+1];

    int count = 0;;
    for (int i=1;i<=leafCount;i++)
    {
        for (int j=i+1;j<=leafCount;j++)
        {
            returnMatrix(count,0) = D[min(leafNumberToNodeNumber[i],\
                leafNumberToNodeNumber[j])][max(leafNumberToNodeNumber[i],\
                    leafNumberToNodeNumber[j])];
            count++;
        }
    }

    return returnMatrix;
}

void NewickToDistanceMatrix (string &treeString, int beg, double D[MAX_NUM_NODES][MAX_NUM_NODES],\
        int &nodeCount, int &leafCount, int transition[MAX_NUM_NODES],\
        int leafNumberToNodeNumber[MAX_NUM_NODES],int distanceOne, double &treeLengthFactor, list <set <string> > &splits)
{
    recurseLevel++;
    // Increase the node count.
    nodeCount++;
    // Use this to remember this nodenumber
    int thisNodeNumber = nodeCount;
    // This is a new node. We now update the distance matrix.

    if (DEBUG_OUTPUT >= 3){
        DEBUGPAD
        cout << "NewickToDistanceMatrix" << endl;
        DEBUGPAD
        cout << "Node(" << nodeCount << ")" << endl;
        //cout << "treeString(beg): " << treeString.substr(beg, treeString.size() - beg) << endl;
        DEBUGPAD
        cout << "beg: = " << beg << endl;
    }
    // Find the node in the string treeString starting at beg.
    // The node is indicated by '(' and ')'.
    // We will also find the distance to the parent.

    int leftP = beg; // Start looking at beg for '('.

    while (treeString[leftP] != '(' && treeString[leftP] != 0){
        leftP++;
    }

    // Now find matching ')', taking into account '(' ')' pairs.
    int nestLevel = 0;   // Track number of nested '(' ')' seen.
    int rightP = leftP + 1;   // Start where we left off.
    while (treeString[rightP] != ')' || nestLevel != 0)
    {
        if (treeString[rightP] == '('){
            nestLevel++;
        }
        if (treeString[rightP] == ')' && nestLevel != 0){
            nestLevel--;
        }
        rightP++;
    }
    // Is there a way to read between the matching ( and ) and list the nodes.
    // This would be a split.
    // I think all the node names will be the characters between [(),] and the next [(),].
    // Omit the distance of course.
    int nodeFinderCount = leftP+1;
    if (DEBUG_OUTPUT >= 3){
        DEBUGPAD
        cout << "Potential split: " << endl;
    }
    set <string> thisSplit;
    while (nodeFinderCount < rightP)
    {
        if (treeString[nodeFinderCount] == ':')
        {
            nodeFinderCount++;
            while(treeString[nodeFinderCount] != '(' && treeString[nodeFinderCount] != ')' && treeString[nodeFinderCount] != ',')
            {
                nodeFinderCount++;
            }
        }
        if (DEBUG_OUTPUT >= 3){
            DEBUGPAD
        }
        // This should be a taxa name. 
        int begTaxaName = nodeFinderCount;
        int endTaxaName = nodeFinderCount;
        if(treeString[nodeFinderCount] != '(' && treeString[nodeFinderCount] != ')' && treeString[nodeFinderCount] != ',')
        {
            while(treeString[nodeFinderCount] != '(' && treeString[nodeFinderCount] != ')' && treeString[nodeFinderCount] != ',')
            {
                if (DEBUG_OUTPUT >= 3){
                    cout << treeString[nodeFinderCount];
                }
                nodeFinderCount++;
                endTaxaName++;
                if (treeString[nodeFinderCount] == ':')
                {
                    nodeFinderCount++;
                    while(treeString[nodeFinderCount] != '(' && treeString[nodeFinderCount] != ')' && treeString[nodeFinderCount] != ',')
                    {
                        nodeFinderCount++;
                    }
                }
            }
            thisSplit.insert(treeString.substr(begTaxaName, endTaxaName - begTaxaName));
            if (DEBUG_OUTPUT >= 3){
                cout << endl;
                cout << "treeString.substr: " << treeString.substr(begTaxaName, endTaxaName - begTaxaName) << endl;
            }
        }
        nodeFinderCount++;
    }

    splits.push_back(thisSplit);

    


    if (DEBUG_OUTPUT >= 3){
        DEBUGPAD
        cout << "treeString(leftP,rightP): " << treeString.substr(leftP, rightP - leftP + 1) << endl;
    }
    // Now look for distance from this node to parent.
    if (treeString[rightP+1] == ':'){
        int tempPos = rightP+2;
        double sign = 1;
        // Do some error check. We should never encounter negative numbers
        if (treeString[rightP+2] == '-'){
            if (DEBUG_OUTPUT >= 3){
                cout << "treeString[rightP+2] = "  << treeString[rightP+2] << endl;
            }
            sign = 0;
            tempPos++; // skip '-'
        }
        // Read digits and '.'
        while((treeString[tempPos] >= '0' && treeString[tempPos] <= '9') || treeString[tempPos] == '.'){
            tempPos++;
        }
        int E_tempPos = tempPos;
        double factor = 1; // This will be changed if 'E' is found in string
        string E_string = "";
        // There is an E to process.
        if (treeString[E_tempPos] == 'E'){
            E_tempPos++; // Skip E
            while((treeString[E_tempPos] >= '0' && treeString[E_tempPos] <= '9') || treeString[E_tempPos] == '.' || treeString[E_tempPos] == '-'){
                E_tempPos++;
            }
            E_string = treeString.substr(tempPos+1,E_tempPos-tempPos);
        }
        string tempString = treeString.substr(rightP+2,tempPos - rightP - 2);
        if (DEBUG_OUTPUT >= 3){
            DEBUGPAD
            cout << "Distance to parent: " << sign << "*" << tempString;
            if (E_tempPos != tempPos){
                cout << "E" << E_string;
            }
            cout << endl;

            DEBUGPAD
            cout << "thisNodeNumber: " << thisNodeNumber << endl;
            DEBUGPAD
            cout << "transition[thisNodeNumber-1]: " << transition[thisNodeNumber-1] << endl;
            //DEBUGPAD
            //cout << "atof: " << atof(tempString.c_str()) << endl;
        }
        if (E_tempPos != tempPos){
            factor = pow(10,atof(E_string.c_str()));
        }
        if (DEBUG_OUTPUT >= 3){
            DEBUGPAD
            cout << "Factor = " << factor << endl;
            DEBUGPAD
            cout << "atof(tempString.c_str()) = " << atof(tempString.c_str()) << endl;
        }

        double distanceToParent;
        if (distanceOne == 0){
            distanceToParent = atof(tempString.c_str())*factor*sign;
            if (sign == 0){
                distanceToParent = 0;
            }
        }
        else if(distanceOne == 1){
            distanceToParent = 1;
        }
        if (treeLengthFactor != 0){
            // Scale by this factor
            distanceToParent*= treeLengthFactor;

        }
        // Update distance from this new node to parent
        if (transition[thisNodeNumber-1] != -1){
            D[min(transition[thisNodeNumber-1],thisNodeNumber-1)]\
                [max(transition[thisNodeNumber-1],thisNodeNumber-1)] = distanceToParent;
            //D[transition[thisNodeNumber-1]][thisNodeNumber-1] = distanceToParent;
            //D[thisNodeNumber-1][transition[thisNodeNumber-1]] = distanceToParent;
        }


        // Use the distance to update D
        for (int l=0;l<MAX_NUM_NODES;l++)
        {
            if (transition[thisNodeNumber-1] != -1 && transition[thisNodeNumber-1] != l &&\
                    thisNodeNumber-1 != l){
                if (D[min(transition[thisNodeNumber-1],l)][max(transition[thisNodeNumber-1],l)] != -1){
                    D[min(thisNodeNumber-1,l)][max(thisNodeNumber-1,l)]\
                        = D[min(transition[thisNodeNumber-1],l)][max(transition[thisNodeNumber-1],l)]\
                        + distanceToParent;
                    //D[thisNodeNumber-1][l] = D[transition[thisNodeNumber-1]][l] + distanceToParent;
                    //D[l][thisNodeNumber-1] = D[transition[thisNodeNumber-1]][l] + distanceToParent;
                }
            }
        }
    }
    else { // There is no distance given
        if (DEBUG_OUTPUT >= 3){
            DEBUGPAD
            cout << "No distance." << endl;
        }
    }
    if (DEBUG_OUTPUT >= 3){
        cout << endl;
    }

    // Now we will parse this nodes data and 
    // make recursive calls if necessary.

    // Scan until ')' or ','
    int nextNodeBegPos = leftP + 1;
    int nextNodeEndPos = leftP + 1;
    while (nextNodeBegPos < rightP){
        nextNodeEndPos = nextNodeBegPos + 1;
        //If the next character isn't '(' then easy to handle. Leaf.
        if (treeString[nextNodeBegPos] != '('){ 
            while (treeString[nextNodeEndPos] != ')' && treeString[nextNodeEndPos] != ','){
                nextNodeEndPos++;
            }
            // Now parse the node

            // This call of NewickToDistanceMatrix is the parent of this node.
            transition[nodeCount] = thisNodeNumber-1;
            nodeCount++;
            leafCount++;
            if(DEBUG_OUTPUT >= 3){
                cout << endl;
                DEBUGPAD
                cout <<  "LEAF:" << endl;
                DEBUGPAD
                cout << "Node(" << nodeCount << "): " << treeString.substr(nextNodeBegPos,\
                        nextNodeEndPos - nextNodeBegPos) << endl; 
                DEBUGPAD
                cout << "nodeCount: " << nodeCount << endl;
                DEBUGPAD
                cout << "leafCount: " << leafCount << endl;
                DEBUGPAD
                cout << "transition[nodeCount-1]: " << transition[nodeCount-1] << endl;
            }
            int tempPos;
            //See if there is a name.
            if (treeString[nextNodeBegPos] != ':'){
                //Everything up until ')', ',' or ':' is the name. We will assume it is a number.
                tempPos = nextNodeBegPos;
                while(treeString[tempPos] != ':' && treeString[tempPos] != ')' &&\
                        treeString[tempPos] != ','){
                    tempPos++;
                }
                //Now everything from nextNodeBegPos to tempPos-1 is the name
                set <string> thisLeafSplit;
                thisLeafSplit.insert(treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos ));
                splits.push_back(thisLeafSplit);
                if(DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "Node name: " << treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos )\
                    << endl;
                    DEBUGPAD
                    cout << "isdigit(((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str())[0]): "\
                        << isdigit(((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str())[0]) << endl;
                    DEBUGPAD
                    const char *tmptmpStr = ((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str());
                    DEBUGPAD
                    cout << "tmptmpStr = " << tmptmpStr << endl;
                    DEBUGPAD
                    cout << "isdigit(tmptmpStr[0]) = " << isdigit(tmptmpStr[0]) << endl;
                    DEBUGPAD
                    cout << "((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str() = " << ((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str()) << endl;
                    DEBUGPAD
                    cout << "atoi: " \
                    << atoi((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str()) << endl;
                }
                // This will reference the position in D corresponding to this nodes 'name'
                if (isdigit(((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str())[0]) != 0){
                    leafNumberToNodeNumber[atoi((treeString.substr(nextNodeBegPos,tempPos\
                                    - nextNodeBegPos)).c_str())] = nodeCount-1;
                }
                else {
                    if(DEBUG_OUTPUT >= 3){
                        DEBUGPAD
                        cout << "Name \"" << isdigit((treeString.substr(nextNodeBegPos,tempPos - nextNodeBegPos)).c_str()[0]) << "\" is not a digit." << endl;
                    }
                }
            }
            else {
                if(DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "No name" << endl;
                }
            }
            
            //Try to find ':' for distance
            tempPos = nextNodeBegPos;
            while(treeString[tempPos] != ':' && treeString[tempPos] != ')' && treeString[tempPos] != ','){
                tempPos++;
            }
            double distance;
            if (treeString[tempPos] == ':'){
                // Do some error check. We should never encounter negative numbers
                double sign = 1;
                if (treeString[tempPos+1] == '-'){
                    if (DEBUG_OUTPUT >= 3){
                        cout << "treeString[tempPos+1] = "  << treeString[tempPos+1] << endl;
                    }
                    sign = 0;
                    tempPos++; // skip '-'
                }
                if(DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "':' found." << endl;
                    DEBUGPAD
                    cout << " treeString.substr(tempPos+1,nextNodeEndPos-tempPos-1): "\
                        << treeString.substr(tempPos+1,nextNodeEndPos-tempPos-1) << endl;

                }
                // There might be 'E' in the distance. Get the digits first.
                int tempTempPos = tempPos + 1;
                while((treeString[tempTempPos] >= '0' && treeString[tempTempPos] <= '9') || treeString[tempTempPos] == '.'){
                    tempTempPos++;
                }
                int E_tempPos = tempTempPos;
                double factor = 1; // This will be changed if 'E' is found in string
                string E_string = "";
                // There is an E to process.
                if (treeString[E_tempPos] == 'E'){
                    E_tempPos++; // Skip E
                    while((treeString[E_tempPos] >= '0' && treeString[E_tempPos] <= '9') || treeString[E_tempPos] == '.' || treeString[E_tempPos] == '-'){
                        E_tempPos++;
                    }
                    E_string = treeString.substr(tempTempPos+1,E_tempPos-tempTempPos);
                }
                // This should hold all the digits
                string tempString = treeString.substr(tempPos+1,tempTempPos - tempPos);
                if (DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "Distance to parent: " << sign << "*" << tempString;
                    if (E_tempPos != tempTempPos){
                        cout << "E" << E_string;
                    }
                    cout << endl;

                    DEBUGPAD
                    cout << "thisNodeNumber: " << thisNodeNumber << endl;
                    DEBUGPAD
                    cout << "transition[thisNodeNumber-1]: " << transition[thisNodeNumber-1] << endl;
                    //DEBUGPAD
                    //cout << "atof: " << atof(tempString.c_str()) << endl;
                }
                if (E_tempPos != tempTempPos){
                    factor = pow(10,atof(E_string.c_str()));
                }
                if (DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "Factor = " << factor << endl;
                    DEBUGPAD
                    cout << "distanceOne = " << distanceOne << endl;
                    DEBUGPAD
                    cout << "treeLengthFactor = " << treeLengthFactor << endl;
                    DEBUGPAD
                    cout << "atof(tempString.c_str()) = " << atof(tempString.c_str()) << endl;
                }

                if (distanceOne == 0){
                    distance = atof(tempString.c_str())*factor*sign; 
                    if (sign == 0){
                        distance = 0;
                    }
                } else if (distanceOne == 1){
                    distance = 1;
                }
                if (treeLengthFactor != 0){
                    // Scale distance by this factor
                    distance *= treeLengthFactor;
                }
                if(DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "Distance: " << distance << endl;
                }
            }
            else {
                if(DEBUG_OUTPUT >= 3){
                    DEBUGPAD
                    cout << "No distance." << endl;
                }
            }

            // This is a new node. We now update the distance matrix. Easy to do for leaf.
            D[min(transition[nodeCount-1],nodeCount-1)][max(transition[nodeCount-1],nodeCount-1)]\
                = distance;
            //D[transition[nodeCount-1]][nodeCount-1] = distance;
            //D[nodeCount-1][transition[nodeCount-1]] = distance;
            // Use the distance to update D
            for (int l=0;l<MAX_NUM_NODES;l++)
            {
                if (transition[nodeCount-1] != -1 && transition[nodeCount-1] != l && nodeCount-1 != l){
                    if (D[min(transition[nodeCount-1],l)][max(transition[nodeCount-1],l)] != -1){
                        D[min(nodeCount-1,l)][max(nodeCount-1,l)] = D[min(transition[nodeCount-1],l)]\
                                                                    [max(transition[nodeCount-1],l)]\
                                                                    + distance;
                        //D[nodeCount-1][l] = D[transition[nodeCount-1]][l] + distance;
                        //D[l][nodeCount-1] = D[transition[nodeCount-1]][l] + distance;
                    }
                }
            }
        }
        else { // Else the next character is '('
            if(DEBUG_OUTPUT >= 3){
                DEBUGPAD
                cout << "Recurse:" << endl;
            }
            // We are about to recurse, so setup the transition matrix.
            // nodeCount will be incrased once NewickToDistanceMatrix is called
            // so nodeCount+1 will be the next node.
            transition[nodeCount] = thisNodeNumber-1;
            NewickToDistanceMatrix(treeString,nextNodeBegPos,D,nodeCount,leafCount,transition,\
                    leafNumberToNodeNumber, distanceOne, treeLengthFactor, splits);

            // We need to find the matching closing ')' to continue parsing.
            nestLevel = 0;   // Track number of nested '(' ')' seen.
            nextNodeEndPos = nextNodeBegPos + 1;   // Start where we left off.
            while (treeString[nextNodeEndPos] != ')' || nestLevel != 0)
            {
                if (treeString[nextNodeEndPos] == '('){
                    nestLevel++;
                }
                if (treeString[nextNodeEndPos] == ')' && nestLevel != 0){
                    nestLevel--;
                }
                nextNodeEndPos++;
            }
            // If the we have a distance, scan past it.
            if (treeString[nextNodeEndPos+1] == ':'){
                nextNodeEndPos++;
                while(treeString[nextNodeEndPos] != ',' && treeString[nextNodeEndPos] != ')'){
                    nextNodeEndPos++;
                }
            }
        }
        if (DEBUG_OUTPUT >= 3){
            //cout << "treeString[nextNodeEndPos] = " << treeString[nextNodeEndPos] << endl;
        }
        nextNodeBegPos = nextNodeEndPos + 1;
    }


    if (DEBUG_OUTPUT >= 3){
        DEBUGPAD
        cout << "leftP = " << leftP << "   " << treeString[leftP] << endl;
        DEBUGPAD
        cout << "rightP = " << rightP << "   " << treeString[rightP]  << endl;
        //string tempString;
        //tempString = treeString.substr(leftP, rightP-leftP+1);
        //cout << "tempString: " << tempString << endl;
    }
    recurseLevel--;
}

void permute_svm_node(svm_node *tempNode)
{
    int largestIndex=0;
    while (tempNode[largestIndex].index != -1){
        largestIndex++;
    }
    largestIndex--;

    double tempVal;
    for(int i=0;i<30;i++){
        int randIndex1 = (rand() % (largestIndex + 1));
        int randIndex2 = (rand() % (largestIndex + 1));

        tempVal = tempNode[randIndex1].value;
        tempNode[randIndex1].value = tempNode[randIndex2].value;
        tempNode[randIndex2].value = tempVal;
    }
}

list <svm_node *> *readNexusTreesUniform(list <string> treeFileNames, SampleParameters &tsp)
{
    if (DEBUG_OUTPUT >= 0){
        cout << "readNexusTreesUniform: tsp.SVM_sampleSize: " << tsp.SVM_sampleSize;
        cout << "  numTrees: " << tsp.numTreesPerFile << endl;
    }
    int totalTrees; // This is the total number of trees considering
                    // burnin and number of runs
    int treesPerRun; // This is the number of valid trees to consider
                    // for each file considering the burnin.
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    string TreeInput;
    string nexInput;
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];

    //svm_node **return_svm_node = new svm_node*[tsp.SVM_sampleSize];
    list <svm_node *> *returnNodes = new list <svm_node *>;

    if (tsp.burninFormat == -1){ // No burnin
        totalTrees = tsp.numTreesPerFile*treeFileNames.size();
        treesPerRun = tsp.numTreesPerFile;
    }
    else if (tsp.burninFormat == 0){ // Percentage burnin
        totalTrees = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent))*treeFileNames.size();
        treesPerRun = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent));
    }
    else if (tsp.burninFormat == 1){ // Fixed number burnin
        totalTrees = (tsp.numTreesPerFile - tsp.burninNumber)*treeFileNames.size();
        treesPerRun = (tsp.numTreesPerFile - tsp.burninNumber);
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "totalTrees: " << totalTrees << endl;
        cout << "treesPerRun: " << treesPerRun << endl;
    }

    // Pick uniformly what trees to sample from the first file.
    set <int> sampleTreeIndex;
    while (sampleTreeIndex.size() != (unsigned)tsp.SVM_sampleSize){
        sampleTreeIndex.insert((int)floor( ((double)rand()/(double)RAND_MAX)*(double)totalTrees));
    }
    if (DEBUG_OUTPUT == 1){
        cout << "sampleTreeIndex.size(): " << sampleTreeIndex.size() << endl;
    }

    // Open all the run files
    list <fstream *> allFiles;
    if (DEBUG_OUTPUT >= 0){
        cout << "Number of tree Files: " << treeFileNames.size() << endl;
    }
    for (list<string>::const_iterator lit=treeFileNames.begin();lit!=treeFileNames.end();lit++){
        if (DEBUG_OUTPUT >= 0){
            cout << "Processing file name: " << *lit << endl;
        }

        list <fstream *> tempFileList; // List of run files for this first file
        fstream *tempFstream = new fstream;
        tempFstream->open((*lit).c_str(),fstream::in);
        if (!tempFstream->good()){
            cout << " ##### Error trying to open " << *lit << " #####" << endl;
            exit (1);
        }
        allFiles.push_back(tempFstream);
    }


    // Now go through the sampleTreeIndex and pick from the run
    // files according to its value for the first filename

    list <fstream *>::iterator firstFile = (allFiles.begin());
    string fileInput;
    fstream *tempFstream = *firstFile;
    //while (getline(*tempFstream,fileInput) && fileInput.find("tree rep") == string::npos){
    // Search for '=' since this should indicate end of the tree
    seekToFirstTree(*tempFstream);
    //while (getline(*tempFstream,fileInput) && ( fileInputUpper.find("TREE") != string::npos && fileInputUpper.find("=") != string::npos)){
    //    //cout << "#" << endl;
    //}
    //tempFstream->seekg(-fileInput.size(),ios_base::cur);
    // Now do the burnin for this file.
    for (int i=0;i<=(tsp.numTreesPerFile-treesPerRun);i++){
        getline(*tempFstream,fileInput);
    }
    if (tempFstream->eof()){
        cout << "Reached end of file before sampling!" << endl;
        for(set<int>::const_iterator temp_sit=sampleTreeIndex.begin();temp_sit!= sampleTreeIndex.end();temp_sit++){
            cout << *temp_sit << endl;
        }
        cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "  treesPerRun: " << treesPerRun << endl;
        cout << "totalTrees: " << totalTrees << endl;
        exit(1);
    }
    // Indicate the current tree of the file.
    int currentTree = 0;
    int currentRunNumber = 1;

    if (DEBUG_OUTPUT >= 3){
        cout << "First tree rep line: " << endl;
        cout << fileInput << endl;
        //getline(*tempFstream,fileInput);
        //cout << "After seek, next line: " << endl;
        //cout << fileInput << endl;

    }
    int sampleCount = 1;
    for(set<int>::const_iterator sit=sampleTreeIndex.begin();sit!= sampleTreeIndex.end();sit++){
        int runNumber = (int)floor((double)*sit/treesPerRun)+1;
        int filePos = *sit % treesPerRun;

        // We should read from <basefilename>_OUT.run<runNumber>.t
        // We should read line (*sit % treesPerRun) + (tsp.numTreesPerFile - treesPerRun)
        // and convert it to a distance vector.
        
        if (DEBUG_OUTPUT >= 3){
            cout << endl;
            cout << "*sit: " << *sit << "  runNumber: " << runNumber;
            cout << "  filePos: " << filePos << "  sampleCount: " << sampleCount << endl;
            cout << "   OLD: currentTree: " << currentTree << "  currentRunNumber: " <<\
                currentRunNumber << endl;
        }
      
        if (runNumber > currentRunNumber) { //We need to move on to the next file

            while (runNumber > currentRunNumber){
                firstFile++;
                currentRunNumber++;
            }
            tempFstream = *firstFile;
            if (DEBUG_OUTPUT >= 1){
                cout << "   Moving on to new file." << endl;
                //cout << "   " << *tempFstream << endl;
            }
            // Why are we searching for tree rep? Not good enough
            //while (getline(*tempFstream,fileInput) && fileInput.find("tree rep") == string::npos){
            // Search for '=' since this should indicate end of the tree
            seekToFirstTree(*tempFstream);
            //while (getline(*tempFstream,fileInput) && fileInput.find("=") == string::npos){
            //    //cout << "##" << endl;
            //}
            //tempFstream->seekg(-fileInput.size(),ios_base::cur);
            // Now do the burnin for this file.
            for (int i=0;i<=(tsp.numTreesPerFile-treesPerRun);i++){
                getline(*tempFstream,fileInput);
            }
            // Indicate the current tree of the file.
            currentTree = 1;
            //currentRunNumber = runNumber;
        }

        while (currentTree < filePos) {
            getline(*tempFstream,fileInput);
            currentTree++;
        }
        if (tempFstream->eof()){
            cout << "Reached end of file!" << endl;
            cout << "filePos: " << filePos << "    runNumber: " << runNumber << "   currentRunNumber: " << currentRunNumber << endl;
            cout << "*sit: " << *sit << endl;

            for(set<int>::const_iterator temp_sit=sampleTreeIndex.begin();temp_sit!= sampleTreeIndex.end();temp_sit++){
                cout << *temp_sit << endl;
            }
            cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "  treesPerRun: " << treesPerRun << endl;
            cout << "totalTrees: " << totalTrees << endl;

            //totalTrees = floor(tsp.numTreesPerFile*(1-tsp.burninPercent))*treeFileNames.size();
            cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "   treeFileNames.size(): " << treeFileNames.size() << endl;
            for (list<string>::const_iterator lit=treeFileNames.begin();lit!=treeFileNames.end();lit++){
                cout << "*lit: " << *lit << endl;
            }
            exit(1);
        }

        // Now convert to distance vector

        if (DEBUG_OUTPUT >= 3) {
            cout << "fileInput: " << fileInput << endl;
        }
        string TreeInput = fileInput.substr(fileInput.find("("),fileInput.size()-fileInput.find("(")-1);
        if (DEBUG_OUTPUT >= 3){
            cout << "String" << endl << TreeInput << endl; // << "NewickToDistanceMatrix:"
        }
        //cout << "found." << endl;
        for (int i=0;i<MAX_NUM_NODES;i++){
            for (int j=0;j<MAX_NUM_NODES;j++){
                distanceMatrix[i][j] = -1; // -1 means infinity or no connection
            }
        }

        
        nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
        leafCount = 0;
        treeLength = 0;
        treeLengthFactor = 0;

        // This is a matrix where transition[i] is the parent of node i
        for (int i=0;i<MAX_NUM_NODES;i++)
        {
            transition[i] = -1; //This signifies no parent
            leafNumberToNodeNumber[i] = -1; //This signifies no parent
        }
        // Now if told, we can scale treeLength to 1 by dividing all distances by treeLength
        if (tsp.scaleToOne == 1 && tsp.distanceOne != 1){
            //Go through the distanceMatrix
            // Before we call NewickToDistanceMatrix we need to compute the treeLength
            // and the treeLengthFactor
            int curIndex = 0; // Current index in TreeInput
            int distanceCount = 0; // Number of distances encountered.
            double curTreeLength = 0;
            // Go through TreeInput and when we encounter ':', read until no more 
            // digits are read. This is a distance.

            //cout << "TreeInput: " << TreeInput << endl;
            while ((unsigned)curIndex <= TreeInput.size()){
                // Find next ':'
                while (TreeInput[curIndex] != ':' && (unsigned)curIndex <= TreeInput.size()){
                    curIndex++;
                }
                //cout << TreeInput[curIndex]; 
                curIndex++;
                distanceCount++;
                string tmpDistance = "";

                // while we have a digit
                while (((TreeInput[curIndex] >= '0' && TreeInput[curIndex] <= '9') || TreeInput[curIndex] == '.') && (unsigned)curIndex <= TreeInput.size()){ 
                    //cout << TreeInput[curIndex]; 
                    tmpDistance.append(TreeInput.substr(curIndex,1));
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
            

        }
        list <set <string> > splits;
        NewickToDistanceMatrix(TreeInput,0,distanceMatrix,nodeCount,leafCount,transition,\
                leafNumberToNodeNumber,tsp.distanceOne, treeLengthFactor,splits); 
        // and pass it to NewickToDistanceMatrix so it can scale properly.
        svm_node    *tempReturnValue = MatrixTo_svm_node_Vector(distanceMatrix, nodeCount,\
                leafNumberToNodeNumber, leafCount);
        
        // Test a theory 
        //cout << "Permuting svm_nodes." << endl;
        //permute_svm_node(tempReturnValue);

        returnNodes->push_back(tempReturnValue);
        //cout << "returnNodes->size(): " << returnNodes->size() << endl;;
        if (DEBUG_OUTPUT >= 3){
                cout << "svm_node vector:    ";
                svm_node *x = tempReturnValue;
                //svm_node *x = *(returnNodes->begin());
                //svm_node *x = MatrixTo_svm_node_Vector(distanceMatrix, nodeCount, leafNumberToNodeNumber, leafCount);
                cout.setf(ios::fixed,ios::floatfield);
                cout.precision(8);
                int j = 0;
                while(x[j].index != -1){
                    //cout << "I: " << x[j].index << " V: ";
                    //cout << setw(8) << x[j].value << " "; 
                    cout << j+1 << ":" << setw(8) << x[j].value << " "; 
                    //cout << x[j].value << " "; 
                    j++;
                }
                cout << ";" << endl;
        }
        // Test that our distance vector is correct with neighbor and treedist
        if (1 == 0){
            fstream distMatrix;
            distMatrix.open("tempdist",fstream::out);
            if (!distMatrix.good()){
                cout << "Could not open tempdist for output." << endl;
                exit (1);
            }
            distMatrix << leafCount << endl;
            for (int j=1;j<=leafCount;j++){
                distMatrix << j << setw(10) << " ";
                for (int k=1;k<=leafCount;k++){
                    if (j==k){
                        distMatrix << "0.00000000 ";
                    }
                    else {
                        distMatrix << setw(8) << distanceMatrix[min(leafNumberToNodeNumber[j],leafNumberToNodeNumber[k])][max(leafNumberToNodeNumber[j],leafNumberToNodeNumber[k])] << " ";
                    }
                }
                distMatrix << endl;
            }
            distMatrix.close();
            system("./neighbor < runneighbor.txt");
            fstream treeFile;
            treeFile.open("outtree",fstream::app | fstream::out);
            if (!treeFile.good()){
                cout << "treeFile not good." << endl;
                exit (1);
            }
            treeFile << TreeInput << ";";
            treeFile.close();
            string copyCommand = "cp outtree outtree";
            char tempstr[60];
            sprintf(tempstr,"%d",sampleCount);
            copyCommand.append(tempstr);
            system(copyCommand.c_str());

            system("./treedist < runtreedist.txt");
            copyCommand = "cp outfile outfile";
            copyCommand.append(tempstr);
            system(copyCommand.c_str());
        }



        if (DEBUG_OUTPUT >= 3){
            cout << "nodeCount: " << nodeCount << endl;
            cout << "leafCount: " << leafCount << endl;
            cout << "Transition matrix" << endl;
            for (int i=0;i<MAX_NUM_NODES;i++)
            {
                if (transition[i] != -1){
                    cout << transition[i] << " "; //This signifies no parent
                }
            }
            cout << endl;
            cout << "leafNumberToNodeNumber matrix:" << endl;
            for (int i=0;i<MAX_NUM_NODES;i++)
            {
                if (leafNumberToNodeNumber[i] != -1){
                    cout << leafNumberToNodeNumber[i] << " "; //This signifies no parent
                }
            }
            cout << endl;
            //cout << "Distance matrix:" << endl;
            //for (int i=0;i<nodeCount;i++){
            //    for (int j=0;j<nodeCount;j++){
            //        cout << setw(3) << distanceMatrix[i][j] << " ";
            //    }
            //    cout << endl;
            //}
            cout << endl;
        }
        if (DEBUG_OUTPUT >= 3){
            cout << "Tree: " << fileInput << endl;
            cout << "   New: currentTree: " << currentTree << "  currentRunNumber: " <<\
                currentRunNumber << endl;
        }

        // Used to control output tempo. 
        if ( 1 == 0) {
            string blah;
            getline(cin,blah);
        }
        sampleCount++;
    }
    // Close all files

    for(list <fstream *>::iterator fit=allFiles.begin();fit != allFiles.end();fit++){
        fstream *tempFstream = *fit;
        tempFstream->close();
        delete *fit;
    }

    return returnNodes;
}

void readNexusTreesUniform(list <string> treeFileNames, SampleParameters &tsp,\
        set <unsigned> treesToRead, NewickWorker &someNewickWorker)
{
    if (DEBUG_OUTPUT >= 1){
        cout << "readNexusTreesUniform (Single NewickWorker): " << tsp.SVM_sampleSize;
        cout << "  numTrees: " << tsp.numTreesPerFile << endl;
        cout << "  treesToRead.size(): " << treesToRead.size() << endl;
    }
    int totalTrees; // This is the total number of trees considering
                    // burnin and number of runs
    int treesPerRun; // This is the number of valid trees to consider
                    // for each file considering the burnin.
    double distanceMatrix[MAX_NUM_NODES][MAX_NUM_NODES];
    string TreeInput;
    string nexInput;
    int nodeCount = 0; // Need to specify node count to properly index distanceMatrix.
    int leafCount = 0; 
    double treeLength = 0;
    double treeLengthFactor = 0;
    int transition[MAX_NUM_NODES];
    int leafNumberToNodeNumber[MAX_NUM_NODES];

    // If the tree files are to contain a different number of trees, then we cant pre-compute
    if (tsp.burninFormat == -1){ // No burnin
        totalTrees = tsp.numTreesPerFile*treeFileNames.size();
        treesPerRun = tsp.numTreesPerFile;
    }
    else if (tsp.burninFormat == 0){ // Percentage burnin
        totalTrees = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent))*treeFileNames.size();
        treesPerRun = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent));
    }
    else if (tsp.burninFormat == 1){ // Fixed number burnin
        totalTrees = (tsp.numTreesPerFile - tsp.burninNumber)*treeFileNames.size();
        treesPerRun = (tsp.numTreesPerFile - tsp.burninNumber);
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "totalTrees: " << totalTrees << endl;
        cout << "treesPerRun: " << treesPerRun << endl;
    }

    if (DEBUG_OUTPUT >= 1){
        cout << "treesToRead.size(): " << treesToRead.size() << endl;
    }

    // Open all the run files
    list <fstream *> allFiles;
    if (DEBUG_OUTPUT >= 1){
        cout << "Number of tree Files: " << treeFileNames.size() << endl;
    }
    for (list<string>::const_iterator lit=treeFileNames.begin();lit!=treeFileNames.end();lit++){
        if (DEBUG_OUTPUT >= 1){
            cout << "Processing file name: " << *lit << endl;
        }

        list <fstream *> tempFileList; // List of run files for this first file
        fstream *tempFstream = new fstream;
        tempFstream->open((*lit).c_str(),fstream::in);
        if (!tempFstream->good()){
            cout << " ##### Error trying to open " << *lit << " #####" << endl;
            exit (1);
        }
        allFiles.push_back(tempFstream);
    }


    // Now go through the treesToRead and pick from the run
    // files according to its value for the first filename

    list <fstream *>::iterator firstFile = (allFiles.begin());
    string fileInput;
    fstream *tempFstream = *firstFile;
    //while (getline(*tempFstream,fileInput) && fileInput.find("tree rep") == string::npos){
    // Search for '=' since this should indicate the first tree
    seekToFirstTree(*tempFstream);
    //while (getline(*tempFstream,fileInput) && fileInput.find("=") == string::npos){
    //    //cout << "#" << endl;
    //}
    //tempFstream->seekg(-fileInput.size(),ios_base::cur);
    // Now do the burnin for this file.
    for (int i=0;i<=(tsp.numTreesPerFile-treesPerRun);i++){
        getline(*tempFstream,fileInput);
    }
    if (tempFstream->eof()){
        cout << "Reached end of file before sampling!" << endl;
        for(set<unsigned>::const_iterator temp_sit=treesToRead.begin();temp_sit!= treesToRead.end();temp_sit++){
            cout << *temp_sit << endl;
        }
        cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "  treesPerRun: " << treesPerRun << endl;
        cout << "totalTrees: " << totalTrees << endl;
        exit(1);
    }
    // Indicate the current tree of the file.
    int currentTree = 0;
    int currentRunNumber = 1;

    if (DEBUG_OUTPUT >= 3){
        cout << "First tree rep line: " << endl;
        cout << fileInput << endl;
        //getline(*tempFstream,fileInput);
        //cout << "After seek, next line: " << endl;
        //cout << fileInput << endl;

    }
    int sampleCount = 1;
    for(set<unsigned>::const_iterator sit=treesToRead.begin();sit!= treesToRead.end();sit++){
        int runNumber = (int)floor((double)*sit/treesPerRun)+1;
        int filePos = *sit % treesPerRun; // This is considering after burn in

        // We should read from <basefilename>_OUT.run<runNumber>.t
        // We should read line (*sit % treesPerRun) + (tsp.numTreesPerFile - treesPerRun)
        // and convert it to a distance vector.
        
        if (DEBUG_OUTPUT >= 2){
            cout << endl;
            cout << "*sit: " << *sit << "  runNumber: " << runNumber;
            cout << "  filePos: " << filePos << "  sampleCount: " << sampleCount << endl;
            cout << "   OLD: currentTree: " << currentTree << "  currentRunNumber: " <<\
                currentRunNumber << endl;
        }
      
        if (runNumber > currentRunNumber) { //We need to move on to the next file
            if (DEBUG_OUTPUT >= 1){
                cout << "   Moving on to new file." << endl;
                //cout << "   " << *tempFstream << endl;
            }
            while (runNumber > currentRunNumber){
                firstFile++;
                currentRunNumber++;
            }
            tempFstream = *firstFile;
            seekToFirstTree(*tempFstream);
            for (int i=0;i<=(tsp.numTreesPerFile-treesPerRun);i++){
                getline(*tempFstream,fileInput);
            }
            // Indicate the current tree of the file.
            currentTree = 0;
            //currentRunNumber = runNumber;
        }
        if (tempFstream->eof()){
            cout << "Reached end of file! Pre getLine." << endl;
        }

        while (currentTree < filePos) {
            if (DEBUG_OUTPUT >= 3){
                cout << "       getline(*tempFstream,fileInput);" << endl;
                cout << "       currentTree++; " << endl;
            }
            getline(*tempFstream,fileInput);
            currentTree++;
            if (DEBUG_OUTPUT >= 3){
                cout << "           done." << endl;
            }
        }
        if (tempFstream->eof()){
            cout << "Reached end of file!" << endl;
            cout << "filePos: " << filePos << "    runNumber: " << runNumber << "   currentRunNumber: " << currentRunNumber << endl;
            cout << "*sit: " << *sit << endl;

            for(set<unsigned>::const_iterator temp_sit=treesToRead.begin();temp_sit!= treesToRead.end();temp_sit++){
                cout << *temp_sit << endl;
            }
            cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "  treesPerRun: " << treesPerRun << endl;
            cout << "totalTrees: " << totalTrees << endl;

            //totalTrees = floor(tsp.numTreesPerFile*(1-tsp.burninPercent))*treeFileNames.size();
            cout << "tsp.numTreesPerFile: " << tsp.numTreesPerFile << "   treeFileNames.size(): " << treeFileNames.size() << endl;
            for (list<string>::const_iterator lit=treeFileNames.begin();lit!=treeFileNames.end();lit++){
                cout << "*lit: " << *lit << endl;
            }
            exit(1);
        }

        if (DEBUG_OUTPUT >= 2) {
            cout << "string TreeInput = fileInput.substr(fileInput.find(\"(\"),fileInput.size()-fileInput.find(\"(\")-1); " << endl;
        }
        string TreeInput = fileInput.substr(fileInput.find("("),fileInput.size()-fileInput.find("(")-1);
        if (DEBUG_OUTPUT >= 2) {
            cout << "fileInput: " << fileInput << endl;
        }
        if (DEBUG_OUTPUT >= 3){
            cout << "String" << endl << TreeInput << endl; // << "NewickToDistanceMatrix:"
        }

        someNewickWorker.work(TreeInput);

        if (DEBUG_OUTPUT >= 3){
            cout << "Tree: " << fileInput << endl;
            cout << "   New: currentTree: " << currentTree << "  currentRunNumber: " <<\
                currentRunNumber << endl;
        }

        sampleCount++;
    }
    // Close all files

    for(list <fstream *>::iterator fit=allFiles.begin();fit != allFiles.end();fit++){
        fstream *tempFstream = *fit;
        tempFstream->close();
        delete *fit;
    }

}

void readNexusTreesUniformAllPairs(list <string> treeFileNames, SampleParameters &tsp,\
        set <unsigned> treesToReadOne, set <unsigned> treesToReadTwo, NewickWorker &someNewickWorker)
{
    if (DEBUG_OUTPUT >= 1){
        cout << "readNexusTreesUniform:  " << tsp.SVM_sampleSize;
        cout << "  numTrees: " << tsp.numTreesPerFile << endl;
    }
    int totalTrees; // This is the total number of trees considering
                    // burnin and number of runs
    int treesPerRun; // This is the number of valid trees to consider
                    // for each file considering the burnin.
    int burninPerTreeFile; // This is the burnin number for each file.

    if (tsp.burninFormat == -1){ // No burnin
        totalTrees = tsp.numTreesPerFile*treeFileNames.size();
        treesPerRun = tsp.numTreesPerFile;
        burninPerTreeFile = 0;
    }
    else if (tsp.burninFormat == 0){ // Percentage burnin
        totalTrees = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent))*treeFileNames.size();
        treesPerRun = (int)floor(tsp.numTreesPerFile*(1-tsp.burninPercent));
        burninPerTreeFile = (int)floor((double)tsp.numTreesPerFile*(tsp.burninPercent));
    }
    else if (tsp.burninFormat == 1){ // Fixed number burnin
        totalTrees = (tsp.numTreesPerFile - tsp.burninNumber)*treeFileNames.size();
        treesPerRun = (tsp.numTreesPerFile - tsp.burninNumber);
        burninPerTreeFile = tsp.burninNumber;
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "totalTrees: " << totalTrees << endl;
        cout << "treesPerRun: " << treesPerRun << endl;
    }
    if (DEBUG_OUTPUT == 1){
        cout << "treesToReadOne.size(): " << treesToReadOne.size() << endl;
        cout << "treesToReadTwo.size(): " << treesToReadTwo.size() << endl;
    }

    // Go through and open nexusTreeFile objects for each file.
    list <nexusTreeFile *> nexusTreeFileList;
    for (list <string>::iterator sit=treeFileNames.begin();sit!=treeFileNames.end();sit++){
        nexusTreeFile *tempNexusTreeFile;
        tempNexusTreeFile = new nexusTreeFile(*sit);
        nexusTreeFileList.push_back(tempNexusTreeFile);
    }

    list <nexusTreeFile *>::iterator fileOne = nexusTreeFileList.begin();

    // Indicate the current tree of the file.
    int currentTreeOne = 0;
    int currentRunNumberOne = 1;

    for (set <unsigned>::const_iterator sitOne=treesToReadOne.begin();sitOne!=treesToReadOne.end();sitOne++){
        int runNumberOne = (int)floor((double)*sitOne/treesPerRun)+1;
        int filePosOne = *sitOne % treesPerRun; // This is considering after burn in
        if (DEBUG_OUTPUT >= 1){
            cout << "runNumberOne: " << runNumberOne << endl;
            cout << "filePosOne: " << filePosOne << endl;
            cout << "currentTreeOne: " << currentTreeOne << endl;
            cout << "currentRunNumberOne: " << currentRunNumberOne << endl;
            //cout << "*fileOne: " << *fileOne << endl;
        }

        if (runNumberOne > currentRunNumberOne) { //We need to move on to the next file

            while (runNumberOne > currentRunNumberOne){
                fileOne++;
                currentRunNumberOne++;
            }
        }
        if (DEBUG_OUTPUT >= 3){
            cout << "currentRunNumberOne: " << currentRunNumberOne << endl;
        }

        // Indicate the current tree of the file.
        int currentTreeTwo = 0;
        int currentRunNumberTwo = 1;
        list <nexusTreeFile *>::iterator fileTwo = nexusTreeFileList.begin();

        string A = (*fileOne)->getTree(filePosOne + burninPerTreeFile);

        for (set <unsigned>::const_iterator sitTwo=treesToReadTwo.begin();sitTwo!=treesToReadTwo.end();sitTwo++){
            // For *sitOne and *sitTwo, determine which file they are
            // referring to and also the proper tree index for each
            // file.
            int runNumberTwo = (int)floor((double)*sitTwo/treesPerRun)+1;
            int filePosTwo = *sitTwo % treesPerRun; // This is considering after burn in
            if (DEBUG_OUTPUT >= 1){
                cout << "runNumberTwo: " << runNumberTwo << endl;
                cout << "filePosTwo: " << filePosTwo << endl;
                cout << "currentTreeTwo: " << currentTreeTwo << endl;
                cout << "currentRunNumberTwo: " << currentRunNumberTwo << endl;
                //cout << "*fileTwo: " << *fileTwo << endl;
            }

            if (runNumberTwo > currentRunNumberTwo) { //We need to move on to the next file

                while (runNumberTwo > currentRunNumberTwo){
                    fileTwo++;
                    currentRunNumberTwo++;
                }
            }
            if (DEBUG_OUTPUT >= 1){
                cout << "currentRunNumberTwo: " << currentRunNumberTwo << endl;
            }
            string B = (*fileTwo)->getTree(filePosTwo + burninPerTreeFile);
            
            someNewickWorker.work(A,B);
        }
    }

    list <nexusTreeFile *>::iterator nit = nexusTreeFileList.begin();
    for (;nit!=nexusTreeFileList.end();nit++){
        delete *nit;
    }
}

void Sample_SVM(list <string> treeFileNamesGroupOne, list <string> treeFileNamesGroupTwo, SampleParameters &tsp, list <svm_node *> &groupOne, list <svm_node *> &groupTwo)
{
    if (DEBUG_OUTPUT >= 1){
        cout << "numTrees: " << tsp.numTreesPerFile << endl;
    }
    if (treeFileNamesGroupOne.empty() || treeFileNamesGroupTwo.empty()){
        cout << "Sample_SVM requires at least two filenames." << endl;
        return;
    }
    // Here we should check the total number of trees vs the tsp.sampleSize to determine
    // which readNexus function to call
    int treesInGroupOne = 0;
    for (list <string>::const_iterator lsit=treeFileNamesGroupOne.begin();lsit!=treeFileNamesGroupOne.end();lsit++){
        if (DEBUG_OUTPUT >= 2) {
            cout << "Checking number of trees in " << *lsit << endl;
        }
        treesInGroupOne += numTreesNexus(*lsit);
    }
    int treesInGroupTwo = 0;
    for (list <string>::const_iterator lsit=treeFileNamesGroupTwo.begin();lsit!=treeFileNamesGroupTwo.end();lsit++){
        if (DEBUG_OUTPUT >= 2) {
            cout << "Checking number of trees in " << *lsit << endl;
        }
        treesInGroupTwo += numTreesNexus(*lsit);
    }
    if (DEBUG_OUTPUT >= 2) {
        cout << "Sample_SVM: treesInGroupOne = " << treesInGroupOne << endl;
        cout << "Sample_SVM: treesInGroupTwo = " << treesInGroupTwo << endl;
        cout << "Sample_SVM: tsp.sampleSize = " << tsp.SVM_sampleSize << endl;
    }
    if (treesInGroupOne < tsp.SVM_sampleSize || treesInGroupTwo < tsp.SVM_sampleSize ){
        cout << "Sample_SVM: not enough trees in either group one or two." << endl;
        cout << "Sample_SVM: treesInGroupOne = " << treesInGroupOne << endl;
        cout << "Sample_SVM: treesInGroupTwo = " << treesInGroupTwo << endl;
        cout << "Sample_SVM: tsp.sampleSize = " << tsp.SVM_sampleSize << endl;
        exit (0);
    }

    // Sample from treeFileNamesGroupOne
    
    if (DEBUG_OUTPUT >= 0){
        cout << "   Sampling " << tsp.SVM_sampleSize << " trees from group one." << endl;
    }
    int sampleStartTime = time(0);
    if (treesInGroupOne == tsp.SVM_sampleSize || tsp.numTreesPerFile == -1) {
        NewickWorker_getSVMList *myNewickWorker;
        if (tsp.distanceOne == 1 && tsp.scaleToOne == 0) {
            myNewickWorker = new NewickWorker_getSVMList_topology;
        }
        else if (tsp.distanceOne == 0 && tsp.scaleToOne == 1) {
            myNewickWorker = new NewickWorker_getSVMList_scaleToOne;
        }
        else {
            myNewickWorker = new NewickWorker_getSVMList;
        }
        set <unsigned> treesToRead;
        fillUnsignedSet(treesToRead,0,treesInGroupOne-1);
        SampleParameters new_tsp = tsp;
        new_tsp.numTreesPerFile = tsp.SVM_sampleSize;
        if (DEBUG_OUTPUT >= 3) {
            cout << "Sample_SVM: Calling readNexusTreesUniform(treeFileNamesGroupOne,new_tsp,treesToRead,*myNewickWorker)" << endl;
        }
        readNexusTreesUniform(treeFileNamesGroupOne,new_tsp,treesToRead,*myNewickWorker);
        if (DEBUG_OUTPUT >= 3) {
            cout << "done" << endl;
        }
        groupOne = myNewickWorker->getNodeList();
        delete myNewickWorker;
    }
    else {
        groupOne = *readNexusTreesUniform(treeFileNamesGroupOne, tsp);
    }
    int sampleEndTime = time(0);
    if (DEBUG_OUTPUT >= 0){
        cout << "Done sampling " << tsp.SVM_sampleSize << " trees from group one. " << sampleEndTime-sampleStartTime << " seconds." << endl << endl;
    }

    // For remaining treeFileNames;

    if (DEBUG_OUTPUT >= 0){
        cout << "   Sampling " << tsp.SVM_sampleSize << " trees from group two." << endl;
    }
    sampleStartTime = time(0);
    if (treesInGroupTwo == tsp.SVM_sampleSize || tsp.numTreesPerFile == -1) {
        NewickWorker_getSVMList *myNewickWorker;
        if (tsp.distanceOne == 1 && tsp.scaleToOne == 0) {
            myNewickWorker = new NewickWorker_getSVMList_topology;
        }
        else if (tsp.distanceOne == 0 && tsp.scaleToOne == 1) {
            myNewickWorker = new NewickWorker_getSVMList_scaleToOne;
        }
        else {
            myNewickWorker = new NewickWorker_getSVMList;
        }
        set <unsigned> treesToRead;
        fillUnsignedSet(treesToRead,0,treesInGroupTwo-1);
        SampleParameters new_tsp = tsp;
        new_tsp.numTreesPerFile = tsp.SVM_sampleSize;
        if (DEBUG_OUTPUT >= 3) {
            cout << "Sample_SVM: Calling readNexusTreesUniform(treeFileNamesGroupOne,new_tsp,treesToRead,*myNewickWorker)" << endl;
        }
        readNexusTreesUniform(treeFileNamesGroupTwo,new_tsp,treesToRead,*myNewickWorker);
        if (DEBUG_OUTPUT >= 3) {
            cout << "done" << endl;
        }
        groupTwo = myNewickWorker->getNodeList();
        delete myNewickWorker;
    }
    else {
        groupTwo = *readNexusTreesUniform(treeFileNamesGroupTwo, tsp);
    }
    sampleEndTime = time(0);
    if (DEBUG_OUTPUT >= 0){
        cout << "Done sampling " << tsp.SVM_sampleSize << " trees from group two. " << sampleEndTime-sampleStartTime << " seconds."  << endl;
    }
}

list <svm_node *> *project_SVM_nodesMatrix (list <svm_node *> &svmVecs, double cutOffValue, int numSingValues)
{
    if (DEBUG_OUTPUT >= 0){
        cout << " ***** project_SVM_nodesMatrix called *****" << endl;
        cout << "svmVecs.size() = " << svmVecs.size() << endl;
    }


    // Convert list of svm_nodes into LAPACK compatible matrix


    // Find longest vector in list
    int vectorLength = 0;
    svm_node *tempNode;
    for (list <svm_node *>::iterator svmit=svmVecs.begin();svmit != svmVecs.end();svmit++){
       tempNode = *svmit; 
       int tempIndex=0;
       while(tempNode[tempIndex].index != -1){
           tempIndex++;
       }
       if (tempIndex > vectorLength){
           vectorLength = tempIndex;
       }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "Largest vector length: " << vectorLength << endl;
    }

    int rows = svmVecs.size();
    int cols = vectorLength;

    // Have to allocate new. I think if I dont, it will try to allocate on the 
    // stack, and there is not enough room. Must free up after done;
    double *S = new double[min(rows,cols)];
    double *AT = new double[rows*cols];
    double *U = new double[rows*rows];
    double *VT = new double[cols*cols];
    double *WORK = new double[2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4)];

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
    
    //cout.setf(ios::fixed,ios::floatfield);
    //cout.precision(3);
    int rowCount = 0;
    for (list <svm_node *>::iterator svmit=svmVecs.begin();svmit != svmVecs.end();svmit++){
       tempNode = *svmit; 
       int tempIndex=0;
       //cout << "row " << rowCount+1 << ": ";
       while(tempNode[tempIndex].index != -1){
           //cout << setw(3) << tempNode[tempIndex].value << " ";
           AT[rowCount + tempIndex*rows] = tempNode[tempIndex].value;
           tempIndex++;
       }
       //cout << ";" << endl;
       rowCount++;
    }
    if (DEBUG_OUTPUT >= 2){
        cout << "A=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for(unsigned j=0;j<(unsigned)rows;j++)
        {
            for(unsigned i=0;i<(unsigned)cols;i++)
            {
                cout << setw(8) << AT[i*rows + j] << " ";
            }
            cout << endl;
        }
    }

    dgesvd_(&JOBU, &JOBVT, &M, &N, AT, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);

    if (INFO != 0)
    {
        cout << "INFO = " << INFO << endl;
        cout << "Exiting!" << endl;
        exit (1);
    }
    if (DEBUG_OUTPUT >= 2){
        cout << "U=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for(unsigned j=0;j<(unsigned)M;j++)
        {
            for(unsigned i=0;i<(unsigned)M;i++)
            {
                cout << setw(8) << U[i*M + j] << " ";
            }
            cout << endl;
        }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "S = ";
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for(unsigned j=0;j<(unsigned)min(rows,cols);j++)
        {
            cout << setw(8) << S[j] << " ";
        }
        cout << endl;
    }
    if (DEBUG_OUTPUT >= 2){
        cout << "VT=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for(unsigned j=0;j<(unsigned)N;j++)
        {
            for(unsigned i=0;i<(unsigned)N;i++)
            {
                cout << setw(8) << VT[i*N + j] << " ";
            }
            cout << endl;
        }
    }
    
    // Sum up the squares of the singular values
    double singularSquaredSum = 0;
    for (int i=0;i< min(rows,cols);i++){
        singularSquaredSum += S[i]*S[i];
    }
    

    int numNewCols; 
    if (cutOffValue == -1) {
        numNewCols = numSingValues;
        if (DEBUG_OUTPUT >= 0){
            cout << "Taking " << numNewCols << " singular vectors." << endl;
            double currentSum = 0;
            for (int i=0;i<numNewCols;i++) {
                currentSum += S[i]*S[i];
            }
            cout << "Constitutes " << currentSum/singularSquaredSum << "% of the singular vectors" << endl;
        }
        if (numNewCols < 1) {
            cout << "numNewCols < 1" << endl;
            exit (0);
        }
    }
    else {
        // Now determine the singular values that sum under cutOffValue*singularSquaredSum.
        int largestIndex = 0;
        double currentSum = 0;
        while (currentSum < singularSquaredSum*cutOffValue && largestIndex < min(rows,cols)){
            currentSum += S[largestIndex]*S[largestIndex];
            largestIndex++;
        }
        // Go back a sum since if we exited the while, then we added too much
        currentSum -= S[largestIndex]*S[largestIndex];
        largestIndex--;

        numNewCols = largestIndex; 
        if (numNewCols <= 0) {
            numNewCols = 1;
        }
        
        
        // Old method
        // Count the number of sinular values > cutOffValue
        //for (int i=0;i< min(rows,cols);i++){
        //    if (S[i] > cutOffValue){
        //        numNewCols++;
        //    }
        //}
        if (DEBUG_OUTPUT >= 0){
            cout << "singularSquaredSum = " << singularSquaredSum << endl;
            cout << "currentSum = " << currentSum << endl;
            cout << cutOffValue*100 << "% of singular value yields " << numNewCols << " singular vectors" << endl;
        }
    }

    list <svm_node *> *returnNodes = new list <svm_node *>;
    svm_node *new_svm_node;

    for(unsigned j=0;j<(unsigned)cols;j++)
    {
        new_svm_node = new svm_node[numNewCols+1];
        for(unsigned i=0;i<(unsigned)numNewCols;i++)
        {
            new_svm_node[i].index = i;
            new_svm_node[i].value = VT[j*cols + i];

        }
        new_svm_node[numNewCols].index = -1;
        returnNodes->push_back(new_svm_node);
    }

    // Free up dynamically allocated data
    delete S;
    delete AT;
    delete U;
    delete VT;
    delete WORK;
    return returnNodes;
}

list <svm_node *> *project_SVM_nodes (list <svm_node *> &svmVecs, double cutOffValue, int numSingValues)
{
    if (DEBUG_OUTPUT >= 0){
        cout << " ***** project_SVM_nodes called. *****" << endl;
        cout << "svmVecs.size() = " << svmVecs.size() << endl;
    }


    // Convert list of svm_nodes into LAPACK compatible matrix


    // Find longest vector in list
    int vectorLength = 0;
    svm_node *tempNode;
    for (list <svm_node *>::iterator svmit=svmVecs.begin();svmit != svmVecs.end();svmit++){
       tempNode = *svmit; 
       int tempIndex=0;
       while(tempNode[tempIndex].index != -1){
           tempIndex++;
       }
       if (tempIndex > vectorLength){
           vectorLength = tempIndex;
       }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "Largest vector length: " << vectorLength << endl;
    }

    int rows = svmVecs.size();
    int cols = vectorLength;

    double AT[rows*cols];
    double U[rows*rows];
    double VT[cols*cols];
    double S[min(rows,cols)];
    double WORK[2*max(3*min(rows,cols) + max(rows,cols),5*min(rows,cols)-4)];

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
    
    //cout.setf(ios::fixed,ios::floatfield);
    //cout.precision(3);
    int rowCount = 0;
    for (list <svm_node *>::iterator svmit=svmVecs.begin();svmit != svmVecs.end();svmit++){
       tempNode = *svmit; 
       int tempIndex=0;
       //cout << "row " << rowCount+1 << ": ";
       while(tempNode[tempIndex].index != -1){
           //cout << setw(3) << tempNode[tempIndex].value << " ";
           AT[rowCount + tempIndex*rows] = tempNode[tempIndex].value;
           tempIndex++;
       }
       //cout << ";" << endl;
       rowCount++;
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "A=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for(unsigned j=0;j<(unsigned)rows;j++)
        {
            for(unsigned i=0;i<(unsigned)cols;i++)
            {
                cout << setw(8) << AT[i*rows + j] << " ";
            }
            cout << endl;
        }
    }

    dgesvd_(&JOBU, &JOBVT, &M, &N, AT, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);

    if (INFO != 0)
    {
        cout << "INFO = " << INFO << endl;
        cout << "Exiting!" << endl;
        exit (1);
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "U=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(3);
        for(unsigned j=0;j<(unsigned)M;j++)
        {
            for(unsigned i=0;i<(unsigned)M;i++)
            {
                cout << setw(3) << U[i*M + j] << " ";
            }
            cout << endl;
        }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "S = ";
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(6);
        for(unsigned j=0;j<(unsigned)min(rows,cols);j++)
        {
            cout << setw(6) << S[j] << " ";
        }
        cout << endl;
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "VT=" << endl;
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(3);
        for(unsigned j=0;j<(unsigned)N;j++)
        {
            for(unsigned i=0;i<(unsigned)N;i++)
            {
                cout << setw(3) << VT[i*N + j] << " ";
            }
            cout << endl;
        }
    }
    list <svm_node *> *returnNodes = new list <svm_node *>;

    
    int numNewCols = 0; // Count the number of sinular values > cutOffValue
    for (int i=0;i< min(rows,cols);i++){
        if (S[i] > cutOffValue){
            numNewCols++;
        }
    }
    if (DEBUG_OUTPUT >= 0){
        cout << "Number of singular values < cutOffValue = " << cutOffValue << ": " << numNewCols << endl;
    }
    svm_node *new_svm_node;

    for(unsigned j=0;j<(unsigned)M;j++)
    {
        new_svm_node = new svm_node[numNewCols+1];
        for(unsigned i=0;i<(unsigned)numNewCols;i++)
        {
            new_svm_node[i].index = i;
            new_svm_node[i].value = U[i*M + j]*S[i];

        }
        new_svm_node[numNewCols].index = -1;
        returnNodes->push_back(new_svm_node);
    }

    return returnNodes;
}

void getTreesMrBayes(MrBayesParameters &tmbp, list <list <string> > &treeFileNames, MrBayesResults &MB_results)
{
    MB_results.parametersCopy = tmbp;
    char tempStr[300];
    // Run the mrbayes data
    if (DEBUG_OUTPUT >= 0){
        cout << " ***** Forming Mr Bayes batch files *****" << endl;
    }
    int         fileCount = 0;
    string      tempMrbayesBatchFileName;
    fstream     mrbayesBatchFile;
    char        tempBuff[50]; 
    list        <string> tempMrbayesBatchFileNames;
    string      fileNamePostfix = "_MB_OUT";
    for (list<string>::const_iterator lit=tmbp.inputNexFiles.begin();lit!=tmbp.inputNexFiles.end();\
            lit++){
        if (DEBUG_OUTPUT >= 0){
            cout << "   MrBayes file: " << fileCount << endl;
            cout << "   File name: " << *lit << endl;
            cout << "   Forming MrBayes batch file:" << endl;
        }
        tempMrbayesBatchFileName = tempPrefix;
        tempMrbayesBatchFileName.append("temp_mrbayes_batch");
        sprintf(tempBuff,"_PID%d_",(int)getpid());
        tempMrbayesBatchFileName.append(tempBuff);
        sprintf(tempBuff,"_%d_",fileCount);
        tempMrbayesBatchFileName.append(tempBuff);
        //tempMrbayesBatchFileName.append(*lit);
        tempMrbayesBatchFileName.append("_.txt");
        tempMrbayesBatchFileNames.push_back(tempMrbayesBatchFileName);
        if (DEBUG_OUTPUT >= 0){
            cout << "       tempMrbayesBatchFileName: " << tempMrbayesBatchFileName << endl;
        }

        mrbayesBatchFile.open(tempMrbayesBatchFileName.c_str(),fstream::out);
        if (!mrbayesBatchFile.good()){
            cout << "Error opening mrbayesBatchFile for output." << endl;
            exit (1);
        }

        mrbayesBatchFile << "begin mrbayes;" << endl;
        mrbayesBatchFile << "   execute " << *lit << ";" << endl;
        if (tmbp.MBP_nst == 1) { // JC model
            mrbayesBatchFile << "   lset nst=" << tmbp.MBP_nst << " rates=equal" << ";" << endl;
            mrbayesBatchFile << "prset statefreqpr=fixed(equal);" << endl;
        }
        else {
            mrbayesBatchFile << "   lset nst=" << tmbp.MBP_nst << " rates=" << tmbp.MBP_rates << ";" << endl;
        }
        mrbayesBatchFile << "   mcmc nruns=" << tmbp.MBP_nruns << " ngen=" << tmbp.MBP_ngen\
            << " samplefreq=" << tmbp.MBP_sampleFreq << " file=" << *lit << fileNamePostfix << ";" << endl;
        mrbayesBatchFile << "end;" << endl;
        mrbayesBatchFile.close();

        // Now for each file fill treeFileNames
        list <string> thisinputNexFilesFileNames;
        thisinputNexFilesFileNames.clear(); // Clear out any previous filenames
        string oneRunFileName;
        for (int i=0;i<tmbp.MBP_nruns;i++){
               oneRunFileName = *lit;
               oneRunFileName.append(fileNamePostfix);
               oneRunFileName.append(".run");
               sprintf(tempStr,"%d",i+1);
               oneRunFileName.append(tempStr);
               oneRunFileName.append(".t");
               thisinputNexFilesFileNames.push_back(oneRunFileName);
        }
        treeFileNames.push_back(thisinputNexFilesFileNames);
       
        fileCount++;
    }
    if (DEBUG_OUTPUT >= 0){
        cout << endl;
    }
    // All of our tempMrbayesBatch files are created. We can either run
    // the programs in series or we can attempt to run in parallel.
    

    // Run parallel mb runs using homemade script.
    if (1 == 0){
        string mrbayesExecCommand = "./runmultmrbayes";
        for (list<string>::const_iterator lit=tempMrbayesBatchFileNames.begin();\
                lit!=tempMrbayesBatchFileNames.end();lit++){
            mrbayesExecCommand.append(" ");
            mrbayesExecCommand.append(*lit);
        }
        system(mrbayesExecCommand.c_str());
    }


    // Run parallel mb using mpi IF MPI_np > 1
    if (1 == 1){
        if (tmbp.MPI_np > 1){
            if (DEBUG_OUTPUT >= 0){
                cout << " ***** Executing Mr Bayes batch files using MPI *****" << endl;
            }
            string mrbayesExecCommand;
            for (list<string>::const_iterator lit=tempMrbayesBatchFileNames.begin();\
                    lit!=tempMrbayesBatchFileNames.end();lit++){
                sprintf(tempStr,"%d ",tmbp.MPI_np);
                mrbayesExecCommand = "mpirun -np ";
                mrbayesExecCommand.append(tempStr);
                mrbayesExecCommand.append(tmbp.MBP_pathToMB);
                mrbayesExecCommand.append(" ");
                mrbayesExecCommand.append(*lit);
                mrbayesExecCommand.append(" > ");
                if (tmbp.MBP_saveOutput == 1){
                    mrbayesExecCommand.append(tempPrefix);
                    mrbayesExecCommand.append(*lit);
                    mrbayesExecCommand.append("_mrbayesoutput.txt");
                }
                else {
                    mrbayesExecCommand.append(" /dev/null");
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "   Executing " << mrbayesExecCommand << endl;
                }
                int startTime = time(0);
                system(mrbayesExecCommand.c_str());
                int endTime = time(0);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       " << endTime - startTime << " wall clock seconds." << endl << endl;
                }
            }
        }
        else if (tmbp.MPI_np == 1){
            if (DEBUG_OUTPUT >= 0){
                cout << " ***** Executing Mr Bayes batch files *****" << endl;
            }
            string mrbayesExecCommand;
            for (list<string>::const_iterator lit=tempMrbayesBatchFileNames.begin();\
                    lit!=tempMrbayesBatchFileNames.end();lit++){
                mrbayesExecCommand = tmbp.MBP_pathToMB;
                mrbayesExecCommand.append(" ");
                mrbayesExecCommand.append(*lit);
                mrbayesExecCommand.append(" > ");
                if (tmbp.MBP_saveOutput == 1){
                    mrbayesExecCommand.append(tempPrefix);
                    mrbayesExecCommand.append(*lit);
                    mrbayesExecCommand.append("_mrbayesoutput.txt");
                }
                else {
                    mrbayesExecCommand.append(" /dev/null");
                }
                if (DEBUG_OUTPUT >= 0){
                    cout << "   Executing " << mrbayesExecCommand << endl;
                }
                int startTime = time(0);
                system(mrbayesExecCommand.c_str());
                int endTime = time(0);
                if (DEBUG_OUTPUT >= 0){
                    cout << "       " << endTime - startTime << " wall clock seconds." << endl << endl;
                }
            }

        }
    }

    // Look for split frequency here.
    // Look in *lit + fileNamePostfix ".mcmc" for the last line and last argument.
    // That is the split frequency.
    for (list<string>::const_iterator lit=tmbp.inputNexFiles.begin();lit!=tmbp.inputNexFiles.end();lit++){

        // This is for using a script, but we should do it by hand.
        //string tempSplitFreqFile = tempPrefix;
        //tempSplitFreqFile.append("temp_splitfreq");
        //sprintf(tempBuff,"_PID%d_ID%d%d.txt",(int)getpid(),(int)rand() % 32768);
        //tempSplitFreqFile.append(tempBuff);

        //// tempSplitFreqFile is the file name that we want to write the split frequencies to.

        string tempFileName = *lit;
        tempFileName.append(fileNamePostfix);
        tempFileName.append(".mcmc");
        fstream tempFile;
        tempFile.open(tempFileName.c_str(),fstream::in);
        if (tempFile.bad())
        {
            cout << tempFileName << " could not be opened." << endl;
            exit (0);
        }
        // Seek to the last line and get the last field
        string fileInput, tempFileInput;
        // This reads lines until it can no longer
        while (getline(tempFile,tempFileInput))
        {
            fileInput = tempFileInput;
        }
        if (fileInput.length() == 0){
            cout << "fileInput.size()=0 in getTreesMrBayes." << endl;
            cout << "   Dumping first and last 40 lines:" << endl;
            tempFile.seekg (0, ios::beg);
            if (tempFile.bad())
            {
                cout << tempFileName << " could not be seekg'd." << endl;
                exit (0);
            }
            int nextOpenIndex=0;
            string lastFourtyLines[40];
            int numFirstLines=0;
            while (getline(tempFile,tempFileInput))
            {
                if (numFirstLines < 40) {
                    cout << tempFileInput << endl;
                    numFirstLines++;
                }
                lastFourtyLines[nextOpenIndex] = tempFileInput;
                nextOpenIndex = (nextOpenIndex + 1) % 40;
            }
            for (int i=0;i<40;i++){
                cout << lastFourtyLines[(nextOpenIndex+i+1) % 40];
            }
            if (tmbp.MBP_saveOutput == 1){
                cout << "Listing mb output file names." << endl;

                for (list<string>::const_iterator lit=tempMrbayesBatchFileNames.begin();\
                        lit!=tempMrbayesBatchFileNames.end();lit++){
                    string tempMBOutputFileName = tempPrefix;
                    tempMBOutputFileName.append(*lit);
                    tempMBOutputFileName.append("_mrbayesoutput.txt");
                    cout << tempMBOutputFileName << endl;

                    //fstream tempMBOutputFile;
                    //tempMBOutputFile.open(tempMBOutputFileName.c_str(),fstream::in);
                    //if (tempMBOutputFile.bad()){
                    //    cout << "Could not read " << tempMBOutputFileName << endl;
                    //}
                    //while (getline(tempMBOutputFile,tempFileInput)){
                    //    cout << tempFileInput << endl;
                    //}
                    //tempMBOutputFile.close();
                }
            }
            
            exit(0);
        }

        if (DEBUG_OUTPUT >= 0){
            cout << "Reading split freq from " << tempFileName << endl;
            cout << "   fileInput: " << fileInput << endl;
        }
        // Now the last field of fileInput has the split frequency.
        int lastSpacePosFromEnd=1;
        while (lastSpacePosFromEnd < fileInput.length() && ((fileInput[fileInput.length()-lastSpacePosFromEnd] >= '0' && fileInput[fileInput.length()-lastSpacePosFromEnd] <= '9') || fileInput[fileInput.length()-lastSpacePosFromEnd] == '.')) {
            lastSpacePosFromEnd++;
        }
        if (DEBUG_OUTPUT >= 0){
            cout << "   lastSpacePosFromEnd = " << lastSpacePosFromEnd << endl;
            cout << "   substr = " << (fileInput.substr(fileInput.length()-lastSpacePosFromEnd + 1, lastSpacePosFromEnd)) << endl;
        }
        double tempSplitFreq = atof((fileInput.substr(fileInput.length()-lastSpacePosFromEnd + 1, lastSpacePosFromEnd)).c_str());
        if (DEBUG_OUTPUT >= 0){
            cout << "   splitfreq: " << tempSplitFreq << endl;
        }
        MB_results.splitFreqs.push_back(tempSplitFreq);

        tempFile.close ();
    }

    string delCommand = "rm -f ";
    for (list<string>::const_iterator lit=tempMrbayesBatchFileNames.begin();\
            lit!=tempMrbayesBatchFileNames.end();lit++){
        delCommand.append(*lit);
        delCommand.append(" ");
    }
    system(delCommand.c_str());

    // We should clean up all the non-tree mb generated files. We know files will have MB_OUT to the names.
    // look for (*lit)_MB_OUT.mcmc and (*lit)_MB_OUT.run?.p
    list <string> filesToDelete;
    for (list<string>::const_iterator lit=tmbp.inputNexFiles.begin();lit!=tmbp.inputNexFiles.end();lit++){
        string tmpFileName = (*lit);
        tmpFileName.append(fileNamePostfix);
        tmpFileName.append(".mcmc");
        filesToDelete.push_back(tmpFileName);
        tmpFileName = (*lit);
        tmpFileName.append(fileNamePostfix);
        tmpFileName.append(".run*.p");
        filesToDelete.push_back(tmpFileName);
    }
    delCommand = "rm -f ";
    for (list<string>::const_iterator lit=filesToDelete.begin();lit!=filesToDelete.end();lit++){
        delCommand.append(*lit);
        delCommand.append(" ");
    }
    system(delCommand.c_str());

}

void getTreesMrBayes(MrBayesParameters &tmbp,list <Alignments *> &inputAlignments, list <list <string> > &treeFileNames, MrBayesResults &MB_results)
{
    // Write the alignments from inputAlignments down as temporary files.
    list <string> newFileNames;
    string tempAlignmentsBaseName = tempPrefix;
    tempAlignmentsBaseName.append("tempalignments_ID");
    char tmpstr[80];
    sprintf (tmpstr,"%d%d%d%d",(int)getpid(),rand(),rand(),rand());
    tempAlignmentsBaseName.append(tmpstr);

    int tempFileCount = 0;
    for (list <Alignments *>::iterator lait = inputAlignments.begin();lait != inputAlignments.end();lait++){
        string tempString = tempAlignmentsBaseName;
        sprintf(tmpstr,"%d",tempFileCount);
        tempString.append("_file");
        tempString.append(tmpstr);
        tempString.append(".nex");
        
        fstream tempFile;
        tempFile.open(tempString.c_str(),fstream::out);
        if (tempFile.bad()){
            cout << "Could not open " << tempString << endl;
            exit(0);
        }
        (*lait)->setOutputFormat(OF_NEXUS);
        tempFile << *(*lait);
        tempFile.close();
        newFileNames.push_back(tempString);
        tempFileCount++;
    }

    MrBayesParameters new_tmbp = tmbp;
    new_tmbp.inputNexFiles = newFileNames;

    getTreesMrBayes(new_tmbp, treeFileNames, MB_results);

    // We should delete these temporary alignments
    for (list <string>::iterator lsit = newFileNames.begin();lsit!=newFileNames.end();lsit++){
        string delCommand = "rm -f ";
        delCommand.append(*lsit);
        system(delCommand.c_str());
    }
}

void getTreesJackknife(JackknifeParameters &tjp,list <Alignments *> &inputAlignments, list <string> outputFileNames)
{
    if (inputAlignments.size() != outputFileNames.size()){
        cout << "getTreesJackknife: inputAlignments.size() != outputFileNames.size()" << endl;
        exit (0);
    }
    char randIdString[80];
    sprintf(randIdString,"%d%d%d%d",(int)getpid(),rand() % 100000,rand() % 100000,rand() % 100000);
    string tempJackknifeFileName = tempPrefix;
    tempJackknifeFileName.append("tempjackknifefile_ID");
    tempJackknifeFileName.append(randIdString);
    tempJackknifeFileName.append("_.nex");

    fstream tempJackknifeFile;
    //string runalignToTreeCommand = "./run_dnadist_neighbor_mult tempjackknifefile.nex ";
    //string runalignToTreeCommand = "./run_dnadist_neighbor_mult ";
    string runalignToTreeCommand = tjp.alignToTreeCommand;
    runalignToTreeCommand.append(" ");
    runalignToTreeCommand.append(tempJackknifeFileName);
    runalignToTreeCommand.append(" ");

    char tmpStr[200];
    sprintf(tmpStr,"%d",tjp.jackknifeCount);
    runalignToTreeCommand.append(tmpStr);    
    runalignToTreeCommand.append(" >> ");
    string tempCommandStr;

    list <Alignments *>::iterator lait=inputAlignments.begin();
    list <string>::const_iterator lsit=outputFileNames.begin();
    int i = 0;
    while (lait!=inputAlignments.end()){
        if (DEBUG_OUTPUT >= 0){
            cout << "   Processing " << i << " with output file name: " << *lsit << endl; 
        }
        i++;
        //cout << *(*lait);
        tempCommandStr.clear();
        tempCommandStr = runalignToTreeCommand;
        tempCommandStr.append(*lsit);
        // Write preproc stuff for tempCommandStr;
        //fstream runFile;
        //runFile.open((*lsit).c_str(),fstream::out);
        //if (!runFile.good()){
        //    cout << "Can not open file " << *lsit << endl;
        //    exit(0);
        //}

        //tempJackknifeFile.open("tempjackknifefile.nex",fstream::out);
        tempJackknifeFile.open(tempJackknifeFileName.c_str(),fstream::out);
        if (!tempJackknifeFile.good ()){
            cout << "Can not open " << tempJackknifeFileName << endl;
            exit(0);
        }
        
        //runFile << "#NEXUS\n";
        //// Why? ID is not required.
        ////runFile << "[ID: " << rand() % 1000 << "]\n";
        //runFile << "begin trees;\n";
        //runFile.close();

        int tenPercent=(int)floor((double)tjp.jackknifeCount/10);
        int jackknifeStartTime = time(0);
        int jackknifeInnerStartTime = time(0);
        for (int l=0;l<tjp.jackknifeCount;l++){
            if (DEBUG_OUTPUT >= 2){
                if (((l+1) % tenPercent) == 0){
                    cout << "   " << setw(3) << 10*(l+1)/tenPercent << "%";
                    int tempEndTime = time(0);
                    cout << "   " << tempEndTime - jackknifeInnerStartTime << " seconds." << endl;
                    jackknifeInnerStartTime = time(0);
                }
            }
            Alignments jackAlignments = (*lait)->getJackknife(tjp.jackknifeColSize);
            jackAlignments.setOutputFormat(OF_PHYLIP);
            tempJackknifeFile << jackAlignments;
        }
        tempJackknifeFile.close();
        //Now run program to generate tree and append to tempCommandStr
        if (DEBUG_OUTPUT >= 0){
            cout << "Running: \"" << tempCommandStr << "\"" << endl;
        }
        // Perhaps with MAC system, can not append to this file.
        system(tempCommandStr.c_str());

        int jackknifeEndTime = time(0);
        if (DEBUG_OUTPUT >= 0){
            cout << "Total seconds: " << jackknifeEndTime - jackknifeStartTime << endl;
        }
        
        // Place the end at the end of the filename
        // MAC is having problems here
        //runFile.open((*lsit).c_str(),fstream::app);
        //if (!runFile.good()){
        //    cout << "Can not open file (41303) \"" << *lsit << "\"" << endl;
        //    cout << "Trying to read state." << endl;
        //    if ((runFile.rdstate() & ifstream::eofbit) != 0){
        //        cout << "       eofbit set" << endl;
        //    }
        //    if ((runFile.rdstate() & ifstream::failbit) != 0){
        //        cout << "       failbit set" << endl;
        //    }
        //    if ((runFile.rdstate() & ifstream::badbit) != 0){
        //        cout << "       badbit set" << endl;
        //    }
        //    //cout << "Ommitting \"end;\" from file" << endl;
        //    //runFile.close();
        //    //cout << "Trying to open for normal output. " << endl;
        //    //runFile.open((*lsit).c_str(),fstream::out);
        //    //cout << "Seeking to end of file." << endl;
        //    //runFile.seekp(ios_base::end);
        //    //if (!runFile.good()){
        //    //    cout << "Not Good, exiting." << endl;
        //    //}
        //    //else {
        //    //    cout << "Good, continuing!" << endl;
        //    //}
        //    //exit(0);
        //    cout << "Using system command to append \"end;\"." << endl;
        //    string tmptmpcommand = "echo 'end;' >> ";
        //    tmptmpcommand.append(*lsit);
        //    system(tmptmpcommand.c_str());
        //}
        //else {
        //    runFile << "end;" << endl;
        //    runFile.close ();
        //}
        // Lets delete the tempfile, so we don't unkownling use it again
        string delTempJackknifeFile = "rm -f ";
        delTempJackknifeFile.append(tempJackknifeFileName);
        system(delTempJackknifeFile.c_str()); 
        lait++;
        lsit++;
    }
}

void getTreesJackknife(JackknifeParameters &tjp, SampleParameters &tsp, list <Alignments *> &inputAlignments, string &outputFileName)
{
    char randIdString[80];
    sprintf(randIdString,"%d%d%d%d",(int)getpid(),rand() % 100000,rand() % 100000,rand() % 100000);
    string tempJackknifeFileName = tempPrefix;
    tempJackknifeFileName.append("tempjackknifefile_ID");
    tempJackknifeFileName.append(randIdString);
    tempJackknifeFileName.append("_.nex");

    outputFileName = tempPrefix;
    outputFileName.append("temp_calcSVMseparationJackknife_ID");
    sprintf(randIdString,"%d%d%d%d",(int)getpid(),rand() % 100000,rand() % 100000,rand() % 100000);
    outputFileName.append(randIdString);
    outputFileName.append("_.nex");

    fstream tempJackknifeFile;
    //string runalignToTreeCommand = "./run_dnadist_neighbor_mult tempjackknifefile.nex ";
    //string runalignToTreeCommand = "./run_dnadist_neighbor_mult ";
    string runalignToTreeCommand = tjp.alignToTreeCommand;
    runalignToTreeCommand.append(" ");
    runalignToTreeCommand.append(tempJackknifeFileName);
    runalignToTreeCommand.append(" ");

    char tmpStr[200];
    sprintf(tmpStr,"%d",tsp.SVM_sampleSize);
    runalignToTreeCommand.append(tmpStr);    
    runalignToTreeCommand.append(" >> ");
    runalignToTreeCommand.append(outputFileName);
    if (DEBUG_OUTPUT >= 2) {
        cout << "getTreesJackknife(2): runalignToTreeCommand = " << runalignToTreeCommand << endl;
    }

    // REPLACE WITH UNIFORM -DAVE 9/22/10
    // Maybe we want to be able to specify the svm sample size for each group
    int samplesPerAlignments[inputAlignments.size()];
    if (tsp.SVM_sampleUniform == 1 ) {
        if ( tsp.SVM_sampleSize % inputAlignments.size() != 0)
        {
            cerr << "SVM_sampleUniform \% SVM_sampleSize != 0" << endl;
            exit (0);
        }
        for (int i=0;i<inputAlignments.size();i++){
            samplesPerAlignments[i]=tsp.SVM_sampleSize/inputAlignments.size();
        }
    }
    if (tsp.SVM_sampleUniform == 0 ) {
        for (int i=0;i<inputAlignments.size();i++){
            samplesPerAlignments[i]=0;
        }
        for (int i=0;i<tsp.SVM_sampleSize;i++){
            samplesPerAlignments[(rand() % inputAlignments.size())]++;
        }
    }
    if (DEBUG_OUTPUT >= 2) {
        cout << "tsp.SVM_sampleSize: " << tsp.SVM_sampleSize << endl;
        cout << "inputAlignments.size(): " << inputAlignments.size() << endl;
        for (int i=0;i<inputAlignments.size();i++){
            cout << "samplesPerAlignments[" << i << "]: " << samplesPerAlignments[i] << endl;
        }
    }

    //fstream runFile;
    //runFile.open(outputFileName.c_str(),fstream::out);
    //if (!runFile.good()){
    //    cout << "Can not open file " << outputFileName << endl;
    //    exit(0);
    //}

    tempJackknifeFile.open(tempJackknifeFileName.c_str(),fstream::out);
    if (!tempJackknifeFile.good ()){
        cout << "Can not open " << tempJackknifeFileName << endl;
        exit(0);
    }
    
    //runFile << "#NEXUS\n";
    //runFile << "begin trees;\n";
    //runFile.close();


    list <Alignments *>::iterator lait=inputAlignments.begin();
    int i = 0;
    int jackknifeStartTime = time(0);
    while (lait!=inputAlignments.end()){
        int tenPercent=(int)floor((double)tsp.SVM_sampleSize/10);
        int jackknifeInnerStartTime = time(0);
        for (int l=0;l<samplesPerAlignments[i];l++){
            if (DEBUG_OUTPUT >= 2){
                if (((l+1) % tenPercent) == 0){
                    cout << "   " << setw(3) << 10*(l+1)/tenPercent << "%";
                    int tempEndTime = time(0);
                    cout << "   " << tempEndTime - jackknifeInnerStartTime << " seconds." << endl;
                    jackknifeInnerStartTime = time(0);
                }
            }
            Alignments jackAlignments = (*lait)->getJackknife(tjp.jackknifeColSize);
            jackAlignments.setOutputFormat(OF_PHYLIP);
            tempJackknifeFile << jackAlignments;
        }

        
        lait++;
        i++;
    }
    tempJackknifeFile.close();
    //Now run program to generate tree and append to tempCommandStr
    if (DEBUG_OUTPUT >= 0){
        cout << "Running: \"" << runalignToTreeCommand << "\"" << endl;
    }
    // Perhaps with MAC system, can not append to this file.
    system(runalignToTreeCommand.c_str());
    int jackknifeEndTime = time(0);
    if (DEBUG_OUTPUT >= 0){
        cout << "Total seconds: " << jackknifeEndTime - jackknifeStartTime << endl;
    }

    // Place the end at the end of the filename
    // MAC is having problems here. So is the HPC
    //runFile.open(outputFileName.c_str(),fstream::app);
    //if (!runFile.good()){
    //    cout << "Can not open file (41303) \"" << outputFileName << "\"" << endl;
    //    cout << "Trying to read state." << endl;
    //    if ((runFile.rdstate() & ifstream::eofbit) != 0){
    //        cout << "       eofbit set" << endl;
    //    }
    //    if ((runFile.rdstate() & ifstream::failbit) != 0){
    //        cout << "       failbit set" << endl;
    //    }
    //    if ((runFile.rdstate() & ifstream::badbit) != 0){
    //        cout << "       badbit set" << endl;
    //    }
    //    //cout << "Ommitting \"end;\" from file" << endl;
    //    //runFile.close();
    //    //cout << "Trying to open for normal output. " << endl;
    //    //runFile.open((*lsit).c_str(),fstream::out);
    //    //cout << "Seeking to end of file." << endl;
    //    //runFile.seekp(ios_base::end);
    //    //if (!runFile.good()){
    //    //    cout << "Not Good, exiting." << endl;
    //    //}
    //    //else {
    //    //    cout << "Good, continuing!" << endl;
    //    //}
    //    //exit(0);
    //    cout << "Using system command to append \"end;\"." << endl;
    //    string tmptmpcommand = "echo 'end;' >> ";
    //    tmptmpcommand.append(outputFileName);
    //    system(tmptmpcommand.c_str());
    //}
    //else {
    //    runFile << "end;" << endl;
    //    runFile.close ();
    //}
    // Lets delete the tempfile, so we don't unknowing use it again
    string delTempJackknifeFile = "rm -f ";
    delTempJackknifeFile.append(tempJackknifeFileName);
    system(delTempJackknifeFile.c_str()); 
}


void getTreesJackknife(JackknifeParameters &tjp, list <list <string> > &treeFileNames)
{
    // Now we should just open all the alignments, and pass the appropriate parameters to
    // void getTreesJackknife(JackknifeParameters &tjp,list <Alignments> &inputAlignments, list <string> outputFileNames)
    list <Alignments *> myAlignments;
    list <string> outputFileNames;
    string      fileNamePostfix = "_DNADIST_and_NEIGHBOR_OUT";
    for (list <string>::const_iterator lsit=tjp.inputNexFiles.begin();lsit!=tjp.inputNexFiles.end();lsit++){
        list <string> thisinputNexFilesFileNames;
        thisinputNexFilesFileNames.clear(); // Clear out any previous filenames
        string oneRunFileName = tempPrefix;
        oneRunFileName.append(*lsit);
        oneRunFileName.append(fileNamePostfix);
        oneRunFileName.append(".run1.t");
        outputFileNames.push_back(oneRunFileName);

        thisinputNexFilesFileNames.push_back(oneRunFileName);
        treeFileNames.push_back(thisinputNexFilesFileNames);
        Alignments *newAlignments = new Alignments(*lsit);
        myAlignments.push_back(newAlignments);
        //cout << "newAlignments" << endl << newAlignments;
    }
    getTreesJackknife(tjp,myAlignments,outputFileNames);
    //Lets clean up myAlignments. 
    for (list <Alignments *>::iterator lait = myAlignments.begin();lait!= myAlignments.end();lait++){
        delete *lait;
    }
    myAlignments.clear();
}

void getTreesJackknife(JackknifeParameters &tjp, list <Alignments *> &inputAlignments, list <list <string> > &treeFileNames)
{
    string tempAlignmentsBaseName = tempPrefix;
    tempAlignmentsBaseName.append("tempalignments_ID");
    list <string> outputFileNames;
    char tmpstr[80];
    sprintf (tmpstr,"%d%d%d%d",(int)getpid(),rand(),rand(),rand());
    tempAlignmentsBaseName.append(tmpstr);

    int tempFileCount = 0;
    for (list <Alignments *>::iterator lait = inputAlignments.begin();lait != inputAlignments.end();lait++){
        list <string> thisinputNexFilesFileNames;
        thisinputNexFilesFileNames.clear(); // Clear out any previous filenames
        string tempString = tempAlignmentsBaseName;
        sprintf(tmpstr,"%d",tempFileCount);
        tempString.append("_file");
        tempString.append(tmpstr);
        tempString.append(".run1.t");
        thisinputNexFilesFileNames.push_back(tempString);
        treeFileNames.push_back(thisinputNexFilesFileNames);
        outputFileNames.push_back(tempString);
        tempFileCount++;
    }

    getTreesJackknife(tjp,inputAlignments,outputFileNames);
}

void getTreesMultInd(MultIndParameters &tmip,list <Alignments *> &inputAlignments, list <string> outputFileNames)
{
    if (inputAlignments.size() != outputFileNames.size()){
        cout << "getTreesMultInd: inputAlignments.size() != outputFileNames.size()" << endl;
        exit (0);
    }
    char randIdString[80];
    sprintf(randIdString,"%d%d%d%d",(int)getpid(),rand() % 100000,rand() % 100000,rand() % 100000);
    string tempMultIndFileName = tempPrefix;
    tempMultIndFileName.append("temp_MultInt_ID");
    tempMultIndFileName.append(randIdString);
    tempMultIndFileName.append("_.nex");

    fstream tempMultIndFile;
    string runalignToTreeCommand = tmip.alignToTreeCommand;
    runalignToTreeCommand.append(" ");
    runalignToTreeCommand.append(tempMultIndFileName);
    runalignToTreeCommand.append(" ");

    char tmpStr[200];
    sprintf(tmpStr,"%d",tmip.numTreesReconstruct);
    runalignToTreeCommand.append(tmpStr);    
    runalignToTreeCommand.append(" >> ");
    string tempCommandStr;

    list <Alignments *>::iterator lait=inputAlignments.begin();
    list <string>::const_iterator lsit=outputFileNames.begin();
    int i = 0;
    while (lait!=inputAlignments.end()){
        if (DEBUG_OUTPUT >= 0){
            cout << "   Processing " << i << " with output file name: " << *lsit << endl; 
        }
        i++;
        //cout << *(*lait);
        tempCommandStr.clear();
        tempCommandStr = runalignToTreeCommand;
        tempCommandStr.append(*lsit);
        // Write preproc stuff for tempCommandStr;
        fstream runFile;
        runFile.open((*lsit).c_str(),fstream::out);
        if (!runFile.good()){
            cout << "Can not open file " << *lsit << endl;
            exit(0);
        }

        tempMultIndFile.open(tempMultIndFileName.c_str(),fstream::out);
        if (!tempMultIndFile.good ()){
            cout << "Can not open " << tempMultIndFileName << endl;
            exit(0);
        }
        
        runFile << "#NEXUS\n";
        // Why? ID is not required.
        //runFile << "[ID: " << rand() % 1000 << "]\n";
        runFile << "begin trees;\n";
        runFile.close();

        int tenPercent=(int)floor((double)tmip.numTreesReconstruct/10);
        int multIndStartTime = time(0);
        int multIndInnerStartTime = time(0);
        for (int l=0;l<tmip.numTreesReconstruct;l++){
            if (DEBUG_OUTPUT >= 2){
                if (((l+1) % tenPercent) == 0){
                    cout << "   " << setw(3) << 10*(l+1)/tenPercent << "%";
                    int tempEndTime = time(0);
                    cout << "   " << tempEndTime - multIndInnerStartTime << " seconds." << endl;
                    multIndInnerStartTime = time(0);
                }
            }
            // Go through each species in tmip.speciesGroupings and pick one individual randomly
            set <unsigned> someTaxa;
            someTaxa.clear();
            for (list <set <unsigned> >::iterator lsit=tmip.speciesGroupings.begin();lsit!=tmip.speciesGroupings.end();lsit++){
                int randIndex = rand() % (*lsit).size();
                int setIndex = 0;
                for (set <unsigned>::iterator sit=(*lsit).begin();sit!=(*lsit).end();sit++){
                    if (setIndex == randIndex){
                        someTaxa.insert(*sit);
                    }
                    setIndex++;
                }
            }
            Alignments jackAlignments = (*lait)->getTaxaSubset(someTaxa);
            jackAlignments.setOutputFormat(OF_PHYLIP);
            tempMultIndFile << jackAlignments;
        }
        tempMultIndFile.close();
        //Now run program to generate tree and append to tempCommandStr
        if (DEBUG_OUTPUT >= 0){
            cout << "Running: \"" << tempCommandStr << "\"" << endl;
        }
        // Perhaps with MAC system, can not append to this file.
        system(tempCommandStr.c_str());

        int multIndEndTime = time(0);
        if (DEBUG_OUTPUT >= 0){
            cout << "Total seconds: " << multIndEndTime - multIndStartTime << endl;
        }
        
        // Place the end at the end of the filename
        // MAC is having problems here. HPC also having problems here.
        runFile.open((*lsit).c_str(),fstream::app);
        if (!runFile.good()){
            cout << "Can not open file getTreesMultInd \"" << *lsit << "\"" << endl;
            cout << "Trying to read state." << endl;
            if ((runFile.rdstate() & ifstream::eofbit) != 0){
                cout << "       eofbit set" << endl;
            }
            if ((runFile.rdstate() & ifstream::failbit) != 0){
                cout << "       failbit set" << endl;
            }
            if ((runFile.rdstate() & ifstream::badbit) != 0){
                cout << "       badbit set" << endl;
            }
            //cout << "Ommitting \"end;\" from file" << endl;
            //runFile.close();
            //cout << "Trying to open for normal output. " << endl;
            //runFile.open((*lsit).c_str(),fstream::out);
            //cout << "Seeking to end of file." << endl;
            //runFile.seekp(ios_base::end);
            //if (!runFile.good()){
            //    cout << "Not Good, exiting." << endl;
            //}
            //else {
            //    cout << "Good, continuing!" << endl;
            //}
            //exit(0);
            cout << "Using system command to append \"end;\"." << endl;
            string tmptmpcommand = "echo 'end;' >> ";
            tmptmpcommand.append(*lsit);
            system(tmptmpcommand.c_str());
        }
        else {
            runFile << "end;" << endl;
            runFile.close ();
        }
        // Lets delete the tempfile, so we don't unknownling use it again.
        string delTempJackknifeFile = "rm -f ";
        delTempJackknifeFile.append(tempMultIndFileName);
        system(delTempJackknifeFile.c_str()); 
        lait++;
        lsit++;
    }
}

void doPCA(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, list <svm_node *> *projectionMatrix, svm_node *empiricalMean)
{
    if (DEBUG_OUTPUT >= 2){
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        cout << "Matrix before PCA projection:" << endl;
        for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
            int i=0;
            while ((*svmit)[i].index != -1){
                cout << setw(8) << (*svmit)[i].value << " ";
                i++;
            }
            cout << ";" << endl;
        }
        for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
            int i=0;
            while ((*svmit)[i].index != -1){
                cout << setw(8) << (*svmit)[i].value << " ";
                i++;
            }
            cout << ";" << endl;
        }
    }

    //cout << "Size: " << projectionMatrix->size() << endl;
    int empiricalMeanLength=0;

    // Trust that the first vector in groupOne is the same length as all vectors
    // in groupOne and groupTwo.
    while ((*(groupOne.begin()))[empiricalMeanLength].index != -1){
        empiricalMeanLength++;
    }
    empiricalMeanLength++;

    if (DEBUG_OUTPUT >= 1){
        cout << "Empirical mean: ";
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for (int i=0;i<empiricalMeanLength-1;i++){
            cout << setw(8) << empiricalMean[i].value << " ";
        }
        cout << endl;
    }

    // Now subtract from groupOne and groupTwo
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value -= empiricalMean[i].value;
            }
        }
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value -= empiricalMean[i].value;
            }
        }
    }
    // Done subtracting empirical mean

    list <svm_node *> *combined = new list <svm_node *>;
    svm_node *new_svm_node;
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        int i=0;
        while ((*svmit)[i].index != -1){
            i++;
        }
        new_svm_node = new svm_node[i+1];
        for (int j=0;j<=i;j++){
            new_svm_node[j].index = (*svmit)[j].index;
            new_svm_node[j].value = (*svmit)[j].value;
        }
        combined->push_back(new_svm_node); 
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        int i=0;
        while ((*svmit)[i].index != -1){
            i++;
        }
        new_svm_node = new svm_node[i+1];
        for (int j=0;j<=i;j++){
            new_svm_node[j].index = (*svmit)[j].index;
            new_svm_node[j].value = (*svmit)[j].value;
        }
        combined->push_back(new_svm_node); 
    }
    int GroupOneSize = groupOne.size();
    int GroupTwoSize = groupTwo.size();
    // Need to actually delete the svm nodes
    for (list <svm_node *>::iterator lsvmit=groupOne.begin();lsvmit!=groupOne.end();lsvmit++){
        delete [] (*lsvmit);
    }
    for (list <svm_node *>::iterator lsvmit=groupTwo.begin();lsvmit!=groupTwo.end();lsvmit++){
        delete [] (*lsvmit);
    }
    groupOne.clear();
    groupTwo.clear();

    if (DEBUG_OUTPUT >= 1){
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        cout << "Projection matrix:" << endl;
        for (list <svm_node *>::iterator svmit=projectionMatrix->begin();svmit != projectionMatrix->end();svmit++){
            int i=0;
            while ((*svmit)[i].index != -1){
                cout << setw(8) << (*svmit)[i].value << " ";
                i++;
            }
            cout << ";"  << endl;
        }
    }

    int newColCount=0;
    while( (*(projectionMatrix->begin()))[newColCount].index != -1){
        newColCount++;
    }

    // Multiply combined * projectionMatrix and place in groupOne and groupTwo
    list <svm_node *>::iterator combit=combined->begin();
    for (int j=0;j<GroupOneSize;j++){
        // New row of groupOne
        svm_node *new_svm_node = new svm_node[newColCount+1];
        for (int k=0;k<newColCount;k++){
            new_svm_node[k].index=k;
            new_svm_node[k].value=0;
            int l=0;
            for (list <svm_node *>::iterator svmit=projectionMatrix->begin();svmit!=projectionMatrix->end();svmit++){
                new_svm_node[k].value += (*combit)[l].value*(*svmit)[k].value;
                l++;
            }
        }
        new_svm_node[newColCount].index = -1;
        groupOne.push_back(new_svm_node);
        combit++;
    }

    for (int j=0;j<GroupTwoSize;j++){
        // New row of groupOne
        svm_node *new_svm_node = new svm_node[newColCount+1];
        for (int k=0;k<newColCount;k++){
            new_svm_node[k].index = k;
            new_svm_node[k].value=0;
            int l=0;
            for (list <svm_node *>::iterator svmit=projectionMatrix->begin();svmit!=projectionMatrix->end();svmit++){
                new_svm_node[k].value += (*combit)[l].value*(*svmit)[k].value;
                l++;
            }
        }
        new_svm_node[newColCount].index = -1;
        groupTwo.push_back(new_svm_node);
        combit++;
    }

    if (DEBUG_OUTPUT >= 2){
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        cout << "Matrix after PCA projection:" << endl;
        for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
            int i=0;
            while ((*svmit)[i].index != -1){
                cout << setw(8) << (*svmit)[i].value << " ";
                i++;
            }
            cout << ";"  << endl;
        }
        for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
            int i=0;
            while ((*svmit)[i].index != -1){
                cout << setw(8) << (*svmit)[i].value << " ";
                i++;
            }
            cout << ";"  << endl;
        }
    }
    // Delete things we created!
    for (list <svm_node *>::iterator svmit=combined->begin();svmit != combined->end();svmit++){
        delete *svmit;
    }
    delete combined;
}

void getPCA_data(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, double projectSVDcutOff, list <svm_node *> **projectionMatrix, svm_node **empiricalMean)
{
    if (projectSVDcutOff != -1)
    {
        getPCA_data(groupOne, groupTwo, projectSVDcutOff, projectionMatrix, empiricalMean);
    }
    else {
        cout << "getPCA_data called with projectSVDcutOff == -1" << endl;
        exit (0);
    }
}

void getPCA_data(list <svm_node * > &groupOne, list <svm_node * > &groupTwo, double projectSVDcutOff, list <svm_node *> **projectionMatrix, svm_node **empiricalMean, int numSingValues)
{
    int empiricalMeanLength=0;

    // Subtract the mean to make zero empirical mean.
    if (DEBUG_OUTPUT >= 0){
        cout << "Calculating mean." << endl;
    }
    // Trust that the first vector in groupOne is the same length as all vectors
    // in groupOne and groupTwo.
    while ((*(groupOne.begin()))[empiricalMeanLength].index != -1){
        empiricalMeanLength++;
    }
    empiricalMeanLength++;
    svm_node *tempempiricalMean = new svm_node[empiricalMeanLength];
    for (int i=0;i<empiricalMeanLength-1;i++){
        tempempiricalMean[i].value = 0;
        tempempiricalMean[i].index = i+1;

    }
    tempempiricalMean[empiricalMeanLength-1].value = 0;
    tempempiricalMean[empiricalMeanLength-1].index = -1;

    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                tempempiricalMean[i].value += (*svmit)[i].value;
            }
        }
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                tempempiricalMean[i].value += (*svmit)[i].value;
            }
        }
    }
    // Now divide by the total number of vectors to get the mean.
    for (int i=0;i<empiricalMeanLength-1;i++){
        tempempiricalMean[i].value /= (groupOne.size() + groupTwo.size());
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "Empirical mean: ";
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for (int i=0;i<empiricalMeanLength-1;i++){
            cout << setw(8) << tempempiricalMean[i].value << " ";
        }
        cout << endl;
    }
    
    // PROBLEM HERE
    //list <svm_node * > groupOne = groupOne; 
    //list <svm_node * > groupTwo = groupTwo;
    // Now subtract from groupOne and groupTwo
    // WE HAVE TO UNDO THIS!!!!
    // TODO incorperate this into combined calculation so that we don't have to undo
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value -= tempempiricalMean[i].value;
            }
        }
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value -= tempempiricalMean[i].value;
            }
        }
    }
    *empiricalMean = tempempiricalMean;
    // Done subtracting empirical mean

    list <svm_node *> *combined = new list <svm_node *>;
    svm_node *new_svm_node;
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        int i=0;
        while ((*svmit)[i].index != -1){
            i++;
        }
        new_svm_node = new svm_node[i+1];
        for (int j=0;j<=i;j++){
            new_svm_node[j].index = (*svmit)[j].index;
            new_svm_node[j].value = (*svmit)[j].value;
        }
        combined->push_back(new_svm_node); 
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        int i=0;
        while ((*svmit)[i].index != -1){
            i++;
        }
        new_svm_node = new svm_node[i+1];
        for (int j=0;j<=i;j++){
            new_svm_node[j].index = (*svmit)[j].index;
            new_svm_node[j].value = (*svmit)[j].value;
        }
        combined->push_back(new_svm_node); 
    }

    list <svm_node *> *tempprojectionMatrix = project_SVM_nodesMatrix(*combined,projectSVDcutOff,numSingValues);
    *projectionMatrix = tempprojectionMatrix;
    // UNDO THE SUBTRACTION OF THE MEAN
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value += tempempiricalMean[i].value;
            }
        }
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                (*svmit)[i].value += tempempiricalMean[i].value;
            }
        }
    }
    // Delete things we created!
    for (list <svm_node *>::iterator svmit=combined->begin();svmit != combined->end();svmit++){
        delete *svmit;
    }
    delete combined;
}

svm_problem form_svm_problem(list <svm_node * > &groupOne, list <svm_node * > &groupTwo)
{
    svm_problem SVM_returnProblem;
    svm_node **trainingData;
    trainingData = new svm_node*[groupOne.size() + groupTwo.size()];
    int trainingDataCount=0;
    double *y = new double[groupOne.size() + groupTwo.size()];
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        trainingData[trainingDataCount] = *svmit;
        y[trainingDataCount] = 1; // Put these points in class 1.
        trainingDataCount++;
    }
    for (list <svm_node *>::iterator svmit=groupTwo.begin();svmit != groupTwo.end();svmit++){
        trainingData[trainingDataCount] = *svmit;
        y[trainingDataCount] = 2; // Put these points in class 2.
        trainingDataCount++;
    }
    SVM_returnProblem.l = groupOne.size() + groupTwo.size();
    SVM_returnProblem.y = y;
    SVM_returnProblem.x = trainingData; 
    return SVM_returnProblem;
}

svm_parameter getLinear_svm_parameter()
{
    svm_parameter return_svm_parameter;
    return_svm_parameter.svm_type = C_SVC;
    return_svm_parameter.kernel_type = LINEAR;
    return_svm_parameter.cache_size = 100; //default
    return_svm_parameter.C = 1; //default
    return_svm_parameter.eps = 0.001; //default
    return_svm_parameter.nr_weight = 0; // Do not change weight
    return_svm_parameter.shrinking = 1; //default
    return_svm_parameter.probability = 0; //default
    return return_svm_parameter;
}

void delete_svm_problem(svm_problem &prob)
{
    delete [] prob.x;
    delete [] prob.y;
}

void get_svm_predictions(svm_problem &myProblem, svm_model &myModel, int &groupOnePredictCorrectCount, int &groupTwoPredictCorrectCount)
{
    for (int i=0;i<myProblem.l;i++){
        if (DEBUG_OUTPUT >= 2)
        {
            int j = 0;
            cout << "i=" << i << "  ";
            while(myProblem.x[i][j].index != -1){
                //cout << "Index: " << trainingData[i][j].index << " Value: ";
                //cout << setw(5) << trainingData[i][j].value << " ";
                cout << setw(3) << j+1 << ":" << setw(5) << myProblem.x[i][j].value << " ";
                //cout << trainingData[i][j].value << " ";
                j++;
            }
            cout << ";" <<  endl;
        }
        int predictValue = (int)svm_predict(&myModel,myProblem.x[i]);
        if (myProblem.y[i] == 1 && myProblem.y[i] == floor((double)predictValue)){
            groupOnePredictCorrectCount++;
        }
        if (myProblem.y[i] == 2 && myProblem.y[i] == floor((double)predictValue)){
            groupTwoPredictCorrectCount++;
        }

        if (DEBUG_OUTPUT >= 2){
            cout << setw(4) << "i=" << i << "   svm_predict: " << predictValue << endl;
        }
    }
}

void calcSVMseparationMrBayes(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, MrBayesParameters &tmbp, MrBayesResults &MB_results, SampleParameters &tsvmp, SVM_separationResults &results)
{
    // Write out these alignments and then read them into mr bayes

    list <string> newFileNames;
    string tempAlignmentsBaseName = tempPrefix;
    tempAlignmentsBaseName.append("tempalignments_ID");
    char tmpstr[80];
    sprintf (tmpstr,"%d%d%d%d",(int)getpid(),rand(),rand(),rand());
    tempAlignmentsBaseName.append(tmpstr);

    int tempFileCount = 0;
    for (list <Alignments *>::iterator lait = AlignmentsOne.begin();lait != AlignmentsOne.end();lait++){
        string tempString = tempAlignmentsBaseName;
        sprintf(tmpstr,"%d",tempFileCount);
        tempString.append("_file");
        tempString.append(tmpstr);
        tempString.append(".nex");
        
        fstream tempFile;
        tempFile.open(tempString.c_str(),fstream::out);
        if (tempFile.bad()){
            cout << "Could not open " << tempString << endl;
            exit(0);
        }
        (*lait)->setOutputFormat(OF_NEXUS);
        tempFile << *(*lait);
        tempFile.close();
        newFileNames.push_back(tempString);
        tempFileCount++;
    }
    for (list <Alignments *>::iterator lait = AlignmentsTwo.begin();lait != AlignmentsTwo.end();lait++){
        string tempString = tempAlignmentsBaseName;
        sprintf(tmpstr,"%d",tempFileCount);
        tempString.append("_file");
        tempString.append(tmpstr);
        tempString.append(".nex");
        
        fstream tempFile;
        tempFile.open(tempString.c_str(),fstream::out);
        if (tempFile.bad()){
            cout << "Could not open " << tempString << endl;
            exit(0);
        }
        (*lait)->setOutputFormat(OF_NEXUS);
        tempFile << *(*lait);
        tempFile.close();
        newFileNames.push_back(tempString);
        tempFileCount++;
    }

    MrBayesParameters new_tmbp = tmbp;
    new_tmbp.inputNexFiles = newFileNames;

    calcSVMseparationMrBayes(new_tmbp,MB_results,tsvmp,results,AlignmentsOne.size());
    // We should delete these temporary alignments
    for (list <string>::iterator lsit = newFileNames.begin();lsit!=newFileNames.end();lsit++){
        string delCommand = "rm -f ";
        delCommand.append(*lsit);
        system(delCommand.c_str());
    }
}

void calcSVMseparationMrBayes(MrBayesParameters &tmbp, MrBayesResults &MB_results, SampleParameters &tsvmp, SVM_separationResults &results, int numGroupOne)
{
    list <list <string> > TreeListListFileNames; // For getTreeMrBayes

    getTreesMrBayes(tmbp,TreeListListFileNames,MB_results);

    SampleParameters new_tsvmp = tsvmp;

    new_tsvmp.numTreesPerFile = (tmbp.MBP_ngen/tmbp.MBP_sampleFreq);

    list <string> treeFileNamesGroupOne;
    list <string> treeFileNamesGroupTwo;

    // Default is to compare first input filename with the rest
    // Lets use numGroupOne
    int fileGroupOneCount = 0;
    list <list <string> >::iterator llsit=TreeListListFileNames.begin();
    while (fileGroupOneCount < numGroupOne){
        list <string>::iterator lsit = (*llsit).begin();
        while(lsit != (*llsit).end()){
            treeFileNamesGroupOne.push_back(*lsit);
            lsit++;
        }
        llsit++;
        //treeFileNamesGroupOne = *llsit;
        //llsit++;
        fileGroupOneCount++;
    }
    // The rest go in treeFileNamesGroupTwo
    while (llsit != TreeListListFileNames.end()){
        list <string>::iterator lsit = (*llsit).begin();
        while(lsit != (*llsit).end()){
            treeFileNamesGroupTwo.push_back(*lsit);
            lsit++;
        }
        llsit++;
    }

    calcSVMseparation(treeFileNamesGroupOne,treeFileNamesGroupTwo,new_tsvmp,results);

    // We should clean up the tree files here.
    string delOutputFiles = "rm -f ";
    for (list <string>::iterator lsit=treeFileNamesGroupOne.begin();lsit!= treeFileNamesGroupOne.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    for (list <string>::iterator lsit=treeFileNamesGroupTwo.begin();lsit!= treeFileNamesGroupTwo.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    //cout << "Tree deletion command: " << endl << delOutputFiles << endl;
    system(delOutputFiles.c_str());
}

void calcSVMseparationJackknife(list <string> inputFileNames, int numGroupOne, JackknifeParameters &tjp, SampleParameters &tsvmp, SVM_separationResults &results)
{
    // Read in alignments then call calcSVMseparationJackknife using alignments
    list <Alignments *> AlignmentsOne, AlignmentsTwo;

    int fileGroupOneCount = 0;
    Alignments *tempAlignments; 
    for (list <string>::iterator lsit = inputFileNames.begin();lsit!= inputFileNames.end();lsit++){
        tempAlignments = new Alignments(*lsit);
        if (fileGroupOneCount < numGroupOne){
            AlignmentsOne.push_back(tempAlignments);
        }
        else {
            AlignmentsTwo.push_back(tempAlignments);
        }
        fileGroupOneCount++;
    }
    calcSVMseparationJackknife(AlignmentsOne,AlignmentsTwo,tjp,tsvmp,results);
}

void calcSVMseparationJackknife(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, JackknifeParameters &tjp, SampleParameters &tsp, SVM_separationResults &results)
{
    list <Alignments *> allAlignments;
    for (list <Alignments *>::iterator lait = AlignmentsOne.begin();lait!=AlignmentsOne.end();lait++){
        allAlignments.push_back(*lait);
    }
    for (list <Alignments *>::iterator lait = AlignmentsTwo.begin();lait!=AlignmentsTwo.end();lait++){
        allAlignments.push_back(*lait);
    }
    // First get the trees. Pick outputFileNames first
    list <string> outputFileNames; 
    list <string> outputFileNameGroupOne; 
    list <string> outputFileNameGroupTwo; 
    list <string> outputFileNameGroupOneResample; 
    list <string> outputFileNameGroupTwoResample; 
    string tempFileName;
    int ii=0;
    char randIdString[80];
    sprintf(randIdString,"%d%d%d",rand() % 100000,rand() % 100000,rand() % 100000);
    
    //// When are these files going to be deleted? We need to delete them at the end of this function
    //for (list <Alignments *>::iterator lait = allAlignments.begin();lait!=allAlignments.end();lait++){
    //    string newFileName = tempPrefix;
    //    newFileName.append("temp_calcSVMseparationJackknife_ID");
    //    newFileName.append(randIdString);
    //    newFileName.append("_");
    //    char tempChar[80];
    //    sprintf(tempChar,"%d.jk.t",ii);
    //    newFileName.append(tempChar);
    //    outputFileNames.push_back(newFileName);
    //    ii++;
    //}
    
    // Why should we calculate jackkknifeCount number of reconstructions for each alignment?
    // For efficiency we should calculate jackknifeCount for each group?  
    // WARNING: This might not give us a uniform sampling of the distribution of each group!
    // On the other hand, by taking discrete samples, what does uniform mean? With respect
    // to the continuous distribution?
    // We can essentially pass to getTrees the number of trees to sample from each distribution


    //getTreesJackknife(tjp,allAlignments,outputFileNames);

    //list <string> treeFileNamesGroupOne;
    //list <string> treeFileNamesGroupTwo;

    // Place the output filenames into their groups according
    // to the sizes of AlignmentsOne and AlignmentsTwo
    //list <string>::const_iterator lsit = outputFileNames.begin();
    //int insertOutputFileNamesCount = 0;
    //while (insertOutputFileNamesCount < AlignmentsOne.size()){
    //    treeFileNamesGroupOne.push_back(*lsit);
    //    lsit++;
    //    insertOutputFileNamesCount++;
    //}
    //while (insertOutputFileNamesCount < AlignmentsOne.size() + AlignmentsTwo.size())
    //{
    //    treeFileNamesGroupTwo.push_back(*lsit);
    //    lsit++;
    //    insertOutputFileNamesCount++;
    //}

    getTreesJackknife(tjp,tsp,AlignmentsOne,tempFileName);
    outputFileNameGroupOne.push_back(tempFileName);
    getTreesJackknife(tjp,tsp,AlignmentsTwo,tempFileName);
    outputFileNameGroupTwo.push_back(tempFileName);
    SampleParameters    new_tsp = tsp;
    new_tsp.SVM_sampleSize = tsp.SVM_resampleSize;
    getTreesJackknife(tjp,new_tsp,AlignmentsOne,tempFileName);
    outputFileNameGroupOneResample.push_back(tempFileName);
    getTreesJackknife(tjp,new_tsp,AlignmentsTwo,tempFileName);
    outputFileNameGroupTwoResample.push_back(tempFileName); 

    tsp.numTreesPerFile = -1; // Signifies to sample all trees

    calcSVMseparation(outputFileNameGroupOne,outputFileNameGroupTwo,outputFileNameGroupOneResample,outputFileNameGroupTwoResample,tsp,results);

    // Now we need to delete the output file names
    string delOutputFiles = "rm -f ";
    for (list <string>::iterator lsit=outputFileNameGroupOne.begin();lsit!= outputFileNameGroupOne.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    for (list <string>::iterator lsit=outputFileNameGroupTwo.begin();lsit!= outputFileNameGroupTwo.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    for (list <string>::iterator lsit=outputFileNameGroupOneResample.begin();lsit!= outputFileNameGroupOneResample.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    for (list <string>::iterator lsit=outputFileNameGroupTwoResample.begin();lsit!= outputFileNameGroupTwoResample.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    system(delOutputFiles.c_str());
}

void calcSVMseparationMultInd(list <Alignments *> AlignmentsOne, list <Alignments *> AlignmentsTwo, MultIndParameters &tmip, SampleParameters &tsp, SVM_separationResults &results)
{
    list <Alignments *> allAlignments;
    for (list <Alignments *>::iterator lait = AlignmentsOne.begin();lait!=AlignmentsOne.end();lait++){
        allAlignments.push_back(*lait);
    }
    for (list <Alignments *>::iterator lait = AlignmentsTwo.begin();lait!=AlignmentsTwo.end();lait++){
        allAlignments.push_back(*lait);
    }
    // First get the trees. Pick outputFileNames first
    list <string> outputFileNames; 
    int ii=0;
    char randIdString[80];
    sprintf(randIdString,"%d%d%d%d",(int)getpid(),rand() % 100000,rand() % 100000,rand() % 100000);
    
    // When are these files going to be deleted? We need to delete them at the end of this function
    for (list <Alignments *>::iterator lait = allAlignments.begin();lait!=allAlignments.end();lait++){
        string newFileName = tempPrefix;
        newFileName.append("temp_calcSVMseparationMultInd_ID");
        newFileName.append(randIdString);
        newFileName.append("_");
        char tempChar[80];
        sprintf(tempChar,"%d.run1.t",ii);
        newFileName.append(tempChar);
        outputFileNames.push_back(newFileName);
        ii++;
    }
    
    getTreesMultInd(tmip,allAlignments,outputFileNames);

    list <string> treeFileNamesGroupOne;
    list <string> treeFileNamesGroupTwo;

    // Place the output filenames into their groups according
    // to the sizes of AlignmentsOne and AlignmentsTwo
    list <string>::const_iterator lsit = outputFileNames.begin();
    int insertOutputFileNamesCount = 0;
    while (insertOutputFileNamesCount < AlignmentsOne.size()){
        treeFileNamesGroupOne.push_back(*lsit);
        lsit++;
        insertOutputFileNamesCount++;
    }
    while (insertOutputFileNamesCount < AlignmentsOne.size() + AlignmentsTwo.size())
    {
        treeFileNamesGroupTwo.push_back(*lsit);
        lsit++;
        insertOutputFileNamesCount++;
    }


    tsp.numTreesPerFile = tmip.numTreesReconstruct;

    calcSVMseparation(treeFileNamesGroupOne,treeFileNamesGroupTwo,tsp,results);

    // Now we need to delete the output file names
    string delOutputFiles = "rm -f ";
    for (list <string>::iterator lsit=outputFileNames.begin();lsit!= outputFileNames.end();lsit++){
        delOutputFiles.append(*lsit);
        delOutputFiles.append(" ");
    }
    system(delOutputFiles.c_str());

}

void calcSVMseparationMultInd(list <string> inputFileNames, MultIndParameters &tmip, SampleParameters &tsvmp, SVM_separationResults &results, int numGroupOne)
{
    // Read in alignments then call calcSVMseparationJackknife using alignments
    list <Alignments *> AlignmentsOne, AlignmentsTwo;

    int fileGroupOneCount = 0;
    Alignments *tempAlignments; 
    for (list <string>::iterator lsit = inputFileNames.begin();lsit!= inputFileNames.end();lsit++){
        tempAlignments = new Alignments(*lsit);
        if (fileGroupOneCount < numGroupOne){
            AlignmentsOne.push_back(tempAlignments);
        }
        else {
            AlignmentsTwo.push_back(tempAlignments);
        }
        fileGroupOneCount++;
    }
    calcSVMseparationMultInd(AlignmentsOne,AlignmentsTwo,tmip,tsvmp,results);
}

void calcSVMseparation(list <string> &treeFileNamesGroupOne, list <string> &treeFileNamesGroupTwo, SampleParameters &tsp, SVM_separationResults &results)
{
    calcSVMseparation(treeFileNamesGroupOne, treeFileNamesGroupTwo, treeFileNamesGroupOne,treeFileNamesGroupTwo,tsp,results);
}

void calcSVMseparation(list <string> &treeFileNamesGroupOne, list <string> &treeFileNamesGroupTwo,list <string> &treeFileNamesGroupOneResample, list <string> &treeFileNamesGroupTwoResample, SampleParameters &tsp, SVM_separationResults &results)
{
    time_t startTime = time(0);
    list <svm_node *> groupOne;
    list <svm_node *> groupTwo;
    // If SVM_sampleSize == -1 or SVM_resampleSize == -1, then autocalculate based off one of the trees
    if (tsp.SVM_sampleSize == -1 || tsp.SVM_resampleSize == -1) {
        int numLeaves = countNumLeaves(getTree(*(treeFileNamesGroupOne.begin()),0));
        tsp.SVM_sampleSize = 3*numLeaves*(numLeaves-1);
        tsp.SVM_resampleSize = 6*numLeaves*(numLeaves-1);
    }

    // Now all the trees are saved in outputFileNames
    Sample_SVM(treeFileNamesGroupOne,treeFileNamesGroupTwo,tsp,groupOne,groupTwo);

    // Here we do difference of means if we are supposed to.
    // Do this before PCA?? Ask Rudy/Peter
    if (tsp.doDiffMeans == 1) {
        results.varianceGroupOne = SampleVariance(groupOne);
        results.traceVarianceGroupOne = results.varianceGroupOne.trace();
        results.detVarianceGroupOne = results.varianceGroupOne.det();
        results.logDetVarianceGroupOne = results.varianceGroupOne.log_abs_det();

        results.varianceGroupTwo = SampleVariance(groupTwo);
        results.traceVarianceGroupTwo = results.varianceGroupTwo.trace();
        results.detVarianceGroupTwo = results.varianceGroupTwo.det();
        results.logDetVarianceGroupTwo = results.varianceGroupTwo.log_abs_det();


        results.diffMeans = DifferenceOfMeansMatrix(groupOne,groupTwo);
        results.diffMeans_l2_norm = results.diffMeans.twoNorm();

        if (DEBUG_OUTPUT >= 1) {
            cout << "Variance of groupOne: ";
            cout << "Trace: " << results.traceVarianceGroupOne << "  ";
            cout << "Det: " << results.detVarianceGroupOne << "  ";
            cout << "Log Abs Det: " << results.logDetVarianceGroupOne << endl;
        }

        if (DEBUG_OUTPUT >= 1) {
            cout << "Variance of groupTwo: ";
            cout << "Trace: " << results.traceVarianceGroupTwo << "  ";
            cout << "Det: " << results.detVarianceGroupTwo << "  ";
            cout << "Log Abs Det: " << results.logDetVarianceGroupTwo << endl;
        }

        if (DEBUG_OUTPUT >= 1) {
            cout << "Difference of means: ";
            cout << results.diffMeans;
            cout << "Difference of means l2 norm: ";
            cout << results.diffMeans_l2_norm << endl;;
        }
    }


    list <svm_node *> *projectionMatrix; 
    svm_node    *empiricalMean;
    results.numSingularValuesTaken = -1; // Set to -1 to signify PCA not done. Might change below.
    if (tsp.projectSVD == 1){
        if (DEBUG_OUTPUT >= 0){
            cout << " ***** Projecting data using PCA ***** " << endl;
        }
        if (DEBUG_OUTPUT >= 1){
            cout << " groupOne.size() = " << groupOne.size() << "   groupTwo.size() = " << groupTwo.size() << endl;
        }
        getPCA_data(groupOne,groupTwo,tsp.projectSVDcutOff,&projectionMatrix,&empiricalMean);
        if (DEBUG_OUTPUT >= 1){
            cout << " groupOne.size() = " << groupOne.size() << "   groupTwo.size() = " << groupTwo.size() << endl;
        }
        // projectionMatrix and empiricalMean need to be deleted!
        doPCA(groupOne,groupTwo,projectionMatrix,empiricalMean);
        if (DEBUG_OUTPUT >= 1){
            cout << " groupOne.size() = " << groupOne.size() << "   groupTwo.size() = " << groupTwo.size() << endl;
        }
        svm_node *tempNode = *(groupOne.begin());
        results.numSingularValuesTaken = 0;
        while (tempNode[results.numSingularValuesTaken].index != -1){
            results.numSingularValuesTaken++;
        }


        // Here we do difference of means if we are supposed to.
        // Do this before PCA?? Ask Rudy/Peter
        if (tsp.doDiffMeans == 1) {
            if (DEBUG_OUTPUT >= 1) {
                cout << " **** PCA ****" << endl;
                cout << "Variance of groupOne: ";
            }
            results.varianceGroupOnePCA = SampleVariance(groupOne);
            results.traceVarianceGroupOnePCA = results.varianceGroupOne.trace();
            results.detVarianceGroupOnePCA = results.varianceGroupOne.det();
            results.logDetVarianceGroupOnePCA = results.varianceGroupOne.log_abs_det();

            results.varianceGroupTwoPCA = SampleVariance(groupTwo);
            results.traceVarianceGroupTwoPCA = results.varianceGroupTwo.trace();
            results.detVarianceGroupTwoPCA = results.varianceGroupTwo.det();
            results.logDetVarianceGroupTwoPCA = results.varianceGroupTwo.log_abs_det();


            results.diffMeansPCA = DifferenceOfMeansMatrix(groupOne,groupTwo);
            results.diffMeans_l2_normPCA = results.diffMeans.twoNorm();

            if (DEBUG_OUTPUT >= 1) {
                cout << "Variance of groupOne: ";
                cout << "Trace: " << results.traceVarianceGroupOnePCA << "  ";
                cout << "Det: " << results.detVarianceGroupOnePCA << "  ";
                cout << "Log Abs Det: " << results.logDetVarianceGroupOnePCA << endl;
            }

            if (DEBUG_OUTPUT >= 1) {
                cout << "Variance of groupTwo: ";
                cout << "Trace: " << results.traceVarianceGroupTwoPCA << "  ";
                cout << "Det: " << results.detVarianceGroupTwoPCA << "  ";
                cout << "Log Abs Det: " << results.logDetVarianceGroupTwoPCA << endl;
            }

            if (DEBUG_OUTPUT >= 0) {
                cout << "Difference of means: ";
                cout << results.diffMeansPCA;
                cout << "Difference of means l2 norm: ";
                cout << results.diffMeans_l2_normPCA << endl;;
            }
        }
    }

    // Combine data set one and Two
    svm_problem SVM_problem = form_svm_problem(groupOne,groupTwo);
    svm_parameter firstParameter;
    if (tsp.modelType == LINEAR){
        firstParameter = getLinear_svm_parameter();
    }
    else{
        cout << "Unknown model type for svm" << endl;
        exit(0);
    }
    svm_check_parameter(&SVM_problem, &firstParameter);
    svm_model *myFirstModel;

    if (DEBUG_OUTPUT >= 0){
        cout << "Calling svm_train." << endl;
    }

    int svmTrainTimeStart = time(0);
    myFirstModel = svm_train(&SVM_problem,&firstParameter);
    int svmTrainTimeEnd = time(0);

    if (DEBUG_OUTPUT >= 0){
        cout << "   " << svmTrainTimeEnd-svmTrainTimeStart << " seconds (wall clock)." << endl;
        cout << endl << "Now running svm_predict on all svm sample points." << endl;
    }

    int svmPredictStartTime = time(0);
    int groupOnePredictCorrectCount = 0;
    int groupTwoPredictCorrectCount = 0;

    get_svm_predictions(SVM_problem, *myFirstModel, groupOnePredictCorrectCount,groupTwoPredictCorrectCount);
    time_t svmPredictEndTime;
    if (DEBUG_OUTPUT >= 0){
        cout << "   groupOnePredictCorrectCount: " << groupOnePredictCorrectCount << " out of " << tsp.SVM_sampleSize << endl;
        cout << "   groupTwoPredictCorrectCount: " << groupTwoPredictCorrectCount << " out of " << tsp.SVM_sampleSize  << endl;
        cout << "       " << groupOnePredictCorrectCount + groupTwoPredictCorrectCount << " out of " << 2*tsp.SVM_sampleSize << endl;
        cout << "       Percent: " << 100*(double)((double)groupOnePredictCorrectCount + (double)groupTwoPredictCorrectCount)/((double)groupOne.size() + (double)groupTwo.size()) << "% correct" << endl;;
        svmPredictEndTime = time(0);
        cout << svmPredictEndTime - svmPredictStartTime << " seconds (wall clock)." << endl;
    }
    results.groupOneCount = groupOnePredictCorrectCount;
    results.groupTwoCount = groupTwoPredictCorrectCount;
    results.separationPercentage = 100*(double)((double)groupOnePredictCorrectCount + (double)groupTwoPredictCorrectCount)/((double)groupOne.size() + (double)groupTwo.size());

    // Now do resampling
    if (DEBUG_OUTPUT >= 0){
        cout << endl << " ***** New resampling ***** " << endl;
    }

    list <svm_node * > resamplingGroupOne, resamplingGroupTwo;
    int oldSVM_sampleSize = tsp.SVM_sampleSize;
    tsp.SVM_sampleSize = tsp.SVM_resampleSize;

    // Use the resample tree files.
    Sample_SVM(treeFileNamesGroupOneResample, treeFileNamesGroupTwoResample, tsp, resamplingGroupOne, resamplingGroupTwo);
    tsp.SVM_sampleSize = oldSVM_sampleSize;
    //cout << resamplingGroupOne.size() << " " << resamplingGroupTwo.size() << endl;

    if (tsp.projectSVD == 1){
        if (DEBUG_OUTPUT >= 0){
            cout << " ***** Projecting data using previous PCA matrix ***** " << endl;
        }
        doPCA(resamplingGroupOne,resamplingGroupTwo,projectionMatrix,empiricalMean);
    }


    // Combine data set one and Two
    svm_problem SVM_resamplingProblem = form_svm_problem(resamplingGroupOne, resamplingGroupTwo);

    if (DEBUG_OUTPUT >= 0){
        cout << endl << "Now running svm_predict on all svm sample points." << endl;
    }
    svmPredictStartTime = time(0);
    int resamplingGroupOnePredictCorrectCount = 0;
    int resamplingGroupTwoPredictCorrectCount = 0;
    get_svm_predictions(SVM_resamplingProblem, *myFirstModel, resamplingGroupOnePredictCorrectCount,resamplingGroupTwoPredictCorrectCount);

    if (DEBUG_OUTPUT >= 0){
        cout << "   resamplingGroupOnePredictCorrectCount: " << resamplingGroupOnePredictCorrectCount << " out of " << tsp.SVM_resampleSize  << endl;
        cout << "   resamplingGroupTwoPredictCorrectCount: " << resamplingGroupTwoPredictCorrectCount << " out of " << tsp.SVM_resampleSize  << endl;
        cout << "       " << resamplingGroupOnePredictCorrectCount + resamplingGroupTwoPredictCorrectCount << " out of " << 2*tsp.SVM_resampleSize << endl;
        cout << "       Percent: " << 100*(double)((double)resamplingGroupOnePredictCorrectCount + (double)resamplingGroupTwoPredictCorrectCount)/((double)resamplingGroupOne.size() + (double)resamplingGroupTwo.size()) << "% correct. Resample" << endl;;
        svmPredictEndTime = time(0);
        cout << svmPredictEndTime - svmPredictStartTime << " seconds (wall clock)." << endl;
    }
    results.resampleGroupOneCount = groupOnePredictCorrectCount;
    results.resampleGroupTwoCount = groupTwoPredictCorrectCount;
    results.resampleSeparationPercentage = 100*(double)((double)groupOnePredictCorrectCount + (double)groupTwoPredictCorrectCount)/((double)groupOne.size() + (double)groupTwo.size());

    //groupOneCount = resamplingGroupOnePredictCorrectCount;
    //groupTwoCount = resamplingGroupTwoPredictCorrectCount;

    // groupOne and groupTwo and resamplingGroupOne and resamplingGroupTwo need to be deleted
    for (list <svm_node *>::iterator lsvmit=groupOne.begin();lsvmit!=groupOne.end();lsvmit++){
        delete [] (*lsvmit);
    }
    for (list <svm_node *>::iterator lsvmit=groupTwo.begin();lsvmit!=groupTwo.end();lsvmit++){
        delete [] (*lsvmit);
    }

    // resamplingGroupOne and resamplingGroupTwo and resamplingGroupOne and resamplingGroupTwo need to be deleted
    for (list <svm_node *>::iterator lsvmit=resamplingGroupOne.begin();lsvmit!=resamplingGroupOne.end();lsvmit++){
        delete [] (*lsvmit);
    }
    for (list <svm_node *>::iterator lsvmit=resamplingGroupTwo.begin();lsvmit!=resamplingGroupTwo.end();lsvmit++){
        delete [] (*lsvmit);
    }

    if (tsp.projectSVD == 1){
        // Delete projectionMatrix and empiricalMean
        for (list <svm_node *>::iterator lsvmit=projectionMatrix->begin();lsvmit!=projectionMatrix->end();lsvmit++){
            delete [] (*lsvmit);
        }
        delete [] empiricalMean;
    }

    time_t endTime = time(0);
    results.secondsToCompute = endTime-startTime;
}

void getSomeJackknifeAlignments(Alignments &origAlignments, list <Alignments *> &AlignmentsOne, list <int> &AlignmentsLengths, int colSize)
{
    AlignmentsOne.clear();

    Alignments *tempAlignments;
    list <int>::const_iterator liit = AlignmentsLengths.begin();
    for(;liit!=AlignmentsLengths.end();liit++)
    {
        tempAlignments = new Alignments;
        *tempAlignments = origAlignments.getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsOne.push_back(tempAlignments);
    }
}

void getSomeJackknifeAlignments(Alignments &origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize)
{
    AlignmentsOne.clear();
    AlignmentsTwo.clear();
    list <int>::const_iterator liit = AlignmentsLengths.begin();
    
    if (AlignmentsLengths.size() < 2){
        cout << "getSomeJackknifeAlignments called with AlignmentsLengths.size() < 2" << endl;
        exit(0);
    }

    int groupOneCount = 0;
    Alignments *tempAlignments;

    // The first groupOneSize go in AlignmentsOne
    while (groupOneCount < groupOneSize){
        tempAlignments = new Alignments;
        *tempAlignments = origAlignments.getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsOne.push_back(tempAlignments);
        liit++;
        groupOneCount++;
    }

    // The rest go in AlignmentsTwo
    while (liit != AlignmentsLengths.end()){
        tempAlignments = new Alignments;
        *tempAlignments = origAlignments.getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsTwo.push_back(tempAlignments);
        liit++;
    } 
}

void getSomeJackknifeAlignments(list <Alignments *> &origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize)
{ 
    AlignmentsOne.clear();
    AlignmentsTwo.clear();
    list <int>::const_iterator liit = AlignmentsLengths.begin();
    list <Alignments *>::iterator lait;
    
    if (AlignmentsLengths.size() < 2){
        cout << "getSomeJackknifeAlignments called with AlignmentsLengths.size() < 2" << endl;
        exit(0);
    }

    int groupOneCount = 0;
    Alignments *tempAlignments;

    // The first groupOneSize go in AlignmentsOne
    while (groupOneCount < groupOneSize){
        tempAlignments = new Alignments;
        // Get a random alignment from origAlignments.
        lait = origAlignments.begin();
        int stopNum = (rand() % origAlignments.size());
        for (int i=0;i<stopNum;i++) {
            lait++;
        }
        *tempAlignments = (*lait)->getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsOne.push_back(tempAlignments);
        liit++;
        groupOneCount++;
    }

    // The rest go in AlignmentsTwo
    while (liit != AlignmentsLengths.end()){
        tempAlignments = new Alignments;
        lait = origAlignments.begin();
        int stopNum = (rand() % origAlignments.size());
        for (int i=0;i<stopNum;i++) {
            lait++;
        }
        *tempAlignments = (*lait)->getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsTwo.push_back(tempAlignments);
        liit++;
    } 

}


// Should make option here to not allow any permutation where two groups are as the original
void getSomeJackknifeAlignmentsPermute(list <Alignments *> origAlignments, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, list <int> &AlignmentsLengths, int colSize, int groupOneSize, int allowAnyPermutation)
{ 
    AlignmentsOne.clear();
    AlignmentsTwo.clear();
    list <int>::const_iterator liit = AlignmentsLengths.begin();
    list <Alignments *>::iterator lait;
    
    if (AlignmentsLengths.size() < 2){
        cout << "getSomeJackknifeAlignments called with AlignmentsLengths.size() < 2" << endl;
        exit(0);
    }

    int groupOneCount = 0;
    Alignments *tempAlignments;

    int firstRunAlignmentsOne = 0; // This indicates if we have used an alignment from groupTwo
    // The first groupOneSize go in AlignmentsOne
    while (groupOneCount < groupOneSize){
        tempAlignments = new Alignments;
        // Get a random alignment from origAlignments.
        lait = origAlignments.begin();
        int stopNum = (rand() % origAlignments.size());
        // Make the first alignment from groupTwo at least
        while (stopNum < groupOneSize && firstRunAlignmentsOne == 0 && allowAnyPermutation != 1) {
            stopNum = (rand() % origAlignments.size());
        }
        firstRunAlignmentsOne=1;
        for (int i=0;i<stopNum;i++) {
            lait++;
        }
        *tempAlignments = (*lait)->getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsOne.push_back(tempAlignments);

        // Now that we have used alignment *lait, we should remove it from the list
        // We can safely remove them from origAlignments since it is a local variable now
        origAlignments.erase(lait);

        liit++; // Increase the AlignmentsLengths counter
        groupOneCount++;
    }

    // The rest go in AlignmentsTwo
    while (liit != AlignmentsLengths.end()){
        tempAlignments = new Alignments;
        lait = origAlignments.begin();
        int stopNum = (rand() % origAlignments.size());
        for (int i=0;i<stopNum;i++) {
            lait++;
        }
        *tempAlignments = (*lait)->getJackknife(colSize, (int)floor((double)*liit/colSize));
        AlignmentsTwo.push_back(tempAlignments);

        // Now that we have used alignment *lait, we should remove it from the list
        // We can safely remove them from origAlignments since it is a local variable now
        if (DEBUG_OUTPUT >= 1){
            cout << "origAlignments.erase(lait)" << endl;
        }
        origAlignments.erase(lait);

        liit++;
    } 

}

void userGenerateNewAlignmentsUniform(string genCommand, list <Alignments *> &AlignmentsOne, list <Alignments *> &AlignmentsTwo, int numNewAlignments, int groupOneSize)
{
    if (numNewAlignments < 2){
        cout << "in userGenerateNewAlignmentsUniform, numNewAlignments < 2" << endl;
    }
    Alignments *tempAlignments;

    char tmptmpStr[80];
    sprintf(tmptmpStr,"%d%d%d",rand(),rand(),rand());
    string tempFileName = "tempuserGenerateNewAlignmentsUniform";
    tempFileName.append(tmptmpStr);
    tempFileName.append(".nex");

    string sysCommand = genCommand;
    sysCommand.append(" ");
    sysCommand.append(tempFileName);

    string deleteCommand = "rm -f ";
    deleteCommand.append(tempFileName);

    int groupOneCount = 0;
    // The first groupOneSize go in AlignmentsOne
    while (groupOneCount < groupOneSize){
        system(sysCommand.c_str()); // Actually generate the file
        tempAlignments = new Alignments(tempFileName); // Read in the new file
        AlignmentsOne.push_back(tempAlignments);
        system(deleteCommand.c_str()); // delete the file
        groupOneCount++;
    }

    // Rest go into AlignmentsTwo
    for (int i=0;i<numNewAlignments-groupOneSize;i++){
        system(sysCommand.c_str()); // Actually generate the file

        tempAlignments = new Alignments(tempFileName); // Read in the new file
        AlignmentsTwo.push_back(tempAlignments);
         system(deleteCommand.c_str()); // delete the file
    }
}

svm_node *Mean(list <svm_node *> &groupOne)
{
    int empiricalMeanLength=0;

    // Subtract the mean to make zero empirical mean.
    if (DEBUG_OUTPUT >= 1){
        cout << "Calculating mean." << endl;
    }
    // Trust that the first vector in groupOne is the same length as all vectors
    // in groupOne.
    while ((*(groupOne.begin()))[empiricalMeanLength].index != -1){
        empiricalMeanLength++;
    }
    empiricalMeanLength++;
    if (DEBUG_OUTPUT >= 2){
        cout << "empiricalMeanLength: " << empiricalMeanLength << endl;
    }
    svm_node *tempempiricalMean = new svm_node[empiricalMeanLength];
    for (int i=0;i<empiricalMeanLength-1;i++){
        tempempiricalMean[i].value = 0;
        tempempiricalMean[i].index = i+1;

    }
    tempempiricalMean[empiricalMeanLength-1].value = -1;
    tempempiricalMean[empiricalMeanLength-1].index = -1;

    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength-1;i++){
            if ((*svmit)[i].index != -1){
                tempempiricalMean[i].value += (*svmit)[i].value;
            }
        }
    }

    // Now divide by the total number of vectors to get the mean.
    for (int i=0;i<empiricalMeanLength-1;i++){
        tempempiricalMean[i].value /= (groupOne.size());
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "Empirical mean: ";
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(8);
        for (int i=0;i<empiricalMeanLength-1;i++){
            cout << setw(8) << tempempiricalMean[i].value << " ";
        }
        cout << endl;
    }

    return tempempiricalMean;
}

double DifferenceOfMeans_l2(list <svm_node *> &groupOne, list <svm_node *> &groupTwo)
{
    svm_node *nodeOne, *nodeTwo;

    nodeOne = Mean(groupOne); 
    nodeTwo = Mean(groupTwo); 

    int i = 0;
    while(nodeOne[i].index != -1 && nodeTwo[i].index != -1)
    {
        nodeOne[i].value -= nodeTwo[i].value;
        i++;
    }
    Matrix tempM = svm_node_to_Matrix(nodeOne);
    delete [] nodeOne;
    delete [] nodeTwo;
    return tempM.twoNorm();

}

Matrix DifferenceOfMeansMatrix(list <svm_node *> &groupOne, list <svm_node *> &groupTwo)
{
    svm_node *nodeOne, *nodeTwo;

    nodeOne = Mean(groupOne); 
    nodeTwo = Mean(groupTwo); 

    int i = 0;
    while(nodeOne[i].index != -1 && nodeTwo[i].index != -1)
    {
        nodeOne[i].value -= nodeTwo[i].value;
        i++;
    }
    Matrix returnMatrix;
    returnMatrix = svm_node_to_Matrix(nodeOne);

    delete [] nodeOne;
    delete [] nodeTwo;
    return returnMatrix;
}

svm_node *DifferenceOfMeans(list <svm_node *> &groupOne, list <svm_node *> &groupTwo)
{
    svm_node *nodeOne, *nodeTwo;

    nodeOne = Mean(groupOne); 
    nodeTwo = Mean(groupTwo); 

    int i = 0;
    while(nodeOne[i].index != -1 && nodeTwo[i].index != -1)
    {
        nodeOne[i].value -= nodeTwo[i].value;
        i++;
    }

    delete [] nodeTwo;
    return nodeOne;
}

svm_node *Variance(list <svm_node *> &groupOne)
{
    
    svm_node *mean = Mean(groupOne);
    svm_node *variance;

    int empiricalMeanLength=0;

    while (mean[empiricalMeanLength].index != -1){
        empiricalMeanLength++;
    }

    variance = new svm_node[empiricalMeanLength+1];
    for (int i=0;i<empiricalMeanLength;i++){
        variance[i].value = 0;
        variance[i].index = i+1;
    }
    variance[empiricalMeanLength].value = 0;
    variance[empiricalMeanLength].index = -1;

    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        for (int i=0;i<empiricalMeanLength;i++){
            if ((*svmit)[i].index != -1){
                variance[i].value += ((*svmit)[i].value - mean[i].value)*((*svmit)[i].value - mean[i].value);
            }
        }
    }

    for (int i=0;i<empiricalMeanLength;i++){
        variance[i].value /= groupOne.size();
    }
    return variance;
}

Matrix SampleVariance(list <svm_node *> &groupOne)
{
    svm_node *meanGroupOne_svm_node = Mean(groupOne);

    Matrix meanGroupOne = svm_node_to_Matrix(meanGroupOne_svm_node);
    delete [] meanGroupOne_svm_node;
    if (DEBUG_OUTPUT >= 3) {
        cout << "meanGroupOne = " << endl << meanGroupOne;
    }

    Matrix sampleVariance(meanGroupOne.rows,meanGroupOne.rows);
    Matrix tempM;
    if (DEBUG_OUTPUT >= 3) {
        cout << "sampleVariance" << endl << sampleVariance << endl;
    }


    if (DEBUG_OUTPUT >= 3) {
        cout << "Starting loop." << endl;
    }
    for (list <svm_node *>::iterator svmit=groupOne.begin();svmit != groupOne.end();svmit++){
        tempM = svm_node_to_Matrix(*svmit);
        if (DEBUG_OUTPUT >= 3) {
            cout << "tempM" << endl << tempM;
            cout << "(tempM - meanGroupOne)" << endl << (tempM - meanGroupOne);
            cout << "(tempM - meanGroupOne).transpose()" << endl << (tempM - meanGroupOne).transpose();
            cout << "(tempM - meanGroupOne)*((tempM - meanGroupOne).transpose())" << endl << (tempM - meanGroupOne)*((tempM - meanGroupOne).transpose());
        }
        sampleVariance = sampleVariance + (tempM - meanGroupOne)*((tempM - meanGroupOne).transpose());
        if (DEBUG_OUTPUT >= 3) {
            cout << "sampleVariance" << endl << sampleVariance << endl;
        }
    }
    if (DEBUG_OUTPUT >= 3) {
        cout << "Done." << endl;
        cout << "sampleVariance" << endl << sampleVariance << endl;
    }
    for (int i=0;i<meanGroupOne.rows;i++){
        for (int j=0;j<meanGroupOne.rows;j++){
            sampleVariance(i,j) = sampleVariance(i,j)/(meanGroupOne.rows - 1);
        }
    }
    return sampleVariance;
}

Matrix svm_node_to_Matrix(svm_node *someNode)
{
    int nodeLength=0;

    while (someNode[nodeLength].index != -1){
        nodeLength++;
    }

    Matrix returnMatrix(nodeLength,1);
    for (int i=0;i<nodeLength;i++){
        returnMatrix(i,0) = someNode[i].value;
    }

    return returnMatrix;
}

svm_node *Matrix_to_svm_node(Matrix &someMatrix)
{
    svm_node *returnNode = new svm_node[someMatrix.rows + 1];

    for(int i=0;i<someMatrix.rows;i++){
        returnNode[i].value = someMatrix(i,0);
        returnNode[i].index = i+1;
    }
    
    returnNode[someMatrix.rows].index = -1;

    return returnNode;
}

void fillUnsignedSet(set <unsigned> &someSet, unsigned low, unsigned high)
{
    for (int i=low;i<=high;i++){
        someSet.insert(i);
    }
}

void fillUnsignedSetExclude(set <unsigned> &someSet, unsigned low, unsigned high, set <unsigned> &excludeSet)
{
    for (int i=low;i<=high;i++){
        if (excludeSet.find(i) == excludeSet.end()){
            someSet.insert(i);
        }
    }
}

list <string> listListStringToListString(list <list <string> > &TreeListListFileNames)
{
    list <string> treeFileNames;
    for (list <list <string> >::iterator llit=TreeListListFileNames.begin();llit!=TreeListListFileNames.end();llit++){
        for (list <string>::iterator lit=(*llit).begin();lit!=(*llit).end();lit++){
            treeFileNames.push_back(*lit);
        }
    }
    return treeFileNames;
}


double pvalueMin(list <double> groupOne, list <double> groupTwo)
{
    double currMin = *(groupOne.begin());
    for (list <double>::const_iterator ldit = groupOne.begin();ldit!=groupOne.end();ldit++){
        if (*ldit < currMin) {
            currMin = *ldit;
        }
    }

    double pvalue = 0;
    for (list <double>::const_iterator ldit = groupTwo.begin();ldit!=groupTwo.end();ldit++){
        if (*ldit > currMin) {
            pvalue++;
        }
    }
    pvalue /= groupTwo.size();

    return pvalue;
}


int GCD(int x,int y)
{
  if(y==0)  // base case, the programs stops if y reaches 0.
     return x;     //it returns the GCD
  else 
    return GCD(y,x%y); //if y doesn't reach 0 then recursion continues
}

int lcm(int a,int b)
{
    int n;
    for(n=1;;n++)
    {
        if(n%a == 0 && n%b == 0)
        {
          return n;
        }
    }
}

