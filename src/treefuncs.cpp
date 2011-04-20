// $Rev: 774 $ $Date: 2011-02-09 17:53:15 -0500 (Wed, 09 Feb 2011) $

/** \file treefuncs.cpp */
#include "treefuncs.h"


string getTree(string treeFileName, int treeNumber)
{
    if (treeNumber < 0) {
        cout << "getTree: treeFileName = "  << treeFileName << "        treeNumber = " << treeNumber << endl;
        cout << "   treeNumber < 0" << endl;
        exit(1);
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "getTree: treeFileName = "  << treeFileName << "        treeNumber = " << treeNumber << endl;
    }
    fstream treeFile;
    treeFile.open(treeFileName.c_str(),fstream::in);
    if (!treeFile.good()){
        cout << " ##### Error trying to open " << treeFileName << " #####" << endl;
        exit (1);
    }

    string fileInput;
    // Search for '=' since this should indicate the first tree
    seekToFirstTree(treeFile);
    //while (getline(treeFile,fileInput) && fileInput.find("=") == string::npos){
    //    //cout << "#" << endl;
    //}
    //// Seek to first tree: we read past it.
    //treeFile.seekg(-fileInput.size(),ios_base::cur);

    getline(treeFile,fileInput);
    for (int i=1;i<treeNumber;i++){
        getline(treeFile,fileInput);
    }
    if (treeFile.eof()){
        cout << "getTree: Reached end of file before sampling!" << endl;
        exit(1);
    }

    string TreeInput = fileInput.substr(fileInput.find("("),fileInput.size()-fileInput.find("(")-1);
    treeFile.close();
    if (DEBUG_OUTPUT >= 1){
        cout << "getTree: " << TreeInput << endl;
    }
    return TreeInput;
}

void seekToFirstTree(istream &inputFile)
{
    if (DEBUG_OUTPUT >= 3){
        cout << "seekToFirstTree called." << endl;
    }
    string fileInputUpper;
    getline(inputFile,fileInputUpper);
    if (inputFile.bad()) {
        cout << "Can not read file given to seekToFirstTree." << endl;
        exit (0);
    }
    transform(fileInputUpper.begin(), fileInputUpper.end(),fileInputUpper.begin(), ::toupper);
    while (!( fileInputUpper.find("TREE") != string::npos && fileInputUpper.find("=") != string::npos))
    {
        if (DEBUG_OUTPUT >= 3){
            cout << fileInputUpper.substr(0,40) << endl;
        }
        getline(inputFile,fileInputUpper);
        if (inputFile.bad()) {
            cout << "Can not read file given to seekToFirstTree." << endl;
            exit (0);
        }
        transform(fileInputUpper.begin(), fileInputUpper.end(),fileInputUpper.begin(), ::toupper);
    }
    if (DEBUG_OUTPUT >= 3){
        cout << fileInputUpper.substr(0,40) << endl;
        cout << "inputFile.tellg(): "   << inputFile.tellg() << endl;
        cout << "inputFile.bad(): "     << inputFile.bad() << endl;
        cout << "inputFile.good(): "    << inputFile.good() << endl;
        cout << "inputFile.eof(): "     << inputFile.eof() << endl;
        cout << "inputFile.rdstate(): " << inputFile.rdstate() << endl;
    }
    // Something is happening to Peter's computer here!
    int curPosition = inputFile.tellg();
    // This might be causing problems to seek negative int from cur position
    // inputFile.seekg(-(fileInputUpper.size()+1),ios_base::cur);
    inputFile.seekg( curPosition - ((fileInputUpper.size()+1)), ios_base::beg);
    if (DEBUG_OUTPUT >= 3){
        cout << "inputFile.tellg(): " << inputFile.tellg() << endl;
        cout << "inputFile.bad(): " << inputFile.bad() << endl;
        cout << "inputFile.good(): " << inputFile.good() << endl;
        cout << "inputFile.eof(): " << inputFile.eof() << endl;
        cout << "inputFile.rdstate(): " << inputFile.rdstate() << endl;
        cout << "seekToFirstTree done." << endl;
    }
}

int numTreesNexus(string fileName)
{
    if (DEBUG_OUTPUT >= 3) {
        cout << "numTreesNexus: Opening " << fileName << endl;
    }
    fstream inputFile;
    inputFile.open(fileName.c_str(),fstream::in);
    if (inputFile.bad()){
        cout << "numTreesNexus: Could not open " << fileName << endl;
    }
    int returnValue = numTreesNexus(inputFile);
    inputFile.close();
    return returnValue;
}

int numTreesNexus(istream &in)
{
    int treeCount = 0;
    seekToFirstTree(in);
    if (in.eof()){
        cout << "numTreesNexus: Reached end of file before counting!" << endl;
    }
    string fileInputUpper;
    getline(in,fileInputUpper);
    if (DEBUG_OUTPUT >= 3) {
        cout << "numTreesNexus: fileInputUpper = " << fileInputUpper << endl;
    }
    transform(fileInputUpper.begin(), fileInputUpper.end(),fileInputUpper.begin(), ::toupper);
    while (( fileInputUpper.find("TREE") != string::npos && fileInputUpper.find("=") != string::npos)){
        //if (DEBUG_OUTPUT >= 3) {
        //    cout << "numTreesNexus: fileInputUpper = " << fileInputUpper << endl;
        //}
        //if (DEBUG_OUTPUT >= 3) {
        //    cout << "numTreesNexus: Found tree." << endl;
        //}
        getline(in,fileInputUpper);
        transform(fileInputUpper.begin(), fileInputUpper.end(),fileInputUpper.begin(), ::toupper);
        treeCount++;
    }

    return treeCount;
}

int countNumLeaves(string someTree)
{
    string TreeInput = someTree.substr(someTree.find("("),someTree.size()-someTree.find("(")-1);
    int leafCount = 0;
    if (DEBUG_OUTPUT >= 1){
        cout << "someTree = " << someTree << endl;
    }

    int currColon = TreeInput.find(":");
    while (currColon != string::npos) {
        if (DEBUG_OUTPUT >= 1){
            cout << "TreeInput = " << TreeInput << endl;
            cout << "currColon = " << currColon << endl;
            cout << "TreeInput[currColon-1] = " << TreeInput[currColon-1] << endl;
            cout << "TreeInput[currColon+1] = " << TreeInput[currColon+1] << endl;
        }
        // Check if there is a alphanum, non '(' or ')' before ':' and number after.
        if (isalpha(TreeInput[currColon-1]) || isdigit(TreeInput[currColon-1]) && isdigit(TreeInput[currColon+1])) {
            leafCount++;
        }
        TreeInput = TreeInput.substr(currColon+1,TreeInput.size()-currColon);
        currColon = TreeInput.find(":");
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "leafCount = " << leafCount << endl;
    }

    return leafCount;
}
