// $Rev: 359 $ $Date: 2009-11-15 23:57:12 -0500 (Sun, 15 Nov 2009) $

/** \file nexustreefile.cpp */
#include "nexustreefile.h"

using namespace::std;

nexusTreeFile::nexusTreeFile()
{
}

nexusTreeFile::nexusTreeFile(string fileName)
{
    init(fileName);
}

nexusTreeFile::~nexusTreeFile()
{
    if (DEBUG_OUTPUT >= 1) {
        cout << "nexusTreeFile::~nexusTreeFile     Closing " << thisFileName << endl;
    }
    inputFile.close();
}

void nexusTreeFile::init(string fileName)
{
    thisFileName = fileName;
    currentTree = 0;
    inputFile.open(fileName.c_str(),fstream::in);
    if(inputFile.bad()){
        cout << "nexusTreeFile::init   can not open " << fileName << endl;
    }
    seekToFirstTree(inputFile);
}

string nexusTreeFile::getTree(int treeNumber)
{
    if (treeNumber < 0) {
        cout << "nexusTreeFile::getTree: treeFileName = "  << thisFileName << "     treeNumber = " << treeNumber << "   currentTree = " << currentTree << endl;      
        cout << "   treeNumber < 0" << endl;
        exit(1);
    }
    if (DEBUG_OUTPUT >= 1){
        cout << "nexusTreeFile::getTree: treeFileName = "  << thisFileName << "     treeNumber = " << treeNumber << "   currentTree = " << currentTree << endl;      
    }
    if (inputFile.eof()){
        cout << "nexusTreeFile::getTree: Reached end of file before sampling (pre-read)!" << endl;
        exit(1);
    }
    // If it is less than treeNumber, seek to begining.
    // This could later be improved to seek backwards in a smart way
    if (treeNumber < currentTree){
        if (DEBUG_OUTPUT >= 1){
            cout << "   Seeking to begining. " << endl;
        }
        inputFile.seekg(0,ios_base::beg);
        seekToFirstTree(inputFile);
        currentTree = 0;
    }
    string fileInput;

    for (int i=currentTree;i<=treeNumber;i++){
        getline(inputFile,fileInput);
        if (DEBUG_OUTPUT >= 3){
            cout << fileInput.substr(0,40) << endl;
        }
    }
    currentTree=treeNumber+1;
    if (inputFile.eof()){
        cout << "nexusTreeFile::getTree: Reached end of file before sampling!" << endl;
        exit(1);
    }

    string TreeInput = fileInput.substr(fileInput.find("("),fileInput.size()-fileInput.find("(")-1);
    if (DEBUG_OUTPUT >= 1){
        cout << "nexusTreeFile::getTree: " << TreeInput << endl;
    }
    return TreeInput;
}

