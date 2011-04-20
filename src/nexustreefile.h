// $Rev: 359 $ $Date: 2009-11-15 23:57:12 -0500 (Sun, 15 Nov 2009) $

/** \file nexustreefile.h */
#ifndef NEXUSTREEFILE_H
#define NEXUSTREEFILE_H 1

#include <stdlib.h>
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include "debugoutput.h"
#include "treefuncs.h"

using namespace::std;

/// This class will be attached to some nexus file containing trees
/// Its main purpose is to retrieve user specified trees in the file.
/// Perhaps later can add functionality to dump entire trees.
class nexusTreeFile {
public:
    nexusTreeFile();
    nexusTreeFile(string fileName);
    ~nexusTreeFile();
    void init(string fileName);
    string getTree(int treeNumber);
protected:
    string  thisFileName;
    fstream inputFile;
    int currentTree; ///< Starts at 0. Store the tree which is next to be read.
};

#endif
