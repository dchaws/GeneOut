
/** \file treefuncs.h */
#ifndef TREEFUNCS_H
#define TREEFUNCS_H 1

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "debugoutput.h"
#include <algorithm>

using namespace::std;

/// It will attempt to return the number of trees contained in the file fileName.
int numTreesNexus(string fileName);

/// It will attempt to return the number of trees contained in the file &in.
int numTreesNexus(istream &in);

///
/// This will open the file treeFile, and return the treeNumber tree in the nexus file.
/// The first tree in the file will have index 0
///
string getTree(string treeFile, int treeNumber);

/// This takes in a newick tree and returns the number of leaves.
int countNumLeaves(string someTree);

/// This function will seek to the first tree of the nexus file inputFile.
void seekToFirstTree(istream &inputFile);

#endif
