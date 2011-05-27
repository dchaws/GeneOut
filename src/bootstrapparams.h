/** \file bootstrapparams.h */
#ifndef BOOTSTRAPPARAMS_H
#define BOOTSTRAPPARAMS_H 1

#include <string>
#include <stdlib.h>
#include <list>

using namespace::std;

///  Data structure to hold parameters for bootstraping with multiple columns
class BootstrapParameters {
public:
    list        <string> inputNexFiles;
    int         bootstrapColSize;// = 10; // Default
    int         bootstrapCount;// = 1000; // Default
    string      alignToTreeCommand;
private:
};

#endif
