// $Rev: 359 $ $Date: 2009-11-15 23:57:12 -0500 (Sun, 15 Nov 2009) $

/** \file multindparams.h */
#ifndef MULTINDPARAMS_H
#define MULTINDPARAMS_H 1

#include <stdlib.h>
#include <string>
#include <list>
#include <set>

using namespace::std;

///  Data structure to hold parameters for multiple individual posterior distribution sampling.
class MultIndParameters {
public:
    list        <string> inputNexFiles;
    int         numTreesReconstruct;// = 1000; // Default
    string      alignToTreeCommand;
    list <set <unsigned> > speciesGroupings;
private:
};

#endif
