// $Rev: 359 $ $Date: 2009-11-15 23:57:12 -0500 (Sun, 15 Nov 2009) $

/** \file jkparams.h */
#ifndef JKPARAMS_H
#define JKPARAMS_H 1

#include <string>
#include <stdlib.h>
#include <list>

using namespace::std;

///  Data structure to hold parameters for jackknifing(bootstraping with multiple columns)
class JackknifeParameters {
public:
    list        <string> inputNexFiles;
    int         jackknifeColSize;// = 10; // Default
    int         jackknifeCount;// = 1000; // Default
    string      alignToTreeCommand;
private:
};

#endif
