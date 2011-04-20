// $Rev: 369 $ $Date: 2009-11-17 14:57:48 -0500 (Tue, 17 Nov 2009) $

/** \file nexus.cpp */
#ifndef NEXUS_H
#define NEXUS_H 1

#include <string>
#include <list>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "debugoutput.h"


using namespace::std;

// This command eats up comments in the istream that start
// with '[' and end with ']'.
// We should handle nested ']', although it is  
// bad practice that any program would use nested comments...
void eatUpComments(std::istream &in,char &c);

int isPunctuation(char c);

// This function reads in the next Nexus token from the input
// A Nexux token is defined as a sequence of alphanum characters or punctions,
// e.g. "( ) , . ; : ! ? = "
//Ignore whitespaces,tabs,newlines,carraige returns. Also treat ';' as a token.
// This should ignore any characters between '[' and ']'
string nextNexusToken(std::istream &in);

string nextNexusTokenUpper(std::istream &in);

/// This will return everything starting at the current position betweeen '"' and '"'.
/// This will not handle [ ], i.e. it will return them. If this reads a ';' before
/// the first '"' it will return ';'
string nextQuoteBlock(std::istream &in);

/// This will return the next sequence of digits and digit-related alphanum characters
/// such as {-,+,.,E}.
string nextNexusNumber(std::istream &in);

// This function seeks/reads from the input until the next ; is found
// The input stream will be placed at the next character after ;
// This should ignore any characters between '[' and ']'
void readUntilSemiColonNexus(std::istream &in);

int errorCheckNexusToken(string someString,const char *str);

#endif

