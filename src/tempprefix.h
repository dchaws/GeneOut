// $Rev: 342 $ $Date: 2009-11-14 21:34:26 -0500 (Sat, 14 Nov 2009) $

/** \file tempprefix.h */
#ifndef TEMPPREFIX_H
#define TEMPPREFIX_H 1

#include <string>
/// This is the prefix that all functions should use when creating temporary files.
/// If this is non-empty, the last character should be '/'
extern string tempPrefix;
#endif

