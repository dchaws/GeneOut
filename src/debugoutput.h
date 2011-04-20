// $Rev: 342 $ $Date: 2009-11-14 21:34:26 -0500 (Sat, 14 Nov 2009) $
#ifndef DEBUGOUTPUT_H
#define DEBUGOUTPUT_H 1

// Used to control debug output. 
//int DEBUG_OUTPUT = 0;
extern int DEBUG_OUTPUT;
// Used to give nice looking output.
#define DEBUGPAD for (int tempAGIEJSS=0;tempAGIEJSS<recurseLevel*3;tempAGIEJSS++){ cout << " "; }
#endif
