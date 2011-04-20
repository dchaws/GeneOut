// $Rev:  $ $Date:  $
#include <iostream>
#include "alignments.h"
#include <string>

int DEBUG_OUTPUT = -1;
int main (int argc, char **argv)
{
    if (argc > 1 ) {
        string tmpStr = argv[1];
        if (tmpStr == "-v") {
            DEBUG_OUTPUT = 3;
        }
    }
    Alignments myalign;

    cin >> myalign;
    myalign.setOutputFormat(OF_PHYLIP);
    cout << myalign;
}

