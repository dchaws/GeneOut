// $Rev: 304 $ $Date: 2009-11-04 21:28:04 -0500 (Wed, 04 Nov 2009) $
#include <iostream>
#include "alignments.h"

int DEBUG_OUTPUT = 0;
int main (int argc, char **argv)
{
    int inputFormat = IF_NEXUS;
    if (argc > 1)
    {
        for (int i=1;i<argc;i++)
        {
            //cout << argv[i];
            string tempString = argv[i];
            if (tempString == "-p")
            {
                inputFormat = IF_PHYLIP;
            }
        }
    }
    Alignments myalign;
    myalign.setInputFormat(inputFormat);

    cin >> myalign;
    //cout << myalign;
    //myalign.setOutputFormat(OF_NEXUS);
    //cout << myalign;
    //myalign.setOutputFormat(OF_PHYLIP);
    //cout << myalign;
    cout.precision(10);
    cout << "Min: " << myalign.sequenceDivergenceMin() << endl;
    cout << "Avg: " << myalign.sequenceDivergenceAvg() << endl;
    myalign.printSequenceDivergencePairs();

}

