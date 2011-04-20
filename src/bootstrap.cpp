// $Rev: 304 $ $Date: 2009-11-04 21:28:04 -0500 (Wed, 04 Nov 2009) $
#include <iostream>
#include "alignments.h"

int DEBUG_OUTPUT = 0;
int main (int argc, char **argv)
{
    int inputFormat = IF_NEXUS;
    int outputFormat = OF_STANDARD;
    if (argc > 2)
    {
        for (int i=2;i<argc;i++)
        {
            //cout << argv[i];
            string tempString = argv[i];
            if (tempString == "-p")
            {
                inputFormat = IF_PHYLIP;
                outputFormat = OF_PHYLIP;
            }
        }
    }
    else {
        cout << "Need size of bootstrap columns. " << endl;
        exit(0);
    }

    // First argument is the bootstrap size.
    int numCols = atoi(argv[1]);

    int randomSeed = time(0);

    // Second argument is the random seed.
    if (argc > 3)
    {
        randomSeed = atoi(argv[2]);
    }

    srand(randomSeed);

    Alignments myalign;
    myalign.setInputFormat(inputFormat);

    cin >> myalign;
    //cout << myalign;
    //myalign.setOutputFormat(OF_NEXUS);
    //cout << myalign;
    //myalign.setOutputFormat(OF_PHYLIP);
    //cout << myalign;
    //cout.precision(10);
    //cout << "Min: " << myalign.sequenceDivergenceMin() << endl;
    //cout << "Avg: " << myalign.sequenceDivergenceAvg() << endl;
    //myalign.printSequenceDivergencePairs();


    Alignments bootstrapAlign = myalign.getJackknife(1, numCols);

    bootstrapAlign.setOutputFormat(outputFormat);
    cout << bootstrapAlign;


}

