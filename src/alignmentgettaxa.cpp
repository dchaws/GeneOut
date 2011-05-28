#include <iostream>
#include "alignment.h"

int DEBUG_OUTPUT = 0;
int main (int argc, char **argv)
{
    //cout << argc << endl;
    //cout << argv[0] << endl;
    //cout << argv[1] << endl;
    //cout << argv[2] << endl;
    int numRepeats = 1;
    if (argc < 3) {
        printf ("Usage: alignmentsgettaxa <number of randomly selected taxa> <number of repeats>\n");
        exit(0);
    }
    if (argc > 2){
        numRepeats = atoi(argv[2]);
        if (numRepeats < 1) {
            printf("numRepeats < 1.\n");
            exit(0);
        }
    }
        
    Alignment myalign;

    cin >> myalign;
    //myalign.setOutputFormat(OF_NEXUS);
    //cout << myalign;
    int numTaxa = atoi (argv[1]);
    if (numTaxa < 1 || numTaxa >= myalign.get_ntax()) {
        printf(" numTaxa < 1 ||  numTaxa >= myalign.get_ntax()\n");
        exit(0);
    }
    //cout << "numTaxa = " << numTaxa << ".   numRepeats = " << numRepeats << endl;

    set <unsigned> taxaList;
    srand(time(0));

    for (int i=0;i<numRepeats;i++){
        taxaList.clear();
        while (taxaList.size() != numTaxa) {
            taxaList.insert(rand() % myalign.get_ntax());
        }
        Alignment subalign = myalign.getTaxaSubset(taxaList);
        subalign.setOutputFormat(OF_NEXUS);
        cout << subalign;
    }

}

