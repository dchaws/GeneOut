#include <iostream>
#include "alignment.h"

int DEBUG_OUTPUT = 3;
int main (int argc, char **argv)
{
    Alignment myalign;

    cin >> myalign;
    cout << myalign;
    myalign.setOutputFormat(OF_NEXUS);
    cout << myalign;
    myalign.setOutputFormat(OF_PHYLIP);
    cout << myalign;
    cout << "Taking subcolumns" << endl;
    cout << myalign.getContiguousColumns(495,10);

    cout << "Deleting all non-ATCG columns" << endl;
    cout << myalign.delColumnsWithMissingData();

    cout << "Concatenating alignments" << endl;
    Alignment A = myalign + myalign;
    cout << A;
    //myalign += myalign;
    //cout << myalign;

    //cout << "Testing myalign again" << endl;
    //cout << myalign;
    cout << "Doing bootstrap test." << endl;
    srand(time(0));
    cout << myalign.getBootstrap(10,50);
    cout << myalign.getBootstrap(10,100);

    cout << "Doing subtaxa test." << endl;
    cout << "myalign has " << myalign.get_ntax() << " taxa." << endl;
    //cout << "Number of taxa to randomly select? ";
    int numTaxaSelect = 8;
    //cin >> numTaxaSelect;
    cout << "Selecting " << numTaxaSelect << " taxa." << endl;
    set <unsigned> someTaxa;
    srand((unsigned)time(0));
    while (someTaxa.size() != numTaxaSelect){
        int tmpint = rand() % myalign.get_ntax();
        someTaxa.insert(tmpint);
    }
    cout << endl;
    for (set <unsigned>::const_iterator sit=someTaxa.begin();sit!=someTaxa.end();sit++){
        cout << setw(3) << *sit << " ";
    }
    cout << endl;
    cout << myalign.getTaxaSubset(someTaxa);

}
