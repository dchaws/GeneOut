# GeneOut
GeneOut software for the paper "A support vector machine based test for incongruence between sets of trees in tree space", Haws, David C and Huggins, Peter and O’Neill, Eric M and Weisrock, David W and Yoshida, Ruriko, BMC bioinformatics, 2012. http://www.biomedcentral.com/1471-2105/13/210

# INSTALLATION:

Requires libsvm, blas, lapack, and a program to perform phylogenetic reconstruction (dnadist, phyml, mrbayes).

## libsvm
See: [libsvm](http://www.csie.ntu.edu.tw/~cjlin/libsvm/)

Download libsvm version 2.9 from (libsvm2.9)[http://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/]

    tar zxf libsvm-2.9.tar.gz
    cd libsvm-2.9
    make  # or make lib

    # Move libsvm.so.* and svm.h to a location such as ~/local/lib and ~/local/include respectively.
    mkdir -p ~/local/lib
    mkdir -p ~/local/include
    cp libsvm.so.1 ~/local/lib
    cp svm.h ~/local/include
    # create symbolic link
    ln -s ~/local/lib/libsvm.so.1 ~/local/lib/libsvm.so


## Compiling


Specify libsvm directory.

    ./configure --with-svm-prefix=/Users/haws/local/
    make


This will compile the program 'geneout' which will be in the ./src/ directory.


# Usage

# geneout

The geneout program performs the SVM separation test described in the paper above. Roughly,
it takes in two sets of alignments, A and B. It then performs multiple bootstraps of the 
alignment data, does phylogenetic reconstruction, and returns a p-value for the null-hypothesis
that the two sets of alignments were generated from the same distribution of (true) phylogenetic
histories. Thus, for geneout, one needs to specify the alignments (or trees directly), the
bootstrapping parameters, the number of permutation tests to perform, and the phylogenetic
reconstruction method (and parameters needed).


The geneout program is controlled by a .nex file which determines the input files and
the parameters.


## Example

    begin multsvmsep;
        method jk;
        distmethod jk;
        jkp colsize=10 count=1000;
        svmp samplesize=168 resamplesize=336;
        testparameters numinitcalc=100 numstattests=100;
        tempprefix "/tmp/";
        grouponefiles "GeneFromSpecies1Set1_Mesq/GeneFromSpecies1Set1_Mesq_0.nex";
        grouptwofiles "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_0.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_1.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_2.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_3.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_4.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_5.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_6.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_7.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_8.nex" "GeneFromSpecies2Set1_Mesq/GeneFromSpecies2Set1_Mesq_9.nex";
    end;


