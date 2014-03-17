# GeneOut
GeneOut software for the paper "A support vector machine based test for incongruence between sets of trees in tree space", Haws, David C and Huggins, Peter and O’Neill, Eric M and Weisrock, David W and Yoshida, Ruriko, BMC bioinformatics, 2012. http://www.biomedcentral.com/1471-2105/13/210

# INSTALLATION:

Requires libsvm, blas, lapack, and a program to perform phylogenetic reconstruction (phylip (dnadist, neighbor), phyml, mrbayes).

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

The general nexus format is described in many places:
* http://sysbio.oxfordjournals.org/content/46/4/590
* http://en.wikipedia.org/wiki/Nexus_file
* http://wiki.christophchamp.com/index.php/NEXUS_file_format

The geneout program is controlled by a .nex file which determines the input files and
the parameters. Specifically, geneout is controlled by a new nexus 'multsvmsep' block.
All other parameters are given by lines starting with a command then with arguments, and
finally terminated by a semicolon.


    begin multsvmsep;
        command1 arg1 arg2 ... argN;
        command2 arg1 arg2 ... argN;
        command3 arg1 arg2 ... argN;
        ...
        commandM arg1 arg2 ... argN;
    end;

Note: Currently, geneout should be run from the ./src directory for simplicity.
This is because by default geneout looks for some specific scripts to perform
phylogenetic reconstruction.

## Control Block Commands

### testparameters

    testparameters numinitcalx=X numstattests=Y (setinitcalc=Z);

Controls the parameters of the overall permuation test being performed. The
command "numinitcalc=X" specifies that X svm separation tests should be
performed between the two sets of alignments. The average is then used as the
initial seperation percentage. The command "numstattests=Y" specifies that Y
permutation tests should be performed, that is the data is permunted and Y svm
seperations are performed. The optional command "setinitcal=Z" can be used to
instead set the initial svm seperation as Z.

### grouponefiles

    grouponefiles "file1" "file2" ... "fileN";

Specifies the N input alignment files associated with the first group of
alignements in the test.  This can be one or many files.

### grouptwofiles

    grouptwofiles "file1" "file2" ... "fileN";

Specifies the N input alignment files associated with the second group of
alignements in the test.  This can be one or many files.

### svmp 

    svmp samplesize=X resamplesize=Y (pathdifference=1) (scaletoone=1) (burninpercent=X) (burninnumber=X);

Controls the parameters of the svm tests. The command "samplesize=X" specifies
that X trees should be sampled from the space of possible trees to train the
svm seperating hyperplane. The command "resamplesize=Y" specifies that Y trees
should be sampled from the space of possible trees to compute the svm
seperation percentage. The optional command "pathdifference=1" specifies that
all the tree edge lengths should be set to 1. The optional command
"scaletoone=1" specifies that the vector of tree edge lengths should be scaled
such that the sum of all tree lengths is 1.

If Mr. Bayes is being used to reconstruct the trees, then the optional command
"burninpercent=X" specifies that the first X% of trees should be ignored, i.e.
treated as burnin. Alternatively, the optional command "burninnumber=X"
specifies that the first X lines of the tree file should be ignored.

### method

    method <bs|nobs|tree>;

Controls if bootstrapping is performed on the alignment (bs), if no
boostrapping is performed (nobs) or if trees are to be read in (tree).

### distmethod

    distmethod <bs|mb|trees>;

Controls which method is used to sample from the distrubution of phylogenetic
histories of the alignements.  Bootstrap the alignments and do single tree
reconstruction (bs), use Mr. Bayes. (mb), or read in a set of trees (trees).
For each option, other commands are used to give paramaters.

### bsp
    
        bsp colsize=X count=Y (aligntotreecommand="somecommand");

If "distmethod bs;" is selected, then the "bsp" command controls how
bootstrapping is performed to sample from the distrubition of phylogenetic
histories of the alignments. The parameter "colsize=X" controls how many
columns to contiguously sample during bootstrapping. The parameter "count=Y"
controls how many sets of bootstrap columns are taken. Thus, if "colsize=10"
and "count=1000" then 10 columns are taken at a time, and this is repeated 1000
times to create a new alignment of size 10000.

The command 'aligntotreecommand="somecommand"' is optional and is used to
specify the program that will take in an alignment and output a tree. [More
details on this needed. -Dave].

### mbp

    mbp mbp_nst=X mbp_ngen=X mbp_nruns=X mbp_samplefreq=X mpi_np=X (mbp_rates=X) (mbp_pathtomb="pathtomb") (mbp_paramfile="mbparamfile");

If "distmethod mb" is selected, then the "mbp" command controls how mr bayes is
called and used. The command "mbp_nst=X" specifies that X processors should be used,
default is 1. The command "mbp_nruns=X" specifies that X chains should be run, default is 2.
The command "mbp_samplefreq=X" specifies that trees should be sampled every X trees. The
optional command "mbp_rates=X" passes X along to Mr. Bayes. The optional command
'mbp_pathtomb="pathtomb"' specifies the specific path to Mr. Bayes program. The
optional command 'mbp_paramfile="mbparamfile"' specifies that Mr. Bayes should read
all of its parameters from the file "mbparamfile".

### tempprefix 

    tempprefix "/tmp/"

This command is optional and specifies the path to use to store temporary
files.  For example, if running on a HPC cluster, one should specify a
temporary directory that is local and not on the network.


### randseed

    randseed X;

This command is optional and specifies the random seed. Otherwise, the clock is
used.

## Example Nexus Files


Example using bootstrap to sample phylogenetic space.

    begin multsvmsep;
        method bs;
        distmethod bs;
        bsp colsize=10 count=1000;
        svmp samplesize=168 resamplesize=336;
        testparameters numinitcalc=10 numstattests=10;
        tempprefix "/tmp/";
        grouponefiles "../Examples/Alignments_SD_0.1Ne/Species1/GeneFromSpecies1Set1_Mesq_0.nex";
        grouptwofiles "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_0.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_1.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_2.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_3.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_4.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_5.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_6.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_7.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_8.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_9.nex";
    end;

Example using mrbayes to sample phylogenetic space.

    begin multsvmsep;
        method bs;
        distmethod mb;
        mbp mbp_nst=1 mbp_ngen=100000 mbp_nruns=2 mbp_samplefreq=100 mpi_np=4;
        svmp samplesize=168 resamplesize=336 burninpercent=0.25;
        testparameters numinitcalc=100 numstattests=100;
        tempprefix "/tmp/";
        grouponefiles "../Examples/Alignments_SD_0.1Ne/Species1/GeneFromSpecies1Set1_Mesq_0.nex";
        grouptwofiles "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_0.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_1.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_2.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_3.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_4.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_5.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_6.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_7.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_8.nex" "../Examples/Alignments_SD_0.1Ne/Species2/GeneFromSpecies2Set1_Mesq_9.nex";
    end;

# Example 
Example run:

    cd src
    ./geneout -nex ../Examples/test_bs_0.1Ne.nex

# Auxillory Scripts

The geneout program uses auxillory shell scripts in ./src directory in order to do phylogenetic
reconstruction using phlyip(defualt), phyml, or mr bayes.

The scripts files

* run_dnadist_neighbor_mult
* run_dnadist_F84-GR4_neighbor_mult
* run_dnadist_HKY85-GR4_neighbor_mult
* run_dnadist_JC_neighbor_mult

each take as input a single argument, which is a phylip ready alignment file. They each
in turn call 'dnadist' followed by 'neighbor' to perform phylogenetic reconstruction.
Each script in turn uses a template file to control dnadist to perform a specific type
of analysis. By default the script 'run_dnadist_neighbor_mult' is called, which in turn
just calls 'run_dnadist_F84-GR4_mult'. 

The template input files for 'dnadist' are

* dnadist_F84-GR4_input_mult_TEMPLATE
* dnadist_HKY85-GR4_input_mult_TEMPLAT
* dnadist_JC_input_mult_TEMPLATE	

For example, to use run_dnadist_JC_neighbor_mult for phylogenetic analysis, add the command
aligntotreecommand to the bsp command.

    bsp colsize=10 count=1000 aligntotreecommand="./run_dnadist_JC_neighbor_mult";

The script run_phyml_mult will use phyml to perform the phylogenetic analysis.

For example, to use run_phyml_mult for phylogenetic analysis, add the command
aligntotreecommand to the bsp command.

    bsp colsize=10 count=1000 aligntotreecommand="./run_phyml_mult";

to the nexus control file.


## run_phyml_mult

By default this script tries to run the command "PhyML_3.0". If another version is installed
edit the file run_phyml_mult and change 6th line:

    $PhyML="PhyML_3.0"

to the appropriate command.
