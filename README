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



`

