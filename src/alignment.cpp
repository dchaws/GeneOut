#include "alignment.h"
#include "debugoutput.h"
#include "nexus.h"


using namespace::std;

void Alignment::init ()
{
    //cout << "Alignment " << this << " allocated." << endl;
    nchar = 0;
    ntax = 0;
    //if (alignmentAllocated == 1){
    //    deleteAlignment();
    //}
    alignmentAllocated = 0;
    alignment = 0; // Set alignment to 0
    AlignmentName.clear();
    taxaNames.clear();
    outputFormat = OF_STANDARD;
    inputFormat = IF_NEXUS;
}
 
//Alignment::Alignment (Alignment &someAlignment)
//{
//    init();
//    *this = someAlignment;
//}

Alignment::Alignment (int tmp_ntax, int tmp_nchar, string tmp_AlignmentName, list <string> tmp_taxaNames, char **tmp_alignment, int tmp_alignmentAllocated )
{
    init (); // Some other stuff might need to be set
    ntax = tmp_ntax;
    nchar = tmp_nchar;
    AlignmentName = tmp_AlignmentName;
    taxaNames = tmp_taxaNames;
    alignment = tmp_alignment;
    alignmentAllocated = tmp_alignmentAllocated;
}

Alignment::Alignment (string fileName, int setInputFormat)
{
    init();
    inputFormat = setInputFormat;
    readFile(fileName);
}

Alignment::Alignment (string fileName)
{
    init ();
    readFile(fileName);
}

Alignment::Alignment ()
{
    init ();
}

Alignment::~Alignment ()
{
    //cout << "Alignment " << this << " deleted." << endl;
    // Delete alignment
    deleteAlignment();

}

void Alignment::readFile (string fileName)
{
    fstream inputFile;
    inputFile.open(fileName.c_str());
    if (inputFile.good()){
        inputFile >> *this;
    }
    else{
        cout << "Alignment::Alignment (string fileName): could not open " << fileName << endl;
        exit(0);
    }

}

void Alignment::deleteAlignment ()
{
    if (alignmentAllocated == 1)
    {
        for (int i=0;i<nchar;i++)
        {
            delete [] alignment[i];
        }
        delete [] alignment;
        alignment = 0;
    }
    alignmentAllocated = 0;
}

void Alignment::allocateAlignment()
{
    allocateAlignment(ntax,nchar);

}

void Alignment::allocateAlignment(int tmp_ntax,int tmp_nchar)
{
    if (alignmentAllocated == 0){
        if (tmp_ntax >= 1 && tmp_nchar >= 1)
        {
            alignment = new char*[tmp_nchar];
            for (int i=0;i<tmp_nchar;i++){
                alignment[i] = new char[tmp_ntax];
                for (int j=0;j<tmp_ntax;j++){
                    alignment[i][j] = 0;
                }
            }
            alignmentAllocated = 1;
        }
    }
    else {
        cout << "Alignment::allocateAlignment called with already allocated alignment." << endl;
        exit(0);
    }
}

Alignment::Alignment (std::istream &in)
{
    in >> *this;
}

std::ostream& operator << (std::ostream &out, const Alignment &someAlignment)
{
    if (someAlignment.outputFormat == OF_STANDARD){ // standard format

        if (!someAlignment.AlignmentName.empty())
        {
            out << "Alignment: " << someAlignment.AlignmentName << endl;

        }
        else {
            out << "Alignment:" << endl;
        }
        out << "Number of taxa: " << someAlignment.ntax << "    Number of characters: " << someAlignment.nchar << endl;
        if (someAlignment.alignmentAllocated == 1)
        {
            list <string>::const_iterator lsit;
            if (!someAlignment.taxaNames.empty())
            {
                lsit = someAlignment.taxaNames.begin();
            }
            for(int j=0;j<someAlignment.ntax;j++){
                if (!someAlignment.taxaNames.empty())
                {
                    out << *lsit << "          ";
                }
                for(int i=0;i<someAlignment.nchar;i++){
                    out << someAlignment.alignment[i][j];
                }
                out << endl;
                if (!someAlignment.taxaNames.empty())
                {
                    lsit++;
                }
            }
        }
        else {
            cout << "   No alignment." << endl;
        }
    }
    else if (someAlignment.outputFormat == OF_NEXUS){
        out << "#NEXUS" << endl;
        out << "begin data;" << endl;
        out << "dimensions ntax=" << someAlignment.ntax << " nchar=" << someAlignment.nchar << ";" << endl;
        out << "format datatype=dna missing=? gap=-;" << endl;
        out << "matrix" << endl;
        list <string>::const_iterator lsit;
        if (!someAlignment.taxaNames.empty())
        {
            lsit = someAlignment.taxaNames.begin();
        }
        for (int i=0;i<someAlignment.ntax;i++){
            if (!someAlignment.taxaNames.empty()) {
                out << "   " << *lsit << "     ";
            }
            else {
                out << "   " << i+1 << "     ";
            }
            for (int j=0;j<someAlignment.nchar;j++){
                out << someAlignment.alignment[j][i];
            }
            out << endl;
            if (!someAlignment.taxaNames.empty()) {
                lsit++;
            }
        }
        out << ";" << endl;
        out << "end;" << endl;
    }
    else if (someAlignment.outputFormat == OF_PHYLIP){
        out << someAlignment.ntax << " " << someAlignment.nchar << endl;

        list <string>::const_iterator lsit;
        if (!someAlignment.taxaNames.empty())
        {
            lsit = someAlignment.taxaNames.begin();
        }
        for (int i=0;i<someAlignment.ntax;i++){
            if (!someAlignment.taxaNames.empty() && 1 == 0) { //phylip can not handle these names
                if (DEBUG_OUTPUT >= 2){
                    cout << "Output taxa name: \"" << *lsit << setw(10) << "\"" << endl;
                }
                out << *lsit << setw(10);
                lsit++;
            }
            else {
                if (DEBUG_OUTPUT >= 2){
                    cout << "Output generic taxa name: \"" << "   " << i+1 << setw(10) << "\"" << endl;
                }
                //out << "   " << i+1 << setw(10);
                out << i+1 << setw(10); // For some reason phylip likes this best!!! So annoying
            }
            //out << i+1 << setw(10);
            for (int j=0;j<someAlignment.nchar;j++){
                out << someAlignment.alignment[j][i];
            }
            out << endl;
        }
    }

    return out;
}

std::istream& operator >> (std::istream &in, Alignment &someAlignment)
{
    // If we are going to read in a new alignment, delete previous alignment
    // and clear everything out, except the inputFormat;
    int tempFormat = someAlignment.inputFormat;
    someAlignment.init ();
    someAlignment.inputFormat = tempFormat;
    if (in.bad()){
        cout << "Alignment friend operator >>, input bad bit set" << endl;
        exit (0);
    }

    if (someAlignment.inputFormat == IF_NEXUS)
    {
        string input;
        // First look for '#NEXUS'
        input = nextNexusToken(in);
        if (input != "#NEXUS"){
            if (DEBUG_OUTPUT >= 0){
                cout << "Algorithms friend operator >>, '#NEXUS' not first line" << endl;
            }
            // Lets not be so strict about it.
            // Seek to begining
            in.seekg(0,ios::beg);
            //exit(0);
        }
        if (DEBUG_OUTPUT >= 1) {
            cout << "Algorithms friend operator >>, '#NEXUS' found" << endl;
        }

        // Now look for "TAXA", "CHARACTERS", or "DATA"
        int taxaRead = 0; // Change to 1 if "TAXA" block parsed.

        while (in.good()){
            input = nextNexusToken(in);
            transform(input.begin(), input.end(),input.begin(), ::toupper);
            if (input == "BEGIN"){
                input = nextNexusTokenUpper(in);
                if (input == "TAXA"){
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "BEGIN TAXA;" << endl;
                    }
                    input = nextNexusTokenUpper(in);
                    while (input != "END" && input != "ENDBLOCK"){
                        //cout << "input = " << input << endl;
                        if (input == "DIMENSIONS"){
                            input = nextNexusTokenUpper(in);
                            while (input != ";"){
                                // Look for ntax
                                if (input == "NTAX"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of taxa
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignment.ntax = " << input << endl;
                                    }
                                    someAlignment.ntax = atoi(input.c_str());
                                } 
                                input = nextNexusTokenUpper(in);
                            }
                            //readUntilSemiColonNexus(in);
                        }
                        else if (input == "TAXLABELS"){
                            input = nextNexusToken(in);
                            while (input != ";"){
                                someAlignment.taxaNames.push_back(input);
                                if (DEBUG_OUTPUT >= 1) {
                                    cout << "Adding label: " << input << endl;
                                }
                                input = nextNexusToken(in);
                            }
                        }
                        else {
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Skipping command" << endl;
                            }
                            readUntilSemiColonNexus(in);
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Done Skipping" << endl;
                            }
                        }
                        // Skip this command
                        input = nextNexusTokenUpper(in);
                    }
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "END;" << endl;
                    }
                }
                if (input == "CHARACTERS"){
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "BEGIN CHARACTERS;" << endl;
                    }
                    input = nextNexusTokenUpper(in);
                    while (input != "END" && input != "ENDBLOCK"){
                        //cout << "input = " << input << endl;
                        if (input == "DIMENSIONS"){
                            // Look for ntax
                            input = nextNexusTokenUpper(in);
                            while (input != ";"){
                                if (input == "NTAX"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of taxa
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignment.ntax = " << input << endl;
                                    }
                                    someAlignment.ntax = atoi(input.c_str());
                                } 
                                else if (input == "NCHAR"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of characters 
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignment.nchar = " << input << endl;
                                    }
                                    someAlignment.nchar = atoi(input.c_str());
                                }
                                input = nextNexusTokenUpper(in);
                            }
                        }
                        else if (input == "FORMAT"){
                            input = nextNexusToken(in);
                            while (input != ";"){
                                if (input == "DATATYPE"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the data type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Data Type is " << input << endl;
                                    }
                                }
                                else if (input == "MISSING"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the missing character type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Missing character " << input << endl;
                                    }
                                }
                                else if (input == "GAP"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the missing character type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Gap character " << input << endl;
                                    }
                                }
                                input = nextNexusTokenUpper(in);
                            }
                        }
                        else if (input == "MATRIX"){
                            // Allocate the alignment
                            someAlignment.allocateAlignment();

                            unsigned curAlignment = 0;
                            input = nextNexusToken(in); // This should be the name of the first alignment
                            if (input != ";" ) {
                                if (DEBUG_OUTPUT > 0 ) {
                                    cout << "Taxa name: " << input << endl;
                                }
                                someAlignment.taxaNames.push_back(input);
                            }
                            while (input != ";"){
                                string alignment;
                                unsigned curAlignmentPos = 0;
                                unsigned nexStringPos = 0;
                                string nexString;
                                getline(in,nexString);
                                if (DEBUG_OUTPUT >= 1) {
                                    cout << "nexString = " << nexString << endl;
                                }
                                while (nexString[nexStringPos] != 0 && curAlignmentPos < someAlignment.nchar){
                                    if (nexString[nexStringPos] == '['){
                                        nexStringPos++;
                                        while (nexString[nexStringPos] != ']' && nexString[nexStringPos] != 0){
                                            nexStringPos++;
                                        }
                                    }

                                    if (nexString[nexStringPos] != ' '){
                                        someAlignment.alignment[curAlignmentPos][curAlignment] = nexString[nexStringPos]; 
                                        curAlignmentPos++;
                                    }
                                    nexStringPos++;
                                }
                                curAlignment++;
                                input = nextNexusToken(in); // This should be the name of the next alignment
                                if (input != ";" ) {
                                    if (DEBUG_OUTPUT > 0 ) {
                                        cout << "Taxa name: " << input << endl;
                                    }
                                    someAlignment.taxaNames.push_back(input);
                                }
                            }
                        }
                        else {
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Skipping command" << endl;
                            }
                            readUntilSemiColonNexus(in);
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Done Skipping" << endl;
                            }
                        }
                        // Skip this command
                        input = nextNexusTokenUpper(in);
                    }
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "END;" << endl;
                    }
                }
                if (input == "DATA"){
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "BEGIN DATA;" << endl;
                    }
                    input = nextNexusTokenUpper(in);
                    while (input != "END" && input != "ENDBLOCK"){
                        //cout << "input = " << input << endl;
                        if (input == "DIMENSIONS"){
                            // Look for ntax
                            input = nextNexusTokenUpper(in);
                            while (input != ";"){
                                if (input == "NTAX"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of taxa
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignment.ntax = " << input << endl;
                                    }
                                    someAlignment.ntax = atoi(input.c_str());
                                } 
                                else if (input == "NCHAR"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of characters 
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignment.nchar = " << input << endl;
                                    }
                                    someAlignment.nchar = atoi(input.c_str());
                                }
                                input = nextNexusTokenUpper(in);
                            }
                        }
                        else if (input == "FORMAT"){
                            input = nextNexusToken(in);
                            while (input != ";"){
                                if (input == "DATATYPE"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the data type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Data Type is " << input << endl;
                                    }
                                }
                                else if (input == "MISSING"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the missing character type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Missing character " << input << endl;
                                    }
                                }
                                else if (input == "GAP"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the missing character type
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Gap character " << input << endl;
                                    }
                                }
                                input = nextNexusTokenUpper(in);
                            }
                        }
                        else if (input == "MATRIX"){
                            // Allocate the alignment
                            someAlignment.allocateAlignment();

                            unsigned curAlignment = 0;
                            input = nextNexusToken(in); // This should be the name of the first alignment
                            if (input != ";" ) {
                                if (DEBUG_OUTPUT > 0 ) {
                                    cout << "Taxa name: " << input << endl;
                                }
                                someAlignment.taxaNames.push_back(input);
                            }
                            while (input != ";"){
                                string alignment;
                                unsigned curAlignmentPos = 0;
                                unsigned nexStringPos = 0;
                                string nexString;
                                getline(in,nexString);
                                if (DEBUG_OUTPUT >= 1) {
                                    cout << "nexString = " << nexString << endl;
                                }
                                while (nexString[nexStringPos] != 0 && curAlignmentPos < someAlignment.nchar){
                                    if (nexString[nexStringPos] == '['){
                                        nexStringPos++;
                                        while (nexString[nexStringPos] != ']' && nexString[nexStringPos] != 0){
                                            nexStringPos++;
                                        }
                                    }
                                    if (nexString[nexStringPos] != ' '){
                                        someAlignment.alignment[curAlignmentPos][curAlignment] = nexString[nexStringPos]; 
                                        curAlignmentPos++;
                                    }
                                    nexStringPos++;
                                }
                                curAlignment++;
                                input = nextNexusToken(in);
                                if (input != ";" ) {
                                    if (DEBUG_OUTPUT > 0 ) {
                                        cout << "Taxa name: " << input << endl;
                                    }
                                    someAlignment.taxaNames.push_back(input);
                                }

                            }
                        }
                        else {
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Skipping command" << endl;
                            }
                            readUntilSemiColonNexus(in);
                            if (DEBUG_OUTPUT >= 1) {
                                cout << "Done Skipping" << endl;
                            }
                        }
                        // Skip this command
                        input = nextNexusTokenUpper(in);
                    }
                    readUntilSemiColonNexus(in);
                    if (DEBUG_OUTPUT >= 1) {
                        cout << "END;" << endl;
                    }
                }
                if (input == "DATA"){
                    input = nextNexusTokenUpper(in);
                    if (input == ";") {
                        if (DEBUG_OUTPUT >= 1) {
                            cout << "BEGIN DATA;" << endl;
                        }
                    }
                }
            }
        }
    }
    if (someAlignment.inputFormat == IF_PHYLIP)
    {
        // First item should be taxa
        in >> someAlignment.ntax;
        in >> someAlignment.nchar;
        someAlignment.allocateAlignment(someAlignment.ntax,someAlignment.nchar);

        // Now read in taxa names, then string of characters.
        string tempInput;
        for (int i=0;i<someAlignment.ntax;i++)
        {
            // Read in name
            in >> tempInput;
            someAlignment.taxaNames.push_back(tempInput);
            // Read in DNA
            in >> tempInput;
            for (int j=0;j<someAlignment.nchar;j++)
            {
                someAlignment.alignment[j][i] = tempInput[j];
            }
            
        }
    }

    return in;
}


Alignment & Alignment::operator = (const Alignment &someAlignment)
{
    this->deleteAlignment();
    this->init();
    this->ntax = someAlignment.ntax;
    this->nchar = someAlignment.nchar;
    this->taxaNames = someAlignment.taxaNames;
    this->AlignmentName = someAlignment.AlignmentName;
    this->allocateAlignment(); //Allocate the alignment with ntax and nchar
    for (int i=0;i<this->ntax;i++){
        for (int j=0;j<this->nchar;j++){
            this->alignment[j][i] = someAlignment.alignment[j][i];
        }
    }
    this->alignmentAllocated = 1;
    this->outputFormat = someAlignment.outputFormat;

    return *this;
}

bool Alignment::operator == (const Alignment &someAlignment)
{
    int allSame = 1;

    for (int i=0;i<this->nchar;i++){
        for (int j=0;j<this->ntax;j++){
            if (this->alignment[i][j] != someAlignment.alignment[i][j]) {
                allSame = 0;
            }
        }
    }
    if (allSame == 1){
        return 1;
    }
    else{
        return 0;
    }
}

Alignment Alignment::operator + (const Alignment &someAlignment)
{
    if (this->ntax != someAlignment.ntax){
        cout << "Alignment::operator +: this->ntax != sA->ntax" << endl;
    }

    if (this->ntax != someAlignment.ntax){
        cout << "Alignment::operator +=: this->ntax != someAlignment.get_ntax()" << endl;
        exit(0);
    }
    char **new_alignment;
    new_alignment = new char*[this->nchar + someAlignment.nchar];
    for (int i=0;i<this->nchar + someAlignment.nchar;i++){
        new_alignment[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            if (i < this->nchar){
                new_alignment[i][j] = this->alignment[i][j];
            }
            else {
                new_alignment[i][j] = someAlignment(i-this->nchar,j);
            }
        }
    }
    Alignment newAlignment(this->ntax,this->nchar + someAlignment.nchar, this->AlignmentName, this->taxaNames,new_alignment,1);

    return newAlignment;
}

Alignment & Alignment::operator += (const Alignment &someAlignment)
{
    //if (this->ntax != someAlignment.get_ntax()){
    if (this->ntax != someAlignment.ntax){
        cout << "Alignment::operator +=: this->ntax != someAlignment.get_ntax()" << endl;
        exit(0);
    }
    char **new_alignment;
    new_alignment = new char*[this->nchar + someAlignment.nchar];
    //for (int i=0;i<this->nchar + someAlignment.get_nchar();i++){
    for (int i=0;i<this->nchar + someAlignment.nchar;i++){
        new_alignment[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            if (i < this->nchar){
                new_alignment[i][j] = this->alignment[i][j];
            }
            else {
                new_alignment[i][j] = someAlignment(i-this->nchar,j);
            }
        }
    }

    this->deleteAlignment();
    this->alignment = new_alignment;
    this->alignmentAllocated = 1;
    this->nchar = this->nchar + someAlignment.nchar;

    // We need to concat someAlignents to this object
    return *this;
}

char Alignment::operator() (unsigned row, unsigned col) const
{
    return alignment[row][col];
}

char & Alignment::operator() (unsigned row, unsigned col)
{
    return alignment[row][col];
}

const int Alignment::get_ntax()
{
    return ntax;
}

const int Alignment::get_nchar()
{
    return nchar;
}

const list <string> Alignment::get_taxaNames()
{
    return taxaNames;
}

const string Alignment::get_AlignmentName()
{   
    return AlignmentName;
}

void Alignment::setOutputFormat(int format)
{
    outputFormat = format;
}

void Alignment::setInputFormat(int format)
{
    inputFormat = format;
}

Alignment Alignment::getContiguousColumns(int pos, int k)
{
    if (pos >= nchar || k > nchar) {
        cout << "Alignment::getContiguousColumns: pos >= nchar || k > nchar" << endl;
        exit(0);
    }
    char **new_alignment;
    new_alignment = new char*[k];
    for (int i=0;i<k;i++){
        new_alignment[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            new_alignment[i][j] = this->alignment[(pos + i) % nchar][j];
        }
    }
    Alignment newAlignment(this->ntax,k, this->AlignmentName, this->taxaNames,new_alignment,1);

    return newAlignment;
}

Alignment Alignment::delColumnsWithMissingData()
{ 
    char **newAlignment;
    int new_nchar = 0;
    int delColumn[nchar];

    for(int i=0;i<nchar;i++){
        delColumn[i]=0;
    }

    // Go through every column
    for(int i=0;i<nchar;i++){
        //cout << "Checking column " << i << endl;
        for (int j=0;j<ntax;j++){
            if (alignment[i][j] != 'A' && alignment[i][j] != 'T' && alignment[i][j] != 'C' && alignment[i][j] != 'G'){
                //cout << "       Bad character alignment[" << i << "][" << j << "] = " << alignment[i][j] << endl;
                delColumn[i]=1; // Delete this column
            }
        }
        if (delColumn[i] == 0){
            //cout << "   Column GOOD" << endl;
            new_nchar++;
        }
        else {
            //cout << "   Column BAD" << endl;
        }
    }

    newAlignment = new char*[new_nchar];
    for(int i=0;i<new_nchar;i++){
        newAlignment[i] = new char[ntax];
    }
    int tempCount = 0;

    for(int i=0;i<nchar;i++){
        if (delColumn[i] == 0){
            for (int j=0;j<ntax;j++){
                newAlignment[tempCount][j] = alignment[i][j];
            }

            tempCount++;
        }
    }
    
    Alignment returnAlignment(this->ntax,new_nchar, this->AlignmentName, this->taxaNames,newAlignment,1);

    return returnAlignment;

}

Alignment Alignment::getBootstrap(int k)
{
    return getBootstrap(k,(int)floor((double)nchar/k));

}

Alignment Alignment::getBootstrap(int k, int count)
{
    Alignment returnAlignment;

    if (count > 0){
        returnAlignment = this->getContiguousColumns(rand() % this->nchar,k);

        for (int i=1;i<count;i++){
            returnAlignment += this->getContiguousColumns(rand() % this->nchar,k);
        }
    }
    return returnAlignment;

}

Alignment Alignment::getTaxaSubset(set <unsigned> someTaxa)
{
    if (DEBUG_OUTPUT > 1 ){
        cout << "Alignment::getTaxaSubset called." << endl;
        for (set <unsigned>::const_iterator sit=someTaxa.begin();sit != someTaxa.end();sit++) {
            cout << *sit << " ";
        }
        cout << endl;
    }

    char **new_alignment;
    new_alignment = new char*[nchar];
    for (int i=0;i<nchar;i++){
        new_alignment[i] = new char[someTaxa.size()];
        unsigned rowCount = 0;
        for (set <unsigned>::const_iterator sit=someTaxa.begin();sit!=someTaxa.end();sit++){
            new_alignment[i][rowCount] = this->alignment[i][*sit];
            rowCount++;
        }
    }
    list <string> newTaxaNames;
    if (!this->taxaNames.empty()) {
        list <string>::const_iterator lsit=this->taxaNames.begin();
        set <unsigned>::const_iterator sit=someTaxa.begin();
        for (int j=0;j<ntax;j++){
            if (DEBUG_OUTPUT > 0 ){
                cout << "j=" << j << "  *sit=" << *sit << "   *lsit=" << *lsit << "   newTaxaNames.size()=" << newTaxaNames.size() << endl;
            }
            if (j==*sit && sit != someTaxa.end()){
                if (DEBUG_OUTPUT > 0 ){
                    cout << "j==*sit && sit != someTaxa.end()"  << endl;
                }
                newTaxaNames.push_back(*lsit);
                sit++;
            }
            lsit++;
        }
    }
    Alignment newAlignment(someTaxa.size(),nchar, this->AlignmentName, newTaxaNames,new_alignment,1);

    if (DEBUG_OUTPUT > 1 ){
        cout << "Alignment::getTaxaSubset done." << endl;
    }
    return newAlignment;
}

double Alignment::sequenceDivergenceMin()
{
    if (alignmentAllocated == 0) {
        cout << "Alignment::sequenceDivergence called with no alignment allocated." << endl;
        exit (1);
    }
    int minHamming = -1;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignment[][i] and alignment[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignment[k][i] != alignment[k][j]) {
                    currHamming++;
                }
            }
            if (currHamming < minHamming || minHamming == -1) {
                minHamming = currHamming;
            }
        }
    }
    
    return (double)(((double)minHamming)/(double)nchar);
}

double Alignment::sequenceDivergenceAvg()
{
    if (alignmentAllocated == 0) {
        cout << "Alignment::sequenceDivergence called with no alignment allocated." << endl;
        exit (1);
    }
    int totalHamming = 0;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignment[][i] and alignment[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignment[k][i] != alignment[k][j]) {
                    currHamming++;
                }
            }
            totalHamming += currHamming;
        }
    }
    double tmpDouble = (double)(ntax-1)*(double)ntax/2.0;
    
    return (double)(((double)totalHamming)/(double)((double)nchar*tmpDouble));
}


void Alignment::printSequenceDivergencePairs()
{
    if (alignmentAllocated == 0) {
        cout << "Alignment::sequenceDivergence called with no alignment allocated." << endl;
        exit (1);
    }
    int totalHamming = 0;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignment[][i] and alignment[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignment[k][i] != alignment[k][j]) {
                    currHamming++;
                }
            }
            // This is the i+1 and j+1 alignment, but we need the name!
            // cout << "(" << i+1 << "," << j+1 << ") " << (double)currHamming/(double)nchar << endl;
            string ithName, jthName;
            list <string>::const_iterator lsit;
            lsit = taxaNames.begin();
            for (int k=0;k<i && lsit != taxaNames.end();k++)
            {
                lsit++;
            }
            ithName = *lsit;

            lsit = taxaNames.begin();
            for (int k=0;k<j && lsit != taxaNames.end();k++)
            {
                lsit++;
            }
            jthName = *lsit;
            cout << "(" << ithName << "," << jthName << ") " << (double)currHamming/(double)nchar << endl;
        }
    }
    double tmpDouble = (double)(ntax-1)*(double)ntax/2.0;
    
}
