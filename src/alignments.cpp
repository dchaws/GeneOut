// $Rev: 807 $ $Date: 2011-04-10 13:21:55 -0400 (Sun, 10 Apr 2011) $
#include "alignments.h"
#include "debugoutput.h"
#include "nexus.h"


using namespace::std;

void Alignments::init ()
{
    //cout << "Alignments " << this << " allocated." << endl;
    nchar = 0;
    ntax = 0;
    //if (alignmentsAllocated == 1){
    //    deleteAlignments();
    //}
    alignmentsAllocated = 0;
    alignments = 0; // Set alignments to 0
    AlignmentsName.clear();
    taxaNames.clear();
    outputFormat = OF_STANDARD;
    inputFormat = IF_NEXUS;
}
 
//Alignments::Alignments (Alignments &someAlignments)
//{
//    init();
//    *this = someAlignments;
//}

Alignments::Alignments (int tmp_ntax, int tmp_nchar, string tmp_AlignmentsName, list <string> tmp_taxaNames, char **tmp_alignments, int tmp_alignmentsAllocated )
{
    init (); // Some other stuff might need to be set
    ntax = tmp_ntax;
    nchar = tmp_nchar;
    AlignmentsName = tmp_AlignmentsName;
    taxaNames = tmp_taxaNames;
    alignments = tmp_alignments;
    alignmentsAllocated = tmp_alignmentsAllocated;
}

Alignments::Alignments (string fileName, int setInputFormat)
{
    init();
    inputFormat = setInputFormat;
    readFile(fileName);
}

Alignments::Alignments (string fileName)
{
    init ();
    readFile(fileName);
}

Alignments::Alignments ()
{
    init ();
}

Alignments::~Alignments ()
{
    //cout << "Alignments " << this << " deleted." << endl;
    // Delete alignments
    deleteAlignments();

}

void Alignments::readFile (string fileName)
{
    fstream inputFile;
    inputFile.open(fileName.c_str());
    if (inputFile.good()){
        inputFile >> *this;
    }
    else{
        cout << "Alignments::Alignments (string fileName): could not open " << fileName << endl;
        exit(0);
    }

}

void Alignments::deleteAlignments ()
{
    if (alignmentsAllocated == 1)
    {
        for (int i=0;i<nchar;i++)
        {
            delete [] alignments[i];
        }
        delete [] alignments;
        alignments = 0;
    }
    alignmentsAllocated = 0;
}

void Alignments::allocateAlignments()
{
    allocateAlignments(ntax,nchar);

}

void Alignments::allocateAlignments(int tmp_ntax,int tmp_nchar)
{
    if (alignmentsAllocated == 0){
        if (tmp_ntax >= 1 && tmp_nchar >= 1)
        {
            alignments = new char*[tmp_nchar];
            for (int i=0;i<tmp_nchar;i++){
                alignments[i] = new char[tmp_ntax];
                for (int j=0;j<tmp_ntax;j++){
                    alignments[i][j] = 0;
                }
            }
            alignmentsAllocated = 1;
        }
    }
    else {
        cout << "Alignments::allocateAlignments called with already allocated alignments." << endl;
        exit(0);
    }
}

Alignments::Alignments (std::istream &in)
{
    in >> *this;
}

std::ostream& operator << (std::ostream &out, const Alignments &someAlignments)
{
    if (someAlignments.outputFormat == OF_STANDARD){ // standard format

        if (!someAlignments.AlignmentsName.empty())
        {
            out << "Alignments: " << someAlignments.AlignmentsName << endl;

        }
        else {
            out << "Alignments:" << endl;
        }
        out << "Number of taxa: " << someAlignments.ntax << "    Number of characters: " << someAlignments.nchar << endl;
        if (someAlignments.alignmentsAllocated == 1)
        {
            list <string>::const_iterator lsit;
            if (!someAlignments.taxaNames.empty())
            {
                lsit = someAlignments.taxaNames.begin();
            }
            for(int j=0;j<someAlignments.ntax;j++){
                if (!someAlignments.taxaNames.empty())
                {
                    out << *lsit << "          ";
                }
                for(int i=0;i<someAlignments.nchar;i++){
                    out << someAlignments.alignments[i][j];
                }
                out << endl;
                if (!someAlignments.taxaNames.empty())
                {
                    lsit++;
                }
            }
        }
        else {
            cout << "   No alignments." << endl;
        }
    }
    else if (someAlignments.outputFormat == OF_NEXUS){
        out << "#NEXUS" << endl;
        out << "begin data;" << endl;
        out << "dimensions ntax=" << someAlignments.ntax << " nchar=" << someAlignments.nchar << ";" << endl;
        out << "format datatype=dna missing=? gap=-;" << endl;
        out << "matrix" << endl;
        list <string>::const_iterator lsit;
        if (!someAlignments.taxaNames.empty())
        {
            lsit = someAlignments.taxaNames.begin();
        }
        for (int i=0;i<someAlignments.ntax;i++){
            if (!someAlignments.taxaNames.empty()) {
                out << "   " << *lsit << "     ";
            }
            else {
                out << "   " << i+1 << "     ";
            }
            for (int j=0;j<someAlignments.nchar;j++){
                out << someAlignments.alignments[j][i];
            }
            out << endl;
            if (!someAlignments.taxaNames.empty()) {
                lsit++;
            }
        }
        out << ";" << endl;
        out << "end;" << endl;
    }
    else if (someAlignments.outputFormat == OF_PHYLIP){
        out << someAlignments.ntax << " " << someAlignments.nchar << endl;

        list <string>::const_iterator lsit;
        if (!someAlignments.taxaNames.empty())
        {
            lsit = someAlignments.taxaNames.begin();
        }
        for (int i=0;i<someAlignments.ntax;i++){
            if (!someAlignments.taxaNames.empty() && 1 == 0) { //phylip can not handle these names
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
            for (int j=0;j<someAlignments.nchar;j++){
                out << someAlignments.alignments[j][i];
            }
            out << endl;
        }
    }

    return out;
}

std::istream& operator >> (std::istream &in, Alignments &someAlignments)
{
    // If we are going to read in a new alignment, delete previous alignments
    // and clear everything out, except the inputFormat;
    int tempFormat = someAlignments.inputFormat;
    someAlignments.init ();
    someAlignments.inputFormat = tempFormat;
    if (in.bad()){
        cout << "Alignments friend operator >>, input bad bit set" << endl;
        exit (0);
    }

    if (someAlignments.inputFormat == IF_NEXUS)
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
                                        cout << "Setting someAlignments.ntax = " << input << endl;
                                    }
                                    someAlignments.ntax = atoi(input.c_str());
                                } 
                                input = nextNexusTokenUpper(in);
                            }
                            //readUntilSemiColonNexus(in);
                        }
                        else if (input == "TAXLABELS"){
                            input = nextNexusToken(in);
                            while (input != ";"){
                                someAlignments.taxaNames.push_back(input);
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
                                        cout << "Setting someAlignments.ntax = " << input << endl;
                                    }
                                    someAlignments.ntax = atoi(input.c_str());
                                } 
                                else if (input == "NCHAR"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of characters 
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignments.nchar = " << input << endl;
                                    }
                                    someAlignments.nchar = atoi(input.c_str());
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
                            // Allocate the alignments
                            someAlignments.allocateAlignments();

                            unsigned curAlignment = 0;
                            input = nextNexusToken(in); // This should be the name of the first alignment
                            if (input != ";" ) {
                                if (DEBUG_OUTPUT > 0 ) {
                                    cout << "Taxa name: " << input << endl;
                                }
                                someAlignments.taxaNames.push_back(input);
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
                                while (nexString[nexStringPos] != 0 && curAlignmentPos < someAlignments.nchar){
                                    if (nexString[nexStringPos] == '['){
                                        nexStringPos++;
                                        while (nexString[nexStringPos] != ']' && nexString[nexStringPos] != 0){
                                            nexStringPos++;
                                        }
                                    }

                                    if (nexString[nexStringPos] != ' '){
                                        someAlignments.alignments[curAlignmentPos][curAlignment] = nexString[nexStringPos]; 
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
                                    someAlignments.taxaNames.push_back(input);
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
                                        cout << "Setting someAlignments.ntax = " << input << endl;
                                    }
                                    someAlignments.ntax = atoi(input.c_str());
                                } 
                                else if (input == "NCHAR"){
                                    input = nextNexusToken(in);
                                    errorCheckNexusToken(input,"=");
                                    input = nextNexusToken(in);
                                    // input should now be the number of characters 
                                    if (DEBUG_OUTPUT >= 1) {
                                        cout << "Setting someAlignments.nchar = " << input << endl;
                                    }
                                    someAlignments.nchar = atoi(input.c_str());
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
                            // Allocate the alignments
                            someAlignments.allocateAlignments();

                            unsigned curAlignment = 0;
                            input = nextNexusToken(in); // This should be the name of the first alignment
                            if (input != ";" ) {
                                if (DEBUG_OUTPUT > 0 ) {
                                    cout << "Taxa name: " << input << endl;
                                }
                                someAlignments.taxaNames.push_back(input);
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
                                while (nexString[nexStringPos] != 0 && curAlignmentPos < someAlignments.nchar){
                                    if (nexString[nexStringPos] == '['){
                                        nexStringPos++;
                                        while (nexString[nexStringPos] != ']' && nexString[nexStringPos] != 0){
                                            nexStringPos++;
                                        }
                                    }
                                    if (nexString[nexStringPos] != ' '){
                                        someAlignments.alignments[curAlignmentPos][curAlignment] = nexString[nexStringPos]; 
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
                                    someAlignments.taxaNames.push_back(input);
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
    if (someAlignments.inputFormat == IF_PHYLIP)
    {
        // First item should be taxa
        in >> someAlignments.ntax;
        in >> someAlignments.nchar;
        someAlignments.allocateAlignments(someAlignments.ntax,someAlignments.nchar);

        // Now read in taxa names, then string of characters.
        string tempInput;
        for (int i=0;i<someAlignments.ntax;i++)
        {
            // Read in name
            in >> tempInput;
            someAlignments.taxaNames.push_back(tempInput);
            // Read in DNA
            in >> tempInput;
            for (int j=0;j<someAlignments.nchar;j++)
            {
                someAlignments.alignments[j][i] = tempInput[j];
            }
            
        }
    }

    return in;
}


Alignments & Alignments::operator = (const Alignments &someAlignments)
{
    this->deleteAlignments();
    this->init();
    this->ntax = someAlignments.ntax;
    this->nchar = someAlignments.nchar;
    this->taxaNames = someAlignments.taxaNames;
    this->AlignmentsName = someAlignments.AlignmentsName;
    this->allocateAlignments(); //Allocate the alignments with ntax and nchar
    for (int i=0;i<this->ntax;i++){
        for (int j=0;j<this->nchar;j++){
            this->alignments[j][i] = someAlignments.alignments[j][i];
        }
    }
    this->alignmentsAllocated = 1;
    this->outputFormat = someAlignments.outputFormat;

    return *this;
}

bool Alignments::operator == (const Alignments &someAlignments)
{
    int allSame = 1;

    for (int i=0;i<this->nchar;i++){
        for (int j=0;j<this->ntax;j++){
            if (this->alignments[i][j] != someAlignments.alignments[i][j]) {
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

Alignments Alignments::operator + (const Alignments &someAlignments)
{
    if (this->ntax != someAlignments.ntax){
        cout << "Alignments::operator +: this->ntax != sA->ntax" << endl;
    }

    if (this->ntax != someAlignments.ntax){
        cout << "Alignments::operator +=: this->ntax != someAlignments.get_ntax()" << endl;
        exit(0);
    }
    char **new_alignments;
    new_alignments = new char*[this->nchar + someAlignments.nchar];
    for (int i=0;i<this->nchar + someAlignments.nchar;i++){
        new_alignments[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            if (i < this->nchar){
                new_alignments[i][j] = this->alignments[i][j];
            }
            else {
                new_alignments[i][j] = someAlignments(i-this->nchar,j);
            }
        }
    }
    Alignments newAlignments(this->ntax,this->nchar + someAlignments.nchar, this->AlignmentsName, this->taxaNames,new_alignments,1);

    return newAlignments;
}

Alignments & Alignments::operator += (const Alignments &someAlignments)
{
    //if (this->ntax != someAlignments.get_ntax()){
    if (this->ntax != someAlignments.ntax){
        cout << "Alignments::operator +=: this->ntax != someAlignments.get_ntax()" << endl;
        exit(0);
    }
    char **new_alignments;
    new_alignments = new char*[this->nchar + someAlignments.nchar];
    //for (int i=0;i<this->nchar + someAlignments.get_nchar();i++){
    for (int i=0;i<this->nchar + someAlignments.nchar;i++){
        new_alignments[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            if (i < this->nchar){
                new_alignments[i][j] = this->alignments[i][j];
            }
            else {
                new_alignments[i][j] = someAlignments(i-this->nchar,j);
            }
        }
    }

    this->deleteAlignments();
    this->alignments = new_alignments;
    this->alignmentsAllocated = 1;
    this->nchar = this->nchar + someAlignments.nchar;

    // We need to concat someAlignents to this object
    return *this;
}

char Alignments::operator() (unsigned row, unsigned col) const
{
    return alignments[row][col];
}

char & Alignments::operator() (unsigned row, unsigned col)
{
    return alignments[row][col];
}

const int Alignments::get_ntax()
{
    return ntax;
}

const int Alignments::get_nchar()
{
    return nchar;
}

const list <string> Alignments::get_taxaNames()
{
    return taxaNames;
}

const string Alignments::get_AlignmentsName()
{   
    return AlignmentsName;
}

void Alignments::setOutputFormat(int format)
{
    outputFormat = format;
}

void Alignments::setInputFormat(int format)
{
    inputFormat = format;
}

Alignments Alignments::getContiguousColumns(int pos, int k)
{
    if (pos >= nchar || k > nchar) {
        cout << "Alignments::getContiguousColumns: pos >= nchar || k > nchar" << endl;
        exit(0);
    }
    char **new_alignments;
    new_alignments = new char*[k];
    for (int i=0;i<k;i++){
        new_alignments[i] = new char[ntax];
        for (int j=0;j<ntax;j++){
            new_alignments[i][j] = this->alignments[(pos + i) % nchar][j];
        }
    }
    Alignments newAlignments(this->ntax,k, this->AlignmentsName, this->taxaNames,new_alignments,1);

    return newAlignments;
}

Alignments Alignments::delColumnsWithMissingData()
{ 
    char **newAlignments;
    int new_nchar = 0;
    int delColumn[nchar];

    for(int i=0;i<nchar;i++){
        delColumn[i]=0;
    }

    // Go through every column
    for(int i=0;i<nchar;i++){
        //cout << "Checking column " << i << endl;
        for (int j=0;j<ntax;j++){
            if (alignments[i][j] != 'A' && alignments[i][j] != 'T' && alignments[i][j] != 'C' && alignments[i][j] != 'G'){
                //cout << "       Bad character alignments[" << i << "][" << j << "] = " << alignments[i][j] << endl;
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

    newAlignments = new char*[new_nchar];
    for(int i=0;i<new_nchar;i++){
        newAlignments[i] = new char[ntax];
    }
    int tempCount = 0;

    for(int i=0;i<nchar;i++){
        if (delColumn[i] == 0){
            for (int j=0;j<ntax;j++){
                newAlignments[tempCount][j] = alignments[i][j];
            }

            tempCount++;
        }
    }
    
    Alignments returnAlignments(this->ntax,new_nchar, this->AlignmentsName, this->taxaNames,newAlignments,1);

    return returnAlignments;

}

Alignments Alignments::getJackknife(int k)
{
    return getJackknife(k,(int)floor((double)nchar/k));

}

Alignments Alignments::getJackknife(int k, int count)
{
    Alignments returnAlignments;

    if (count > 0){
        returnAlignments = this->getContiguousColumns(rand() % this->nchar,k);

        for (int i=1;i<count;i++){
            returnAlignments += this->getContiguousColumns(rand() % this->nchar,k);
        }
    }
    return returnAlignments;

}

Alignments Alignments::getTaxaSubset(set <unsigned> someTaxa)
{
    if (DEBUG_OUTPUT > 1 ){
        cout << "Alignments::getTaxaSubset called." << endl;
        for (set <unsigned>::const_iterator sit=someTaxa.begin();sit != someTaxa.end();sit++) {
            cout << *sit << " ";
        }
        cout << endl;
    }

    char **new_alignments;
    new_alignments = new char*[nchar];
    for (int i=0;i<nchar;i++){
        new_alignments[i] = new char[someTaxa.size()];
        unsigned rowCount = 0;
        for (set <unsigned>::const_iterator sit=someTaxa.begin();sit!=someTaxa.end();sit++){
            new_alignments[i][rowCount] = this->alignments[i][*sit];
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
    Alignments newAlignments(someTaxa.size(),nchar, this->AlignmentsName, newTaxaNames,new_alignments,1);

    if (DEBUG_OUTPUT > 1 ){
        cout << "Alignments::getTaxaSubset done." << endl;
    }
    return newAlignments;
}

double Alignments::sequenceDivergenceMin()
{
    if (alignmentsAllocated == 0) {
        cout << "Alignments::sequenceDivergence called with no alignments allocated." << endl;
        exit (1);
    }
    int minHamming = -1;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignments[][i] and alignments[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignments[k][i] != alignments[k][j]) {
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

double Alignments::sequenceDivergenceAvg()
{
    if (alignmentsAllocated == 0) {
        cout << "Alignments::sequenceDivergence called with no alignments allocated." << endl;
        exit (1);
    }
    int totalHamming = 0;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignments[][i] and alignments[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignments[k][i] != alignments[k][j]) {
                    currHamming++;
                }
            }
            totalHamming += currHamming;
        }
    }
    double tmpDouble = (double)(ntax-1)*(double)ntax/2.0;
    
    return (double)(((double)totalHamming)/(double)((double)nchar*tmpDouble));
}


void Alignments::printSequenceDivergencePairs()
{
    if (alignmentsAllocated == 0) {
        cout << "Alignments::sequenceDivergence called with no alignments allocated." << endl;
        exit (1);
    }
    int totalHamming = 0;
    int currHamming = -1;
    for(int i=0;i<ntax;i++){
        for (int j=i+1;j<ntax;j++){
            //Compute the hamming distance of alignments[][i] and alignments[][j]
            currHamming=0;
            for (int k=0;k<nchar;k++){
                if (alignments[k][i] != alignments[k][j]) {
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
