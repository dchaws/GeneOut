#include "nexus.h"

// This command eats up comments in the istream that start
// with '[' and end with ']'.
// We should handle nested ']', although it is completly 
// stupid that any program would use nested comments
void eatUpComments(std::istream &in,char &c)
{
    int nestCount = 0;
    if (c == '['){ // Check for comment which starts with '['
        if (DEBUG_OUTPUT >= 1) {
            cout << "Comment found: " << c;
        }
        // Read until ']' seen
        if (in.bad()) { return;}
        c = (char)in.get();
        if (c == '['){
            nestCount++;
        }
        if (DEBUG_OUTPUT >= 1) {
            cout << c;
        }
        while (c != ']' || nestCount != 0){
            if (c == ']'){
                nestCount--;
            }
            if (in.bad()) { return;}
            c = (char)in.get();
            if (c == '['){
                nestCount++;
            }
            if (DEBUG_OUTPUT >= 1) {
                cout << c;
            }
        }
        // c == ']', so scan past it
        if (in.bad()) { return;}
        c = (char)in.get();
        if (DEBUG_OUTPUT >= 1) {
            cout << c;
        }
    }
}

int isPunctuation(char c)
{
    if (c == ',' || c == '.' || c == '!' || c == '?' || c == '@' || c == '*' || c == '(' || c == ')' || c == ';' || c == ':' || c == '=' ){ 
        return 1;
    }
    return 0;
}

// This function reads in the next Nexus token from the input
// A Nexux token is defined as a sequence of alphanum characters or punctions,
// e.g. "( ) , . ; : ! ? = "
//Ignore preceeding whitespaces,tabs,newlines,carraige returns. Also treat ';' as a token.
// This should ignore any characters between '[' and ']'
string nextNexusToken(std::istream &in)
{
    char c; 
    string returnString;

    // Read in and ignore whitespaces,tabs,newlines,carraige returns
    if (in.bad()) { return returnString;}
    c = (char)in.get();
    if (c == EOF){
        return returnString;
    }

    while (c == ' ' || c == '\n' || c == '\r' || c == '\t'){
        if (in.bad()) { return returnString;}
        c = (char)in.get();
        eatUpComments(in,c);
    }
    eatUpComments(in,c);
    if (c == EOF){
        return returnString;
    }

    // If the next character is a puncuation return it.
    if (isPunctuation(c)){
        returnString.append(1,c);
        if (DEBUG_OUTPUT >= 1) {
            cout << "nextNexusToken(isPunctuation): " << returnString << endl;
        }
        return returnString;
    }
    
    // Else we must be on a legitimate token
    returnString.append(1,c);
    
    if (in.bad()) { return returnString;}
    c = (char)in.get();
    eatUpComments(in,c);
    // Read until whitespace, tab, newline, carraige return OR ';'
    while (c != ' ' && c != '\n' && c != '\r' && c != '\t' && !isPunctuation(c) && c != EOF){
        eatUpComments(in,c);
        returnString.append(1,c);
        if (in.bad()) { return returnString;}
        c = (char)in.get();
    }
    // We read to far, so put on character back
    in.putback(c);
    eatUpComments(in,c);

    if (DEBUG_OUTPUT >= 1) {
        cout << "nextNexusToken               : " << returnString << endl;
    }

    return returnString;
}

string nextNexusNumber(std::istream &in)
{
    string returnString;
    // Read until first digit seen.
    char c = (char)in.get();
    if (in.bad()) { return returnString;}
    while (!isdigit(c) && c != '.' && c != 'E' && c != '-' && c == '+'){
        if (c == ';') {
            returnString = ';';
            if (DEBUG_OUTPUT >= 1) {
                cout << "nextNexusNumber: " << returnString << endl;
            }
            return returnString;
        }
        if (in.bad()) { return returnString;}
        c = (char)in.get();
    }
    // Now read everything that is a digit, or digit-related alphanum characters.
    while (isdigit(c) || c == '.' || c == 'E' || c == '-' || c == '+'){
        returnString.append(1,c);
        if (in.bad()) { return returnString;}
        c = (char)in.get();
    }
    if (DEBUG_OUTPUT >= 1) {
        cout << "nextNexusNumber: " << returnString << endl;
    }
    // We read to far, so put on character back
    in.putback(c);
    return returnString;

}

string nextNexusTokenUpper(std::istream &in)
{
    string input = nextNexusToken(in);
    transform(input.begin(), input.end(),input.begin(), ::toupper);
    return input;
}

string nextQuoteBlock(std::istream &in)
{
    string returnString;

    // Read until first '"' seen.
    char c = (char)in.get();
    if (in.bad()) { return returnString;}
    while (c != '"'){
        if (c == ';') {
            returnString = ';';
            if (DEBUG_OUTPUT >= 1) {
                cout << "nextQuoteBlock: " << returnString << endl;
            }
            return returnString;
        }
        if (in.bad()) { return returnString;}
        c = (char)in.get();
    }
    // Now read everything into returnString until next '"' is seen.
    c = (char)in.get();
    while (c != '"'){
        returnString.append(1,c);
        if (in.bad()) { return returnString;}
        c = (char)in.get();
    }
    if (DEBUG_OUTPUT >= 1) {
        cout << "nextQuoteBlock: " << returnString << endl;
    }
    return returnString;
}


// This function seeks/reads from the input until the next ; is found
// The input stream will be placed at the next character after ;
// This should ignore any characters between '[' and ']'
void readUntilSemiColonNexus(std::istream &in)
{
    string input;
    input = nextNexusToken(in);
    while (input != ";"){
        input = nextNexusToken(in);
    }

}

int errorCheckNexusToken(string someString,const char *str)
{
    if (someString != str){
        cout << someString << " != " << str << endl;
        exit(0);
    }
    else return 1; // Everything ok
}

