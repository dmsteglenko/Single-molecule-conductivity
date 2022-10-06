#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED

#include "dataTypes.h"

void  stringParse(string sbuffer, vector<string> &tokens, const char * delim);

void inputFileParser(const char * inputFileName, inputs & variables);

#endif // PARSER_H_INCLUDED
