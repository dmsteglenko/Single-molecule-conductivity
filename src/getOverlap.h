#ifndef GETOVERLAP_H_INCLUDED
#define GETOVERLAP_H_INCLUDED

#include "dataTypes.h"

size_t getNBasis(const char * inFileName);

void getOverlapMat(const char * inFileName, dmat & A);

#endif // GETOVERLAP_H_INCLUDED
