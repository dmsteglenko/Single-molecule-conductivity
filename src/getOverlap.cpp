#include "getOverlap.h"
#include "dataTypes.h"



size_t getNBasis(const char * inFileName)
{
    ifstream input;
    input.open(inFileName, ios::in);

    char buffer[256];
    size_t NBasis;
    do
    {
        input.getline(buffer, 256);
    }
    while(strstr(buffer, "NBasis =") == NULL);
    sscanf(buffer, "%*s = %i", &NBasis);

    input.close();

    return NBasis;
}



size_t getNAtoms(const char * inFileName)
{
    ifstream input;
    input.open(inFileName, ios::in);

    char buffer[256];
    size_t NAtoms;
    do
    {
        input.getline(buffer, 256);
    }
    while(strstr(buffer, "NAtoms=") == NULL);
    sscanf(buffer, "%*s %i", &NAtoms);

    input.close();

    return NAtoms;
}



void getOverlapMat(const char * inFileName, dmat & A)
{
    ifstream input;
    input.open(inFileName, ios::in);

    size_t m, col1, col2, col3, col4, col5;
    string sbuffer;
    do
    {
        getline(input, sbuffer);
    }
    while(strstr(sbuffer.c_str(), "*** Overlap ***") == NULL);

    getline(input, sbuffer);

    do
    {
        replace(sbuffer.begin(), sbuffer.end(), 'D', 'E');
        if(sbuffer[6] == ' ')
        {
            // reads header line numbering of overlap matrix
            switch (sbuffer.length())
            {
            case 17:
                sscanf(sbuffer.c_str(), "%i", &col1);
                break;
            case 31:
                sscanf(sbuffer.c_str(), "%i %i", &col1, &col2);
                break;
            case 45:
                sscanf(sbuffer.c_str(), "%i %i %i", &col1, &col2, &col3);
                break;
            case 59:
                sscanf(sbuffer.c_str(), "%i %i %i %i", &col1, &col2, &col3, &col4);
                break;
            default :
                sscanf(sbuffer.c_str(), "%i %i %i %i %i", &col1, &col2, &col3, &col4, &col5);
                break;
            }
        }
        else
        {
            sscanf(sbuffer.c_str(), "%i %*s", &m);
            // reads overlap matrix columns
            switch (sbuffer.length())
            {
            case 21:
                sscanf(sbuffer.c_str(), "%*i %lf", &A(m-1,col1-1));
                break;
            case 35:
                sscanf(sbuffer.c_str(), "%*i %lf %lf", &A(m-1,col1-1), &A(m-1,col2-1));
                break;
            case 49:
                sscanf(sbuffer.c_str(), "%*i %lf %lf %lf", &A(m-1,col1-1), &A(m-1,col2-1), &A(m-1,col3-1));
                break;
            case 63:
                sscanf(sbuffer.c_str(), "%*i %lf %lf %lf %lf", &A(m-1,col1-1), &A(m-1,col2-1), &A(m-1,col3-1), &A(m-1,col4-1));
                break;
            default :
                sscanf(sbuffer.c_str(), "%*i %lf %lf %lf %lf %lf", &A(m-1,col1-1), &A(m-1,col2-1), &A(m-1,col3-1), &A(m-1,col4-1), &A(m-1,col5-1));
                break;
            }
        }
        getline(input, sbuffer);
    }
    while(strstr(sbuffer.c_str(), "*** Kinetic Energy ***") == NULL);

    A = symmatl(A);

    input.close();
}
