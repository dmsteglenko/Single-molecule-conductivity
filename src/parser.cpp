#include "parser.h"
#include "dataTypes.h"



void  stringParse(string sbuffer, vector<string> &tokens, const char * delim = " ")
{
    tokens.clear();
    char * cbuffer = new char [sbuffer.length() + 1];
    strcpy(cbuffer, sbuffer.c_str());
    char * token = strtok(cbuffer, delim);
    while (token != NULL)
    {
        tokens.push_back(token);
        token = strtok(NULL, delim);
    }
    delete [] cbuffer;
}




void inputFileParser(const char * inputFileName, inputs & variables)
{
    ifstream inputFile;
    inputFile.open(inputFileName, ios::in);

    string sbuffer;
    vector<string> tokens;

    do
    {
        getline(inputFile, sbuffer);

        //------------------------------ trim string whitespace and comment - # ------------------------------//
        sbuffer.erase(sbuffer.begin(), find_if(sbuffer.begin(), sbuffer.end(), [](char s)
        {
            return s != ' ';
        }));
        sbuffer.erase(find_if(sbuffer.begin(), sbuffer.end(), [](char s)
        {
            return s == '#';
        }), sbuffer.end());
        //------------------------------------  ************************  -----------------------------------//

        //------------------------------    convert key string to upper case    ------------------------------//
        for_each(sbuffer.begin(), find_if(sbuffer.begin(), sbuffer.end(), [] (char s)
        {
            return s == '=';
        }), [] (char &c)
        {
            c = toupper(c);
        });
        //------------------------------------  ************************  -----------------------------------//


        if (sbuffer.length())
        {
            stringParse(sbuffer, tokens, " =");


            //-------------------------- set values --------------------------//
            if (tokens[0] == "INPUT_FILE_COMPLEX")
            {
                variables.INPUT_FILE_COMPLEX = tokens[1];
            }

            if (tokens[0] == "WF")
            {
                variables.WF = tokens[1];
            }

            if(tokens[0] == "INPUT_FILE_ISOLATED")
            {
                variables.INPUT_FILE_ISOLATED = tokens[1];
            }

            if(tokens[0] == "OUT_FILE")
            {
                variables.OUT_FILE = tokens[1];
            }

            if(tokens[0] == "ORB_WINDOW")
            {
                variables.ORB_WINDOW = tokens[1];
            }

            if(tokens[0] == "HOMO_RES")
            {
                variables.HOMO_RES = strtod(tokens[1].c_str(), nullptr);
            }

            if(tokens[0] == "U_LWR")
            {
                variables.U_LWR = strtod(tokens[1].c_str(), nullptr);
            }

            if(tokens[0] == "U_UPR")
            {
                variables.U_UPR = strtod(tokens[1].c_str(), nullptr);
            }

            if(tokens[0] == "DEG_THRH")
            {
                variables.DEG_THRH = strtod(tokens[1].c_str(), nullptr);
            }

            if(tokens[0] == "T")
            {
                variables.T = strtod(tokens[1].c_str(), nullptr);
            }

            if(tokens[0] == "SN_1")
            {
                variables.SN_1 = atoi(tokens[1].c_str());
            }

            if(tokens[0] == "SN_N")
            {
                variables.SN_N = atoi(tokens[1].c_str());
            }

            if(tokens[0] == "INTEGRATOR")
            {
                variables.INTEGRATOR = atoi(tokens[1].c_str());
            }

            if(tokens[0] == "N_PTS")
            {
                variables.N_PTS = atoi(tokens[1].c_str());
            }
            //-------------------------- ***************** --------------------------//

        }

    }
    while(!inputFile.eof());

    inputFile.close();
}

