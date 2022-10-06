#include "dataTypes.h"
#include "parser.h"

void getMolecOrbitals(const char * inFileName, vector <MolecOrb> & Orbitals, const size_t NBasis, string spinPolarization = "no")
{
    ifstream input;
    input.open(inFileName, ios::in);

    string startKeyword, endKeyword1, endKeyword2;
    if(spinPolarization == "no")
    {
        startKeyword = "Molecular Orbital Coefficients";
        endKeyword1 = "DENSITY MATRIX";
        endKeyword2 = "Density Matrix";
    }
    else if (spinPolarization == "alfa")
    {
        startKeyword = "Alpha Molecular Orbital Coefficients";
        endKeyword1 = "Beta Molecular Orbital Coefficients";
        endKeyword2 = "Beta Molecular Orbital Coefficients";
    }
    else if (spinPolarization == "beta")
    {
        startKeyword = "Beta Molecular Orbital Coefficients";
        endKeyword1 = "Alpha Density Matrix";
        endKeyword2 = "Alpha Density Matrix";
    }



    string sbuffer;
    do
    {
        getline(input, sbuffer);
    }
    while(strstr(sbuffer.c_str(), startKeyword.c_str()) == NULL);


    int MoPos = input.tellg();


    vector <Atom> AtomsInfo;
    Atom atom;
    char s1[12], s2[12], s3[12], s4[12], s5[12];
    size_t m, aNumber, col1, col2, col3, col4, col5;
    m = aNumber = 0;

//-----Skip lines to get first Atom line into sbuffer-----//
    getline(input, sbuffer);
    getline(input, sbuffer);
    getline(input, sbuffer);
    getline(input, sbuffer);
//-----********************************************-----//


//----------------------Fill AtomsInfo----------------------//
    sbuffer.resize(22);
    do
    {
        if(sscanf(sbuffer.c_str(), "%s %s %s %s", s1, s2, s3, s4) == 4)
        {
            atom.aOrbType.clear();
            atom.atomN = atoi(s2);
            atom.atomSymb = s3;
            atom.aOrbType.push_back(s4);
            AtomsInfo.push_back(atom);
        }
        else if (sscanf(sbuffer.c_str(), "%s %6[0-9A-Z\+ \-] %s %s", s1, s2, s3, s4) == 2)
        {
            AtomsInfo[AtomsInfo.size()-1].aOrbType.push_back(s2);
        }
        getline(input, sbuffer);
        sbuffer.resize(22);
    }
    while(sscanf(sbuffer.c_str(), "%s %s %s %s", s1, s2, s3, s4) != -1);
//-----------------------********************************************-----------------------//



//-------------------Initialize vector of molecular orbitals----------------------//
    Orbitals.resize(NBasis);

    struct
    {
        vector <Atom> AtomsInfo;
        void operator () (MolecOrb & MO)
        {
            MO.mOrbEigenval = 0;
            MO.Atoms = AtomsInfo;

            struct
            {
                void operator () (Atom & atom)
                {
                    atom.aOrbCoeffs.zeros(atom.aOrbType.size());
                }
            } initAtom;

            for_each(MO.Atoms.begin(), MO.Atoms.end(), initAtom);
        }
    } initMO;
    initMO.AtomsInfo = AtomsInfo;

    for_each(Orbitals.begin(),Orbitals.end(), initMO);
//-------------------**************************************----------------------//


    input.seekg(MoPos);
    getline(input, sbuffer);

    do
    {
        size_t counter = 1;

        while (counter <= NBasis + 3)
        {
            if(counter == 1)
            {
                sscanf(sbuffer.substr(20).c_str(), "%i %i %i %i %i", &col1, &col2, &col3, &col4, &col5);
            }
            else if (counter == 2)
            {
                switch (sscanf(sbuffer.substr(20).c_str(), "%s %s %s %s %s", s1, s2, s3, s4, s5))
                {
                case 1:
                    Orbitals[col1-1].mOrbSpace = s1;
                    break;
                case 2:
                    Orbitals[col1-1].mOrbSpace = s1;
                    Orbitals[col2-1].mOrbSpace = s2;
                    break;
                case 3:
                    Orbitals[col1-1].mOrbSpace = s1;
                    Orbitals[col2-1].mOrbSpace = s2;
                    Orbitals[col3-1].mOrbSpace = s3;
                    break;
                case 4:
                    Orbitals[col1-1].mOrbSpace = s1;
                    Orbitals[col2-1].mOrbSpace = s2;
                    Orbitals[col3-1].mOrbSpace = s3;
                    Orbitals[col4-1].mOrbSpace = s4;
                    break;
                default:
                    Orbitals[col1-1].mOrbSpace = s1;
                    Orbitals[col2-1].mOrbSpace = s2;
                    Orbitals[col3-1].mOrbSpace = s3;
                    Orbitals[col4-1].mOrbSpace = s4;
                    Orbitals[col5-1].mOrbSpace = s5;
                    break;
                }
            }
            else if (counter == 3)
            {
                switch (sscanf(sbuffer.substr(20).c_str(), "%s %s %s %s %s", s1, s2, s3, s4, s5))
                {
                case 1:
                    Orbitals[col1-1].mOrbEigenval = strtod(s1, NULL);
                    break;
                case 2:
                    Orbitals[col1-1].mOrbEigenval = strtod(s1, NULL);
                    Orbitals[col2-1].mOrbEigenval = strtod(s2, NULL);
                    break;
                case 3:
                    Orbitals[col1-1].mOrbEigenval = strtod(s1, NULL);
                    Orbitals[col2-1].mOrbEigenval = strtod(s2, NULL);
                    Orbitals[col3-1].mOrbEigenval = strtod(s3, NULL);
                    break;
                case 4:
                    Orbitals[col1-1].mOrbEigenval = strtod(s1, NULL);
                    Orbitals[col2-1].mOrbEigenval = strtod(s2, NULL);
                    Orbitals[col3-1].mOrbEigenval = strtod(s3, NULL);
                    Orbitals[col4-1].mOrbEigenval = strtod(s4, NULL);
                    break;
                default:
                    Orbitals[col1-1].mOrbEigenval = strtod(s1, NULL);
                    Orbitals[col2-1].mOrbEigenval = strtod(s2, NULL);
                    Orbitals[col3-1].mOrbEigenval = strtod(s3, NULL);
                    Orbitals[col4-1].mOrbEigenval = strtod(s4, NULL);
                    Orbitals[col5-1].mOrbEigenval = strtod(s5, NULL);
                    break;
                }
            }
            else if (counter > 3)
            {
                if (sscanf(sbuffer.substr(0, 20).c_str(), "%s %s %s %s", s1, s2, s3, s4) == 4)
                {
                    m = 0;
                    aNumber = atoi(s2);
                }

                switch (sscanf(sbuffer.substr(20).c_str(), "%s %s %s %s %s", s1, s2, s3, s4, s5))
                {
                case 1:
                    Orbitals[col1-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s1, NULL);
                    break;
                case 2:
                    Orbitals[col1-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s1, NULL);
                    Orbitals[col2-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s2, NULL);
                    break;
                case 3:
                    Orbitals[col1-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s1, NULL);
                    Orbitals[col2-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s2, NULL);
                    Orbitals[col3-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s3, NULL);
                    break;
                case 4:
                    Orbitals[col1-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s1, NULL);
                    Orbitals[col2-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s2, NULL);
                    Orbitals[col3-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s3, NULL);
                    Orbitals[col4-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s4, NULL);
                    break;
                default:
                    Orbitals[col1-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s1, NULL);
                    Orbitals[col2-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s2, NULL);
                    Orbitals[col3-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s3, NULL);
                    Orbitals[col4-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s4, NULL);
                    Orbitals[col5-1].Atoms[aNumber-1].aOrbCoeffs[m] = strtod(s5, NULL);
                    break;
                }
                // print number AO readed for each atom// cout << "  " << m << "  " << endl;
                m++;
            }

            counter++;
            getline(input, sbuffer);
        }
    }
    while((strstr(sbuffer.c_str(), endKeyword1.c_str()) || strstr(sbuffer.c_str(), endKeyword2.c_str())) == NULL);
    input.close();
}






void printMO(vector <MolecOrb> & MOs, const size_t Nbas, const size_t numberOfColumns = 5, const size_t columnWidth = 20)
{
    size_t MONum = 0;
    size_t MOSize = 0;

    vector<string> sbufferMOs;
    sbufferMOs.insert(sbufferMOs.end(), Nbas + 3, "");

    MOSize = MOs.size();


    for_each(MOs.begin(), MOs.end(), [&sbufferMOs, columnWidth, &MONum, &MOSize, &Nbas, &numberOfColumns] (MolecOrb &MO)
    {

        vector <string> sbuffer;
        char  cbufferMONum[columnWidth];
        char  cbufferMOSpace[columnWidth];
        char  cbufferMOEigenVal[columnWidth];

        MONum++;

        if (MONum % numberOfColumns != 0)
        {
            //---------------------------------- compose and fill first 3 lines of MO ----------------------------------//
            sprintf(cbufferMONum, "%*s%*s", columnWidth/2+to_string(MONum).size()/2, to_string(MONum).c_str(), columnWidth/2-to_string(MONum).size()/2, "");
            sprintf(cbufferMOSpace, "%*s%*s", columnWidth/2+MO.mOrbSpace.size()/2, MO.mOrbSpace.c_str(), columnWidth/2-MO.mOrbSpace.size()/2, "");
            sprintf(cbufferMOEigenVal, "%*s%*s", columnWidth/2+to_string(MO.mOrbEigenval).size()/2, to_string(MO.mOrbEigenval).c_str(), columnWidth/2-to_string(MO.mOrbEigenval).size()/2, "");

            sbuffer.push_back(cbufferMONum);
            sbuffer.push_back(cbufferMOSpace);
            sbuffer.push_back(cbufferMOEigenVal);
            //---------------------------- ************************************************ ----------------------------//


            //---------------------------------- compose and fill coefficients of MO ----------------------------------//
            for_each(MO.Atoms.begin(), MO.Atoms.end(), [&sbuffer, columnWidth] (Atom & atom)
            {
                for_each(atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end(), [&sbuffer, columnWidth] (double val)
                {
                    char  cbuffer[columnWidth];
                    sprintf(cbuffer, "%*s%*s", columnWidth/2+to_string(val).size()/2, to_string(val).c_str(), columnWidth/2-to_string(val).size()/2, "");
                    sbuffer.push_back(cbuffer);
                });
            });
            //---------------------------- ************************************************ ----------------------------//


            //-------------------- transform resulting empty "Table of strings" into resulting ones --------------------//
            transform(sbuffer.begin(), sbuffer.end(), sbufferMOs.begin(), sbufferMOs.begin(), [] (string &s1, string &s2)
            {
                return s2.append(s1);
            });
            //---------------------------- ************************************************ ----------------------------//

        }
        else
        {
            //---------------------------------- compose and fill first 3 lines of MO ----------------------------------//
            sprintf(cbufferMONum, "%*s%*s", columnWidth/2+to_string(MONum).size()/2, to_string(MONum).c_str(), columnWidth/2-to_string(MONum).size()/2, "");
            sprintf(cbufferMOSpace, "%*s%*s", columnWidth/2+MO.mOrbSpace.size()/2, MO.mOrbSpace.c_str(), columnWidth/2-MO.mOrbSpace.size()/2, "");
            sprintf(cbufferMOEigenVal, "%*s%*s", columnWidth/2+to_string(MO.mOrbEigenval).size()/2, to_string(MO.mOrbEigenval).c_str(), columnWidth/2-to_string(MO.mOrbEigenval).size()/2, "");

            sbuffer.push_back(cbufferMONum);
            sbuffer.push_back(cbufferMOSpace);
            sbuffer.push_back(cbufferMOEigenVal);
            //---------------------------- ************************************************ ----------------------------//


            //---------------------------------- compose and fill coefficients of MO ----------------------------------//
            for_each(MO.Atoms.begin(), MO.Atoms.end(), [&sbuffer, columnWidth] (Atom & atom)
            {
                for_each(atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end(), [&sbuffer, columnWidth] (double val)
                {
                    char  cbuffer[columnWidth];
                    sprintf(cbuffer, "%*s%*s", columnWidth/2+to_string(val).size()/2, to_string(val).c_str(), columnWidth/2-to_string(val).size()/2, "");
                    sbuffer.push_back(cbuffer);
                });
            });
            //---------------------------- ************************************************ ----------------------------//


            //-------------------- transform resulting empty "Table of strings" into resulting ones --------------------//
            transform(sbuffer.begin(), sbuffer.end(), sbufferMOs.begin(), sbufferMOs.begin(), [] (string &s1, string &s2)
            {
                return s2.append(s1);
            });
            //---------------------------- ************************************************ ----------------------------//



            //------------------------------------------ print chunk result -------------------------------------------//
            for_each(sbufferMOs.begin(), sbufferMOs.end(), [] (string s)
            {
                printf("%s \n", s.c_str());
            });
            //------------------------------------------ ****************** -------------------------------------------//


            //----------------------------- reset sbuffer Table before next portion of MO -----------------------------//
            sbufferMOs.clear();
            sbufferMOs.insert(sbufferMOs.end(), Nbas + 3, "");
            //------------------------------------------ ****************** -------------------------------------------//
        }



        //------------------------------------------ print last chunk of result -------------------------------------------//
        if( (MOSize == MONum) && (MONum % numberOfColumns) != 0)
        {
            for_each(sbufferMOs.begin(), sbufferMOs.end(), [] (string s)
            {
                printf("%s \n", s.c_str());
            });
        }
        //------------------------------------------ ************************* -------------------------------------------//
    });

}




size_t getMOAbsPosition(vector <MolecOrb> & MOs, MolecOrb & MO)
{
    return distance(MOs.begin(), find(MOs.begin(), MOs.end(), MO)) + 1;
}



string getMORelPosition(vector <MolecOrb> & MOs, MolecOrb & MO)
{
    string result;
    vector <MolecOrb> :: iterator itLUMO, itHOMO;
    itLUMO = find_if(MOs.begin(), MOs.end(), [] (MolecOrb & MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });
    itHOMO = itLUMO - 1;

    if (distance(find(MOs.begin(), MOs.end(), MO), itHOMO) > 0)
    {
        result.append("HOMO - ");
        return result.append(to_string(distance(find(MOs.begin(), MOs.end(), MO), itHOMO)));
    }

    if (distance(find(MOs.begin(), MOs.end(), MO), itHOMO) == 0)
    {
        return result.append("HOMO ");
    }

    if (distance(find(MOs.begin(), MOs.end(), MO), itLUMO) == 0)
    {
        return result.append("LUMO ");
    }

    if (distance(itLUMO, find(MOs.begin(), MOs.end(), MO)) > 0)
    {
        result.append("LUMO + ");
        return result.append(to_string(distance(itLUMO, find(MOs.begin(), MOs.end(), MO))));
    }

//-------------- if none of above conditions does not fulfill return empty string --------------//
    return result;
}




void getOrbitalWindow( vector <MolecOrb> &MOs, vector <WindowOrbital> &windowMOs, string range )
{

//--------------------------- define variables an HOMO LUMO iterators ---------------------------//
    vector <string> rangeValues;
    vector <MolecOrb> chunkMOs;
    vector <MolecOrb> :: iterator itLUMO, itHOMO;

    itLUMO = find_if(MOs.begin(), MOs.end(), [] (MolecOrb & MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });
    itHOMO = itLUMO - 1;
//--------------------------- ************************************** ---------------------------//


    stringParse(range, rangeValues, " ,");


    if (rangeValues[0] == "ALL" || rangeValues[0] == "all")
    {
        for_each(MOs.begin(), MOs.end(), [&MOs, &windowMOs] (MolecOrb & MO)
        {
            WindowOrbital wMO;
            wMO.absPosition = getMOAbsPosition(MOs, MO);
            wMO.relPosition = getMORelPosition(MOs, MO);
            wMO.MO = MO;
            windowMOs.push_back(wMO);
        });
    }

    if ( !(rangeValues[0] == "ALL" || rangeValues[0] == "all") && rangeValues.size() == 1 )
    {
        size_t halfRange = atoi(rangeValues[0].c_str());

        chunkMOs.insert(chunkMOs.begin(), itHOMO - halfRange, itLUMO + halfRange + 1);


        for_each(chunkMOs.begin(), chunkMOs.end(), [&MOs, &windowMOs] (MolecOrb & MO)
        {
            WindowOrbital wMO;
            wMO.absPosition = getMOAbsPosition(MOs, MO);
            wMO.relPosition = getMORelPosition(MOs, MO);
            wMO.MO = MO;
            windowMOs.push_back(wMO);
        });
    }

    if (rangeValues.size() == 2)
    {
        size_t lowerBound = atoi(rangeValues[0].c_str());
        size_t upperBound = atoi(rangeValues[1].c_str());

        chunkMOs.insert(chunkMOs.begin(), itHOMO - lowerBound, itLUMO + upperBound + 1);


        for_each(chunkMOs.begin(), chunkMOs.end(), [&MOs, &windowMOs] (MolecOrb & MO)
        {
            WindowOrbital wMO;
            wMO.absPosition = getMOAbsPosition(MOs, MO);
            wMO.relPosition = getMORelPosition(MOs, MO);
            wMO.MO = MO;
            windowMOs.push_back(wMO);
        });
    }

    if (rangeValues.size() > 2)
    {
        for_each(rangeValues.begin(), rangeValues.end(), [&chunkMOs, &MOs] (string s)
        {
            chunkMOs.push_back(MOs[atoi(s.c_str()) - 1]);
        });


        for_each(chunkMOs.begin(), chunkMOs.end(), [&MOs, &windowMOs] (MolecOrb & MO)
        {
            WindowOrbital wMO;
            wMO.absPosition = getMOAbsPosition(MOs, MO);
            wMO.relPosition = getMORelPosition(MOs, MO);
            wMO.MO = MO;
            windowMOs.push_back(wMO);
        });

    }

//------------------end of procedure-----------------//

}
