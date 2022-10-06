#ifndef DATATYPES_H_INCLUDED
#define DATATYPES_H_INCLUDED

//------- define physical constants and conversion constants -----//
const double au2eV(27.211385);
const double k_B(0.00008617333);
const double q_e(0.00000000000000000016021766);
const double h_(0.0000000000000006582119514);
const double Pi(3.14159265359);
//------- ************************************************* -----//

#include <vector>
#include <armadillo>
#include <gsl/gsl_integration.h>
#include <string>


using namespace std;
using namespace arma;



struct inputs
{
    string INPUT_FILE_COMPLEX = "";
    string INPUT_FILE_ISOLATED = "";
    string OUT_FILE = "";
    string WF = "RHF";
    string ORB_WINDOW = "0";
    double HOMO_RES = 5.45;
    double U_LWR = 0;
    double U_UPR = 5;
    double DEG_THRH = 0.07565;
    double T = 298.15;
    size_t SN_1 = 1;
    size_t SN_N = 1;
    size_t INTEGRATOR = 0;
    size_t N_PTS = 100;
};



struct Atom
{
    size_t atomN = 0;
    string atomSymb = "";
    vector <string> aOrbType;
    dvec aOrbCoeffs;


//------------------- Operation of comparison two Atom entity -------------------//
    bool operator == (Atom atom)
    {
        if (atomN == atom.atomN)
        {
            if (atomSymb == atom.atomSymb)
            {
                if ( equal(aOrbType.begin(), aOrbType.end(), atom.aOrbType.begin()) )
                {
                    if ( equal(aOrbCoeffs.begin(), aOrbCoeffs.end(), atom.aOrbCoeffs.begin()) )
                    {
                        return true;
                    }
                }
            }
        }

        return false;

    }
//------------------- ************************************ -------------------//
};



struct  MolecOrb
{
    string mOrbSpace = "";
    double mOrbEigenval = 0.0;
    vector <Atom> Atoms;


//------------------- Operation of comparison two Molecular orbital entity -------------------//
    bool operator == (MolecOrb MO)
    {
        if (mOrbSpace == MO.mOrbSpace)
        {
            if (mOrbEigenval == MO.mOrbEigenval)
            {
                if (equal(Atoms.begin(), Atoms.end(), MO.Atoms.begin()))
                {
                    return true;
                }
            }
        }

        return false;

    }
//------------------- ************************************************* -------------------//
};



struct WindowOrbital
{
    size_t absPosition;
    string relPosition;
    MolecOrb MO;
};



struct integrParams
{
    double EnergyMO = 0.0;
    double siteMON = 0.0;
    double site1MO = 0.0;
    double gamma1 = 0.0;
    double gammaN = 0.0;
    double Vd = 0.0;
    double FermiE = 0.0;
    double T = 273.15;
};



struct integrResult
{
    double voltage = 0;
    double current = 0;
    double curr_error = 0;
    double trans_prob = 0;
};


struct orbitalCurrentResult
{
    size_t absPosition;
    string relPosition;
    vector <integrResult> viPoints;
};



#endif // DATATYPES_H_INCLUDED
