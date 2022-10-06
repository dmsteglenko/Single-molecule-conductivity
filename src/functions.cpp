#include "dataTypes.h"
#include "getOrbitals.h"




double siteOrbProjection(const size_t siteNum, MolecOrb & Orbital, const dmat & OverlapMat)
{
    vector<double> aVec, mVec;
    size_t k;
    k = 0;
    double result = 0;

//---------------------------Construct entire MO from set of coeff of each Atom object---------------------------//
    for_each(Orbital.Atoms.begin(), Orbital.Atoms.end(), [&mVec] (Atom & atom)
    {
        mVec.insert(mVec.end(), atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end());
    });
//---------------------------**************************************************---------------------------//


//-------------Construct MO represented of atom with siteNum from set of coeff of Atom object-------------//
    for_each(Orbital.Atoms.begin(), Orbital.Atoms.end(), [&aVec, &siteNum,  &k] (Atom & atom)
    {
        if ( k == siteNum - 1 )
        {
            aVec.insert(aVec.end(), atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end());
        }
        else
        {
            aVec.insert(aVec.end(), distance(atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end()), 0.0);
        };
        k++;
    });
//---------------------------**************************************************---------------------------//






//----------------------------------------Calculation of the result---------------------------------------//
    k = 0;
    for_each(aVec.begin(), aVec.end(), [&k, &mVec, &OverlapMat, &result] (double c)
    {
        size_t p = 0;
        for_each(mVec.begin(), mVec.end(), [&k, &p, &c, &OverlapMat, &result] (double q)
        {
            result += c * q * OverlapMat(k,p);
            p++;
        });
        k++;
    });
//---------------------------**************************************************---------------------------//

    return result;
}




double dExpansion(const size_t siteNum, vector<MolecOrb> Orbitals, const double degThrh = 0, const double T = 298.15)
{

    double lumoEnergy = 0;
    double Bksum = 0;
    vector<MolecOrb> degenerateOrbitals;
    vector<double> dk;


    vector<MolecOrb>::iterator itLumo = find_if(Orbitals.begin(), Orbitals.end(), [] (MolecOrb & MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });

    lumoEnergy = (*itLumo).mOrbEigenval;

    while (  ((*itLumo).mOrbEigenval - lumoEnergy)*au2eV <= degThrh  )
    {
        degenerateOrbitals.push_back(*itLumo);
        itLumo++;
    }


    printf(" Number of degenerate orbitals of isolated molecule: %i \n", degenerateOrbitals.size()); //Ok!
    for_each(degenerateOrbitals.begin(), degenerateOrbitals.end(), [ &Orbitals ] (MolecOrb MO)
    {
        printf(" Orbital number: %i, ( %s ) \n", getMOAbsPosition(Orbitals, MO), getMORelPosition(Orbitals, MO).c_str());
    });


    for_each(degenerateOrbitals.begin(), degenerateOrbitals.end(), [& dk, &siteNum] (MolecOrb & MO)
    {
        vector<double> mVec, aVec;
        double sumC2aVec, sumC2mVec;
        sumC2aVec = sumC2mVec = 0;

        for_each(MO.Atoms.begin(), MO.Atoms.end(), [&mVec] (Atom & atom)
        {
            mVec.insert(mVec.end(), atom.aOrbCoeffs.begin(), atom.aOrbCoeffs.end());
        });

        for_each(MO.Atoms[siteNum - 1].aOrbCoeffs.begin(), MO.Atoms[siteNum - 1].aOrbCoeffs.end(), [&sumC2aVec] (double c)
        {
            sumC2aVec += c * c;
        });

        for_each(mVec.begin(), mVec.end(), [&sumC2mVec] (double c)
        {
            sumC2mVec += c * c;
        });

        dk.push_back(sqrt(sumC2aVec/sumC2mVec));
    });


    transform(degenerateOrbitals.begin(), degenerateOrbitals.end(), dk.begin(), dk.begin(), [&Bksum, &lumoEnergy, &T] (MolecOrb & MO, double d)
    {
        Bksum += exp(-(au2eV * (MO.mOrbEigenval - lumoEnergy)) / (k_B * T));
        return exp(-(au2eV * (MO.mOrbEigenval - lumoEnergy)) / (k_B * T)) * d;
    } );


    return ( accumulate(dk.begin(), dk.end(), 0.0) / Bksum );
}




double V (vector <MolecOrb> orbitalsComplex, vector <MolecOrb> orbitalsIsolated, const double homoReserv)
{
    double extHomoEnrg, extLumoEnrg, isoLumoEnrg, dEnHLext, dEnLumoIsol;

    vector<MolecOrb>::iterator itLumoComplex = find_if(orbitalsComplex.begin(), orbitalsComplex.end(), [] (MolecOrb & MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });

    extLumoEnrg = (*itLumoComplex).mOrbEigenval;
    extHomoEnrg = (*(--itLumoComplex)).mOrbEigenval;

    vector<MolecOrb>::iterator itLumoIsolated = find_if(orbitalsIsolated.begin(), orbitalsIsolated.end(), [] (MolecOrb & MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });
    isoLumoEnrg = (*itLumoIsolated).mOrbEigenval;


    //printf("LUMO extended    energy  %f \nHOMO extended energy  %f \n", extLumoEnrg, extHomoEnrg);
    //printf("LUMO isolated    energy  %f \nHOMO reservoir energy  %f \n", isoLumoEnrg, homoReserv);

    dEnHLext = abs(extLumoEnrg - extHomoEnrg);
    dEnLumoIsol = abs(homoReserv - isoLumoEnrg);


    return sqrt(abs( (dEnHLext - dEnLumoIsol) * dEnLumoIsol)/2);
}



double getFermiEnrg(vector <MolecOrb> Orbitals)
{
    double HomoEnrg, LumoEnrg;

    vector<MolecOrb> :: iterator itLumo = find_if(Orbitals.begin(), Orbitals.end(), [] (MolecOrb MO)
    {
        return strstr(MO.mOrbSpace.c_str(), "V");
    });

    LumoEnrg = (*itLumo).mOrbEigenval;
    HomoEnrg = (*(--itLumo)).mOrbEigenval;

    return -(au2eV*(LumoEnrg + HomoEnrg)) / 2;
}



double f0 (double t, void * params)
{
    // Here i defined new integrand function derived by variable substitution t = 1/(x + 1 -Vd) //
    // Such transformation allow us to use standard adaptive integration module QAG or CQUAD onto (0,1] interval //

    integrParams funcVariables = *(integrParams *) params;


    const double EnergyMO = funcVariables.EnergyMO;
    const double FermiE = funcVariables.FermiE;
    const double Vd = funcVariables.Vd;
    const double site1MO = funcVariables.site1MO;
    const double siteMON = funcVariables.siteMON;
    const double gamma1 = funcVariables.gamma1;
    const double gammaN = funcVariables.gammaN;
    const double T = funcVariables.T;



    double numerator_f0_1 = pow(gamma1*gammaN*site1MO*siteMON,2) * ( log(1 + exp((FermiE + Vd - (1 - t + t*Vd)/t)/(k_B*T))) - log(1 + exp((FermiE - (1 - t + t*Vd)/t)/(k_B*T))) );

    double numerator_f0_2 = -(pow((1 - t + t*Vd), 2)/pow(t,2)) - (2*(1 - t + t*Vd)/t) + (2*Vd*(1 - t + t*Vd)/t) - 1 + 2*Vd - pow(Vd,2);


    double denominator_f0 = pow(-EnergyMO - ((1 - t + t*Vd)/t), 2) + pow((gamma1*site1MO + gammaN*siteMON)/2, 2);


    return -(numerator_f0_1 * numerator_f0_2)/denominator_f0;

}



double f1 (double t, void * params)
{
    // Here i defined new integrand function derived by variable substitution x = Vd + (1-t)/t //
    // Such transformation allow us to use standard adaptive integration QAG or CQUAD onto (0,1] interval //
    // More slowly converges than f0 //

    integrParams funcVariables = *(integrParams *) params;


    const double EnergyMO = funcVariables.EnergyMO;
    const double FermiE = funcVariables.FermiE;
    const double Vd = funcVariables.Vd;
    const double site1MO = funcVariables.site1MO;
    const double siteMON = funcVariables.siteMON;
    const double gamma1 = funcVariables.gamma1;
    const double gammaN = funcVariables.gammaN;
    const double T = funcVariables.T;


    double numerator_f1 = 4*pow(gamma1*gammaN*site1MO*siteMON, 2) * (log(1 + exp((FermiE*t - Vd*t - 1 + t)/(k_B * T * t))) - log(1 + exp((FermiE * t - 1 + t)/(k_B * T * t))));

    double denominator_f1 = 4*pow(EnergyMO*t, 2) + 8*EnergyMO*pow(t,2)*Vd + 8*EnergyMO*t - 8*EnergyMO*pow(t,2) + 4*pow(Vd*t,2) + 8*Vd*t - 8*Vd*pow(t,2) + 4 - 8*t + 4*pow(t,2) + pow(t*gamma1*site1MO,2) + 2*pow(t,2)*gamma1*site1MO*gammaN*siteMON + pow(t*gammaN*siteMON,2);



    return -(numerator_f1)/denominator_f1;
}


double f2 (double x, void * params)
{
    integrParams funcVariables = *(integrParams *) params;


    const double EnergyMO = funcVariables.EnergyMO;
    const double FermiE = funcVariables.FermiE;
    const double Vd = funcVariables.Vd;
    const double site1MO = funcVariables.site1MO;
    const double siteMON = funcVariables.siteMON;
    const double gamma1 = funcVariables.gamma1;
    const double gammaN = funcVariables.gammaN;
    const double T = funcVariables.T;

    const double Gamma = (gamma1 * site1MO + gammaN * siteMON) / 2;

    //--------------------------------------------------------------------------------------------------------------------------------------------//
    // f_sigma = (log(1 + exp((FermiE + Vd - x)/(kB * T))) - log(1 + exp((FermiE - x)/(kB * T))))
    // Tau_sigma = (pow(gamma1, 2) * pow(gammaN, 2) * pow(siteMO1, 2) * pow(siteMON, 2)) / (pow((EnergyMO - x), 2) + pow(Gamma, 2))
    //--------------------------------------------------------------------------------------------------------------------------------------------//

    double f_sigma = (log(1 + exp((FermiE + Vd - x)/(k_B * T))) - log(1 + exp((FermiE - x)/(k_B * T))));
    double Tau_sigma_sq = (pow(gamma1, 2) * pow(gammaN, 2) * pow(site1MO, 2) * pow(siteMON, 2)) / (pow((-EnergyMO - x), 2) + pow(Gamma, 2));


    double fn = Tau_sigma_sq * f_sigma;

    return fn;
}


double TransitionProb(double gamma1, double gammaN, double site1MO, double siteMON, double EnergyMO, double E)
{
    double Gamma = (gamma1 * site1MO + gammaN * siteMON) / 2;
    double Tau_sigma_sq = (pow(gamma1, 2) * pow(gammaN, 2) * pow(site1MO, 2) * pow(siteMON, 2)) / (pow((-EnergyMO - E), 2) + pow(Gamma, 2));
    return Tau_sigma_sq;
}


void integrateRoutine(const double E, integrParams & integrVariables, integrResult & viPoint, size_t routine = 0)
{

    if (routine == 0)
    {
        //------------------------ Integrate it!!! ------------------------//
        gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
        gsl_function F;
        F.function = &f0;
        F.params = &integrVariables;
        gsl_integration_qag(&F, 0, 1, 0, 1e-7, 1000, 6, w, &viPoint.current, &viPoint.curr_error);
        gsl_integration_workspace_free(w);
        //------------------------ ************** ------------------------//
    }
    if (routine == 1)
    {
        //------------------------ Integrate it!!! ------------------------//
        gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(1000);
        gsl_function F;
        F.function = &f0;
        F.params = &integrVariables;
        gsl_integration_cquad(&F, 0, 1, 0, 1e-07, w, &viPoint.current, &viPoint.curr_error, NULL);
        gsl_integration_cquad_workspace_free(w);
        //------------------------ ************** ------------------------//
    }
    if (routine == 2)
    {
        //------------------------ Integrate it!!! ------------------------//
        gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
        gsl_function F;
        F.function = &f2;
        F.params = &integrVariables;
        gsl_integration_qagiu(&F, E, 0, 1e-7, 1000, w, &viPoint.current, &viPoint.curr_error);
        gsl_integration_workspace_free(w);
        //------------------------ ************** ------------------------//
    }

}
