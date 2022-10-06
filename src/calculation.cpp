#include "dataTypes.h"
#include "functions.h"
#include "getOrbitals.h"


void calcOrbitalsCurrent(vector <MolecOrb> & orbitalsComplex, vector <MolecOrb> & orbitalsIsolated,
                         const dmat & Overlap, const inputs & variables, vector <orbitalCurrentResult> & orbitalsVCResult)
{
    printf(" ************************* Entering V-I data calculation module ************************* \n \n");


//----------------- Define and set variables and input data ----------------------------//

    double gamma1, gammaN, VLumo, d1, dN,  FermiE, T;
    vector<WindowOrbital> windowOrbitals;

    T = variables.T;
    FermiE = getFermiEnrg(orbitalsComplex);
    printf(" Fermi energy: %f \n \n", FermiE);

    printf(" calculation d1 expansion coefficient: \n");
    d1 = dExpansion(variables.SN_1, orbitalsIsolated, variables.DEG_THRH, T);
    printf(" d1 expansion coefficient = %f \n \n", d1);

    printf(" calculation dN expansion coefficient: \n");
    dN = dExpansion(variables.SN_N, orbitalsIsolated, variables.DEG_THRH, T);
    printf(" dN expansion coefficient = %f \n \n", dN);

    VLumo = V(orbitalsComplex, orbitalsIsolated, variables.HOMO_RES);
    gammaN = VLumo * dN;
    gamma1 = VLumo * d1;
    printf(" V(LUMO): %f \n gamma1 coupling coefficient: %f \n gammaN coupling coefficient: %f \n \n \n", VLumo, gamma1, gammaN);

    getOrbitalWindow(orbitalsComplex, windowOrbitals, variables.ORB_WINDOW);

//----------------- ************************************* ----------------------------//





//----------------- For each orbital in window calculate voltage-current result ----------------------------//

    for_each(windowOrbitals.begin(), windowOrbitals.end(), [&gamma1, &gammaN, &FermiE, &T, &variables, &Overlap, &orbitalsVCResult] (WindowOrbital wMO)
    {
        orbitalCurrentResult orbitalVC;
        orbitalVC.absPosition = wMO.absPosition;
        orbitalVC.relPosition = wMO.relPosition;
        double site1MO = siteOrbProjection(variables.SN_1, wMO.MO, Overlap);
        double siteMON = siteOrbProjection(variables.SN_N, wMO.MO, Overlap);
        double EnergyMO = wMO.MO.mOrbEigenval * au2eV;

        //printf(" <Site|Orbital> coefficient: site1MO = %f, siteMON = %f \n MO energy %f eV \n", site1MO, siteMON, EnergyMO);


        //--------------------------------- generate voltage range ---------------------------------//
        dvec Urange = linspace(variables.U_LWR, variables.U_UPR, variables.N_PTS);
        //--------------------------------- ********************** ---------------------------------//




        //--------------------- calculate V-I characteristics for current MO ----------------------//
        for_each(Urange.begin(), Urange.end(), [&gamma1, &gammaN, &site1MO, &siteMON, &EnergyMO, &FermiE, &T, &orbitalVC, &variables] (const double V)
        {

            //----------------- set variables and data structure -----------------//
            const double E = V;
            integrParams integrVariables;

            integrVariables.Vd = E;
            integrVariables.T = T;
            integrVariables.FermiE = FermiE;
            integrVariables.EnergyMO = EnergyMO;
            integrVariables.site1MO = site1MO;
            integrVariables.siteMON = siteMON;
            integrVariables.gamma1 = gamma1;
            integrVariables.gammaN = gammaN;

            integrResult viPoint;
            viPoint.voltage = V;
            //----------------- ****************************** -----------------//



            //------------------------ Integrate it!!! ------------------------//
            integrateRoutine(E, integrVariables, viPoint, variables.INTEGRATOR);
            //------------------------ ************** ------------------------//



             //----------------- scale resulting value by pre integral factor ---------------//
             viPoint.current *= 0.5 * pow((81/16), (1/3))*((q_e * k_B * T)/(4 * Pi * h_ * FermiE));
             viPoint.curr_error *= 0.5 * pow((81/16), (1/3))*((q_e * k_B * T)/(4 * Pi * h_ * FermiE));
             //---------------------------------------------------------------------------//

            //--------------- calculate transition probability ---------------//

            viPoint.trans_prob = TransitionProb(gamma1, gammaN, site1MO, siteMON, EnergyMO, E);

            //----------------------------------------------------------------//



            orbitalVC.viPoints.push_back(viPoint);
        });
        //--------------------------------- ********************** ---------------------------------//




        printf(" Orbital number:    %i   (  %s  ) \n", orbitalVC.absPosition, orbitalVC.relPosition.c_str());
        printf(" MO energy: %f eV \n <Site|Orbital> coefficients: site1MO = %f, siteMON = %f \n", EnergyMO, site1MO, siteMON);
        printf(" ------------------------------------------------------------------- \n");
        printf(" Voltage (V)   Current (A)   Current error    Transition probability \n");
        printf(" ------------------------------------------------------------------- \n");
        for_each(orbitalVC.viPoints.begin(), orbitalVC.viPoints.end(), [] (integrResult point)
        {
            printf(" %7.3f %15.2e %15.2e %18.2e \n", point.voltage, point.current, point.curr_error, point.trans_prob);
        });
        printf(" ------------------------------------------------------------------- \n \n");




        orbitalsVCResult.push_back(orbitalVC);
    });

//----------------- ********************************************************** ----------------------------//

    printf(" *************************    Out of V-I data calculation module    ************************* \n");
}


