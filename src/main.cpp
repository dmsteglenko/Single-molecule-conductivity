#include <iostream>
#include "dataTypes.h"
#include "parser.h"
#include "getOverlap.h"
#include "getOrbitals.h"
#include "functions.h"
#include "calculation.h"
#include "print.h"





int main(int argc, char* argv[])
{
    printf(" Entering single molecule conductivity program. It will be calculate Voltage-Current characteristics \n");
    printf(" for single molecule attached by terminal atoms i.e. (site) to electrodes - gold clusters (reservoir). \n \n");

//-------------- Define and set variables and data structure --------------//

    inputs inVariables;
    inputFileParser(argv[1], inVariables);

    printf(" Physical constants used: \n Boltzmann constant  %e eV \n Elementary charge %e C \n Planck constant %e eV*s \n \n", k_B, q_e, h_);

    printf(" The following input variables will be used for calculations: \n");
    printf(" Input file for molecule attached to gold electrodes (ref. as complex): %s \n", inVariables.INPUT_FILE_COMPLEX.c_str());
    printf(" Input file for free standing molecule (ref. as isolated): %s \n", inVariables.INPUT_FILE_ISOLATED.c_str());
    printf(" Out file for printing result: %s \n", inVariables.OUT_FILE.c_str());
    printf(" Size of orbital window for which it will calculate (i.e. orbitals around HOMO-LUMO): %s \n",inVariables.ORB_WINDOW.c_str());
    printf(" Site numbers i.e. atoms attached to gold electrodes: %i,  %i \n", inVariables.SN_1, inVariables.SN_N);
    printf(" Energy of HOMO of gold electrodes (reservoir): %8.4f eV \n", inVariables.HOMO_RES);
    printf(" Energy window within LUMO+k orbitals are degenerate and includes for calculation of dk: %8.4f eV \n",  inVariables.DEG_THRH);
    printf(" Number of points V-I curve for calculation: %i \n", inVariables.N_PTS);
    printf(" Lower an Upper bound of Voltage for calculation: %10.4f V, %10.4f V \n",  inVariables.U_LWR, inVariables.U_UPR);
    printf(" Wave function of molecular system attached to electrodes:  %s \n", inVariables.WF.c_str());
    printf(" Temperature under calculations are carried out:  %8.2f K \n", inVariables.T);

    if(inVariables.INTEGRATOR == 0)
    {
        printf(" Integrator type: QAG adaptive integration procedure with integration corresponding 61 point Gauss-Kronrod rules. \n");
        printf(" It's reasonable compromise between speed and robust. \n \n");
    }
    if(inVariables.INTEGRATOR == 1)
    {
        printf(" Integrator type: CQUAD a new doubly-adaptive general-purpose quadrature routine. \n");
        printf(" It may be more slowly than \"INTEGRATOR = 0\" but it very stable and robust routine. \n \n");
    }

    if(inVariables.INTEGRATOR == 2)
    {
        printf(" Integrator type: QAGIU adaptive integration on infinite upper interval. \n");
        printf(" It holds an original integral equation but frequently fails (integral diverges) on large orbital window. \n \n");
    }



    size_t NBasisComplex = getNBasis(inVariables.INPUT_FILE_COMPLEX.c_str());
    size_t NBasisIsolated = getNBasis(inVariables.INPUT_FILE_ISOLATED.c_str());

    printf(" Number of AO for complex molecule: %i \n Number of AO for isolated molecule: %i \n \n", NBasisComplex, NBasisIsolated);


    dmat overlapMat(NBasisComplex, NBasisComplex, fill::zeros);
    getOverlapMat(inVariables.INPUT_FILE_COMPLEX.c_str(), overlapMat);

    if (inVariables.WF == "RHF" || inVariables.WF == "rhf")
    {
        vector<MolecOrb> orbitalsComplex, orbitalsIsolated;
        getMolecOrbitals(inVariables.INPUT_FILE_COMPLEX.c_str(), orbitalsComplex, NBasisComplex);
        getMolecOrbitals(inVariables.INPUT_FILE_ISOLATED.c_str(), orbitalsIsolated, NBasisIsolated);

        vector<orbitalCurrentResult> currentResults;
        printf("                     ----- Start of closed shell calculations ----- \n");
        calcOrbitalsCurrent(orbitalsComplex, orbitalsIsolated, overlapMat, inVariables, currentResults);
        printf("                        ----- End of closed shell calculations -----                        \n \n");

        printf("                    ----- Total result for closed shell calculations -----                     \n");
        printTotalResult(currentResults);

    }
    else if (inVariables.WF == "UHF" || inVariables.WF == "uhf")
    {
        vector<MolecOrb> orbitalsComplex_alfa, orbitalsComplex_beta, orbitalsIsolated_alfa, orbitalsIsolated_beta;
        getMolecOrbitals(inVariables.INPUT_FILE_COMPLEX.c_str(), orbitalsComplex_alfa, NBasisComplex, "alfa");
        getMolecOrbitals(inVariables.INPUT_FILE_COMPLEX.c_str(), orbitalsComplex_beta, NBasisComplex, "beta");

        getMolecOrbitals(inVariables.INPUT_FILE_ISOLATED.c_str(), orbitalsIsolated_alfa, NBasisIsolated, "alfa");
        getMolecOrbitals(inVariables.INPUT_FILE_ISOLATED.c_str(), orbitalsIsolated_beta, NBasisIsolated, "beta");

        vector<orbitalCurrentResult> currentResults_alfa, currentResults_beta;
        printf("                     ----- Start of open shell calculations: alfa density -----                      \n");
        calcOrbitalsCurrent(orbitalsComplex_alfa, orbitalsIsolated_alfa, overlapMat, inVariables, currentResults_alfa);
        printf("                     ----- End of open shell calculations: alfa density -----                     \n \n");
        printf("                     ----- Start of open shell calculations: beta density -----                      \n");
        calcOrbitalsCurrent(orbitalsComplex_beta, orbitalsIsolated_beta, overlapMat, inVariables, currentResults_beta);
        printf("                     ----- End of open shell calculations: beta density -----                     \n \n");

        printf("              ----- Total result for open shell calculations: alfa density -----              \n");
        printTotalResult(currentResults_alfa);
        printf("              ----- Total result for open shell calculations: beta density -----              \n");
        printTotalResult(currentResults_beta);
    }



//-------------- *********************************** --------------//

    return 0;
}
