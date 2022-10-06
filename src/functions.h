#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "dataTypes.h"


double siteOrbProjection(const size_t siteNum, MolecOrb & Orbital, const dmat & OverlapMat);

double dExpansion(const size_t siteNum, vector<MolecOrb> Orbitals, const double degThrh = 0, const double T = 298.15);

double V (vector <MolecOrb> orbitalsComplex, vector <MolecOrb> orbitalsIsolated, const double homoReserv);

double getFermiEnrg(vector <MolecOrb> Orbitals);

double f (double x, void * params);

double TransitionProb(double gamma1, double gammaN, double site1MO, double siteMON, double EnergyMO, double E);

void integrateRoutine(const double E, integrParams & integrVariables, integrResult & viPoint, size_t routine = 0);

#endif // FUNCTIONS_H_INCLUDED
