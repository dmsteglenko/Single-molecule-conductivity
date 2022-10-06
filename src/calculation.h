#ifndef CALCULATION_H_INCLUDED
#define CALCULATION_H_INCLUDED

void calcOrbitalsCurrent(vector <MolecOrb> & orbitalsComplex, vector <MolecOrb> & orbitalsIsolated,
                         const dmat & Overlap, const inputs & variables, vector <orbitalCurrentResult> & orbitalsVCResult);

#endif // CALCULATION_H_INCLUDED
