#ifndef GETORBITALS_H_INCLUDED
#define GETORBITALS_H_INCLUDED

#include "dataTypes.h"

void getMolecOrbitals(const char * inFileName, vector <MolecOrb> & Orbitals, const size_t NBasis, string spinPolarization = "no");

void printMO(vector <MolecOrb> & MOs, const size_t Nbas, const size_t numberOfColumns = 5, const size_t columnWidth = 20);

size_t getMOAbsPosition(vector <MolecOrb> & MOs, MolecOrb & MO);

string getMORelPosition(vector <MolecOrb> & MOs, MolecOrb & MO);

void getOrbitalWindow(vector <MolecOrb> &MOs, vector <WindowOrbital> &windowMOs, string range);



#endif // GETORBITALS_H_INCLUDED
