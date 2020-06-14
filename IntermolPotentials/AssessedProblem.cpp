#include "AssessedProblem.h"
#include <iostream>
using namespace std;

// Print out the current charges and dipole moments
void printSystem(int nCharges, pointCharge charges[], int nDipoles, dipole dipoles[]){
    cout << "\nPoint charges";
    cout << "\nPosition z / Angstrom    Charge q / e";
    for(int i=0; i<nCharges; i++){
        cout << "\n" << charges[i].getPosition() << "    " << charges[i].getCharge();
    }
    cout << "\n\nDipoles";
    cout << "\nPosition z / Angstrom    Dipole moment mu / D";
    for(int i=0; i<nDipoles; i++){
        cout << "\n" << dipoles[i].getPosition() << "    " << dipoles[i].getMoment();
    }
}

int iterate(pointCharge charges, dipole dipoles){
    pointCharge oldCharges = charges;
    dipole oldDipoles = dipoles;
    return 0;
}

// Default parameters for 1d lattice of identical charges and 1d lattice of identical dipoles
int main(int argc, char *argv[]){
    int nCharges = 0;
    int nDipoles = 0;
    double chargeLatticeSeparation = 0.;
    double dipoleLatticeSeparation = 0.;
    double charge = 0.;
    double polarisability = 0.;
    double moment = 0.;

    // Read command line arguments to update parameters
    for (int i=1; i<argc; i+=2 ){
        string flag = string(argv[i]);
        string value = string(argv[i+1]);
        if (flag == "-nc"){
            nCharges = stoi(value);
        } else if (flag == "-nd"){
            nDipoles = stoi(value);
        } else if (flag == "-cl"){
            chargeLatticeSeparation = stod(value);
        } else if (flag == "-dl"){
            dipoleLatticeSeparation = stod(value);
        } else if (flag == "-c"){
            charge = stod(value);
        } else if (flag == "-p"){
            polarisability = stod(value);
        } else if (flag == "-m"){
            moment = stod(value);
        }
    }

    // Place charges
    pointCharge charges[nCharges];
    for(int i=0; i<nCharges; i++){
        charges[i].setCharge(i*chargeLatticeSeparation, charge);
    }

    // Place dipoles
    dipole dipoles[nDipoles];
    for(int i=0; i<nDipoles; i++){
        dipoles[i].setDipole((i+1)*dipoleLatticeSeparation, polarisability, moment);
    }

    printSystem(nCharges, charges, nDipoles, dipoles);
    return 0;
}