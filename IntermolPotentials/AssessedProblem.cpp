#include "AssessedProblem.h"
#include <cmath>
#include <iostream>
#include <vector>

// Print out the current charges and dipole moments
void printSystem(int nCharges, std::vector<pointCharge> charges, int nDipoles, std::vector<dipole> dipoles){
    std::cout << "\nPoint charges";
    std::cout << "\nPosition z / Angstrom    Charge q / e";
    for(int i=0; i<nCharges; i++){
        std::cout << "\n" << charges[i].getPosition() << "    " << charges[i].getCharge();
    }
    std::cout << "\n\nDipoles";
    std::cout << "\nPosition z / Angstrom    Dipole moment mu / D";
    for(int i=0; i<nDipoles; i++){
        std::cout << "\n" << dipoles[i].getPosition() << "    " << dipoles[i].getMoment();
    }
    std::cout << "\n";
}

// Displacement between two position vectors (in one dimension, this is just the difference)
double displacement(double position1, double position2){
    double disp = position2 - position1;
    return disp;
}

// Distance between two position vectors
double distance(double position1, double position2){
    double disp = displacement(position1, position2);
    double dist = sqrt(pow(disp, 2));
    return dist;
}

// Dot product between two vectors (in one dimension, this is just the product)
double dotproduct(double vector1, double vector2){
    double prod = vector1 * vector2;
    return prod;
}

// Electrostatic potential at a given position
double potential(double position, int nCharges, std::vector<pointCharge> charges, int nDipoles, std::vector<dipole> dipoles){
    double phi = 0.;
    for (int i=0; i<nCharges; i++){
        double dist = distance(charges[i].getPosition(), position);
        phi += charges[i].getCharge() / dist;
    }
    for (int i=0; i<nDipoles; i++){
        double disp = distance(dipoles[i].getPosition(), position);
        double dist = distance(dipoles[i].getPosition(), position);
        phi += dotproduct(dipoles[i].getMoment(), disp) / pow(dist, 3);
    }
    return phi;
}

// Electric field at a given position
double field(double position, int nCharges, std::vector<pointCharge> charges, int nDipoles, std::vector<dipole> dipoles){
    double F = 0;
    for (int i=0; i<nCharges; i++){
        double disp = displacement(charges[i].getPosition(), position);
        double dist = distance(charges[i].getPosition(), position);
        F += charges[i].getCharge() * disp / pow(dist, 2);
    }
    for (int i=0; i<nDipoles; i++){
        double disp = distance(dipoles[i].getPosition(), position);
        double dist = distance(dipoles[i].getPosition(), position);
        F += 2 * dipoles[i].getMoment() / pow(dist, 3);
    }
    return F;
}

std::vector<dipole> newDipoles(int nCharges, std::vector<pointCharge> charges, int nDipoles, std::vector<dipole> dipoles){
    std::vector<dipole> oldDipoles(nDipoles);
    double oldMoment, newMoment, polarisability, F;
    for (int i=0; i<nDipoles; i++){
        oldMoment = oldDipoles[i].getMoment();
        polarisability = oldDipoles[i].getPolarisability();
        F = field(oldDipoles[i].getPosition(), nCharges, charges, nDipoles, oldDipoles);
        dipoles[i].setMoment(oldMoment + polarisability * F);
    }
    return dipoles;
}

int main(int argc, char *argv[]){
    // Default parameters for 1d lattice of identical charges and 1d lattice of identical dipoles
    int nCharges = 0;
    int nDipoles = 0;
    double chargeLatticeSeparation = 0.;
    double dipoleLatticeSeparation = 0.;
    double charge = 0.;
    double polarisability = 0.;
    double moment = 0.;

    // Read command line arguments to update parameters
    for (int i=1; i<argc; i+=2 ){
        std::string flag = std::string(argv[i]);
        std::string value = std::string(argv[i+1]);
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
    std::vector<pointCharge> charges(nCharges);
    for(int i=0; i<nCharges; i++){
        charges[i].setCharge(i*chargeLatticeSeparation, charge);
    }

    // Place dipoles
    std::vector<dipole> dipoles(nDipoles);
    for(int i=0; i<nDipoles; i++){
        dipoles[i].setDipole((i+1)*dipoleLatticeSeparation, polarisability, moment);
    }

    printSystem(nCharges, charges, nDipoles, dipoles);
    dipoles = newDipoles(nCharges, charges, nDipoles, dipoles);
    printSystem(nCharges, charges, nDipoles, dipoles);
    return 0;
}