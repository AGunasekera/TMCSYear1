#include "AssessedProblem.h"
#include <cmath>
#include <iostream>
#include <vector>

// Set parameters, parsing from command line arguments, with appropriate unti conversions.
void setParams(int argC, char *argV[], double &extField, double &extPotential, int &nCharges, int &nDipoles, double &cLattSep, double &dLattSep, double &charge, double &polaris, double &statMom, int &thresh){
    extField  = 0.;
    extPotential = 0.;
    nCharges = 0;
    nDipoles = 0;
    cLattSep = 0.;
    dLattSep = 0.;
    charge = 0.;
    polaris = 0.;
    statMom = 0.;
    thresh = 0;

    for (int i=1; i<argC; i+=2 ){
        std::string flag = std::string(argV[i]);
        std::string value = std::string(argV[i+1]);
        if (flag == "-nc"){
            nCharges = stoi(value);
        } else if (flag == "-nd"){
            nDipoles = stoi(value);
        } else if (flag == "-cl"){
            cLattSep = stod(value);
        } else if (flag == "-dl"){
            dLattSep = stod(value);
        } else if (flag == "-c"){
            charge = stod(value);
        } else if (flag == "-a"){
            polaris = stod(value);
        } else if (flag == "-m"){
            double statMomDebye = stod(value);
            statMom = statMomDebye / (2.99792458 * 1.602176634);
        } else if (flag == "-p"){
            double extPotentialVolt = stod(value);
            extPotential = extPotentialVolt * 4 * 4 * atan(1.) * 8.854187812813 / 1602.176634;
        } else if (flag == "-f"){
            double extFieldVolt = stod(value);
            extField = extFieldVolt * 4 * 4 * atan(1.) * 8.854187812813 / 1602.176634;
        } else if (flag == "-t"){
            thresh = stod(value);
        }
    }
}

// Initialise a point charge
void initialiseCharge(double position, double potential, double charge, pointCharge &pC){
    pC.setPosition(position);
    pC.setPotential(potential);
    pC.setCharge(charge);
}

// Initialise a dipole
void initialiseDipole(double position, double field, double polarisability, double moment, dipole &dp){
    dp.setPosition(position);
    dp.setField(field);
    dp.setPolarisability(polarisability);
    dp.setMoment(moment);
}

// Print out the current charges and dipole moments (converted back to Debye)
void printSystem(std::vector<pointCharge> charges, std::vector<dipole> dipoles){
    std::cout << "\nPoint charges";
    std::cout << "\nPosition z / Angstrom    Charge q / e";
    for(int i=0; i<charges.size(); i++){
        std::cout << "\n" << charges[i].getPosition() << "    " << charges[i].getCharge();
    }
    std::cout << "\n\nDipoles";
    std::cout << "\nPosition z / Angstrom    Dipole moment mu / D";
    for(int i=0; i<dipoles.size(); i++){
        std::cout << "\n" << dipoles[i].getPosition() << "    " << dipoles[i].getMoment() * 2.99792458 * 1.602176634;
    }
    std::cout << "\n";
}

// Displacement vector from position1 to position2 (in one dimension, this is just the difference)
double displacement(double position1, double position2){
    double disp = position2 - position1;
    return disp;
}

// Dot product between two vectors (in one dimension, this is just the product)
double dotproduct(double vector1, double vector2){
    double prod = vector1 * vector2;
    return prod;
}

// Magnitude of a vector
double magnitude(double vector){
    return sqrt(dotproduct(vector, vector));
}

// Electrostatic potential at a given position generated by a given point charge
double pointChargePotential(double position, pointCharge pointCharge_){
    double dist = magnitude(displacement(pointCharge_.getPosition(), position));
    return pointCharge_.getCharge() / dist;
}

// Electrostatic potential at a given position generated by a given dipole
double dipolePotential(double position, dipole dipole_){
    double disp = displacement(dipole_.getPosition(), position);
    double dist = magnitude(disp);
    return dotproduct(dipole_.getMoment(), disp) / pow(dist, 3);
}

// Electric field at a given position generated by a given point charge
double pointChargeField(double position, pointCharge pointCharge_){
    double disp = displacement(pointCharge_.getPosition(), position);
    double dist = magnitude(disp);
    return pointCharge_.getCharge() * disp / pow(dist, 3);
}

// Electric field at a given position generated by a given dipole
double dipoleField(double position, dipole dipole_){
    double disp = displacement(dipole_.getPosition(), position);
    double dist = magnitude(disp);
    return 2 * dipole_.getMoment() / pow(dist, 3);
}

// Calculate the potential experienced by each point charge produced by the rest of the system
void potentialsAtCharges(std::vector<pointCharge> &charges, std::vector<dipole> &dipoles, double externalPotential){
    double position, potential;
    for (int i=0; i<charges.size(); i++){
        position = charges[i].getPosition();
        potential = externalPotential;
        for (int j=0; j<charges.size(); j++){
            if (position != charges[j].getPosition()){
                potential += pointChargePotential(position, charges[j]);
            }
        }
        for (int j=0; j<dipoles.size(); j++){
            if (position != dipoles[j].getPosition()){
                potential += dipolePotential(position, dipoles[j]);
            }
        }
        charges[i].setPotential(potential);
    }
}

// Calculate the field experienced by each dipole produced by the rest of the system
void fieldAtDipoles(std::vector<pointCharge> &charges, std::vector<dipole> &dipoles, double externalField){
    std::vector<dipole> oldDipoles = dipoles;
    double position, field;
    for (int i=0; i<dipoles.size(); i++){
        position = oldDipoles[i].getPosition();
        field = externalField;
        for (int j=0; j<charges.size(); j++){
            if (position != charges[j].getPosition()){
                field += pointChargeField(position, charges[j]);
            }
        }
        for (int j=0; j<dipoles.size(); j++){
            if (position != dipoles[j].getPosition()){
                field += dipoleField(position, oldDipoles[j]);
            }
        }
        dipoles[i].setField(field);
    }
}

// Update dipole moments based on the field they experience, returning the largest absolute change in dipole moment
double updateDipoles(std::vector<dipole> &dipoles, double staticMoment){
    double oldMoment, newMoment, polarisability, field;
    double change = 0;
    for (int i=0; i<dipoles.size(); i++){
        oldMoment = dipoles[i].getMoment();
        polarisability = dipoles[i].getPolarisability();
        field = dipoles[i].getField();
        newMoment = staticMoment + polarisability * field;
        if (magnitude(change) < magnitude(newMoment - oldMoment)){
            change = newMoment - oldMoment;
        }
        dipoles[i].setMoment(newMoment);
    }
    return magnitude(change);
}

// Energy of the system in its given state (note adding on external potential and external field, as only interaction terms have been double counted)
double systemEnergy(std::vector<pointCharge> charges, std::vector<dipole> dipoles, double extPotential, double extField){
    double sumEnergy = 0.;
    for (int i=0; i<charges.size(); i++){
        sumEnergy += charges[i].getCharge() * (charges[i].getPotential() + extPotential);
    }
    for (int i=0; i<dipoles.size(); i++){
        sumEnergy -= dotproduct(dipoles[i].getMoment(), (dipoles[i].getField() + extField));
    }
    return sumEnergy / 2;
}

int main(int argc, char *argv[]){
    // Declare parameters for 1d lattice of identical charges and 1d lattice of identical dipoles
    /* System of units used:
    * Distance in Angstrom
    * Charge in elementary charge e
    * Energy in e^2 / (4pi epsilon0 Angstrom), converted to eV for output
    * Potential in e / (4pi epsilon0 Angstrom), external potential converted from V for input
    * Field in e / (4pi epsilon0 Angstrom^2), external field converted from V / A for input
    * Dipole moment in e * Angstrom = 2.99792458 * 1.602176634 * Debye (converted from Debye and back for input and output)
    * Polarisability in 4pi epsilon0 Angstrom^3
    */
    double externalPotential, externalField;
    int nCharges, nDipoles;
    double chargeLatticeSeparation, dipoleLatticeSeparation;
    double charge, polarisability, staticMoment;

    // Convergence threshold: iteration terminates when all dipoles change by less than 10^(-threshold) Debye
    int threshold;

    // Set parameters
    setParams(argc, argv, externalField, externalPotential, nCharges, nDipoles, chargeLatticeSeparation, dipoleLatticeSeparation, charge, polarisability, staticMoment, threshold);

    // Place charges in 1d lattice, with first charge at 0
    std::vector<pointCharge> charges(nCharges);
    for(int i=0; i<nCharges; i++){
        initialiseCharge(i*chargeLatticeSeparation, externalPotential, charge, charges[i]);
    }

    // Place dipoles, with first dipole at 1 * dipoleLatticeSeparation
    std::vector<dipole> dipoles(nDipoles);
    for(int i=0; i<nDipoles; i++){
        initialiseDipole((i+1)*dipoleLatticeSeparation, externalField, polarisability, staticMoment, dipoles[i]);
    }

    double dipoleChange_D = 0.;
    double energyChange_eV = 0.;
    double totalEnergy_eV, electrostaticEnergy_eV;
    double inductionEnergy_eV = 0.;
    int i = 0, maxIterations = 1000;
    while(i < maxIterations){
        std::cout << "Iteration " << i << "\n";
        potentialsAtCharges(charges, dipoles, externalPotential);
        fieldAtDipoles(charges, dipoles, externalField);
        printSystem(charges, dipoles);
        totalEnergy_eV = systemEnergy(charges, dipoles, externalPotential, externalField) * 1602.176634 / (4 * 4 * atan(1.) * 8.854187812813);
        if (i == 0){
            electrostaticEnergy_eV = totalEnergy_eV;
        }
        energyChange_eV = -inductionEnergy_eV;
        inductionEnergy_eV = totalEnergy_eV - electrostaticEnergy_eV;
        energyChange_eV += inductionEnergy_eV;
        std::cout << "\nEnergy: " << totalEnergy_eV << " eV\n";
        std::cout << "Induction energy: " << inductionEnergy_eV << " eV\n";
        std::cout << "Change in induction energy: " << energyChange_eV << " eV\n";
        dipoleChange_D = updateDipoles(dipoles, staticMoment) * 2.99792458 * 1.602176634;
        std::cout << "\nNew dipoles calculated. Max change: " << dipoleChange_D << " D\n";
        i++;
        if (dipoleChange_D < pow(10, -threshold)){
            std::cout << "\nCalculation converged in " << i << " steps.\n";
            printSystem(charges, dipoles);
            totalEnergy_eV = systemEnergy(charges, dipoles, externalPotential, externalField) * 1602.176634 / (4 * 4 * atan(1.) * 8.854187812813);
            energyChange_eV = -inductionEnergy_eV;
            inductionEnergy_eV = totalEnergy_eV - electrostaticEnergy_eV;
            energyChange_eV += inductionEnergy_eV;
            std::cout << "\nElectrostatic energy: " << electrostaticEnergy_eV << " eV\n";
            std::cout << "Final energy: " << totalEnergy_eV << " eV\n";
            std::cout << "Final induction energy: " << inductionEnergy_eV << " eV\n";
            std::cout << "Last change in induction energy: " << energyChange_eV << " eV\n";
            return 0;
        }
    }
    std::cout << "\nDipoles not converged in " << maxIterations << " steps.\n";
    return 1;
}