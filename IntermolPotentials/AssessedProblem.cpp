#include <iostream>
using namespace std;

class pointCharge {
// Class for a point charge at a given position with a given charge.
    double position;
    double charge;
    public:
        pointCharge(){
            position = 0;
            charge = 0;
        }
        void setCharge(double position_, double charge_){
            position = position_;
            charge = charge_;
        }
        double getPosition(){
            return position;
        }
        double getCharge(){
            return charge;
        }
};

class dipole {
// Class for a dipole at a given position, with a given polarisability, and a corresponding dipole moment.
    double position;
    double polarisability;
    double moment;
    public:
        dipole(){
            position = 0;
            polarisability = 0;
            moment = 0;
        }
        void setDipole(double position_, double polarisability_, double moment_){
            position = position_;
            polarisability = polarisability_;
            moment = moment_;
        }
        double getPosition(){
            return position;
        }
        double getPolarisability(){
            return polarisability;
        }
        double getMoment(){
            return moment;
        }
};

int printSystem(int nCharges, pointCharge charges[], int nDipoles, dipole dipoles[]){
// Print out the current charges and dipole moments
    cout << "Point charges" << endl;
    cout << "Position z / Angstrom    Charge q / e" << endl;
    for(int i=0; i<nCharges; i++){
        cout << charges[i].getPosition() << "    " << charges[i].getCharge() << endl;
    }
    cout << endl;
    cout << "Dipoles" << endl;
    cout << "Position z / Angstrom    Dipole moment mu / D" << endl;
    for(int i=0; i<nDipoles; i++){
        cout << dipoles[i].getPosition() << "    " << dipoles[i].getMoment() << endl;
    }
    return 0;
}

int iterate(pointCharge charges, dipole dipoles){
    pointCharge oldCharges = charges;
    dipole oldDipoles = dipoles;
    return 0;
}

int main(){
// Initial parameters for 1d lattice of identical charges and 1d lattice of identical dipoles
    int nCharges = 0;
    int nDipoles = 5;
    double chargeLatticeSeparation = -1.;
    double dipoleLatticeSeparation = 1.;
    double charge = 1.;
    double polarisability = 1.;
    double moment = 1.;

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