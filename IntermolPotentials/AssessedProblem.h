#ifndef INTERMOLECULAR_POTENTIALS_ASSESSED
#define INTERMOLECULAR_POTENTIALS_ASSESSED

// Class for a point charge at a given position with a given charge.
class pointCharge {
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

// Class for a dipole at a given position, with a given polarisability, and a corresponding dipole moment.
class dipole {
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

#endif