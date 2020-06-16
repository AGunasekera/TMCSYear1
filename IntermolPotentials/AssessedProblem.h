#ifndef INTERMOLECULAR_POTENTIALS_ASSESSED
#define INTERMOLECULAR_POTENTIALS_ASSESSED

// Class for a point charge at a given position with a given charge.
class pointCharge {
    double position;
    double potential;
    double charge;
    public:
        pointCharge(){
            position = 0;
            potential = 0;
            charge = 0;
        }
        void setPosition(double position_){
            position = position_;
        }
        void setPotential(double potential_){
            potential = potential_;
        }
        void setCharge(double charge_){
            charge = charge_;
        }
        double getPosition(){
            return position;
        }
        double getPotential(){
            return potential;
        }
        double getCharge(){
            return charge;
        }
};

// Class for a dipole at a given position, with a given polarisability, and a corresponding dipole moment.
class dipole {
    double position;
    double field;
    double polarisability;
    double moment;
    public:
        dipole(){
            position = 0;
            field = 0;
            polarisability = 0;
            moment = 0;
        }
        void setPosition(double position_){
            position = position_;
        }
        void setField(double field_){
            field = field_;
        }
        void setPolarisability(double polarisability_){
            polarisability = polarisability_;
        }
        void setMoment(double moment_){
            moment = moment_;
        }
        double getPosition(){
            return position;
        }
        double getField(){
            return field;
        }
        double getPolarisability(){
            return polarisability;
        }
        double getMoment(){
            return moment;
        }
};

#endif