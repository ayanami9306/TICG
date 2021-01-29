#ifndef TICG_OBJECT_H
#define TICG_OBJECT_H

#include <vector>

class SimulationBox 
{
    public:
    SimulationBox();
    ~SimulationBox(){};
    //parameter for box
    double Lx, Ly, Lz;  //size of simulation box
    double unitval_phi; //unitval_phi = num_beads / num_grid
    int num_grid_x, num_grid_y, num_grid_z, num_grid_xy, num_grid; //number of grids along each axis / total
    double num_molecules_per_grid;
    double size_grid_x, size_grid_y, size_grid_z; //size of each grids for each axis
    bool is_zwall; //if true, non-periodic condition along the z-axis, also hard-wall (reject if exceed z)

    //coefficient for energy
    double half_kappaN, chiN;
    std::vector<double> ideal_bond_length_square_for_bond_type;
    std::vector<double> ideal_inverse_bond_length_square_for_bond_type;
    
    //paramter for simulation
    int num_beads, num_molecules;
    //sqrt(N_bar) * dV / Re^3
    double const_param_for_Hnb;
    double prob_AtomTranslation, prob_MoleculeTranslation, prob_MoleculeRotation; // probability of action in each MC cycle, summation is 1
    double TranslationMax, RotationMax; // max value for translation distance and rotation angle

    unsigned long long int num_Accepts, num_Rejects;
    int Cur_Cycle_Annealing, Cur_Cycle_Simulation;
    int Max_Cycle_Annealing, Max_Cycle_Simulation;
    int num_Wiggle_Interval, num_Output_Interval;
};

class Bond
{
    public:
    int connected_bead_id;
    int bond_type;
    Bond(int, int);
    ~Bond(){};

    private:
};

class Bead
{
    public:
    Bead(int, int);
    ~Bead(){};
    int grid_id;
    int bead_id;
    int bead_type;
    std::vector<double> position;
    std::vector<Bond> bondlist;

    private:
};

class Molecule
{
    public:
    Molecule();
    ~Molecule(){};
    std::vector<Bead *> beadlist;
};

class Grid
{
    public:
    Grid(int);
    ~Grid(){};
    int grid_id;
    std::vector<double> phi;
};

#endif