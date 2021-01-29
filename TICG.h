#ifndef TICG_H
#define TICG_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <random>
#include <functional>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <ctime>
#include <cmath>
#include <smmintrin.h>
#include <emmintrin.h>
#include <algorithm>
#include <numeric>
#include "TICG_objects.h"
#define MTRAND() rand_distribution(rand_engine)
#define MTRANDDIST() rand_distance(rand_engine)
#define MTRANDBEAD() rand_beadnum(rand_engine)
#define MTRANDGRID() rand_grid(rand_engine)
#define MTRANDAXIS() rand_axis(rand_engine)

class TICG {
    public:
    //public variables

    //public functions
    TICG(char*);
    ~TICG();
    void TICG_Simulation();

    private:
    //private variables
    SimulationBox simbox;
    std::vector<Bead> beads;
    std::vector<Molecule> molecules;
    std::vector<Grid> grids;
    std::mt19937 rand_engine;
    std::uniform_real_distribution<double> rand_distribution;
    std::uniform_real_distribution<double> rand_distance;
    std::uniform_real_distribution<double> rand_grid;
    std::uniform_int_distribution<int> rand_beadnum;
    std::uniform_int_distribution<int> rand_axis;
    std::vector<double> sin_table;
    std::vector<double> cos_table;
    char simname[200];

    //private functions
    void Write_Trigonometric_Table();
    double float_mod(double, double);
    double Calc_Bond_Length2(std::vector<double> *, std::vector<double> *);
    int GridID_from_Position(std::vector<double> *);
    double Calc_Bonded_Potential_for_Bead(std::vector<double> *, std::vector<Bond> *);
    double Calc_Difference_Bonded_Potential_for_Bead(std::vector<double> *, std::vector<double> *, std::vector<Bond> *);
    double Calc_Difference_Nonbonded_Potential_for_Grid(Grid *, double, double);
    double Calc_Nonbonded_Potential_for_Grid(Grid *);
    void MonteCarlo(int);
    void TranslateBead();
    //int TranslateMolecule();
    //int RotateMolecule();
    void WiggleBox();
    void AdjustCM();
    void WriteoutSim(char *);
    void ReadinSim(char *);
    void ConstructSim(char *);
    char READ_CONFIG(char*, char*);
    void AssignMoleculeIndex();
    void AssignBeadtoGrid();
    void Measure_CM(std::vector<Bead *> * , std::vector<double> * );
    void Remove_Bead_from_Grid(Grid *, int bead_type);
    void Add_Bead_to_Grid(Grid *, int bead_type);
    void writePDB(char *, bool);
    void writeEnergyProfile(char *, bool);
    std::vector<std::string> strsplit(std::string str, char delimiter);
};

inline double TICG::float_mod(double a, double n) 
{
    return a - (floor(a/n) * n);
}

inline double TICG::Calc_Bond_Length2(std::vector<double> * pos1, std::vector<double> * pos2)
{
    return ((*pos1)[0]-(*pos2)[0]) * ((*pos1)[0]-(*pos2)[0]) + ((*pos1)[1]-(*pos2)[1]) * ((*pos1)[1]-(*pos2)[1]) + ((*pos1)[2]-(*pos2)[2]) * ((*pos1)[2]-(*pos2)[2]);
}

inline int TICG::GridID_from_Position(std::vector<double> * pos)
{
    //apply periodic boundary, and calc grid index
    return (static_cast <int>(float_mod((*pos)[0], simbox.Lx) / simbox.size_grid_x)) + (static_cast <int>(float_mod((*pos)[1], simbox.Ly) / simbox.size_grid_y) * simbox.num_grid_x) + (static_cast <int>(float_mod((*pos)[2], simbox.Lz) / simbox.size_grid_z) * simbox.num_grid_xy);
}

inline void TICG::Remove_Bead_from_Grid(Grid * target_grid, int bead_type)
{
    target_grid->phi[bead_type] -= simbox.unitval_phi;
}

inline void TICG::Add_Bead_to_Grid(Grid * target_grid, int bead_type)
{
    target_grid->phi[bead_type] += simbox.unitval_phi;
}

#endif