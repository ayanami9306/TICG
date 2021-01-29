#include "TICG_objects.h"

SimulationBox::SimulationBox()
{

}

Bond::Bond(int id, int type)
{
    connected_bead_id = id;
    bond_type = type;
}

Bead::Bead(int id, int type)
{
    bead_id = id;
    bead_type = type;
    grid_id = -1;
}

Molecule::Molecule()
{

}

Grid::Grid(int id) 
{
    grid_id = id;
} 