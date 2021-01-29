#include "TICG.h"


void TICG::WiggleBox()
{
    int axis = MTRANDAXIS();
    double ds = MTRANDGRID();
    int pre_grid_id;
    for (auto &target_bead : beads)
    {
        target_bead.position[axis] += ds;

        pre_grid_id = target_bead.grid_id;
        target_bead.grid_id = GridID_from_Position(&target_bead.position);
        if (pre_grid_id != target_bead.grid_id)
        {
            Remove_Bead_from_Grid(&grids[pre_grid_id], target_bead.bead_type);
            Add_Bead_to_Grid(&grids[target_bead.grid_id], target_bead.bead_type);
        }
    }
}

void TICG::AdjustCM()
{
    for (auto &target_molecule : molecules)
    {
        std::vector<double> CM;
        Measure_CM(&(target_molecule.beadlist), &CM);

        CM[0] -= float_mod(CM[0], simbox.Lx);
        CM[1] -= float_mod(CM[1], simbox.Ly);
        CM[2] -= float_mod(CM[2], simbox.Lz);

        for (auto &target_bead : target_molecule.beadlist)
        {
            target_bead->position[0] -= CM[0];
            target_bead->position[1] -= CM[1];
            target_bead->position[2] -= CM[2];
        }
    }
}