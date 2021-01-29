#include "TICG.h"

void TICG::TranslateBead()
{
    //Hopping Bead
    Bead *hop_bead = &beads[MTRANDBEAD()];
    double dx = MTRANDDIST(), dy = MTRANDDIST(), dz = MTRANDDIST();
            
    std::vector<double> position_after = hop_bead->position;
    std::vector<double> *position_before = &(hop_bead->position);
    position_after[0] += dx;
    position_after[1] += dy;
    position_after[2] += dz;

    //is wall?
    if (simbox.is_zwall) 
    {
        if (position_after[2] < 0 || position_after[2] >= simbox.Lz)
        {
            simbox.num_Rejects ++;
            return;
        }
    }

    //std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
    //std::chrono::duration<double> sec = std::chrono::system_clock::now() - start_time;
    //sec = std::chrono::system_clock::now() - start_time;
    //time_hopping += std::chrono::duration_cast<std::chrono::nanoseconds>(sec).count();
    
    //calc bonded potential -dE = U_pre - U_after
    double minus_delta_E = Calc_Difference_Bonded_Potential_for_Bead(position_before, &position_after, &(hop_bead->bondlist));
    
    //get new grid id for hopped bead
    int new_grid_id = GridID_from_Position(&position_after),
        pre_grid_id = hop_bead->grid_id,
        hop_bead_type = hop_bead->bead_type;
    //if grid id is different : calc nonbonded potential

    //calc difference of nonbonded potential
    if (new_grid_id != pre_grid_id) minus_delta_E += Calc_Difference_Nonbonded_Potential_for_Grid(&grids[new_grid_id], simbox.unitval_phi * !(hop_bead_type ^ 0), simbox.unitval_phi * !(hop_bead_type ^ 1)) + Calc_Difference_Nonbonded_Potential_for_Grid(&grids[pre_grid_id], -simbox.unitval_phi * !(hop_bead_type ^ 0), -simbox.unitval_phi * !(hop_bead_type ^ 1));
    
    //accept or reject
    //accept absolutely
    if (minus_delta_E >= 0)
    {
        (*position_before)[0] = position_after[0];
        (*position_before)[1] = position_after[1];
        (*position_before)[2] = position_after[2];
        if (pre_grid_id != new_grid_id)
        {
            hop_bead->grid_id = new_grid_id;
            Remove_Bead_from_Grid(&grids[pre_grid_id], hop_bead_type);
            Add_Bead_to_Grid(&grids[new_grid_id], hop_bead_type);
        }
        simbox.num_Accepts++;
    }
    //accept hopping
    else if (MTRAND() < exp(minus_delta_E))
    {
        (*position_before)[0] = position_after[0];
        (*position_before)[1] = position_after[1];
        (*position_before)[2] = position_after[2];
        if (pre_grid_id != new_grid_id)
        {
            hop_bead->grid_id = new_grid_id;
            Remove_Bead_from_Grid(&grids[pre_grid_id], hop_bead_type);
            Add_Bead_to_Grid(&grids[new_grid_id], hop_bead_type);
        }
        simbox.num_Accepts++;
    }
    //else reject
    else simbox.num_Rejects++;   
}
/*
int TICG::TranslateMolecule()
{
    //Hopping Molecule
    std::vector<Bead *> *target_beads = &(molecules[static_cast <int>(MTRAND() * simbox.num_molecules)].beadlist);
    double hopping_distance = static_cast <double>(MTRAND() * simbox.TranslationMax);
    //0<=theta<=180
    int theta = static_cast <int>(MTRAND() * 181);
    //0<=phi<360
    int phi = static_cast <int>(MTRAND() * 360);
    double dx = hopping_distance * sin_table[theta] * cos_table[phi];
    double dy = hopping_distance * sin_table[theta] * sin_table[phi];
    double dz = hopping_distance * cos_table[theta];

    //internal coordinate does not affect, so delta H_b is 0.
    //so, measure delta H_nb only.

    std::vector<std::vector <double>> molecule_pos;
    std::vector<int> molecule_new_grid_ids;
    molecule_pos.assign(target_beads->size(), std::vector<double> {0.0, 0.0, 0.0});
    molecule_new_grid_ids.assign(target_beads->size(), 0);
    //pre_grid_id, new_grid_id, bead_type
    std::vector<std::vector<int>> affected_grids;
    std::set<int> affected_grid_ids;
    int bead_count = 0;

    for (auto &target_bead : *target_beads)
    {
        std::vector<double> position_after = target_bead->position;
        position_after[0] += dx;
        position_after[1] += dy;
        position_after[2] += dz;
        //is wall?
        if (simbox.is_zwall) 
        {
            if (position_after[2] < 0 || position_after[2] >= simbox.Lz)
            {
                simbox.num_Rejects ++;
                return target_beads->size()-1;
            }
        }
        //save hopping position
        std::copy(position_after.begin(), position_after.end(), molecule_pos[bead_count].begin());

        //get grid id
        int new_grid_id = GridID_from_Position(&position_after);
        int pre_grid_id = target_bead->grid_id;

        molecule_new_grid_ids[bead_count++] = new_grid_id;

        //if grid needs update, save grid info
        if (pre_grid_id != new_grid_id)
        {
            std::vector<int> grid_info = {pre_grid_id, new_grid_id, target_bead->bead_type};
            affected_grids.push_back(grid_info);
            affected_grid_ids.insert(pre_grid_id);
            affected_grid_ids.insert(new_grid_id);
        }
    }

    double minus_delta_E = 0;
    //measure pre grid potential
    for (auto &target_grid_id : affected_grid_ids) minus_delta_E += Calc_Nonbonded_Potential_for_Grid(&grids[target_grid_id]);
            
    //update phi of grid
    for (auto &grid_info : affected_grids)
    {
        Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
        Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
    }

    //measure after grid potential
    for (auto &target_grid_id : affected_grid_ids) minus_delta_E -= Calc_Nonbonded_Potential_for_Grid(&grids[target_grid_id]);

    //rollback phi of grid
    for (auto &grid_info : affected_grids)
    {
        Remove_Bead_from_Grid(&grids[grid_info[1]], grid_info[2]);
        Add_Bead_to_Grid(&grids[grid_info[0]], grid_info[2]);
    }

    //accept or reject
    //accept absolutely
    if (minus_delta_E >= 0)
    {   
        int pos_counter = 0;
        for (auto &target_bead : *target_beads)
        {
            std::vector<double> position_after = molecule_pos[pos_counter];
            std::copy(position_after.begin(), position_after.end(), target_bead->position.begin());
            target_bead->grid_id = molecule_new_grid_ids[pos_counter++];
        }

        //update phi of grid
        for (auto &grid_info : affected_grids)
        {
            Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
            Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
        }
        simbox.num_Accepts++;
    }
    //accept hopping
    else if (MTRAND() < exp(minus_delta_E))
    {
        int pos_counter = 0;
        for (auto &target_bead : *target_beads)
        {
            std::vector<double> position_after = molecule_pos[pos_counter];
            std::copy(position_after.begin(), position_after.end(), target_bead->position.begin());
            target_bead->grid_id = molecule_new_grid_ids[pos_counter++];
        }

        //update phi of grid
        for (auto &grid_info : affected_grids)
        {
            Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
            Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
        }
        simbox.num_Accepts++;
    }
    //else reject
    else simbox.num_Rejects++;

    return target_beads->size()-1;
}

int TICG::RotateMolecule()
{
    //Rotating Molecule
    std::vector<Bead *> *target_beads = &(molecules[static_cast <int>(MTRAND() * simbox.num_molecules)].beadlist);

    //measure center of mass of molecule
    std::vector<double> CM;
    Measure_CM(target_beads, &CM);

    double x_angle = static_cast <int>(MTRAND() * 360);
    double y_angle = static_cast <int>(MTRAND() * 360);
    double z_angle = static_cast <int>(MTRAND() * 360);

    //calc rotation max
    std::vector<std::vector<double>> Rotation_Matrix { 
        {cos_table[z_angle]*cos_table[y_angle], cos_table[z_angle]*sin_table[y_angle]*sin_table[x_angle] - sin_table[z_angle]*cos_table[x_angle], cos_table[z_angle]*sin_table[y_angle]*cos_table[x_angle]+sin_table[z_angle]*sin_table[x_angle]},
        {sin_table[z_angle]*cos_table[y_angle], sin_table[z_angle]*sin_table[y_angle]*sin_table[x_angle] + cos_table[z_angle]*cos_table[x_angle], sin_table[z_angle]*sin_table[y_angle]*cos_table[x_angle]-cos_table[z_angle]*sin_table[x_angle]},
        {-sin_table[y_angle], cos_table[y_angle]*sin_table[x_angle], cos_table[y_angle]*cos_table[x_angle]}
    };

    //internal coordinate does not affect, so delta H_b is 0.
    //so, measure delta H_nb only.

    std::vector<std::vector <double>> molecule_pos;
    std::vector<int> molecule_new_grid_ids;
    molecule_pos.assign(target_beads->size(), std::vector<double> {0.0, 0.0, 0.0});
    molecule_new_grid_ids.assign(target_beads->size(), 0);
    //pre_grid_id, new_grid_id, bead_type
    std::vector<std::vector<int>> affected_grids;
    std::set<int> affected_grid_ids;

    int bead_count = 0;
    for (auto &target_bead : *target_beads)
    {
        //rotate beads for CM
        std::vector<double> position_after = target_bead->position;

        //translate to origin
        position_after[0] -= CM[0];
        position_after[1] -= CM[1];
        position_after[2] -= CM[2];
        std::vector<double> temp_vector = position_after;

        //rotate...
        position_after[0] = Rotation_Matrix[0][0]*temp_vector[0] + Rotation_Matrix[0][1]*temp_vector[1] + Rotation_Matrix[0][2]*temp_vector[2];
        position_after[1] = Rotation_Matrix[1][0]*temp_vector[0] + Rotation_Matrix[1][1]*temp_vector[1] + Rotation_Matrix[1][2]*temp_vector[2];
        position_after[2] = Rotation_Matrix[2][0]*temp_vector[0] + Rotation_Matrix[2][1]*temp_vector[1] + Rotation_Matrix[2][2]*temp_vector[2];

        //translate again
        position_after[0] += CM[0];
        position_after[1] += CM[1];
        position_after[2] += CM[2];

        //is wall?
        if (simbox.is_zwall) 
        {
            if (position_after[2] < 0 || position_after[2] >= simbox.Lz)
            {
                simbox.num_Rejects ++;
                return target_beads->size()-1;
            }
        }
        //save hopping position
        std::copy(position_after.begin(), position_after.end(), molecule_pos[bead_count].begin());

        //get grid id
        int new_grid_id = GridID_from_Position(&position_after);
        int pre_grid_id = target_bead->grid_id;

        molecule_new_grid_ids[bead_count++] = new_grid_id;

        //if grid needs update, save grid info
        if (pre_grid_id != new_grid_id)
        {
            std::vector<int> grid_info = {pre_grid_id, new_grid_id, target_bead->bead_type};
            affected_grids.push_back(grid_info);
            affected_grid_ids.insert(pre_grid_id);
            affected_grid_ids.insert(new_grid_id);
        }
    }

    double minus_delta_E = 0;
    //measure pre grid potential
    for (auto &target_grid_id : affected_grid_ids) minus_delta_E += Calc_Nonbonded_Potential_for_Grid(&grids[target_grid_id]);
            
    //update phi of grid
    for (auto &grid_info : affected_grids)
    {
        Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
        Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
    }

    //measure after grid potential
    for (auto &target_grid_id : affected_grid_ids) minus_delta_E -= Calc_Nonbonded_Potential_for_Grid(&grids[target_grid_id]);

    //rollback phi of grid
    for (auto &grid_info : affected_grids)
    {
        Remove_Bead_from_Grid(&grids[grid_info[1]], grid_info[2]);
        Add_Bead_to_Grid(&grids[grid_info[0]], grid_info[2]);
    }

    //accept or reject
    //accept absolutely
    if (minus_delta_E >= 0)
    {   
        int pos_counter = 0;
        for (auto &target_bead : *target_beads)
        {
            std::vector<double> position_after = molecule_pos[pos_counter];
            std::copy(position_after.begin(), position_after.end(), target_bead->position.begin());
            target_bead->grid_id = molecule_new_grid_ids[pos_counter++];
        }

        //update phi of grid
        for (auto &grid_info : affected_grids)
        {
            Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
            Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
        }
        simbox.num_Accepts++;
    }
    //accept hopping
    else if (MTRAND() < exp(minus_delta_E))
    {
        int pos_counter = 0;
        for (auto &target_bead : *target_beads)
        {
            std::vector<double> position_after = molecule_pos[pos_counter];
            std::copy(position_after.begin(), position_after.end(), target_bead->position.begin());
            target_bead->grid_id = molecule_new_grid_ids[pos_counter++];
        }

        //update phi of grid
        for (auto &grid_info : affected_grids)
        {
            Remove_Bead_from_Grid(&grids[grid_info[0]], grid_info[2]);
            Add_Bead_to_Grid(&grids[grid_info[1]], grid_info[2]);
        }
        simbox.num_Accepts++;
    }
    //else reject
    else simbox.num_Rejects++;

    return target_beads->size()-1;
}
*/