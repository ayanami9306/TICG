#include "TICG.h"

double TICG::Calc_Bonded_Potential_for_Bead(std::vector<double> * pos, std::vector<Bond> * bondlist)
{
    double bonded_potential = 0;
    //H_b/kbT = (3/2)*sigma((r_ij)^2/b^2) (if connected)
    for (auto &target_bead_bond : *bondlist)
    {
        std::vector<double> *bonded_bead_pos = &(beads[target_bead_bond.connected_bead_id].position);
        bonded_potential += Calc_Bond_Length2(pos, bonded_bead_pos) * simbox.ideal_inverse_bond_length_square_for_bond_type[target_bead_bond.bond_type];
    }
    return 1.5 * bonded_potential;
}

double TICG::Calc_Nonbonded_Potential_for_Grid(Grid * target_grid)
{
    /*
    H_nb/kbT = sqrt(N_bar) * integral_V dr/(Re^3)*(chiN*phi_A*phi_B + 0.5*kappaN*(phi_A+phi_B)^2)
    sqrt(N_bar) = rho(particle number density)*Re^3/N
    sqrt(N_bar) /Re^3 = n(number of molecule)/V
    for discrete grid system, H_nb/kbT = (n/V) * sigma (chiN*phi_A*phi_B + 0.5*kappaN*(phi_A+phi_B)^2)*dV
    dV is same for this system, and dV = V/(number of grids)
    so, H_nb/kbT = (n/(number of grids))*sigma (chiN*phi_A*phi_B + 0.5*kappaN*(phi_A+phi_B)^2)
    for 1 grid box, H'_nb = (n/(number of grids))*(chiN*phi_A*phi_B + 0.5*kappaN*(phi_A+phi_B)^2)
    */
    std::vector<double> * grid_phi = &(target_grid->phi);
    double sum_phi = (*grid_phi)[0] + (*grid_phi)[1] - 1;
    double mul_phi = (*grid_phi)[0] * (*grid_phi)[1];
    //double sum_phi = std::accumulate(grid_phi->begin(), grid_phi->end(), 0);
    //double mul_phi = std::accumulate(grid_phi->begin(), grid_phi->end(), 1.0, std::multiplies<double>());
    return (simbox.chiN * mul_phi + simbox.half_kappaN * sum_phi * sum_phi) * simbox.const_param_for_Hnb;
}

double TICG::Calc_Difference_Nonbonded_Potential_for_Grid(Grid * target_grid, double delta_A, double delta_B)
{
    std::vector<double> *grid_phi = &(target_grid->phi);
    double before_sum_phi = (*grid_phi)[0] + (*grid_phi)[1] - 1,
        before_mul_phi = (*grid_phi)[0] * (*grid_phi)[1];
    double after_sum_phi = before_sum_phi + delta_A + delta_B,
    after_mul_phi = ((*grid_phi)[0] + delta_A) * ((*grid_phi)[1] + delta_B);

    return (simbox.chiN * (before_mul_phi - after_mul_phi) + simbox.half_kappaN * (before_sum_phi * before_sum_phi - after_sum_phi * after_sum_phi)) * simbox.const_param_for_Hnb;
}

double TICG::Calc_Difference_Bonded_Potential_for_Bead(std::vector<double> *before_pos, std::vector<double> *after_pos, std::vector<Bond> * bondlist)
{
    double bonded_potential = 0;
    //H_b/kbT = (3/2)*sigma((r_ij)^2/b^2) (if connected)
    for (auto &target_bead_bond : *bondlist)
    {
        std::vector<double> *bonded_bead_pos = &(beads[target_bead_bond.connected_bead_id].position);
        bonded_potential += (Calc_Bond_Length2(before_pos, bonded_bead_pos) - Calc_Bond_Length2(after_pos, bonded_bead_pos)) * simbox.ideal_inverse_bond_length_square_for_bond_type[target_bead_bond.bond_type];
    }
    return 1.5 * bonded_potential;
}

void TICG::Measure_CM(std::vector<Bead *> * beadlist, std::vector<double> * CM)
{
    *CM = std::vector<double> {0.0, 0.0, 0.0};
    for (auto &target_bead : *beadlist)
    {
        std::vector<double> & pos = target_bead->position;
        (*CM)[0] += pos[0];
        (*CM)[1] += pos[1];
        (*CM)[2] += pos[2];
    }
    (*CM)[0] /= (*beadlist).size();
    (*CM)[1] /= (*beadlist).size();
    (*CM)[2] /= (*beadlist).size();
}

std::vector<std::string> TICG::strsplit(std::string input, char delimiter) {
    std::vector<std::string> strlist;
    std::stringstream iss(input);
    std::string temp;

    //left trim
    input.erase(0, input.find_first_not_of(" \t\n\r\f\v"));
    input.erase(input.find_last_not_of(" \t\n\r\f\v") + 1);
    //right trim
 
    while (getline(iss, temp, delimiter)) {
        strlist.push_back(temp);
    }
 
    return strlist;
}