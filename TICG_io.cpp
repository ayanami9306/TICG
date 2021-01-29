#include "TICG.h"

void TICG::WriteoutSim(char * outname)
{
    std::ofstream fileobj;
    fileobj.open(outname);
    if(fileobj.is_open())
    {
        fileobj << "BOXSIZE " << simbox.Lx << " " << simbox.Ly << " " << simbox.Lz << std::endl;
        fileobj << "NUMGRID " << simbox.num_grid_x << " " << simbox.num_grid_y << " " << simbox.num_grid_z << std::endl;
        //MCTIME Cur_Cycle_Annealing Cur_Cycle_Simulation
        fileobj << "MCTIME " << simbox.Cur_Cycle_Annealing << " " << simbox.Cur_Cycle_Simulation << std::endl;

        //BONDLENGTH b_type1 b_type2 ...
        fileobj << "BONDLENGTH2 ";
        for (auto &bondlength : simbox.ideal_bond_length_square_for_bond_type) fileobj << std::fixed << std::setprecision(4) << bondlength << " ";
        fileobj << std::endl;

        //write ATOMS section
        fileobj << std::endl << "ATOMS" << std::endl;
        //bead_id bead_type X Y Z
        for (auto &bead_i : beads)
        {
            fileobj << std::fixed << std::setprecision(4) << bead_i.bead_id 
            << " " << bead_i.bead_type << " " << bead_i.position[0] << " " 
            << bead_i.position[1] << " " <<bead_i.position[2] << std::endl;
        }
        
        //write BONDS section
        fileobj << std::endl << "BONDS" << std::endl;
        int count = 1;
        //bond_id bond_type bead_id1 bead_id2
        for (auto &bead_i : beads)
        {
            int bead_id1 = bead_i.bead_id;
            for (auto &bond_i : bead_i.bondlist)
            {
                int bead_id2 = bond_i.connected_bead_id;
                if (bead_id1 < bead_id2)
                {
                    fileobj << count ++ << " " << bond_i.bond_type << " " 
                    << bead_id1 << " " << bead_id2 << std::endl;
                }
            }
        }
        fileobj.close();
    }
    else 
    {
        std::cout << "ERROR : write file is used already" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void TICG::ReadinSim(char * inname)
{
    std::ifstream fileobj;
    fileobj.open(inname);
    if (fileobj.is_open())
    {
        int linecount = 1;
        std::string Read_Flag = "";
        std::string temp_str = "";
        simbox.num_beads = 0;
        int num_bondtype = 0;
        while(!fileobj.eof())
        {
            getline(fileobj, temp_str);

            if (temp_str == "") continue;
            else if (temp_str == "ATOMS") Read_Flag = "ATOMS";
            else if (temp_str == "BONDS")
            {
                if (temp_str == "ATOMS") Read_Flag = "BONDS";
                else if (simbox.ideal_bond_length_square_for_bond_type.empty())
                {
                    std::cout << "LINE " << linecount << " : ERROR : BONDLENGTH2 is not defined before." << std::endl;
                    exit(EXIT_FAILURE);
                }
                else 
                {
                    std::cout << "LINE " << linecount << " : ERROR : ATOMS is not defined before." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if (temp_str.find("MCTIME") != std::string::npos)
            {
                std::vector<std::string> tokenlist = strsplit(temp_str, ' ');
                if (tokenlist.size() != 3)
                {
                    std::cout << "LINE " << linecount << " : ERROR : MCTIME value is not valid." << std::endl;
                    exit(EXIT_FAILURE);
                }
                try 
                {
                    simbox.Cur_Cycle_Annealing = std::stoi(tokenlist[1]);
                    simbox.Cur_Cycle_Simulation = std::stoi(tokenlist[2]);
                }
                catch (const std::exception& expn)
                {
                    std::cout << expn.what()  << "LINE " << linecount << " : ERROR : MCTIME value is not valid." << std::endl;
                    exit(EXIT_FAILURE);
                } 
            }
            else if (temp_str.find("BONDLENGTH2") != std::string::npos)
            {
                std::vector<std::string> tokenlist = strsplit(temp_str, ' ');
                if (tokenlist.size() < 2)
                {
                    std::cout  << "LINE " << linecount << " : ERROR : BONDLENGTH2 value is not valid." << std::endl;
                    exit(EXIT_FAILURE);
                }
                try 
                {
                    tokenlist.erase(tokenlist.begin());
                    for (auto strbondlength2 : tokenlist)
                    {
                        double bondlength2 = std::stod(strbondlength2);
                        simbox.ideal_bond_length_square_for_bond_type.push_back(bondlength2);
                        simbox.ideal_inverse_bond_length_square_for_bond_type.push_back(1.0 / bondlength2);
                    }
                    num_bondtype = simbox.ideal_bond_length_square_for_bond_type.size();
                }
                catch (const std::exception& expn)
                {
                    std::cout << expn.what()  << "LINE " << linecount << " : ERROR : BONDLENGTH2 value is not valid." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else 
            {
                //bead_id bead_type X Y Z
                if (Read_Flag == "ATOMS")
                {
                    std::vector<std::string> tokenlist = strsplit(temp_str, ' ');
                    if (tokenlist.size() != 5)
                    {
                        std::cout  << "LINE " << linecount << " : ERROR : ATOM value is not valid." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    try 
                    {
                        Bead temp_atom(std::stoi(tokenlist[1]), std::stoi(tokenlist[2]));
                        if (temp_atom.bead_id != simbox.num_beads) 
                        {
                            std::cout << "LINE " << linecount << " : ERROR : ATOM id is not ordered correctly." << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        temp_atom.position.push_back(std::stod(tokenlist[3]));
                        temp_atom.position.push_back(std::stod(tokenlist[4]));
                        temp_atom.position.push_back(std::stod(tokenlist[5]));

                        beads.push_back(temp_atom);
                    }
                    catch (const std::exception& expn)
                    {
                        std::cout << expn.what()  << "LINE " << linecount << " : ERROR : ATOM value is not valid." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                //bondindex bond_type bead_id1 bead_id2
                else if (Read_Flag == "BONDS")
                {
                    std::vector<std::string> tokenlist = strsplit(temp_str, ' ');
                    if (tokenlist.size() != 5)
                    {
                        std::cout  << "LINE " << linecount << " : ERROR : BOND value is not valid." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    try 
                    {
                        int bead_id1 = std::stoi(tokenlist[3]);
                        int bead_id2 = std::stoi(tokenlist[4]);

                        Bond temp_bond(-1, std::stoi(tokenlist[2]));
                        if (temp_bond.bond_type >= num_bondtype) 
                        {
                            std::cout  << "LINE " << linecount << " : ERROR : BOND type is not defined." << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        temp_bond.connected_bead_id = bead_id2;
                        beads[bead_id1].bondlist.push_back(temp_bond);
                        temp_bond.connected_bead_id = bead_id1;
                        beads[bead_id2].bondlist.push_back(temp_bond);
                    }
                    catch (const std::exception& expn)
                    {
                        std::cout << expn.what()  << "LINE " << linecount << " : ERROR : ATOM value is not valid." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
            linecount++;
        }
        simbox.num_beads = beads.size();
        fileobj.close();
    }
    else 
    {
        std::cout << "ERROR : topo file is not exist." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void TICG::ConstructSim(char * inname)
{
    simbox.Cur_Cycle_Annealing = 0;
    simbox.Cur_Cycle_Simulation = 0;
    std::ifstream fileobj;
    fileobj.open(inname);
    if (fileobj.is_open())
    {
        std::string temp_str;
        int DP = 0;
        double fraction = -1.0;
        double density_chain = 0;
        simbox.is_zwall = false;
        //keyword (value)
        while(!fileobj.eof())
        {
            getline(fileobj, temp_str);
            char propertyname[200];

            //boxdim Lx Ly Lz
            if(temp_str.find("boxdim") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf %lf %lf", propertyname, &simbox.Lx, &simbox.Ly, &simbox.Lz);
            //numgrid numgridx numgridy numgridz
            else if(temp_str.find("numgrid") != std::string::npos) sscanf(temp_str.c_str(), "%s %d %d %d", propertyname, &simbox.num_grid_x, &simbox.num_grid_y, &simbox.num_grid_z); 
            //iszwall
            else if(temp_str.find("iszwall") != std::string::npos) simbox.is_zwall = true;
            //DP valDp
            else if(temp_str.find("dp") != std::string::npos) sscanf(temp_str.c_str(), "%s %d", propertyname, &DP); 
            //fraction valfraction
            else if(temp_str.find("fraction") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf", propertyname, &fraction); 
            //densitychain val
            else if(temp_str.find("densitychain") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf", propertyname, &density_chain); 
        }

        //bondlength - single type
        simbox.ideal_bond_length_square_for_bond_type.push_back(1.0);
        simbox.ideal_inverse_bond_length_square_for_bond_type.push_back(1.0);

        //verify value
        if (simbox.Lx <= 0 || simbox.Ly <= 0 || simbox.Lz <= 0)
        {
            std::cout << "ERROR : simulation box size is not valid." << std::endl;
            exit(EXIT_FAILURE);
        }      
        else if (simbox.num_grid_x <= 0 || simbox.num_grid_y <= 0 || simbox.num_grid_z <= 0)
        {
            std::cout << "ERROR : number of grid is not valid." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if (DP <= 0)
        {
            std::cout << "ERROR : degree of polymerization is lower than 1." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if (fraction < 0 || fraction > 1)
        {
            std::cout << "ERROR : fraction value should be [0, 1]." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if (density_chain <= 0)
        {
            std::cout << "ERROR : average chain density is lower than 1." << std::endl;
            exit(EXIT_FAILURE);
        }

        //calculate grid size from input paramters
        simbox.size_grid_x = simbox.Lx / simbox.num_grid_x;
        simbox.size_grid_y = simbox.Ly / simbox.num_grid_y;
        simbox.size_grid_z = simbox.Lz / simbox.num_grid_z;

        simbox.num_grid_xy = simbox.num_grid_x * simbox.num_grid_y;
        simbox.num_grid = simbox.num_grid_xy * simbox.num_grid_z;

        for (int grid_id=0; grid_id < simbox.num_grid; grid_id++)
        {
            grids.push_back(Grid(grid_id));
            grids[grid_id].phi = std::vector<double> {0.0, 0.0};
        }

        //Re^2 = (N-1)b^2
        double Re = std::sqrt(static_cast <double>(DP - 1) * simbox.ideal_bond_length_square_for_bond_type[0] * simbox.ideal_bond_length_square_for_bond_type[0]);
        double ratio_V_Re = (simbox.Lx/Re) * (simbox.Ly/Re) * (simbox.Lz/Re);
        simbox.num_molecules = ratio_V_Re * density_chain;
        simbox.num_beads = simbox.num_molecules * DP;

        //positioning atoms - diblock copolymer
        double x = 0.5, y = 0.5, z = 0.5;
        double num_atoms_per_grid_axis = pow(density_chain * DP, (1.0/3.0));
        double space = Re / num_atoms_per_grid_axis * 0.8;

        int bead_index = 0;
        for (int chain_i = 0; chain_i < simbox.num_molecules; chain_i++)
        {
            if (x + space * (DP-1) >= simbox.Lx)
            {
                y += space;
                x = 0.5;

                if (y >= simbox.Ly)
                {
                    y = 0.5;
                    z += space;
                }

                if (z >= simbox.Lz) 
                {
                    std::cout << "ERROR : chain density is too high" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            for (int bead_i = 0; bead_i < DP; bead_i++)
            {
                if (bead_i < DP * fraction) beads.push_back(Bead(bead_index, 0));
                else beads.push_back(Bead(bead_index, 1));
                beads[bead_index].position = std::vector<double>{x, y, z};
                if (bead_i != 0)
                    beads[bead_index].bondlist.push_back(Bond (bead_index-1, 0));
                if (bead_i != (DP-1))
                    beads[bead_index].bondlist.push_back(Bond (bead_index+1, 0));
                x += space;
                bead_index++;
            }
        }
        fileobj.close();
    }
    else 
    {
        std::cout << "ERROR : topo parameter file is not exist." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void TICG::writePDB(char * outname, bool is_new_file)
{
    FILE *fp_pdb;
    if(is_new_file)
    {
        fp_pdb = fopen(outname, "w");
    }
    else fp_pdb = fopen(outname, "a");
    fprintf(fp_pdb, "\n\nMODEL%9d\n", simbox.Cur_Cycle_Annealing + simbox.Cur_Cycle_Simulation);
    //print HETATM section
    for (auto &bead : beads)
    {
        std::vector<double> pos = bead.position;
        pos[0] = float_mod(pos[0], simbox.Lx);
        pos[1] = float_mod(pos[1], simbox.Ly);
        pos[2] = float_mod(pos[2], simbox.Lz);
        fprintf(fp_pdb, "HETATM%5x %4s%c%3s %c%4d    %8.2lf%8.2lf%8.2lf%25s\n",
        bead.bead_id+1, "C", 'A'+bead.bead_type, "CAR", 'A'+bead.bead_type, bead.bead_type, 
        pos[0], pos[1], pos[2]," ");
    }
    //print CONECT section
    /*
    for (auto &bead : beads)
    {
        int check = 0;
        for (auto &bond : bead.bondlist)
        {
            if(bead.bead_id < bond.connected_bead_id)
            {
                if(!check)
                {
                    fprintf(fp_pdb, "CONECT%5d%5d", bead.bead_id+1, bond.connected_bead_id+1);
                    check=1;
                }
                else fprintf(fp_pdb, "%5d", bond.connected_bead_id+1);
            }
        }
        if(check) fprintf(fp_pdb, "\n");
    }
    */
    //print ENDMDL
    fprintf(fp_pdb, "ENDMDL\n\n");
    fclose(fp_pdb);
}

void TICG::writeEnergyProfile(char * outname, bool is_new_file)
{
    FILE *fp_csv;
    if(is_new_file)
    {
        fp_csv = fopen(outname, "w");
        fprintf(fp_csv, "TIME, Hb, Hnb, HTotal, AccRatio\n");
    }
    else fp_csv = fopen(outname, "a");
    double H_b = 0, H_nb = 0;
    for (auto &atom : beads) H_b += Calc_Bonded_Potential_for_Bead(&(atom.position), &(atom.bondlist));
    for (auto &grid : grids) H_nb += Calc_Nonbonded_Potential_for_Grid(&grid);
    H_b /= 2.0;
    fprintf(fp_csv, "%d, %lf, %lf, %lf, %lf\n", simbox.Cur_Cycle_Annealing + simbox.Cur_Cycle_Simulation, H_b, H_nb, H_b + H_nb, static_cast <double>(simbox.num_Accepts) / (static_cast <double>(simbox.num_Accepts) + static_cast <double>(simbox.num_Rejects)));
    fclose(fp_csv);
    simbox.num_Accepts = 0;
    simbox.num_Rejects = 0;
}