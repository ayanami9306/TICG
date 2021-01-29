#include "TICG.h"

TICG::TICG(char* input_config)
{
    //initialize random
    std::random_device rd;
    rand_engine = std::mt19937(rd());
    //generate number [0, 1)
    rand_distribution = std::uniform_real_distribution<double> (0, 1);

    char input_topo[200] = "";
    Write_Trigonometric_Table();
    char is_continue = READ_CONFIG(input_config, input_topo);  
    if (is_continue == 'Y') ReadinSim(input_topo);
    else if (is_continue == 'N') ConstructSim(input_topo);
    simbox.unitval_phi = (static_cast <double>(simbox.num_grid)) / (static_cast <double>(simbox.num_beads));
    simbox.const_param_for_Hnb = (static_cast <double>(simbox.num_molecules)) / (static_cast <double>(simbox.num_grid));
    AssignMoleculeIndex();
    AssignBeadtoGrid();
    rand_beadnum = std::uniform_int_distribution<int> (0, simbox.num_beads - 1);
    rand_distance = std::uniform_real_distribution<double> (-simbox.TranslationMax, simbox.TranslationMax);
    rand_grid = std::uniform_real_distribution<double> (-simbox.size_grid_x * 0.5, simbox.size_grid_x * 0.5);
    if (simbox.is_zwall) rand_axis = std::uniform_int_distribution<int> (0, 1);
    else rand_axis = std::uniform_int_distribution<int> (0, 2);
}

TICG::~TICG()
{
    
}

void TICG::Write_Trigonometric_Table()
{
    int angle_range = 360;
    sin_table = std::vector<double> (angle_range);
    cos_table = std::vector<double> (angle_range);
    for (int angle_i=0; angle_i < angle_range; angle_i++)
    {
        sin_table[angle_i] = sin((static_cast <double>(angle_i)) / M_PI);
        cos_table[angle_i] = cos((static_cast <double>(angle_i)) / M_PI);
    }
}

char TICG::READ_CONFIG(char * input_config, char* input_topo)
{
    char is_continue;
    std::ifstream fileobj;
    fileobj.open(input_config);
    if (fileobj.is_open())
    {
        std::string temp_str;
        simbox.num_Accepts = 0;
        simbox.num_Rejects = 0;
        simbox.Cur_Cycle_Annealing = 0;
        simbox.Cur_Cycle_Simulation = 0;
        simbox.Max_Cycle_Annealing = 0;
        simbox.Max_Cycle_Simulation = 0;
        simbox.num_Wiggle_Interval = 10;
        simbox.num_Output_Interval = 1;
        while(!fileobj.eof())
        {
            getline(fileobj, temp_str);
            char propertyname[200];
            
            //keyword (value)
            //run annealstep MCstep
            if(temp_str.find("run") != std::string::npos) sscanf(temp_str.c_str(), "%s %d %d", propertyname, &simbox.Max_Cycle_Annealing, &simbox.Max_Cycle_Simulation);
            //statmc numAccepts, numRejects
            else if(temp_str.find("statMC") != std::string::npos) sscanf(temp_str.c_str(), "%s %llu %llu", propertyname, &simbox.num_Accepts, &simbox.num_Rejects);
            //continue Y/N filename
            else if(temp_str.find("continue") != std::string::npos) sscanf(temp_str.c_str(), "%s %c %s", propertyname, &is_continue, input_topo);
            //outinterval interval
            else if(temp_str.find("outinterval") != std::string::npos) sscanf(temp_str.c_str(), "%s %d", propertyname, &simbox.num_Output_Interval);
            //wiggleinterval interval{...}
            else if(temp_str.find("wiggleinterval") != std::string::npos) sscanf(temp_str.c_str(), "%s %d", propertyname, &simbox.num_Wiggle_Interval);
            //hoppingprob prob_AtomTranslation prob_MoleculeTranslation prob_MoleculeRotation
            else if(temp_str.find("hoppingprob") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf %lf %lf", propertyname, &simbox.prob_AtomTranslation, &simbox.prob_MoleculeTranslation, &simbox.prob_MoleculeRotation);
            //maxvalhopping TranslationMax RotationMax
            else if(temp_str.find("maxhopping") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf %lf", propertyname, &simbox.TranslationMax, &simbox.RotationMax);
            //outputname simname
            else if(temp_str.find("outputname") != std::string::npos) sscanf(temp_str.c_str(), "%s %s", propertyname, simname);
            //ecoeff chiN kappaN
            else if(temp_str.find("ecoeff") != std::string::npos) sscanf(temp_str.c_str(), "%s %lf %lf", propertyname, &simbox.chiN, &simbox.half_kappaN); 
        }
        
        //verify value
        if(simbox.prob_AtomTranslation < 0 || simbox.prob_MoleculeRotation < 0 || simbox.prob_MoleculeTranslation < 0 ||
        simbox.prob_AtomTranslation + simbox.prob_MoleculeRotation + simbox.prob_MoleculeTranslation <= 0)
        {
            std::cout << "ERROR : hopping probability is not valid." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if (is_continue != 'Y' && is_continue != 'N')
        {
            std::cout << "ERROR : topo option is not valid." << std::endl;
            exit(EXIT_FAILURE);
        }
        
        //normalize hopping probability
        double prob_sum = simbox.prob_AtomTranslation + simbox.prob_MoleculeTranslation + simbox.prob_MoleculeRotation;
        simbox.prob_AtomTranslation /= prob_sum;
        simbox.prob_MoleculeTranslation /= prob_sum;
        simbox.prob_MoleculeRotation /= prob_sum;
        
        simbox.prob_MoleculeTranslation += simbox.prob_AtomTranslation;
        simbox.prob_MoleculeRotation += simbox.prob_MoleculeTranslation;

        simbox.half_kappaN /= 2.0;
    }
    else 
    {
        std::cout << "ERROR : config file is not exist." << std::endl;
        exit(EXIT_FAILURE);
    }

    return is_continue;
}

void TICG::AssignMoleculeIndex()
{
    int molecule_index = 0;
    std::vector<int> molecule_id(simbox.num_beads, -1);
    for (int beadindex = 0; beadindex < simbox.num_beads; beadindex++)
    {
        if (molecule_id[beadindex] == -1)
        {
            std::queue<int> queue_floodfill;
            queue_floodfill.push(beadindex);
            Molecule tempmolecule;
            while (!queue_floodfill.empty())
            {
                int cur_bead_index = queue_floodfill.front();
                queue_floodfill.pop();
                
                if (molecule_id[cur_bead_index] == -1)
                {
                    molecule_id[cur_bead_index] = molecule_index;
                    tempmolecule.beadlist.push_back(&beads[cur_bead_index]);
                    
                    for (auto &bonded : beads[cur_bead_index].bondlist)
                    {
                        if (molecule_id[bonded.connected_bead_id] == -1) queue_floodfill.push(bonded.connected_bead_id);
                    }
                }
            }
            molecules.push_back(tempmolecule);
            molecule_index++;
        }
    }
    simbox.num_molecules = molecules.size();
}

void TICG::AssignBeadtoGrid()
{
    for (auto &bead : beads)
    {
        bead.grid_id = GridID_from_Position(&bead.position);
        Add_Bead_to_Grid(&grids[bead.grid_id], bead.bead_type);
    }
}