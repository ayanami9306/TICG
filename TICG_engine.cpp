#include "TICG.h"

void TICG::TICG_Simulation()
{
    char pdbname[200], csvname[200];
    sprintf(pdbname, "%s.pdb", simname);
    sprintf(csvname, "%s.csv", simname);
    writePDB(pdbname, true);
    writeEnergyProfile(csvname, true);
    const double temp_chiN = simbox.chiN;
    simbox.chiN = 0;
    // chi is 0 : temperature is very high
    simbox.Cur_Cycle_Annealing += simbox.num_Output_Interval;
    for (;simbox.Cur_Cycle_Annealing <= simbox.Max_Cycle_Annealing; simbox.Cur_Cycle_Annealing += simbox.num_Output_Interval)
    {
        MonteCarlo(simbox.num_Output_Interval);
        AdjustCM();
        char outname[200];
        sprintf(outname, "%s_Free_%08d.dat", simname, simbox.Cur_Cycle_Annealing);
        WriteoutSim(outname);
        writePDB(pdbname, false);
        writeEnergyProfile(csvname, false);
        std::cout << simbox.Cur_Cycle_Annealing << std::endl;
    }

    simbox.chiN = temp_chiN;
    // chi is not 0 : specific temperature
    simbox.Cur_Cycle_Simulation += simbox.num_Output_Interval;
    for (;simbox.Cur_Cycle_Simulation <= simbox.Max_Cycle_Simulation; simbox.Cur_Cycle_Simulation += simbox.num_Output_Interval)
    {
        MonteCarlo(simbox.num_Output_Interval);
        AdjustCM();
        char outname[200];
        sprintf(outname, "%s_Anneal_%08d.dat", simname, simbox.Cur_Cycle_Simulation);
        WriteoutSim(outname);
        writePDB(pdbname, false);
        writeEnergyProfile(csvname, false);
        std::cout << simbox.Cur_Cycle_Simulation << std::endl;
    }
}

void TICG::MonteCarlo(int iteration)
{
    std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start_time;
    int iter = simbox.num_Wiggle_Interval;
    int iter_hopping = simbox.num_beads * simbox.num_Wiggle_Interval;
    for (;iter <= iteration; iter += simbox.num_Wiggle_Interval)
    {
        for (int num_hopping_in_cycle = 0; num_hopping_in_cycle < iter_hopping; num_hopping_in_cycle++)
        {   
            TranslateBead();
            //double val_hopping_prob = MTRAND();
            //if (val_hopping_prob < simbox.prob_AtomTranslation) TranslateBead();
            //else if (val_hopping_prob < simbox.prob_MoleculeTranslation) num_hopping_in_cycle += TranslateMolecule();
            //else num_hopping_in_cycle += RotateMolecule();
        }
        WiggleBox();
        AdjustCM();
    }
}