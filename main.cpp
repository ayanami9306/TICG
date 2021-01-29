#include "TICG.h"

int main(int argc, char* argv[])
{
    std::chrono::system_clock::time_point sim_start_time = std::chrono::system_clock::now();
    //argv[0] : input_config filename
    //initialize simulation
    TICG *sim = new TICG(argv[1]);
    //perform simulation
    sim->TICG_Simulation();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - sim_start_time;
    std::cout << "Elapsed Time : " << sec.count() << std::endl;
    std::cout << "Simulation End" << std::endl;
    return 0;
}