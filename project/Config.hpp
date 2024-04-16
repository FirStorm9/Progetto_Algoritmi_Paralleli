#include <fstream>
#include <iostream>
#include <thread>
#include "json.hpp"
#include <unistd.h>
#include <iostream>

using json = nlohmann::json;

struct Config {

    //
    // Beginnin of default configuration
    //

    // number of column and row of the matrix
    int n_col = 200;
    int n_row = 200;

    // dimention of the display of Allegro
    int display_dim = 800;

    int mpi_root = 0;

    // number of steps that the program will do
    int steps = 30;

	// number of threads per column
	int mpi_threads = 4;

	// number of threads per row
	int posix_threads = 2;

    // 0: defined matrix, 1: random matrix from file, 2: matrix from file
    int selection = 0;

    //
    // End of default configuration
    //

	Config() {

		std::ifstream config("../config.cfg");
		if(config.is_open()) {
			std::cout << "Picking external file configuration file..." << std::endl;
		}else{
            std::cout << "Picking default configuration" << std::endl;
            return;
        }

		json jsonConfig;
		config >> jsonConfig;

        if(jsonConfig.contains("mpi_threads")) mpi_threads = jsonConfig["mpi_threads"];
 
        unsigned int max_processor_number_in_pc = std::thread::hardware_concurrency();
        if(mpi_threads <= 0 && mpi_threads >= max_processor_number_in_pc){
            std::cout << "In configuration file the number of mpi_threads (" << mpi_threads << ") isn't valid\nThe maximum possible is: " << max_processor_number_in_pc << std::endl;
            exit(1);
        }

        if(jsonConfig.contains("n_col")) n_col = jsonConfig["n_col"];
		if(jsonConfig.contains("n_row")) n_row = jsonConfig["n_row"];

        if(jsonConfig.contains("posix_threads")) posix_threads = jsonConfig["posix_threads"];
        
		if(jsonConfig.contains("mpi_root")) mpi_root = jsonConfig["mpi_root"];

		if(jsonConfig.contains("steps")) steps = jsonConfig["steps"];

		if(jsonConfig.contains("display_dim")) display_dim = jsonConfig["display_dim"];
		
        if(jsonConfig.contains("selection")) selection = jsonConfig["selection"];
        
    }
};
