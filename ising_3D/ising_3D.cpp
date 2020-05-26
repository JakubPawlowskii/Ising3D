#include <iomanip>

#include "ising_framework.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <sstream>
constexpr auto spin_up_probability = 1.0;

int main()
{
	//T* =  4.51152786
	vector<double> T;
	vector<int> lattice_dim = {4, 8, 12, 16, 20};
	auto temp_T = 0.0;
	do
	{
		T.push_back(temp_T);
		temp_T += 0.01;
	}
	while (temp_T <= 6.0);

#pragma omp parallel for
	for (unsigned int index = 0; index < lattice_dim.size(); index++)
		for (auto& temp : T)
		{
			ising_framework ising(lattice_dim[index], spin_up_probability, temp);

			cout << "Thread number: " << omp_get_thread_num() << " running simulation with L =  " << lattice_dim[index]
				<< " and T = " << temp << endl;

			std::stringstream ss;
			ss << std::fixed << std::setprecision(2) << temp;
			auto temp_string = ss.str();
			auto output_name = "mag//magnetization_L_" + to_string(lattice_dim[index]) + "_T_" + temp_string +".txt";
			auto final_conf_output_name = "conf//final_conf_L_" + to_string(lattice_dim[index]) + "_T_" + temp_string +".txt";

			
			ofstream output;
			output.open(output_name);
			if(!output.is_open()) exit(EXIT_FAILURE);
			
			
			for (unsigned int i = 0; i < 30000; i++) ising.simulate_mcs();

			for (unsigned int i = 0; i < 2000; i++) //should be 2000
			{
				for (unsigned int j = 0; j < 100; j++) ising.simulate_mcs(); //should be 100
				output<<ising.magnetization();
			}

			output.close();
			ising.print_lattice_to_file(final_conf_output_name);
		}

	return 0;
}
