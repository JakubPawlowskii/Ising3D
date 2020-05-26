#pragma once
#include <random>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

//Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]

class ising_framework
{
private:
	unsigned int initial_spin_up_ = 0;
	unsigned int lattice_dim_ = 0;
	unsigned int lattice_vol_ = 0;
	int* lattice_ = nullptr;
	int dE_ = 0;
	double T_ = 0.0;
	unsigned int* pbc_forward_; // periodic boundary conditions in forward direction
	unsigned int* pbc_backward_; // periodic boundary conditions in backward direction
	double boltzmann_factor_[7];
public:
	ising_framework(unsigned int n, unsigned int p, double t):
		initial_spin_up_(p),
		lattice_dim_(n),
		lattice_vol_(n * n * n),
		T_(t)
	{
		lattice_ = new int[lattice_vol_];
		initialize_lattice();

		pbc_forward_ = new unsigned int[lattice_dim_];
		pbc_backward_ = new unsigned int[lattice_dim_];

		for (unsigned int i = 0; i < lattice_dim_; i++)
		{
			pbc_backward_[i] = i - 1;
			pbc_forward_[i] = i + 1;
		}
		pbc_backward_[0] = lattice_dim_ - 1;
		pbc_forward_[lattice_dim_ - 1] = 0;

		boltzmann_factor_[0] = exp(static_cast<double>(12) / T_); // dE = -12
		boltzmann_factor_[1] = exp(static_cast<double>(8) / T_); // dE = -8
		boltzmann_factor_[2] = exp(static_cast<double>(4) / T_); // dE = -4
		boltzmann_factor_[3] = 1;									// dE = 0     all possible values of dE
		boltzmann_factor_[4] = exp(static_cast<double>(-4) / T_); // dE = 4
		boltzmann_factor_[5] = exp(static_cast<double>(-8) / T_); // dE = 8
		boltzmann_factor_[6] = exp(static_cast<double>(-12) / T_); // dE = 12
	}

	~ising_framework()
	{
		delete[] lattice_;
		delete[] pbc_backward_;
		delete[] pbc_forward_;
		lattice_ = nullptr;
		pbc_backward_ = nullptr;
		pbc_forward_ = nullptr;
	}

	unsigned int index(unsigned int x, unsigned int y, unsigned int z) const;
	void initialize_lattice() const;
	void print_lattice_to_file(const string file_name) const;
	void spin_swap_energy_difference(unsigned int x, unsigned int y, unsigned int z);
	double swap_probability() const; //TODO:: put possible values in an array
	double magnetization() const;
	void simulate_mcs();
};
