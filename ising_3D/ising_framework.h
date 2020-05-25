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
	int *lattice_ = nullptr;
	int dE_ = 0;
	int si_ = 0, sj_ = 0;
	double T_ = 0.0;
	unsigned int *pbc_forward_; // periodic boundary conditions in forward direction
	unsigned int *pbc_backward_; // periodic boundary conditions in backward direction
public:
	ising_framework(unsigned int n, unsigned int p, double t):
	initial_spin_up_(p),
	lattice_dim_(n),
	lattice_vol_(n*n*n),
	T_(t)
	{
		lattice_ = new int[lattice_vol_];
		initialize_lattice();

		pbc_forward_ = new unsigned int[lattice_dim_];
		pbc_backward_ = new unsigned int[lattice_dim_];
		
		for(unsigned int i = 0; i < lattice_dim_; i++) 
		{
			pbc_backward_[i] = i-1;
			pbc_forward_[i] = i+1;
		}
		pbc_backward_[0] = lattice_dim_-1;
		pbc_forward_[lattice_dim_-1] = 0;
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
