#include "ising_framework.h"

using namespace std;

random_device device;

static inline uint64_t rotl(const uint64_t x, const int k)
{
	return (x << k) | (x >> (64 - k));
}


static uint64_t s[4] = {
	(static_cast<uint64_t>(device()) << 32) | device(), (static_cast<uint64_t>(device()) << 32) | device(),
	(static_cast<uint64_t>(device()) << 32) | device(), (static_cast<uint64_t>(device()) << 32) | device()
};


uint64_t next()
{
	const uint64_t result = s[0] + s[3];

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

static double to_double(uint64_t x)
{
	const union
	{
		uint64_t i;
		double d;
	} u = {UINT64_C(0x3FF) << 52 | x >> 12};
	return u.d - 1.0;
}


unsigned int ising_framework::index(const unsigned int x, const unsigned int y, const unsigned int z) const
{
	return x + lattice_dim_ * (y + lattice_dim_ * z);
}

void ising_framework::initialize_lattice() const
{
	for (unsigned int i = 0; i < lattice_vol_; i++) lattice_[i] = to_double(next()) <= initial_spin_up_ ? 1 : -1;
}

void ising_framework::print_lattice_to_file(const string file_name) const
{
	ofstream lattice_file;
	lattice_file.open(file_name);
	for (unsigned i = 0; i < lattice_dim_; i++)
		for (unsigned j = 0; j < lattice_dim_; j++)
			for (unsigned k = 0; k < lattice_dim_; k++)
				lattice_file << i << ";" << j << ";" << k << ";" << lattice_[index(i, j, k)] << endl;

	lattice_file.close();
}

void ising_framework::spin_swap_energy_difference(const unsigned int x, const unsigned int y, const unsigned int z)
{
	dE_ = 2 * lattice_[index(x, y, z)] * (lattice_[index(pbc_forward_[x], y, z)] + lattice_[index(
			pbc_backward_[x], y, z)] +
		lattice_[index(x, pbc_forward_[y], z)] + lattice_[index(x, pbc_backward_[y], z)] + lattice_[
			index(x, y, pbc_forward_[z])] + lattice_[index(x, y, pbc_backward_[z])]);
}

double ising_framework::swap_probability() const
{
	switch (dE_)
	{
	case -12: return boltzmann_factor_[0];

	case -8: return boltzmann_factor_[1];

	case -4: return boltzmann_factor_[2];

	case 0: return boltzmann_factor_[3];

	case 4: return boltzmann_factor_[4];

	case 8: return boltzmann_factor_[5];

	case 12: return boltzmann_factor_[6];
	default:
		{
			cout<<endl<<"Unhandled dE value."<<endl;
			return 1;
		}
	}
}

double ising_framework::magnetization() const
{
	auto sum = 0.0;
	for (unsigned int i = 0; i < lattice_vol_; i++) sum += lattice_[i];
	return 1 / static_cast<double>(lattice_vol_) * sum;
}

void ising_framework::simulate_mcs()
{
	for (unsigned i = 0; i < lattice_dim_; i++)
		for (unsigned j = 0; j < lattice_dim_; j++)
			for (unsigned k = 0; k < lattice_dim_; k++)
			{
				spin_swap_energy_difference(i, j, k);
				if (dE_ <= 0 || to_double(next()) < swap_probability()) lattice_[index(i, j, k)] *= -1;
			}
}
