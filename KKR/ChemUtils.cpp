
#include "ChemUtils.h"

#include <cassert>

#include <algorithm>

namespace Chemistry {

	const char* ChemUtils::atoms[] = {
		"H", "He",
		"Li", "Be", "B", "C", "N", "O", "F", "Ne",
		"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
		"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
		"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"
	};


	unsigned int ChemUtils::GetZForAtom(const std::string& name)
	{
		for (unsigned int i = 0; i < _countof(atoms); ++i)
			if (0 == name.compare(atoms[i])) return i + 1;

		return 0;
	}


	std::string ChemUtils::GetAtomNameForZ(unsigned int Z)
	{
		assert(Z > 0);

		if (Z > _countof(atoms)) return "Unknown";

		return atoms[Z - 1];
	}

	std::string ChemUtils::GetPseudopotentialFileForZ(unsigned int Z)
	{
		std::string name = GetAtomNameForZ(Z);
		std::transform(name.begin(), name.end(), name.begin(), ::tolower);

		name += ".lda.lps";

		return "data/" + name;
	}


}