#include "BandStructureBasis.h"

#define _USE_MATH_DEFINES
#include <math.h>


#include <algorithm>

namespace KKR
{

	BandStructureBasis::BandStructureBasis(double a, double rmax)
		: 
		m_a(a), m_Rmax(rmax)
	{
		// if the passed value was zero or negative, make them touching spheres
		if (m_Rmax <= 0)
			m_Rmax = sqrt(2.) * m_a / 4.;

		basisVectors.reserve(137);
	}

	bool BandStructureBasis::GenerateBasisVectorsMaxSize(int maxSize, int realSize)
	{
		basisVectors.clear();
		realVectors.clear();

		// the Bravais lattice is a fcc lattice
		// the reciprocal lattice is a bcc lattice 

		// the basis vectors for the reciprocal space (without the 2 pi / a)

		const Vector3D<double> b1(-1, 1, 1), b2(1, -1, 1), b3(1, 1, -1);

		// you can also get them from the Bravais lattice vectors
		// like this:

		// Bravais lattice vectors for fcc (they should be multiplied by the lattice constant):
		const Vector3D<double> a1(0, 0.5, 0.5), a2(0.5, 0, 0.5), a3(0.5, 0.5, 0);

		// the volume of the cell is a1 * (a2 % a3) which gives 1/4 (multiplied with a^3, of course)

		// reciprocal lattice (they should be multiplied by 2 pi / lattice constant):
		//const Vector3D<double> b1 = a2 % a3 / (a1 * (a2 % a3)); 		
		//const Vector3D<double> b2 = a3 % a1 / (a2 * (a3 % a1));
		//const Vector3D<double> b3 = a1 % a2 / (a3 * (a1 % a2));
		// the denominator is the volume, mentioned above

		int maxSize2 = maxSize * maxSize;

		for (int i = -maxSize; i <= maxSize; ++i)
			for (int j = -maxSize; j <= maxSize; ++j)
				for (int k = -maxSize; k <= maxSize; ++k)
				{
					const Vector3D<double> vect = b1 * i + b2 * j + b3 * k; // reciprocal lattice vector

					const double vectSquared = vect * vect;
					if (vectSquared <= maxSize2 + 0.001) // if it's under the cutoff length, add it
						basisVectors.push_back(vect);
				}

		maxSize2 = realSize * realSize;

		for (int i = -2 * realSize; i <= 2 * realSize; ++i)
			for (int j = -2 * realSize; j <= 2 * realSize; ++j)
				for (int k = -2 * realSize; k <= 2 * realSize; ++k)
				{
					if (0 == i && 0 == j && 0 == k) continue; // not needed for D(2)

					const Vector3D<double> vect = a1 * i + a2 * j + a3 * k; // lattice vector in real space

					const double vectSquared = vect * vect;
					if (vectSquared <= maxSize2 + 0.001) // if it's under the cutoff length, add it
						realVectors.push_back(vect);
				}


		return true;
	}



	void BandStructureBasis::Initialize(std::vector<std::string> path, unsigned int nrPoints)
	{
		kpoints.clear();
		kpoints.reserve(nrPoints);

		m_path.swap(path);

		GenerateBasisVectorsMaxSize(3, 1); // 27 vectors in reciprocal space, 18 in the real space

		//GenerateBasisVectorsMaxSize(5, 2); // 137 vectors in reciprocal space, 140 in real space
		//GenerateBasisVectorsMaxSize(8, 3); // 537 vectors in reciprocal space, 458 in real space
		
		// those are really overkill, they are here for tests:

		//GenerateBasisVectorsMaxSize(10, 4); // 1067 vectors in reciprocal space, 1060 in real space
		//GenerateBasisVectorsMaxSize(15, 6); // 3527 vectors in reciprocal space, 3588 in real space

		// Have them sorted by length, eases up an optimization for D(2)
		std::sort(realVectors.begin(), realVectors.end(),
			[](const auto& val1, const auto& val2) -> bool { return val1.Length() < val2.Length();  });


		const double recVectPre = 2. * M_PI / m_a;
		// adjust 'basis' vectors
		size_t size = basisVectors.size();
		for (unsigned int i = 0; i < size; ++i)
			basisVectors[i] *= recVectPre;

		// adjust real space vectors
		size = realVectors.size();
		for (unsigned int i = 0; i < size; ++i)
			realVectors[i] *= m_a;

		kpoints = symmetryPoints.GeneratePoints(m_path, nrPoints, symmetryPointsPositions);

		// adjust kpoints
		size = kpoints.size();
		for (unsigned int i = 0; i < size; ++i)
			kpoints[i] *= recVectPre;
	}

}
