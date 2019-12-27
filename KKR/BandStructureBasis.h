#pragma once

#include <vector>

#include "Vector3D.h"

#include "SymmetryPoints.h"

namespace KKR
{

class BandStructureBasis
{
public:
	// pass rmax as you want it or negative to have it computed as for touching spheres
	// a for Cu: 6.8219117
	// a for Al: 4.046 / 0.5291772106712 (conversion from Angstroms to Bohrs)
	// for Al the muffin should be smaller, not touching as in Cu case. A 0.39 multiplication seems ok.
	// a for Au: 4.065 / 0.5291772106712
	BandStructureBasis(double a = 6.8219117 /*4.046 / 0.5291772106712*/, double rmax = -1./*sqrt(2.) * 4.046 / 0.5291772106712 / 4. * 0.39*/); // commented out, some 'experimental' values for Al

	SymmetryPoints symmetryPoints;

	std::vector<unsigned int> symmetryPointsPositions;

	virtual void Initialize(std::vector<std::string> path, unsigned int nrPoints = 600);

	unsigned int GetPointsNumber() const { return static_cast<unsigned int>(kpoints.size()); }

	const std::vector<std::string>& GetPath() const { return m_path; }

protected:
	std::vector<std::string> m_path;

	std::vector<Vector3D<double>> basisVectors;
	std::vector<Vector3D<double>> realVectors;

	std::vector<Vector3D<double>> kpoints;

	double m_a;
	double m_Rmax;

	bool GenerateBasisVectorsMaxSize(int maxSize, int realMaxSize);
};

}

