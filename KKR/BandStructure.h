#pragma once

#include <vector>
#include <atomic>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vector3D.h"

#include "Options.h"

#include "BandStructureBasis.h"

namespace KKR
{

	class BandStructure : public BandStructureBasis
	{
	public:
		// pass rmax as you want it or negative to have it computed as for touching spheres
		// a for Cu: 6.8219117
		// a for Al: 4.046 / 0.5291772106712 (conversion from Angstroms to Bohrs)
		// for Al the muffin should be smaller, not touching as in Cu case. A 0.39 multiplication seems ok.
		// a for Au: 4.065 / 0.5291772106712
		BandStructure(double a = 6.8219117 /*4.046 / 0.5291772106712*/, double rmax = -1. /*sqrt(2.) * 4.046 / 0.5291772106712 / 4. * 0.39*/) // commented out, some 'experimental' values for Al
			: BandStructureBasis(a, rmax)
		{
		};

		std::vector<std::vector<double>> results;

		virtual void Initialize(std::vector<std::string> path, unsigned int nrPoints = 600) override
		{
			BandStructureBasis::Initialize(path, nrPoints);
			results.clear();
		};


		std::vector<std::vector<double>> Compute(const std::atomic_bool& terminate, const Options& options);
	};

}

