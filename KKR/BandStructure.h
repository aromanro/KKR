#pragma once

#include <vector>
#include <atomic>
#include <future>

#include "Lambda.h"

#include "Vector3D.h"

#include "Options.h"

#include "BandStructureBasis.h"
#include "Numerov.h"
#include "Pseudopotential.h"

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

		void Initialize(std::vector<std::string> path, unsigned int nrPoints = 600) override
		{
			BandStructureBasis::Initialize(path, nrPoints);
			results.clear();
		};


		std::vector<std::vector<double>> Compute(const std::atomic_bool& terminate, const Options& options);

	protected:
		void ComputeSchrodinger(std::vector<std::future<void>>& tasks, Potential& potential, std::vector<std::vector<double>>& ratios, int numIntervals, int numerovGridNodes, int numerovIntervals, double deltaGrid, double minE, double dE, int lMax, const std::atomic_bool& terminate, const Options& options);
		void GetResult(std::vector<std::vector<double>>& res, std::vector<std::vector<double>>& ratios, Lambda& lambda, int k, double E, double posE, double dE, double det, double oldDet, double olderDet, double detLim, double ctgLimit, double smallMinLimit, int lMax);

		static void SetPotential(Potential& potential, int numerovGridNodes, double Rp, double deltaGrid);
		static double LinearInterpolation(double E, double dE, double det, double oldDet);
		static double QuadraticInterpolation(double E, double dE, double det, double oldDet, double olderDet);
	};

}

