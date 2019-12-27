#pragma once

#include <string>
#include <vector>

namespace KKR {

	class Pseudopotential
	{
	public:
		Pseudopotential();
		~Pseudopotential();

		bool Load(const std::string& name);

		double Value(double x) const;

		void Clear();

		unsigned int GetZ() const { return Z; }
		unsigned int GetZion() const { return Zion; }
		unsigned int ElectronsInCore() const { return Z - Zion; }
		double GetMaxRadius() const { return maxRadius; }
		bool IsValid() const { return valid; }

		static double VeffCu(double R)
		{
			const double R2 = R * R;
			const double R3 = R2 * R;
			const double R4 = R2 * R2;

			return 29. * exp(-2.3151241717834 * pow(R, 0.81266614122432) + 2.1984250222603E-2 * pow(R, 4.2246376280056))
				- 0.15595606773483 * R - 3.1350051440417E-3 * R2 + 5.1895222293006e-2 * R3 - 2.8027608685637E-2 * R4;
		}

		// not used but might turn useful with the proper params
		// use it as -StarkloffJoannopoulos(...) / r as the one for Cu
		static double StarkloffJoannopoulos(int Z, double r, double lambda, double rc)
		{
			return Z * (1. - exp(-lambda * r)) / (1. + exp(-lambda * (r - rc)));
		}


		// the dumbest possible
		// constant for R < Rl
		// -Z/r otherwise (Z given by the number of valence electrons)
		// don't add - or /r for this one
		static double VDumb(int Z, double R, double Rl, double C)
		{
			if (R < Rl) return C;

			return -Z / R;
		}


		static double VAl(double R)
		{
			return VDumb(3, R, 2.675, 1.3905 * 0.5);
		}

	protected:
		void ComputeSpline();

		double Interpolate(size_t interval, double x) const;
		size_t GetIndex(double x) const;

		unsigned int Z;
		unsigned int Zion;

		bool valid;

		double maxRadius;

		std::vector<double> position;
		std::vector<double> pseudopotential;

		std::vector<double> h;
		std::vector<double> z;
	};

}