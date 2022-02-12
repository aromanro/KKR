#pragma once

#include <Eigen/eigen>
#include <vector>

#include "SpecialFunctions.h"
#include "Coefficients.h"

#include "Vector3D.h"

namespace KKR
{
	class Lambda
	{
	public:
		Lambda(const std::vector<Vector3D<double>>& basisVectors, const std::vector<Vector3D<double>>& realVectors, double R, double cellVolume, unsigned int lmax = 4)
			: m_basisVectors(basisVectors), m_realVectors(realVectors), m_R(R), m_oneOverR(1. / R), m_cellVolume(cellVolume), m_lMax(lmax)
		{
			unsigned int dim = m_lMax + 1;
			dim *= dim;

			Lmat.resize(dim, dim);
		}

		bool IsCloseToPole(double E, const Vector3D<double>& k, double limit, const std::vector<double>& ratios, double limit2 = 1E-10) const;
		std::complex<double> D(double E, const Vector3D<double>& k, int L, int M, const Coefficients& coeffs) const;
		inline std::map<std::tuple<int, int>, std::complex<double>> ComputeDmap(double E, const Vector3D<double>& k, const Coefficients& coeffs);
		void Compute(double E, const Vector3D<double>& k, const std::vector<double>& ratios, const Coefficients& coeffs);

		std::complex<double> Determinant() const
		{
			return Lmat.determinant();
			//return Lmat.fullPivLu().determinant();
		}

	protected:
		// technically the basis vectors are Ki + k, those here are only Ki
		const std::vector<Vector3D<double>>& m_basisVectors;

		// real space vectors, needed for Ewald summation
		const std::vector<Vector3D<double>>& m_realVectors;

		const double m_R;
		const double m_oneOverR; // used often in computations, so it's here for optimizations
		const double m_cellVolume;
		const int m_lMax;

		Eigen::MatrixXcd Lmat;
	};

}