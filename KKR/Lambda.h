#pragma once


#include <Eigen/eigen>

#include <vector>

#include "Vector3D.h"

#include "SpecialFunctions.h"


#define _USE_MATH_DEFINES
#include <math.h>


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

		bool IsCloseToPole(double E, const Vector3D<double>& k, double limit, const std::vector<double>& ratios, double limit2 = 1E-10) const
		{			
			// singular points of the free Green function:

			for (const auto& Kn : m_basisVectors)
			{
				const Vector3D<double> kn = Kn + k;
				const double kn2 = kn * kn;
				const double Eminuskn2 = 2. * E - kn2;

				if (abs(Eminuskn2) < limit)
					return true;
			}

			// also for cotangent of the phase shift:

			const std::complex<double> kappa((E >= 0 ? sqrt(2. * E) : 0), (E < 0 ? sqrt(-2. * E) : 0));
			const std::complex<double> kappaR = kappa * m_R;
		
			const int lMax = ratios.size();
			for (int l = 0; l <= lMax; ++l)
			{
				const double logDeriv = ratios[l] - m_oneOverR;

				const std::complex<double> check = SpecialFunctions::Bessel::jderiv(l, kappaR) * kappa - SpecialFunctions::Bessel::j(l, kappaR) * logDeriv;

				if (abs(check) < limit2)
					return true;
			}

			return false;
		}


		std::complex<double> D(double E, const Vector3D<double>& k, int L, int M, const Coefficients& coeffs) const
		{
			//const double eta = 1.56;
			const double eta = 4. * M_PI / std::pow(m_cellVolume, 2. / 3.);

			const std::complex<double> kappa((E >= 0 ? sqrt(2. * E) : 0), (E < 0 ? sqrt(-2. * E) : 0));

			const std::complex<double> kappamL = std::pow(kappa, -L);
			const std::complex<double> I(0, 1);
			const std::complex<double> iL = std::pow(I, L);
			const double EpEta = 2. * E / eta;


			/*

			// this is an example of how it's done in Kohn and Rostoker paper
			// without the full Ewald summation
			// you don't need the real space vectors anymore
			// but more vectors are needed in reciprocal space
			// so comment out the GenerateBasisVectorsMaxSize(3, 1) line in BandStructureBasis
			// and uncomment one that generates more vectors, like GenerateBasisVectorsMaxSize(8, 3) or higher
			// the convergence is not very good especially for big energies
			// I didn't play enough with this, so eta is probably far from optimum
			// it also probably interferes with the code that tries to avoid spurious results due of singularities and so on


			std::complex<double> D(0, 0);
			for (const auto& Kn : m_basisVectors)
			{
				const Vector3D<double> kn = Kn + k;
				const double kn2 = kn * kn;
				const double kn_length = sqrt(kn2);
				const double Eminuskn2 = 2. * E - kn2;

				const double theta = kn.getTheta();
				const double phi = kn.getPhi();

				const std::complex<double> Y = SpecialFunctions::Legendre::Y(L, M, theta, phi);

				D += SpecialFunctions::Bessel::j(L, kn_length * m_R) / Eminuskn2 * Y * std::exp(Eminuskn2 / 5000.);
			}
			
			const std::complex<double> kappaR = kappa * m_R;
			
			D *= 4. * M_PI / (m_cellVolume * SpecialFunctions::Bessel::j(L, kappaR));

			if (0 == L)
			{
				assert(0 == M);

				D += kappa / (4. * M_PI * std::tan(kappaR));
			}

			return D;

			*/

			// The three terms for Ewald summation:

			// **************** first term ******************************************************************************************

			std::complex<double> D1(0, 0);
			for (const auto& Kn : m_basisVectors)
			{
				const Vector3D<double> kn = Kn + k;
				const double kn2 = kn * kn;
				const double kn_length = sqrt(kn2);
				const double Eminuskn2 = 2. * E - kn2;

				const double theta = kn.getTheta();
				const double phi = kn.getPhi();

				const std::complex<double> Y = SpecialFunctions::Legendre::Y(L, M, theta, phi);

				D1 += std::pow(kn_length, L) * std::exp(-kn2 / eta) / Eminuskn2 * Y;
			}

			D1 *= 4. * M_PI / m_cellVolume * kappamL * std::exp(EpEta);


			// **************** second term ******************************************************************************************

			std::complex<double> D2(0, 0);

			double integral = 0;
			double oldLength2 = 0;

			for (const auto& Rn : m_realVectors)
			{
				const double rs2 = Rn * Rn;
				const double rs = sqrt(rs2);

				// the vectors are sorted by length
				// calculate the integral only if the length changes
				if (rs2 > oldLength2 + std::numeric_limits<double>::epsilon())
				{
					integral = 0;
					for (int m = 0; m < 16; ++m)
					{
						const double term = std::pow(E * rs2 / 2., m) / coeffs.Factorial(m) * SpecialFunctions::Gamma(0.5 + L - m, rs2 * eta / 4.);
						integral += term;
						if (abs(term) < 1E-13) break;
					}
					integral *= 0.5 / std::pow(rs, 2. * L + 1.);
				}

				const double theta = Rn.getTheta();
				const double phi = Rn.getPhi();

				const std::complex<double> Y = SpecialFunctions::Legendre::Y(L, M, theta, phi);

				D2 += std::pow(rs, L) * std::exp(I * (k * Rn)) * Y * integral;

				oldLength2 = rs2;
			}

			D2 *= 1. / sqrt(M_PI) * std::pow(-2, L + 1) * iL * kappamL;

			// **************** third term ******************************************************************************************

			double D3 = 0;
			if (0 == L)
			{
				assert(0 == M);

				for (int s = 0; s < 16; ++s)
				{
					const double term = std::pow(EpEta, s) / (2. * s - 1.) / coeffs.Factorial(s);
					D3 += term;
					if (abs(term) < 1E-13) break;
				}

				D3 *= -0.5 * sqrt(eta) * M_1_PI;
			}

			return D1 + D2 + D3;
		}

		void Compute(double E, const Vector3D<double>& k, const std::vector<double>& ratios, const Coefficients& coeffs)
		{
			const std::complex<double> kappa((E >= 0 ? sqrt(2. * E) : 0), (E < 0 ? sqrt(-2. * E) : 0));
			const std::complex<double> kappaR = kappa * m_R;

			// precalculate D values
			std::map<std::tuple<int, int>, std::complex<double>> Dmap;
			for (int L = 0; L <= 2 * m_lMax; ++L)
			{
				// first compute for the non negative M
				for (int M = 0; M <= L; ++M)
				{
					const std::complex<double> Dval = D(E, k, L, M, coeffs);
					if (Dval != std::complex<double>(0, 0))
						Dmap[std::make_tuple(L, M)] = Dval;
				}

				// for the negative M, use the non negative value already calculated

				for (int M = -L; M < 0; ++M)
				{
					auto ind = std::make_tuple(L, -M);
					const auto it = Dmap.find(ind);
					if (it != Dmap.end())
						Dmap[std::make_tuple(L, M)] = pow(-1, -M) * std::conj(it->second);
				}
			}

			const std::complex<double> I(0, 1);

			int i = 0; // the index for l, m
			for (int l = 0; l <= m_lMax; ++l)
				for (int m = -l; m <= l; ++m)
				{
					int j = i; // the index for lp, mp

					for (int lp = l; lp <= m_lMax; ++lp)
						for (int mp = (lp == l ? m : -lp); mp <= lp; ++mp)
						{
							std::complex<double> A(0, 0);

							for (int L = abs(l - lp); L <= l + lp; L += 2)
							{
								//const int sumL = L + l + lp;
								//if (sumL % 2) continue; // must be an even integer, no need to check it, the for above ensures the sum is even

								const double C = coeffs.getCoefficient(l, lp, L, m, mp);
								if (C != 0.)
								{
									const auto it = Dmap.find(std::make_tuple(L, m - mp));
									if (it != Dmap.end())
										A += it->second * C;
								}
							}
							A *= 4. * M_PI * std::pow(I, l - lp);

							if (i == j)
							{
								const double logDeriv = ratios[l] - m_oneOverR;

								const std::complex<double> ctgPhaseShift = (SpecialFunctions::Bessel::nderiv(l, kappaR) * kappa - SpecialFunctions::Bessel::n(l, kappaR) * logDeriv) / (SpecialFunctions::Bessel::jderiv(l, kappaR) * kappa - SpecialFunctions::Bessel::j(l, kappaR) * logDeriv);
								Lmat(i, i) = A + kappa * ctgPhaseShift;
							}
							else
							{
								Lmat(i, j) = A;
								Lmat(j, i) = std::conj(A);
							}
						
							++j;
						}

					++i;
				}
		}

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