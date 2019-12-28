#pragma once

#define _USE_MATH_DEFINES
//#include <math.h>

// you can use boost for the same purpose if spherical Bessel functions are not available


#include <cmath>
#include <complex>

namespace SpecialFunctions
{
	// for now I'll let the ones implemented here

	// the recursion formulae are more clear 
	// I first implemented the code directly with them
	// I'll let them here commented out

	class SpecialFunctions
	{
	protected:

		template<typename T> static T GetRecursiveValue(unsigned int l, const T& x, T& v0, T& v1)
		{
			if (0 == l) return v0;
			else if (1 == l) return v1;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const T v2 = (2. * i - 1.) / x * v1 - v0;

				v0 = v1;
				v1 = v2;
			}

			return v1;
		}

	};

	class Bessel : public SpecialFunctions
	{
	public:
		template<typename T> static T j(unsigned int l, const T& x)
		{
			const T sinxx = sin(x) / x;
			const T cosxx = (0 == l) ? 0. : cos(x) / x;

			T j0 = sinxx;
			T j1 = (0 == l) ? 0. : sinxx / x - cosxx;

			return GetRecursiveValue(l, x, j0, j1);
		}

		template<typename T> static T jderiv(unsigned int l, const T& x)
		{
			//return (j(l, x) - j(l, x - .00001)) / 0.00001;
			
			//return T(l) / x * j(l, x) - j(l + 1, x);

			// the above formula is ok, but this gets in one calculation both values needed for the derivative
			const T sinxx = sin(x) / x;
			const T cosxx = cos(x) / x;

			T j0 = sinxx;
			T j1 = sinxx / x - cosxx;

			GetRecursiveValue(l + 1, x, j0, j1);

			return T(l) / x * j0 - j1;
		}

		template<typename T> static T n(unsigned int l, const T& x)
		{
			const T cosxx = cos(x) / x;
			const T sinxx = (0 == l) ? 0. : sin(x) / x;

			T n0 = -cosxx;
			T n1 = (0 == l) ? 0. : -cosxx / x - sinxx;

			return GetRecursiveValue(l, x, n0, n1);
		}

		template<typename T> static T nderiv(unsigned int l, const T& x)
		{
			//return (n(l, x) - n(l, x - .00001)) / 0.00001;
			
			//return T(l) / x * n(l, x) - n(l + 1, x);

			// the above formula is ok, but this gets in one calculation both values needed for the derivative
			const T cosxx = cos(x) / x;
			const T sinxx = sin(x) / x;

			T n0 = -cosxx;
			T n1 = -cosxx / x - sinxx;

			GetRecursiveValue(l + 1, x, n0, n1);

			return T(l) / x * n0 - n1;
		}
	};

	// Legendre polynomials
	class Legendre
	{
	public:
		template<typename T> static T p(unsigned int l, const T& x)
		{
			if (0 == l) return 1.;
			else if (1 == l) return x;

			//return ((2. * l - 1.) * x * p(l - 1, x) - (l - 1) * p(l - 2, x)) / l;

			T p0 = 1.;
			T p1 = x;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const std::complex<double> p2 = ((2. * i - 1.) * x * p1 - (i - 1.) * p0) / T(i);
				p0 = p1;
				p1 = p2;
			}

			return p1;
		}

		template<typename T> static T pderiv(unsigned int l, const T& x)
		{
			return -1 / sqrt(1 - x * x) * p(l, x);
		}

		static std::complex<double> Y(unsigned int l, int m, double theta, double phi)
		{
			if (m < 0) return pow(-1, m) * std::conj(Y(l, -m, theta, phi));

			return std::sph_legendre(l, m, theta) * exp(std::complex<double>(0, m) * phi);
		}

		static double YReal(unsigned int l, int m, double theta, double phi)
		{
			static const double pre = 1. / sqrt(2.);

			if (0 == m) return Y(l, 0, theta, phi).real();
			else if (m > 0)
				return pre * (Y(l, m, theta, phi) + Y(l, -m, theta, phi)).real();

			return pre * (std::complex<double>(0, -1) * (Y(l, -m, theta, phi) - Y(l, m, theta, phi))).real();
		}
	};



#ifndef _GAMMA
#define _GAMMA
	inline double Gamma(double a, double x, int iter = 15)
	{
		double res = (iter + 1.) * (a - iter - 1.) / (2. * iter + 3. + x - a);

		// backwards
		for (int k = iter; k > 0; --k)
			res = k * (a - k) / (2. * k + 1. + x - a + res);

		return exp(-x) * pow(x, a) / (x - a + 1. + res);
	}
#endif
}
