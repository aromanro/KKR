#pragma once

// the referred formulae are from the book
// Computational Physics by J M Thijssen
// isbn: 9780521833462, https://doi.org/10.1017/CBO9781139171397


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


	class Bessel
	{
	public:
		static double j(unsigned int l, double x)
		{
#ifdef USE_BETTER_BESSEL
			return std::sph_bessel(l, x);
#else
			/*
			if (0 == l) return j0(x);
			else if (1 == l) return j1(x);

			return (2. * l - 1.) / x * j(l - 1, x) - j(l - 2, x);
			*/

			if (0 == l) return sin(x) / x;
			else if (1 == l) return  sin(x) / (x * x) - cos(x) / x;

			const double sinxx = sin(x) / x;
			const double cosxx = cos(x) / x;

			double j0 = sinxx;
			double j1 = sinxx / x - cosxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double j2 = (2. * i - 1.) / x * j1 - j0;

				j0 = j1;
				j1 = j2;
			}

			return j1;
#endif
		}

		static std::complex<double> j(unsigned int l, std::complex<double> x)
		{
			if (0 == l) return sin(x) / x;
			else if (1 == l) return  sin(x) / (x * x) - cos(x) / x;

			const std::complex<double> sinxx = sin(x) / x;
			const std::complex<double> cosxx = cos(x) / x;

			std::complex<double> j0 = sinxx;
			std::complex<double> j1 = sinxx / x - cosxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const std::complex<double> j2 = (2. * i - 1.) / x * j1 - j0;

				j0 = j1;
				j1 = j2;
			}

			return j1;
		}

		template<typename T> static T jderiv(unsigned int l, T x)
		{
			return T(l) / x * j(l, x) - j(l + 1, x);
			//return (j(l, x) - j(l, x - .00001)) / 0.00001;
		}



		static double n(unsigned int l, double x)
		{
#ifdef USE_BETTER_BESSEL
			return std::sph_neumann(l, x);
#else
			/*
			if (0 == l) return n0(x);
			else if (1 == l) return n1(x);

			return (2. * l - 1.) / x * n(l - 1, x) - n(l - 2, x);
			*/

			if (0 == l) return -cos(x) / x;
			else if (1 == l) return -cos(x) / (x * x) - sin(x) / x;

			const double sinxx = sin(x) / x;
			const double cosxx = cos(x) / x;

			double n0 = -cosxx;
			double n1 = -cosxx / x - sinxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double n2 = (2. * i - 1.) / x * n1 - n0;

				n0 = n1;
				n1 = n2;
			}

			return n1;
#endif
		}

		static std::complex<double> n(unsigned int l, std::complex<double> x)
		{
			if (0 == l) return -cos(x) / x;
			else if (1 == l) return -cos(x) / (x * x) - sin(x) / x;

			const std::complex<double> sinxx = sin(x) / x;
			const std::complex<double> cosxx = cos(x) / x;

			std::complex<double> n0 = -cosxx;
			std::complex<double> n1 = -cosxx / x - sinxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const std::complex<double> n2 = (2. * i - 1.) / x * n1 - n0;

				n0 = n1;
				n1 = n2;
			}

			return n1;
		}

		template<typename T> static T nderiv(unsigned int l, T x)
		{
			return T(l) / x * n(l, x) - n(l + 1, x);
			//return (n(l, x) - n(l, x - .00001)) / 0.00001;
		}


		/*
		protected:
			inline static double j0(double x) { return sin(x) / x; }
			inline static double j1(double x) { return sin(x) / (x * x) - cos(x) / x; }

			inline static double n0(double x) { return -cos(x) / x; }
			inline static double n1(double x) { return -cos(x) / (x * x) - sin(x) / x; }
		*/
	};

	// Legendre polynomials
	class Legendre
	{
	public:
		static double p(unsigned int l, double x)
		{
#ifdef USE_BETTER_LEGENDRE
			return std::legendre(l, x);
#else
			if (0 == l) return 1.;
			else if (1 == l) return x;

			//return ((2. * l - 1.) * x * p(l - 1, x) - (l - 1) * p(l - 2, x)) / l;

			double p0 = 1.;
			double p1 = x;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double p2 = ((2. * i - 1.) * x * p1 - (i - 1.) * p0) / i;
				p0 = p1;
				p1 = p2;
			}

			return p1;
#endif
		}

		static std::complex<double> p(unsigned int l, std::complex<double> x)
		{
			if (0 == l) return 1.;
			else if (1 == l) return x;

			std::complex<double> p0 = 1.;
			std::complex<double> p1 = x;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const std::complex<double> p2 = ((2. * i - 1.) * x * p1 - (i - 1.) * p0) / static_cast<double>(i);
				p0 = p1;
				p1 = p2;
			}

			return p1;
		}

		template<typename T> static T pderiv(unsigned int l, T x)
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
