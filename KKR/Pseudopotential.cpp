#include "Pseudopotential.h"

#include <fstream>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cassert>

namespace KKR {

	Pseudopotential::Pseudopotential()
		: Z(0), Zion(0), maxRadius(0), valid(false)
	{
	}


	Pseudopotential::~Pseudopotential()
	{
	}

	bool Pseudopotential::Load(const std::string& name)
	{
		std::ifstream file(name);

		if (!file) return false;

		position.resize(0);
		position.reserve(1602);
		pseudopotential.resize(0);
		pseudopotential.reserve(1602);

		std::string line;

		bool failure = false;

		unsigned int lineno = 0;
		unsigned int cnt = 1;
		while (std::getline(file, line))
		{
			if (!line.empty())
			{
				++lineno;

				std::istringstream strstream(line);

				if (2 == lineno)
				{
					// get Z and Zion
					try
					{
						double dZ = 0;
						double dZion = 0;
						strstream >> dZ >> dZion;
						Z = static_cast<unsigned int>(dZ);
						Zion = static_cast<unsigned int>(dZion);
					}
					catch (...)
					{
						failure = true;
						break;
					}

					continue;
				}
				else if (lineno < 8) continue;

				unsigned int counter;
				double radius;
				double potential;

				try {
					strstream >> counter >> radius >> potential;
				}
				catch (...)
				{
					failure = true;
					break;
				}

				if (counter != cnt)
				{
					failure = true;
					break;
				}

				++cnt;

				position.push_back(radius);
				pseudopotential.push_back(potential);

				maxRadius = std::max(radius, maxRadius);
			}
		}

		if (position.size() >= 2)
		{
			const double increment = position.back() - *(position.end() - 2);
			maxRadius += 500*increment;
			position.push_back(maxRadius);
			pseudopotential.push_back(Zion/maxRadius);
		}

		if (!failure) 
		{
			valid = true;
			ComputeSpline();
		}
		else Clear();

		return !failure;
	}


	void Pseudopotential::ComputeSpline()
	{
		if (position.size() < 3) return;

		const unsigned int size = (unsigned int)position.size();

		std::vector<double> b(size - 1ULL);
		std::vector<double> v(size - 2ULL);
		std::vector<double> u(size - 2ULL);

		h.resize(size - 1ULL);
		z.resize(size);

		z[0] = z[size - 1ULL] = 0;

		for (unsigned int i = 0; i < size - 1; ++i)
		{
			h[i] = position[i + 1ULL] - position[i];
			b[i] = (pseudopotential[i+1ULL] - pseudopotential[i]) / h[i];
		
			if (i)
			{
				v[i - 1ULL] = 2. * (h[i] + h[i - 1ULL]);
				u[i - 1ULL] = 6. * (b[i] - b[i - 1ULL]);
			}
		}

		// Gaussian elimination for tridiagonal system

		// elimination
		std::vector<double> vbar(size - 2);

		vbar[0] = v[0];
		for (unsigned int i = 1; i < size - 2; ++i)
		{
			const double m = h[i - 1ULL] / v[i - 1ULL];
			
			vbar[i] = v[i] - m * h[i - 1ULL];
			u[i] -= m * u[i - 1ULL];
		}

		// backward substitution
		z[size - 2] = u[size - 3ULL] / vbar[size - 3ULL];
		for (int i = size - 3; i > 0; --i)
			z[i] = (u[i - 1ULL] - h[i - 1ULL] * z[i + 1ULL]) / vbar[i - 1ULL];
	}


	double Pseudopotential::Interpolate(size_t interval, double x) const
	{
		const double xmt = x - position[interval];
		const double tmx = position[interval + 1ULL] - x;

		return z[interval + 1ULL] / (6. * h[interval]) * xmt * xmt * xmt + z[interval] / (6. * h[interval]) * tmx * tmx * tmx
			+ (pseudopotential[interval + 1ULL] / h[interval] - z[interval + 1ULL] / 6. * h[interval]) * xmt
			+ (pseudopotential[interval] / h[interval] - z[interval] / 6. * h[interval]) * tmx;
	}

	size_t Pseudopotential::GetIndex(double x) const
	{
		auto it = std::upper_bound(position.begin(), position.end(), x);
		if (it != position.begin()) --it;

		return std::distance(position.begin(), it);
	}


	double Pseudopotential::Value(double x) const
	{
		if (x >= maxRadius) return Zion/x;

		return Interpolate(GetIndex(x), x);
	}

	void Pseudopotential::Clear()
	{
		valid = false;
		Zion = 0; 
		Z = 0;
		maxRadius = 0;
		// clean up data, too
		position.resize(0);
		pseudopotential.resize(0);
	}

}