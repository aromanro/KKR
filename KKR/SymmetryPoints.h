#pragma once

#include "Vector3D.h"

#include <string>
#include <vector>

#include <map>

namespace KKR
{

	class SymmetryPoint
	{
	public:
		SymmetryPoint() {}

		SymmetryPoint(const std::string& Name, const Vector3D<double> pos) : name(Name), position(pos) {}

		SymmetryPoint(const SymmetryPoint& sym)
			: name(sym.name), position(sym.position)
		{
		}

		SymmetryPoint& operator=(const SymmetryPoint& sym)
		{
			name = sym.name;
			position = sym.position;

			return *this;
		}

		std::string name;
		Vector3D<double> position;
	};


	class SymmetryPoints
	{
	public:
		SymmetryPoints();

		std::map<std::string, SymmetryPoint> symmetryPoints;


		std::vector<Vector3D<double>> GeneratePoints(const std::vector<std::string>& path, unsigned int nrPoints, std::vector<unsigned int>& symmetryPointsPositions);
	};

}