#ifndef COMPUTE_H
#define COMPUTE_H

#include <functional>

#include "Body.h"
#include "DataOut.h"


// Computation methods
double eulerStep(const double& f, const double& yi, const double& h);

std::vector<double> eulerVectorStep(const std::vector<double>& f, const std::vector<double>& yi, const double& h);

// Motion computers
std::vector<body> comp(const std::function<std::vector<body>(std::vector<body>&, const double&)>& computer, const std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);

void compByLeapFrog(std::vector<body>& bodies,  double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);


// Motion step computers
    // Explicit Euler
body ByEulerBodyStep(body Body, const double& timeStep);

std::vector<body> ByEulerStep(const std::vector<body>& bodies, const double& timeStep);

#endif //COMPUTE_H
