#ifndef COMPUTE_H
#define COMPUTE_H

#include <functional>

#include "Body.h"
#include "DataOut.h"


// Computation methods
double explicitEulerStep(const double& f, double fromValue, double h);

double correctorStep(const double& prediction, const double& f, double fromValue, double h);


// Motion computers
void compBy(const std::function<body(std::vector<body>&, const size_t&, const double&)>& computer, std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);

void compByLeapFrog(std::vector<body>& bodies,  double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);


// Motion step computers
body compByExplicitEulerStep(std::vector<body>& bodies, const size_t& bodyIndex, const double& timeStep);


#endif //COMPUTE_H
