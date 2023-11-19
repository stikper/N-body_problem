#pragma once
#ifndef COMPUTE_H
#define COMPUTE_H

#include <functional>

#include "Body.h"
#include "DataOut.h"


// Computation methods
double eulerStep(const double& f, const double& yi, const double& h);

std::vector<double> eulerVectorStep(const std::vector<double>& f, const std::vector<double>& yi, const double& h);

double correctorStep(const double& f1, const double& f2, const double& yi, const double& h);

std::vector<double> correctorVectorStep(const std::vector<double>& f1, const std::vector<double>& f2, const std::vector<double>& yi, const double& h);

double RK4Step(const double& f1, const double& f2, const double& f3, const double& f4, const double& yi, const double& h);

std::vector<double> RK4VectorStep(const std::vector<double>& f1, const std::vector<double>& f2, const std::vector<double>& f3, const std::vector<double>& f4, const std::vector<double>& yi, const double& h);

// Motion computers
std::vector<body> comp(const std::function<std::vector<body>(std::vector<body>&, const double&)>& computer, const std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);

std::vector<body> compByLF(const std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles);


// Motion step computers
    // Explicit Euler
body ByEulerBody(body Body, const double& timeStep);

std::vector<body> ByEuler(const std::vector<body>& bodies, const double& timeStep);
    // Predictor-corrector
std::vector<body> ByPredictorCorrector(const std::vector<body>& bodies, const double& timeStep);
    // Runge-Kutta 4-th order
body EulerBodyByDataOf(const body& Data, body Body, const double& timeStep);

std::vector<body> EulerByDataOf(const std::vector<body>& data, const std::vector<body>& bodies, const double& timeStep);

std::vector<body> ByRK4(const std::vector<body>& bodies, const double& timeStep);

std::vector<body> ToLF(const std::vector<body>& bodies, const double& timeStep);

std::vector<body> ByLF(const std::vector<body>& bodies, const double& timeStep);

std::vector<body> FromLF(const std::vector<body>& bodies, const double& timeStep);

#endif //COMPUTE_H
