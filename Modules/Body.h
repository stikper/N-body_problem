#pragma once
#ifndef BODY_H
#define BODY_H

#include <string>

#include "MatVec.h"
#include "PhysicalConsts.h"


struct body {
    std::string name;
    double m; //Mass
    std::vector<double> coordinates = {0, 0, 0};
    std::vector<double> velocity = {0, 0, 0};
    std::vector<double> acceleration = {0, 0, 0};
};

// Methods for bodies
std::vector<double> getDistanceVector(const body& body1, const body& body2);

double getDistance(const std::vector<double>& distVec);

std::vector<double> getForce(const std::vector<body>& bodies,const size_t& bodyIndex);

std::vector<double> getAcceleration(const std::vector<body>& bodies,const size_t& bodyIndex);

std::vector<body> getAccelForAll(const std::vector<body>& bodies);

double getGravitationalPotential (const std::vector<body>& bodies,const std::vector<double>& coordinates, const size_t& excludeBodyIndex);

double getKineticEnergy(const body& obj);

double getPotentialEnergy(const std::vector<body>& bodies, const size_t& bodyIndex);

double getTotalEnergy(const std::vector<body>& bodies);


//Methods for changing coordinate system
std::vector<double> coordsToCart(const double& r, const double& theta, const double& phi);

std::vector<double> velocitiesToCart(const double& theta, const double& phi, const double& v_r, const double& v_theta, const double& v_phi);


#endif //BODY_H
