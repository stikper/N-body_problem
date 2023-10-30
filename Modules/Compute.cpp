#include "Compute.h"

// Computation methods
double explicitEulerStep(const double& f, double fromValue, double h) {
    return fromValue + h * f; // y'(x) = f
}

double correctorStep(const double& prediction, const double& f, double fromValue, double h) {
    return fromValue + (1. / 2) * h * (prediction + f);
}


// Motion computers
void compBy(const std::function<body(std::vector<body>&, const size_t&, const double&)>& computer, std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles) {
    while (t < timeEnd) {
        std::vector<body> newBodies;
        for (size_t i = 0; i < bodies.size(); i++) {
            newBodies.push_back(computer(bodies, i, timeStep));
        }

        t += timeStep;

        for (size_t i = 0; i < bodies.size(); i++) {
            bodies[i] = newBodies[i];
        }

        dataOut(bodies, dataFiles);
        dataWriter(dataFiles[bodies.size()], std::vector<double> {t, getTotalEnergy(bodies)});
    }
}

void compByLeapFrog(std::vector<body>& bodies,  double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles) {
    if (t < timeEnd) {
        for (size_t i = 0; i < bodies.size(); i++) {
            bodies[i].acceleration = getAcceleration(bodies, i);
            for (size_t j = 0; j < bodies[i].velocity.size(); j++) {
                bodies[i].velocity[j] = explicitEulerStep(bodies[i].acceleration[j], bodies[i].velocity[j], timeStep / 2);
            }
        } // Get V(i+1/2) (First step)

        while (t < timeEnd) {
            for (size_t i = 0; i < bodies.size(); i++) {
                for (size_t j = 0; j < bodies[i].coordinates.size(); j++) {
                    bodies[i].coordinates[j] = explicitEulerStep(bodies[i].velocity[j], bodies[i].coordinates[j], timeStep);
                } // Get X(i+1)
                bodies[i].acceleration = getAcceleration(bodies, i);
                for (size_t j = 0; j < bodies[i].velocity.size(); j++) {
                    bodies[i].velocity[j] = explicitEulerStep(bodies[i].acceleration[j], bodies[i].velocity[j], timeStep);
                } // Get V(i+1/2)
            }

            t += timeStep;

            dataOut(bodies, dataFiles);
            dataWriter(dataFiles[bodies.size()], std::vector<double> {t, getTotalEnergy(bodies)});
        }

        for (auto & object : bodies) {
            for (size_t j = 0; j < object.velocity.size(); j++) {
                object.velocity[j] = explicitEulerStep(-object.acceleration[j], object.velocity[j], timeStep / 2);
            }
        } // Get V(i)
    }
}


// Motion step computers
body compByExplicitEulerStep(std::vector<body>& bodies, const size_t& bodyIndex, const double& timeStep) {
    bodies[bodyIndex].acceleration = getAcceleration(bodies, bodyIndex);

    body Body = bodies[bodyIndex];

    for (int j = 0; j < Body.velocity.size(); j++) {
        Body.velocity[j] = explicitEulerStep(Body.acceleration[j], Body.velocity[j], timeStep);
    }

    for (int j = 0; j < Body.coordinates.size(); j++) {
        Body.coordinates[j] = explicitEulerStep(Body.velocity[j], Body.coordinates[j], timeStep);
    }

    return Body;
}

