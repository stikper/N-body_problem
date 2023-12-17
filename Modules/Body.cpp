#include "Body.h"


// Methods for bodies
std::vector<double> getDistanceVector(const body& body1, const body& body2) {
    std::vector<double> distances = differenceVectorVector(body2.coordinates, body1.coordinates);
    return distances;
}

double getDistance(const std::vector<double>& distVec) {
    return getVectorMagnitude(distVec);
}

std::vector<double> getForce(const std::vector<body>& bodies, const size_t& bodyIndex) {
    std::vector<double> force(3, 0);
    for (size_t i = 0; i < bodies.size(); i++) {
        if (i == bodyIndex) {
            continue;
        }
        std::vector<double> r = getDistanceVector(bodies[bodyIndex], bodies[i]);
        const double distance = getDistance(r);
        force = sumVectorVector(force, multiplyNumberVector(G * bodies[i].m * bodies[bodyIndex].m / pow(distance, 3), r)); // Law of Gravity
    }
    return force;
}

std::vector<double> getAcceleration(const std::vector<body>& bodies, const size_t& bodyIndex) {
    const std::vector<double> force = getForce(bodies, bodyIndex);
    std::vector<double> acceleration = multiplyNumberVector(1 / bodies[bodyIndex].m, force); // Newton's second law
    return acceleration;
}

std::vector<body> getAccelForAll(const std::vector<body>& bodies) {
    std::vector<body> result;
    for (size_t i = 0; i < bodies.size(); i++) {
        result.push_back(bodies[i]);
        result[i].acceleration = getAcceleration(bodies, i);
    }
    return result;
}

double getGravitationalPotential (const std::vector<body>& bodies,const std::vector<double>& coordinates, const size_t& excludeBodyIndex = -1) {
    double potential = 0;
    body point;
    point.coordinates = coordinates;
    for (size_t i = 0; i < bodies.size(); i++) {
        if (i == excludeBodyIndex) {
            continue;
        }
        potential += -(G * bodies[i].m) / getVectorMagnitude(getDistanceVector(bodies[i], point));
    }
    return potential;
}

double getKineticEnergy(const body& obj) {
    return (obj.m * pow(getVectorMagnitude(obj.velocity), 2)) / 2;
}

double getPotentialEnergy(const std::vector<body>& bodies, const size_t& bodyIndex) {
    return bodies[bodyIndex].m * getGravitationalPotential(bodies, bodies[bodyIndex].coordinates, bodyIndex);
}

double getTotalEnergy(const std::vector<body>& bodies) {
    double energy = 0;
    for(size_t i = 0; i < bodies.size(); i++) {
        energy += getKineticEnergy(bodies[i]); // Kinetic energy
        energy += (1. / 2) * getPotentialEnergy(bodies, i); // Potential Energy
    }

    return energy;
}


//Methods for changing coordinate system
std::vector<double> coordsToCart(const double& r = 0, const double& theta = 0, const double& phi = 0) {
    std::vector<double> coords(3);
    coords[0] = r * sin(theta) * cos(phi);
    coords[1] = r * sin(theta) * sin(phi);
    coords[2] = r * cos(theta);
    return coords;
}

std::vector<double> velocitiesToCart(const double& theta = 0, const double& phi = 0, const double& v_r = 0, const double& v_theta = 0, const double& v_phi = 0) {
    std::vector<double> velocities(3);
    // Axes: r -> z, phi -> y, theta -> x
    velocities[0] = v_theta; // Axis theta
    velocities[1] = v_phi; // Axis phi
    velocities[2] = v_r; // Axis r
    // Rotate to x, y, z axes
    velocities = rotateVectorByY(velocities, theta); // r to z
    velocities = rotateVectorByZ(velocities, phi); // theta to x (phi to y)
    return velocities;
}
