#include <iostream>
#include <vector>
#include "Modules/MatVec.cpp"
#include "Modules/Consts.h"

using namespace std;


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


int main() {
    double mass = 0.000955;
    vector<double> coords = coordsToCart(19.988409, 1.941617, 4.341507);
    vector<double> velocities = velocitiesToCart(1.941617, 4.341507, -1.720970, 0.702991, 0.702991);

    mass *= MASS_OF_SUN;
    coords = multiplyNumberVector(AU, coords);
    velocities = multiplyNumberVector(ASTRONOMICAL_UNIT / YEAR, velocities);

    cout << "\"mass\": " << mass << ",\n";
    cout << "\"coordinates\": [" << coords[0] << ", " << coords[1] <<  ", " << coords[2] << "],\n";
    cout << "\"velocity\": [" << velocities[0] << ", " << velocities[1] << ", " << velocities[2] << "]" << endl;
    return 0;
}