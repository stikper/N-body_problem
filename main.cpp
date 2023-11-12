#define EARTH_VELOCITY 2.978e4
#define TIME_STEP (3600 * 24 * 5)
#define TIME_END (30 * 365.25 * 24 * 3600)

#define METHOD ByPredictorCorrector
#define FILE_NAME "PC"


#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

#include "Modules/Compute.h"


using namespace std;


int main() {
    // Object initialising
        // Sun
    body Sun;
    Sun.name = "Sun";
    Sun.m = MASS_OF_SUN;
        // Earth
    body Earth;
    Earth.name = "Earth";
    Earth.m = MASS_OF_EARTH;
    Earth.coordinates = {ASTRONOMICAL_UNIT, 0, 0};
    Earth.velocity = {0, EARTH_VELOCITY, 0};
        //Body
//    body Body;
//    Body.name = "Body";
//    Body.m = 0.000955 * MASS_OF_SUN;
//    Body.coordinates = coordsToCart(19.988409, 1.941617, 4.341507);
//    Body.velocity = velocitiesToCart(1.941617, 4.341507, -1.720970, 0.702991, 0.702991);
//    Body.coordinates = multiplyNumberVector(ASTRONOMICAL_UNIT, Body.coordinates);
//    Body.velocity = multiplyNumberVector(ASTRONOMICAL_UNIT / (365.25 * 24 * 3600), Body.velocity);
    // All bodies
    vector<body> bodies;
    bodies.push_back(Sun);
    bodies.push_back(Earth);
    bodies = getAccelForAll(bodies);

    vector<ofstream> dataFiles(bodies.size() + 1);

    for(size_t i = 0; i < bodies.size(); i++) {
        dataFiles[i].open(bodies[i].name + FILE_NAME + ".txt");
    }
    dataFiles[bodies.size()].open("total_energy.txt");

    double t = 0;

    system("chcp 65001"); // Fuck Windows

    dataOut(bodies, dataFiles);
    dataWriter(dataFiles[bodies.size()], vector<double> {t, getTotalEnergy(bodies)});
//    compByLeapFrog(bodies, t, TIME_END, TIME_STEP, dataFiles);
     comp(METHOD, bodies, t, TIME_END, TIME_STEP, dataFiles);

    for (size_t i = 0; i < bodies.size(); i++) {
        dataFiles[i].close();
    }

    return 0;
}
