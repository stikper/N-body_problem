#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

#include "Modules/Compute.h"
#include "Modules/json.hpp"
using json = nlohmann::json;


using namespace std;


int main() {
    // The problem initialising
    ifstream problemFile("../data.json");
    json problem = json::parse(problemFile);
    json problemConfig = problem["config"];
    json problemBodies = problem["bodies"];

    // Configurations
    double t = problemConfig["time_start"];
    double TIME_END = problemConfig["time_end"];
    double TIME_STEP = problemConfig["time_step"];
    string METHOD = problemConfig["method"];


    // Bodies
    vector<body> bodies;
    bodies.reserve(problemBodies.size());

    for (auto & problemBody : problemBodies) {
        body Body;
        Body.name = problemBody["name"];
        Body.m = problemBody["mass"];
        for (size_t i = 0; i < problemBody["coordinates"].size(); i++) {
            Body.coordinates[i] = problemBody["coordinates"][i];
        }
        for (size_t i = 0; i < problemBody["velocity"].size(); i++) {
            Body.velocity[i] = problemBody["velocity"][i];
        }

        bodies.push_back(Body);
    }

    bodies = getAccelForAll(bodies);


    vector<ofstream> dataFiles(bodies.size() + 1);

    for(size_t i = 0; i < bodies.size(); i++) {
        dataFiles[i].open(bodies[i].name + METHOD + ".txt");
    }
    dataFiles[bodies.size()].open("total_energy.txt");


    system("chcp 65001"); // Fuck Windows

    dataOut(bodies, dataFiles);
    dataWriter(dataFiles[bodies.size()], vector<double> {t, getTotalEnergy(bodies)});


    if (METHOD == "Eul") {
        comp(ByEuler, bodies, t, TIME_END, TIME_STEP, dataFiles);
    } else if (METHOD == "PC") {
        comp(ByPredictorCorrector, bodies, t, TIME_END, TIME_STEP, dataFiles);
    } else if (METHOD == "LF") {
        compByLeapFrog(bodies, t, TIME_END, TIME_STEP, dataFiles);
    } else if (METHOD == "RK4") {
        comp(ByRK4, bodies, t, TIME_END, TIME_STEP, dataFiles);
    }


    for (size_t i = 0; i < bodies.size(); i++) {
        dataFiles[i].close();
    }

    return 0;
}
