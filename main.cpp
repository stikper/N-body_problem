#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

#include "Modules/Compute.h"


using namespace std;


int main() {
    // The problem initialising
    ifstream problemFile("../Problems/ecOrbit.json");
    Problem problem(problemFile);

    vector<body> bodies = problem.bodies;
    double t = problem.TIME_START;
    double timeStep = problem.TIME_STEP;

    system("chcp 65001"); // Fuck Windows

    dataOut DataOut(bodies, problem.DATA_OUT_TIME_STEP, problem.METHOD);
    DataOut.forceOut(bodies, t);


    if (problem.METHOD == "Eul") {
        comp(ByEuler, bodies, t, problem.TIME_END, timeStep, DataOut);
    } else if (problem.METHOD == "PC") {
        comp(ByPredictorCorrector, bodies, t, problem.TIME_END, timeStep, DataOut);
    } else if (problem.METHOD == "LF") {
        compByLF(bodies, t, problem.TIME_END, timeStep, DataOut);
    } else if (problem.METHOD == "RK4") {
        comp(ByRK4, bodies, t, problem.TIME_END, timeStep, DataOut);
    }


    DataOut.close();
    return 0;
}
