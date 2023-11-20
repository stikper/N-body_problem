#pragma once
#ifndef N_BODY_PROBLEM_PROBLEM_H
#define N_BODY_PROBLEM_PROBLEM_H

#include <vector>
#include <fstream>
#include <string>

#include "Body.h"
#include "json.hpp"
using json = nlohmann::json;


class Problem {
private:
public:
    explicit Problem(std::ifstream& problemFile);

    // Configurations
        // Time
    double TIME_START;
    double TIME_END;
    double TIME_STEP;
        // Dynamic time step
    double DATA_OUT_TIME_STEP;
    double STEP_DECREASE_ACCELERATION_TRIGGER;
    double STEP_DECREASE_FACTOR;
    double MIN_STEP;
    double STEP_INCREASE_FACTOR;
    double STEP_INCREASE_ACCELERATION_TRIGGER;
        // Other
    std::string METHOD;

    // Bodies
    std::vector<body> bodies;
};


#endif //N_BODY_PROBLEM_PROBLEM_H
