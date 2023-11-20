#include "Problem.h"

Problem::Problem(std::ifstream& problemFile) {
    // JSON parsing
    json problem = json::parse(problemFile);
    json problemConfig = problem["config"];
    json problemBodies = problem["bodies"];

    // Configurations
        // Time
    TIME_START = problemConfig["time_start"];
    TIME_END = problemConfig["time_end"];
    TIME_STEP = problemConfig["time_step"];
        // Dynamic time step
    DATA_OUT_TIME_STEP = problemConfig["data_out_time_step"];
    STEP_DECREASE_ACCELERATION_TRIGGER = problemConfig["step_decrease_acceleration_trigger"];
    STEP_DECREASE_FACTOR = problemConfig["step_decrease_factor"];
    MIN_STEP = problemConfig["min_step"];
    STEP_INCREASE_FACTOR = problemConfig["step_increase_factor"];
    STEP_INCREASE_ACCELERATION_TRIGGER = problemConfig["step_increase_acceleration_trigger"];
        // Other
    METHOD = problemConfig["method"];

    // Bodies
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
}