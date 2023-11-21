#include "Compute.h"

// Computation methods
double eulerStep(const double& f, const double& yi, const double& h) {
    return yi + h * f; // y'(x) = f
}

std::vector<double> eulerVectorStep(const std::vector<double>& f, const std::vector<double>& yi, const double& h) {
    std::vector<double> result(yi.size());
    for (size_t i = 0; i < yi.size(); i++) {
        result[i] = eulerStep(f[i], yi[i], h);
    }
    return result;
}

double correctorStep(const double& f1, const double& f2, const double& yi, const double& h) {
    return yi + (h / 2) * (f1 + f2);
}

std::vector<double> correctorVectorStep(const std::vector<double>& f1, const std::vector<double>& f2, const std::vector<double>& yi, const double& h) {
    std::vector<double> result(yi.size());
    for (size_t i = 0; i < yi.size(); i++) {
        result[i] = correctorStep(f1[i], f2[i], yi[i], h);
    }
    return result;
}

double RK4Step(const double& f1, const double& f2, const double& f3, const double& f4, const double& yi, const double& h) {
    return yi + (h / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
}

std::vector<double> RK4VectorStep(const std::vector<double>& f1, const std::vector<double>& f2, const std::vector<double>& f3, const std::vector<double>& f4, const std::vector<double>& yi, const double& h) {
    std::vector<double> result(yi.size());
    for (size_t i = 0; i < yi.size(); i++) {
        result[i] = RK4Step(f1[i], f2[i], f3[i], f4[i], yi[i], h);
    }
    return result;
}


// Motion computers
std::vector<body> comp(const std::function<std::vector<body>(std::vector<body>&, const double&)>& computer, const std::vector<body>& bodies, double& t, const double& timeEnd, const Problem& problem, dataOut& DataOut) {
    std::vector<body> result = bodies;
    std::vector<body> newResult;
    newResult.reserve(result.size());

    double timeStep = problem.TIME_STEP;
    while (t < timeEnd) {
        newResult = computer(result, timeStep);

        switch (checkStep(result, newResult, problem)) {
            case -1:
                if (timeStep == problem.MIN_STEP) {
                    t += timeStep;
                    result = newResult;
                    break;
                }
                decreaseStep(timeStep, problem);
                continue;
            case 1:
                t += timeStep;
                increaseStep(timeStep, problem);
                result = newResult;
                break;
            case 0:
                t += timeStep;
                result = newResult;
                break;
            case default:
                break;
        }

        DataOut.Out(result, t);
    }
    return result;
}



std::vector<body> compByLF(const std::vector<body>& bodies, double& t, const double& timeEnd, const Problem& problem, dataOut& DataOut) {
    std::vector<body> result = bodies;
    std::vector<body> newLFBodies;
    newLFBodies.reserve(bodies.size());

    double timeStep = problem.TIME_STEP;
    if (t < timeEnd) {
        std::vector<body> LFBodies = ToLF(bodies, timeStep); // First step (Get V(i+h/2))
        while (t < timeEnd) {
            newLFBodies = ByLF(LFBodies, timeStep);

            switch (checkStep(LFBodies, newLFBodies, problem)) {
                case -1:
                    if (timeStep == problem.MIN_STEP) {
                        t += timeStep;
                        LFBodies = newLFBodies;
                        result = FromLF(LFBodies, timeStep); // Get V(i)
                        break;
                    }
                    decreaseStep(timeStep, problem);
                    continue;
                case 1:
                    LFBodies = newLFBodies;
                    result = FromLF(LFBodies, timeStep); // Get V(i)
                    t += timeStep;
                    increaseStep(timeStep, problem);
                    break;
                case 0:
                    LFBodies = newLFBodies;
                    result = FromLF(LFBodies, timeStep); // Get V(i)
                    t += timeStep;
                    break;
                case default:
                    break;
            }

            DataOut.Out(result, t);
        }
    }
    return result;
}


// Motion step computers
    // Explicit Euler
body ByEulerBody(body Body, const double& timeStep) {
    Body.velocity = eulerVectorStep(Body.acceleration, Body.velocity, timeStep);
    Body.coordinates = eulerVectorStep(Body.velocity, Body.coordinates, timeStep);

    return Body;
}

std::vector<body> ByEuler(const std::vector<body>& bodies, const double& timeStep) {
    std::vector<body> result;
    result.reserve(bodies.size());

    for (const body& body : bodies) {
        result.push_back(ByEulerBody(body, timeStep));
    }

    result = getAccelForAll(result);
    return result;
}

    // Predictor-corrector
std::vector<body> ByPredictorCorrector(const std::vector<body>& bodies, const double& timeStep) {
    std::vector<body> result;
    result.reserve(bodies.size());

    const std::vector<body>& k1 = bodies;
    const std::vector<body>& k2 = EulerByDataOf(k1, k1, timeStep);

    for (size_t i = 0; i < bodies.size(); i++) {
        body Body = bodies[i];
        const body& k1b = k1[i];
        const body& k2b = k2[i];

        Body.velocity = correctorVectorStep(k1b.acceleration, k2b.acceleration, k1b.velocity, timeStep);
        Body.coordinates = correctorVectorStep(k1b.velocity, k2b.velocity, k1b.coordinates, timeStep);

        result.push_back(Body);
    }

    result = getAccelForAll(result);
    return result;
}

    // Runge-Kutta 4-th order
body EulerBodyByDataOf(const body& Data, body Body, const double& timeStep) {
    Body.velocity = eulerVectorStep(Data.acceleration, Body.velocity, timeStep);
    Body.coordinates = eulerVectorStep(Data.velocity, Body.coordinates, timeStep);

    return Body;
}

std::vector<body> EulerByDataOf(const std::vector<body>& data, const std::vector<body>& bodies, const double& timeStep) {
    std::vector<body> result;
    result.reserve(bodies.size());

    for (size_t i = 0; i < bodies.size(); i++) {
        result.push_back(EulerBodyByDataOf(data[i], bodies[i], timeStep));
    }

    result = getAccelForAll(result);
    return result;
}

std::vector<body> ByRK4(const std::vector<body>& bodies, const double& timeStep) {
    std::vector<body> result;
    result.reserve(bodies.size());

    const std::vector<body>& k1 = bodies;
    const std::vector<body>& k2 = EulerByDataOf(k1, k1, timeStep / 2);
    const std::vector<body>& k3 = EulerByDataOf(k2, k1, timeStep / 2);
    const std::vector<body>& k4 = EulerByDataOf(k3, k1, timeStep);

    for (size_t i = 0; i < bodies.size(); i++) {
        body Body = bodies[i];
        const body& k1b = k1[i];
        const body& k2b = k2[i];
        const body& k3b = k3[i];
        const body& k4b = k4[i];

        Body.velocity = RK4VectorStep(k1b.acceleration, k2b.acceleration, k3b.acceleration, k4b.acceleration, k1b.velocity, timeStep);
        Body.coordinates = RK4VectorStep(k1b.velocity, k2b.velocity, k3b.velocity, k4b.velocity, k1b.coordinates, timeStep);

        result.push_back(Body);
    }

    result = getAccelForAll(result);
    return result;
}

std::vector<body> ToLF(const std::vector<body>& bodies, const double& timeStep) {
    std::vector<body> LFBodies = bodies;

    for (auto & Body : LFBodies) {
        Body.velocity = eulerVectorStep(Body.acceleration, Body.velocity, timeStep / 2);
    }

    return LFBodies;
}

std::vector<body> ByLF(const std::vector<body>& LFBodies, const double& timeStep) {
    std::vector<body> result = LFBodies;

    for (auto & Body : result) {
        Body.coordinates = eulerVectorStep(Body.velocity, Body.coordinates, timeStep);
    }

    result = getAccelForAll(result);

    for (auto & Body : result) {
        Body.velocity = eulerVectorStep(Body.acceleration, Body.velocity, timeStep);
    }

    return result;
}

std::vector<body> FromLF(const std::vector<body>& LFBodies, const double& timeStep) {
    std::vector<body> result = LFBodies;

    for (auto & Body : result) {
        Body.velocity = eulerVectorStep(Body.acceleration, Body.velocity, -(timeStep / 2));
    }

    return result;
}

// Check if step decreasing / increasing is necessary
int checkStep(const std::vector<body>& bodies, const std::vector<body>& result, const Problem& problem) {
    bool incTrig = false;
    bool doNotInc = false;
    bool decTrig = false;
    for (size_t i = 0; i < bodies.size(); i++) {
        double deltaVel = getVectorMagnitude(differenceVectorVector(result[i].velocity, bodies[i].velocity));
        if (deltaVel >= problem.STEP_DECREASE_ACCELERATION_TRIGGER) {
            decTrig = true;
        } else if (deltaVel <= problem.STEP_INCREASE_ACCELERATION_TRIGGER) {
            incTrig = true;
        } else {
            doNotInc = true;
        }
    }
    if(doNotInc) incTrig = false;
    if(decTrig) incTrig = false;

    if (incTrig) return 1;
    if (decTrig) return -1;
    return 0;
}

void increaseStep(double& timeStep, const Problem& problem) {
    double newTimeStep = timeStep * problem.STEP_INCREASE_FACTOR;
    if (newTimeStep > problem.TIME_STEP) {
        newTimeStep = problem.TIME_STEP;
    }
    timeStep = newTimeStep;
}

void decreaseStep(double& timeStep, const Problem& problem) {
    double newTimeStep = timeStep * problem.STEP_DECREASE_FACTOR;
    if (newTimeStep < problem.MIN_STEP) {
        newTimeStep = problem.MIN_STEP;
    }
    timeStep = newTimeStep;
}