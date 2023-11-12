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
std::vector<body> comp(const std::function<std::vector<body>(std::vector<body>&, const double&)>& computer, const std::vector<body>& bodies, double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles) {
    std::vector<body> result = bodies;
    while (t < timeEnd) {
        result = computer(result, timeStep);

        t += timeStep;

        dataOut(result, dataFiles);
        dataWriter(dataFiles[result.size()], std::vector<double> {t, getTotalEnergy(result)});
    }
    return result;
}

void compByLeapFrog(std::vector<body>& bodies,  double& t, const double& timeEnd, const double& timeStep, std::vector<std::ofstream>& dataFiles) {
    if (t < timeEnd) {
        for (size_t i = 0; i < bodies.size(); i++) {
            bodies[i].acceleration = getAcceleration(bodies, i);
            for (size_t j = 0; j < bodies[i].velocity.size(); j++) {
                bodies[i].velocity[j] = eulerStep(bodies[i].acceleration[j], bodies[i].velocity[j], timeStep / 2);
            }
        } // Get V(i+1/2) (First step)

        while (t < timeEnd) {
            for (size_t i = 0; i < bodies.size(); i++) {
                for (size_t j = 0; j < bodies[i].coordinates.size(); j++) {
                    bodies[i].coordinates[j] = eulerStep(bodies[i].velocity[j], bodies[i].coordinates[j], timeStep);
                } // Get X(i+1)
                bodies[i].acceleration = getAcceleration(bodies, i);
                for (size_t j = 0; j < bodies[i].velocity.size(); j++) {
                    bodies[i].velocity[j] = eulerStep(bodies[i].acceleration[j], bodies[i].velocity[j], timeStep);
                } // Get V(i+1/2)
            }

            t += timeStep;

            dataOut(bodies, dataFiles);
            dataWriter(dataFiles[bodies.size()], std::vector<double> {t, getTotalEnergy(bodies)});
        }

        for (auto & object : bodies) {
            for (size_t j = 0; j < object.velocity.size(); j++) {
                object.velocity[j] = eulerStep(-object.acceleration[j], object.velocity[j], timeStep / 2);
            }
        } // Get V(i)
    }
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

    std::vector<body> k1 = bodies;
    std::vector<body> k2 = EulerByDataOf(k1, k1, timeStep);

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

    std::vector<body> k1;
    std::vector<body> k2;
    std::vector<body> k3;
    std::vector<body> k4;

    k1 = bodies;
    k2 = EulerByDataOf(k1, k1, timeStep / 2);
    k3 = EulerByDataOf(k2, k1, timeStep / 2);
    k4 = EulerByDataOf(k3, k1, timeStep);

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
