#define GRAVITATIONAL_CONSTANT 6.67430e-11
#define MASS_OF_SUN 1.98847e30
#define MASS_OF_EARTH 5.9722e24
#define ASTRONOMICAL_UNIT 1.49597870700e11
#define EARTH_VELOCITY 2.978e4
#define TIME_STEP 3600
#define TIME_END (365.25 * 24 * 3600)
#define G GRAVITATIONAL_CONSTANT

#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
#include <fstream>
#include <functional>

using namespace std;

struct object {
    double m = 1; //Mass
    vector<double> coordinates = {0, 0, 0};
    vector<double> velocity = {0, 0, 0};
    vector<double> acceleration = {0, 0, 0};
};

void dataWriter(ofstream& out, const vector<double>& data) {
    for (auto i : data) {
        out << i << " ";
    }
    out << endl;
}

void dataOut(const vector<object>& objects, vector<ofstream>& dataFiles) {
    for (size_t i = 0; i < dataFiles.size(); i++) {
        dataWriter(dataFiles[i], objects[i].coordinates);
    }
}

double degToRad(double degrees) {
    return degrees * (M_PI / 180);
}

double radToDeg(double radians) {
    return radians * (180 / M_PI);
}


// Methods for vectors
double getVectorMagnitude(const vector<double>& a) {
    double sumOfSquares = 0;
    for (auto i : a) {
        sumOfSquares += pow(i, 2);
    }
    return sqrt(sumOfSquares);
}

vector<double> sumVectorVector(const vector<double>& a, const vector<double>& b) {
    vector<double> result(max(a.size(), b.size()), 0);
    for(size_t i = 0; i < a.size(); i++) {
        result[i] += a[i];
    }
    for(size_t i = 0; i < b.size(); i++) {
        result[i] += b[i];
    }
    return result;
}

vector<double> differenceVectorVector(const vector<double>& a, const vector<double>& b) {
    vector<double> result(max(a.size(), b.size()), 0);
    for(size_t i = 0; i < a.size(); i++) {
        result[i] += a[i];
    }
    for(size_t i = 0; i < b.size(); i++) {
        result[i] -= b[i];
    }
    return result;
}

vector<double> multiplyNumberVector(const double& n, const vector<double>& a) {
    vector<double> result(a.size(), 0);
    for (size_t i = 0; i < a.size(); i++) {
        result[i] += a[i] * n;
    }
    return result;
}

double dotProduct(const vector<double>& a, const vector<double>& b) {
    double result = 0;
    for (size_t i = 0; i < min(a.size(), b.size()); i++) {
        result = a[i] * b[i];
    }
    return result;
}

vector<double> crossProduct(vector<double> a, vector<double> b) {
    for (size_t i = a.size(); i < 3; i++) {
        a.push_back(0);
    }
    for (size_t i = b.size(); i < 3; i++) {
        b.push_back(0);
    }
    vector<double> result(3, 0);
    result[0] = (a[1] * b[2] - a[2] * b[1]);
    result[1] = (a[2] * b[0] - a[0] * b[2]);
    result[2] = (a[0] * b[1] - a[1] * b[0]);
    return result;
}

double getProjection(const vector<double>& a, const vector<double>& b) {
    return dotProduct(a, b) / getVectorMagnitude(b);
}


// Methods for bodies
vector<double> getDistanceVector(const object& obj1, const object& obj2) {
    vector<double> distances = differenceVectorVector(obj2.coordinates, obj1.coordinates);
    return distances;
}

double getDistance(const vector<double>& distVec) {
    return getVectorMagnitude(distVec);
}

vector<double> getForce(const vector<object>& objs, const size_t objectIndex) {
    vector<double> force(3, 0);
    for (size_t i = 0; i < objs.size(); i++) {
        if (i == objectIndex) {
            continue;
        }
        vector<double> r = getDistanceVector(objs[objectIndex], objs[i]);
        double distance = getDistance(r);
        force = sumVectorVector(force, multiplyNumberVector(G * objs[i].m * objs[objectIndex].m / pow(distance, 3), r)); // Law of Gravity
    }
    return force;
}

vector<double> getAcceleration(const vector<object>& objects, const size_t objectIndex) {
    vector<double> force = getForce(objects, objectIndex);
    vector<double> acceleration = multiplyNumberVector(1 / objects[objectIndex].m, force); // Newton's second law
    return acceleration;
}


// Computation methods
double explicitEulerStep(const double& f, double fromValue, double h) {
    return fromValue + h * f; // y'(x) = f
}

double correctorStep(const double& prediction, const double& f, double fromValue, double h) {
    return fromValue + (1. / 2) * h * (prediction + f);
}


// Motion step computers
void compByExplicitEulerStep(vector<object>& objects, const size_t& objectIndex, const double& timeStep) {
    objects[objectIndex].acceleration = getAcceleration(objects, objectIndex);

    for (int j = 0; j < objects[objectIndex].velocity.size(); j++) {
        objects[objectIndex].velocity[j] = explicitEulerStep(objects[objectIndex].acceleration[j], objects[objectIndex].velocity[j], timeStep);
    }

    for (int j = 0; j < objects[objectIndex].coordinates.size(); j++) {
        objects[objectIndex].coordinates[j] = explicitEulerStep(objects[objectIndex].velocity[j], objects[objectIndex].coordinates[j], timeStep);
    }
}

void compByPredictorCorrectorStep(vector<object>& objects, const size_t& objectIndex, const double& timeStep) {
    object originObject = objects[objectIndex];

    // Prediction calculation
    compByExplicitEulerStep(objects, objectIndex, timeStep);
    originObject.acceleration = objects[objectIndex].acceleration; // Avoiding double calculation

    // Corrector step
    for (int j = 0; j < objects[objectIndex].velocity.size(); j++) {
        originObject.velocity[j] = correctorStep(objects[objectIndex].acceleration[j], originObject.acceleration[j], originObject.velocity[j], timeStep);
    }

    for (int j = 0; j < objects[objectIndex].coordinates.size(); j++) {
        originObject.coordinates[j] = correctorStep(objects[objectIndex].velocity[j], originObject.velocity[j], originObject.coordinates[j], timeStep);
    }
    objects[objectIndex] = originObject;
}


// Motion computer
void compBy(const function<void(vector<object>&, const size_t&, const double&)>& computer, vector<object>& objects, double& t, const double& timeEnd, const double& timeStep, vector<ofstream>& dataFiles) {
    while (t < timeEnd) {
        for (size_t i = 0; i < objects.size(); i++) {
            computer(objects, i, timeStep);
        }

        t += timeStep;

        dataOut(objects, dataFiles);
    }
}


int main() {
    // Object initialising
    // Sun
    object Sun;
    Sun.m = MASS_OF_SUN;
    // Earth
    object Earth;
    Earth.m = MASS_OF_EARTH;
    Earth.coordinates = {ASTRONOMICAL_UNIT, 0, 0};
    Earth.velocity = {0, EARTH_VELOCITY, 0};
    // All objects
    vector<object> objects;
    objects.push_back(Sun);
    objects.push_back(Earth);

    vector<ofstream> dataFiles(objects.size());

    for(size_t i = 0; i < objects.size(); i++) {
        dataFiles[i].open(to_string(i) + ".txt");
    }


    double t = 0;


    system("chcp 65001"); // Fuck Windows


    compBy(compByPredictorCorrectorStep,objects, t, TIME_END, TIME_STEP, dataFiles);


    for(size_t i = 0; i < objects.size(); i++) {
        dataFiles[i].close();
    }

    return 0;
}
