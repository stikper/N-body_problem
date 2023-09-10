#define GRAVITATIONAL_CONSTANT 6.67430e-11
#define MASS_OF_SUN 1.98847e30
#define MASS_OF_EARTH 5.9722e24
#define ASTRONOMICAL_UNIT 1.49597870700e11
#define EARTH_VELOCITY 0.978e4
#define TIME_STEP 3600
#define TIME_END (320 * 24 * 3600)
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
    vector<double> velocities = {0, 0, 0};
    vector<double> accelerations = {0, 0, 0};
};

double degToRad(double degrees) {
    return degrees * (M_PI / 180);
}

double radToDeg(double radians) {
    return radians * (180 / M_PI);
}

vector<double> sumVectorVector(vector<double> x, vector<double> y) {
    size_t cols = x.size();
    vector<double> res(cols, 0);
    for (int i = 0; i < cols; i++) {
        res[i] = x[i] + y[i];
    }
    return res;
}

double getDistanceBetweenObjects(const object& obj1, const object& obj2) {
    return sqrt(pow(obj2.coordinates[0] - obj1.coordinates[0], 2) + pow(obj2.coordinates[1] - obj1.coordinates[1], 2)  + pow(obj2.coordinates[2] - obj1.coordinates[2], 2));
}

vector<double> getDistanceByCoordinates(const object& obj1, const object& obj2) {
    vector<double> distances= {obj2.coordinates[0] - obj1.coordinates[0], obj2.coordinates[1] - obj1.coordinates[1], obj2.coordinates[2] - obj1.coordinates[2]};
    return distances;
}

vector<double> getAccelerationsByForces(const object& obj, vector<double> forces) {
    vector<double> accelerations(3);
    for (int i = 0; i < 3; i++) {
        accelerations[i] = forces[i] / obj.m; // Newton's second law
    }
    return accelerations;
}

vector<double> getForcesForObject(const vector<object>& objs, const size_t objectIndex) {
    vector<double> forces = {0, 0, 0};
    for (size_t i = 0; i < objs.size(); i++) {
        if (i == objectIndex) {
            continue;
        }
        double distance = getDistanceBetweenObjects(objs[objectIndex], objs[i]);
        double force = (G * objs[i].m * objs[objectIndex].m) / pow(distance, 2);
        vector<double> distances = getDistanceByCoordinates(objs[objectIndex], objs[i]);
        forces[0] += force * (distances[0] / distance);
        forces[1] += force * (distances[1] / distance);
        forces[2] += force * (distances[2] / distance);
    }
    return forces;
}

double explicitEulerForwardStep(const function <double(vector<double>)>& f, vector<double> args, double fromValue, double h) {
    return fromValue + h * f(std::move(args)); // y'(x) = f(x,y)
}

void dataOut(ofstream& out, const vector<double>& data) {
    for (auto i : data) {
        out << i << " ";
    }
    out << endl;
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
    Earth.velocities = {0, EARTH_VELOCITY, 0};
    // All objects
    vector<object> objects;
    objects.push_back(Sun);
    objects.push_back(Earth);

    ofstream EarthFile;
    EarthFile.open("Earth.txt");

    ofstream Forces;
    Forces.open("Forces.txt");

    double t = 0;

    do {
        for (size_t i = 0; i < objects.size(); i++) {
            objects[i].accelerations = getAccelerationsByForces(objects[i], getForcesForObject(objects, i));

            for(int j = 0; j < objects[i].velocities.size(); j++) {
                objects[i].velocities[j] = explicitEulerForwardStep([](vector<double> args) {
                    return args[0];
                }, vector<double> {objects[i].accelerations[j]}, objects[i].velocities[j], TIME_STEP);
            }

            for(int j = 0; j < objects[i].coordinates.size(); j++) {
                objects[i].coordinates[j] = explicitEulerForwardStep([](vector<double> args) {
                    return args[0];
                }, vector<double> {objects[i].velocities[j]}, objects[i].coordinates[j], TIME_STEP);
            }

            dataOut(EarthFile, objects[1].coordinates);
            dataOut(Forces, getForcesForObject(objects, 1));
        }
        t += TIME_STEP;
    } while (t < TIME_END);

    EarthFile.close();
    Forces.close();

    system("chcp 65001"); // Fuck Windows
    return 0;
}
