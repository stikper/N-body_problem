#include "DataOut.h"


// Data output methods
void dataWriter(std::ofstream& out, const std::vector<double>& data) {
    for (auto i : data) {
        out << i << " ";
    }
    out << std::endl;
}

void dataOut(const std::vector<body>& bodies, std::vector<std::ofstream>& dataFiles) {
    for (size_t i = 0; i < bodies.size(); i++) {
        dataWriter(dataFiles[i], bodies[i].coordinates);
    }
}
