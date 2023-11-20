#include "DataOut.h"


dataOut::dataOut(const std::vector<body>& bodies, const double& step, const std::string& method) {
    _n = bodies.size();
    for(size_t i = 0; i < _n; i++) {
        _dataFiles.resize(_n + 1);
        _dataFiles[i].open(bodies[i].name + "_" + method + ".txt");
        _step = step;
    }
}

void dataOut::Out(const std::vector<body>& bodies, const double& t) {
    size_t n = bodies.size();
    for (size_t i = 0; i < n; i++) {
        dataOut::Writer(_dataFiles[i], bodies[i].coordinates);
    }
    dataOut::Writer(_dataFiles[n], std::vector<double>{t, getTotalEnergy(bodies)});
}

void dataOut::Writer(std::ofstream& out, const std::vector<double>& data) {
    for (auto i : data) {
        out << i << " ";
    }
    out << std::endl;
}

void dataOut::close() {
    for (auto & _dataFile : _dataFiles) {
        _dataFile.close();
    }
}