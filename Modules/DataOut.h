#ifndef DATAOUT_H
#define DATAOUT_H

#include <vector>
#include <fstream>

#include "Body.h"


// Data output methods
void dataWriter(std::ofstream& out, const std::vector<double>& data);

void dataOut(const std::vector<body>& bodies, std::vector<std::ofstream>& dataFiles);

#endif //DATAOUT_H
