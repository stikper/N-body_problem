#ifndef MATVEC_H
#define MATVEC_H

#include <cmath>
#include <vector>

// Methods for vectors and matrix
double getVectorMagnitude(const std::vector<double>& a);

std::vector<double> sumVectorVector(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> differenceVectorVector(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> multiplyNumberVector(const double& n, const std::vector<double>& a);

double dotProduct(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b);

double getProjection(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& A,const std::vector<double>& b);

std::vector<double> rotateVectorByX(const std::vector<double>& vec,const double& angle);

std::vector<double> rotateVectorByY(const std::vector<double>& vec,const double& angle);

std::vector<double> rotateVectorByZ(const std::vector<double>& vec,const double& angle);


#endif //MATVEC_H
