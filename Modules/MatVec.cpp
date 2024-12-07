#include "MatVec.h"


// Methods for vectors and matrix
double getVectorMagnitude(const std::vector<double>& a) {
    double sumOfSquares = 0;
    for (const auto i : a) {
        sumOfSquares += pow(i, 2);
    }
    return sqrt(sumOfSquares);
}

std::vector<double> sumVectorVector(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(std::max(a.size(), b.size()), 0);
    for(size_t i = 0; i < a.size(); i++) {
        result[i] += a[i];
    }
    for(size_t i = 0; i < b.size(); i++) {
        result[i] += b[i];
    }
    return result;
}

std::vector<double> differenceVectorVector(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(std::max(a.size(), b.size()), 0);
    for(size_t i = 0; i < a.size(); i++) {
        result[i] += a[i];
    }
    for(size_t i = 0; i < b.size(); i++) {
        result[i] -= b[i];
    }
    return result;
}

std::vector<double> multiplyNumberVector(const double& n, const std::vector<double>& a) {
    std::vector<double> result(a.size(), 0);
    for (size_t i = 0; i < a.size(); i++) {
        result[i] += a[i] * n;
    }
    return result;
}

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (size_t i = 0; i < std::min(a.size(), b.size()); i++) {
        result += a[i] * b[i];
    }
    return result;
}

std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b) {
    for (size_t i = a.size(); i < 3; i++) {
        a.push_back(0);
    }
    for (size_t i = b.size(); i < 3; i++) {
        b.push_back(0);
    }
    std::vector<double> result(3, 0);
    result[0] = (a[1] * b[2] - a[2] * b[1]);
    result[1] = (a[2] * b[0] - a[0] * b[2]);
    result[2] = (a[0] * b[1] - a[1] * b[0]);
    return result;
}

double getProjection(const std::vector<double>& a, const std::vector<double>& b) {
    return dotProduct(a, b) / getVectorMagnitude(b);
}

std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& A,const std::vector<double>& b) {
    const size_t rows = A.size();
    const size_t cols = A[0].size();
    std::vector<double> result(cols, 0.0);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            result[i] += A[i][j] * b[j];
        }
    }

    return result;
}

std::vector<double> rotateVectorByX(const std::vector<double>& vec,const double& angle = 0) {
    const std::vector<std::vector<double>>
    rot = {{1, 0, 0},
           {0, cos(angle), -sin(angle)},
           {0, sin(angle), cos(angle)}}; // Rotation matrix
    return multiplyMatrixVector(rot, vec);
}

std::vector<double> rotateVectorByY(const std::vector<double>& vec,const double& angle = 0) {
    const std::vector<std::vector<double>>
    rot = {{cos(angle), 0, sin(angle)},
           {0, 1, 0},
           {-sin(angle), 0, cos(angle)}}; // Rotation matrix
    return multiplyMatrixVector(rot, vec);
}

std::vector<double> rotateVectorByZ(const std::vector<double>& vec,const double& angle = 0) {
    const std::vector<std::vector<double>>
    rot = {{cos(angle), -sin(angle), 0},
           {sin(angle), cos(angle), 0},
           {0, 0, 1}}; // Rotation matrix
    return multiplyMatrixVector(rot, vec);
}