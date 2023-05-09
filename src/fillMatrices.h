#ifndef FILLMATRICES_H
#define FILLMATRICES_H

#include <random>

void randomNumberGenerator(int matrixSize, int matrixPitch,
    double matrix[], int seed) {
    // Lower bound for random number generation
    const double lowerBound = kRNGLowerBound;  
    // Upper bound for random number generation
    const double upperBound = kRNGUpperBound;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(lowerBound, upperBound);

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrix[i * matrixPitch + j] = dist(rng);
        }
    }
}

void generateMatrices(int martrixSize, double*& leftMatrix,
    double*& rightMatrix) {
    // Generate matrices leftMatrix and  rightMatrix with a different seed
    // each time
    std::random_device rd;
    for (int i = 0; i < 3; i++) {
        randomNumberGenerator(martrixSize, martrixSize, i == 0 ?
            leftMatrix : rightMatrix, rd());
    }
}

#endif // FILLMATRICES_H
