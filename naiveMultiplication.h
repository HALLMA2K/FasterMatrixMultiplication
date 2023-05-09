#ifndef NAIVEMULTIPLICATION_H
#define NAIVEMULTIPLICATION_H

/*
 * Multiplies two square matrices using a conventional algorithm with a
 * time complexity of O(N^3).
 *
 * @param rows: The number of rows (and columns) in the matrices.
 * @param pitch: The pitch of the matrices (the distance between elements
 *               at (i,j) and (i+1,j), in doubles).
 * @param leftMatrix: A pointer to the first matrix to multiply.
 * @param rightMatrix: A pointer to the second matrix to multiply.
 * @param resultMatrix: A pointer to the matrix to store the result of the
 *                      multiplication in. This matrix must be p
 *                      re-allocated with enough space to hold the result.
 */

void naiveMultiplication(const int rows,
    const int pitch, const double* leftMatrix,
    const double* rightMatrix, double* resultMatrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < pitch; j++) {
            double sum = 0.0;
            for (int k = 0; k < rows; k++) {
                sum += leftMatrix[i * pitch + k] * rightMatrix[k *
                    pitch + j];
            }
            resultMatrix[i * pitch + j] = sum;
        }
    }
}

#endif // NAIVEMULTIPLICATION_H
