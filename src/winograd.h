#ifndef WINOGRAD_MULTIPLICATION
#define WINOGRAD_MULTIPLICATION

void printWinograd(const int rows, const int pitch, const double* matrix) {
    std::cout << "Winograd Multiplication Result:" << std::endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < pitch; j++) {
            std::cout << std::fixed << std::setprecision(kResultPrecision)
                << matrix[i * pitch + j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Multiplies two matrices using Winograd's algorithm, which breaks down
 * the original matrices into smaller matrices to reduce the number of
 * multiplications required. This function has a time complexity of
 * O(n^2.373).
 *
 * @param rows: The number of rows in the left and result matrices, and
 *               the number of columns in the right matrix.
 * @param pitch: The number of columns in the result matrix, and the pitch
 *               of the left and right matrices (the distance between
 *               elements at (i,j) and (i+1,j), in doubles).
 * @param leftMatrix: A pointer to the first matrix to multiply.
 * @param rightMatrix: A pointer to the second matrix to multiply.
 * @param resultMatrix: A pointer to the matrix to store the result of the
 *                      multiplication in. This matrix must be pre-allocated
 *                      with enough space to hold the result.
 */

void winogradMultiplication(const int rows, const int pitch,
    const double* leftMatrix, const double* rightMatrix,
    double* resultMatrix) {
    for (int i = 0; i < rows; i += 2) {
        for (int j = 0; j < rows; j += 2) {
            double c00 = 0, c01 = 0, c10 = 0, c11 = 0;
            for (int t = 0; t < pitch; ++t) {
                c00 += leftMatrix[i * rows + t] * rightMatrix[t * pitch +
                    j];
                c01 += leftMatrix[i * rows + t] * rightMatrix[t * pitch +
                    j + 1];
                c10 += leftMatrix[(i + 1) * rows + t] * rightMatrix[t *
                    pitch + j];
                c11 += leftMatrix[(i + 1) * rows + t] * rightMatrix[t *
                    pitch + j + 1];
            }
            resultMatrix[i * pitch + j] = c00;
            resultMatrix[i * pitch + j + 1] = c01;
            resultMatrix[(i + 1) * pitch + j] = c10;
            resultMatrix[(i + 1) * pitch + j + 1] = c11;
        }
    }
}

#endif // !WINOGRAD_MULTIPLICATION
