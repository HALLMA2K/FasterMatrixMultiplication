#ifndef ALMAN_WILLIAMS_MULTIPLICATION
#define ALMAN_WILLIAMS_MULTIPLICATION

void almanWilliamsMultiplication(const int rows, const int pitch,
    const double* leftMatrix, const double* rightMatrix,
    double* resultMatrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            // compute the reordered index for columns of rightMatrix
            int index = (j + pitch * i) % rows;
            double sum = 0;

            // Compute the sum
            for (int k = 0; k < rows; k++) {
                sum += leftMatrix[i * rows + k] * rightMatrix[k * pitch
                    + index];
            }
            // set the value of c_{i,j} in the resultMatrix to the sum
            resultMatrix[i * rows + j] = sum;
        }
    }
}

#endif // !ALMAN_WILLIAMS_MULTIPLICATION
