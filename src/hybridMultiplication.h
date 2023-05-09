#ifndef HYBRID_MULTIPLICATION
#define HYBRID_MULTIPLICATION

void hybridMatrixAddition(const int rows,
    const int pitch, const double leftMatrix[],
    const double rightMatrix[], double resultMatrix[]) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < pitch; j++)
            resultMatrix[i * pitch + j] = leftMatrix[i * pitch + j] +
            rightMatrix[i * pitch + j];
}

void hybridMatrixSubtraction(const int rows,
    const int pitch, const double leftMatrix[],
    const double rightMatrix[], double resultMatrix[]) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < pitch; j++)
            resultMatrix[i * pitch + j] = leftMatrix[i * pitch + j] -
            rightMatrix[i * pitch + j];
}

/*
 * Performs matrix multiplication using a hybrid of three different
 * algorithms: naive, Winograd, and Strassen. The algorithm to use 
 * depends on the size of the input matrices.
 *
 * @param rows Number of rows in the matrices.
 * @param pitch Pitch of the matrices.
 * @param leftMatrix Pointer to the left matrix.
 * @param rightMatrix Pointer to the right matrix.
 * @param resultMatrix Pointer to the result matrix.
 */

void hybridMultiplication(const int rows, const int pitch,
    const double leftMatrix[], const double rightMatrix[],
    double resultMatrix[]) {
    /*
     * Recursive base case. If matrices are kHybridNaiveThreshold or
     * smaller we just use the conventional algorithm. At what size you
     * should switch will vary based on your hardware platform.
     */
    if (rows <= kHybridNaiveThreshold) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < pitch; j++) {
                double sum = 0.0;
                for (int k = 0; k < pitch; k++) {
                    sum += leftMatrix[i * pitch + k] * rightMatrix[k *
                        pitch + j];
                }
                resultMatrix[i * pitch + j] = sum;
            }
        }
        return;
    }
    /*
     * Recursive base case. If matrices are kHybridWinogradThreshold or
     * smaller we just use the Winograd algorithm. At what size you should
     * switch will vary based on your hardware platform.
     */
    else if (rows <= kHybridWinogradThreshold)
    {
        for (int i = 0; i < rows; i += 2) {
            for (int j = 0; j < rows; j += 2) {
                double c00 = 0, c01 = 0, c10 = 0, c11 = 0;
                for (int t = 0; t < rows; ++t) {
                    c00 += leftMatrix[i * rows + t] * rightMatrix[t *
                        pitch + j];
                    c01 += leftMatrix[i * rows + t] * rightMatrix[t *
                        pitch + j + 1];
                    c10 += leftMatrix[(i + 1) * rows + t] *
                        rightMatrix[t * pitch + j];
                    c11 += leftMatrix[(i + 1) * rows + t] *
                        rightMatrix[t * pitch + j + 1];
                }
                resultMatrix[i * pitch + j] = c00;
                resultMatrix[i * pitch + j + 1] = c01;
                resultMatrix[(i + 1) * pitch + j] = c10;
                resultMatrix[(i + 1) * pitch + j + 1] = c11;
            }
        }
    }
    /*
     * Recursive base case. If matrices are larger than kHybridWinogradThreshold
     * we just use the Strassen algorithm. At what size you should
     * switch will vary based on your hardware platform.
     */
    else
    {
        const int n = rows / 2;       // size of sub-matrices

        const double* A = leftMatrix;
        const double* B = leftMatrix + n;
        const double* C = leftMatrix + n * pitch;
        const double* D = C + n;

        const double* E = rightMatrix;
        const double* F = rightMatrix + n;
        const double* G = rightMatrix + n * pitch;
        const double* H = G + n;

        double* P[7];
        const int sz = n * n * sizeof(double);
        for (int i = 0; i < 7; i++) {
            P[i] = (double*)malloc(sz);
        }
        double* T = (double*)malloc(sz);
        double* U = (double*)malloc(sz);
        

        // P0
        hybridMatrixSubtraction(n, pitch, F, H, T);
        hybridMultiplication(n, pitch, A, T, P[0]);

        // P1
        hybridMatrixAddition(n, pitch, A, B, T);
        hybridMultiplication(n, n, T, H, P[1]);

        // P2
        hybridMatrixAddition(n, pitch, C, D, T);
        hybridMultiplication(n, n, T, E, P[2]);

        // P3
        hybridMatrixSubtraction(n, pitch, G, E, T);
        hybridMultiplication(n, pitch, D, T, P[3]);

        // P4
        hybridMatrixAddition(n, pitch, A, D, T);
        hybridMatrixAddition(n, pitch, E, H, U);
        hybridMultiplication(n, n, T, U, P[4]);

        // P5
        hybridMatrixSubtraction(n, pitch, B, D, T);
        hybridMatrixAddition(n, pitch, G, H, U);
        hybridMultiplication(n, n, T, U, P[5]);

        // P6
        hybridMatrixSubtraction(n, pitch, A, C, T);
        hybridMatrixAddition(n, pitch, E, F, U);
        hybridMultiplication(n, n, T, U, P[6]);

        // resultMatrix upper left
        hybridMatrixAddition(n, n, P[4], P[3], T);
        hybridMatrixSubtraction(n, n, P[5], P[1], U);
        hybridMatrixAddition(n, n, T, U, resultMatrix);

        // resultMatrix lower left
        hybridMatrixAddition(n, n, P[2], P[3], resultMatrix + n *
            pitch);

        // resultMatrix upper right
        hybridMatrixAddition(n, n, P[0], P[1], resultMatrix + n);

        // resultMatrix lower right
        hybridMatrixAddition(n, n, P[0], P[4], T);
        hybridMatrixAddition(n, n, P[2], P[6], U);
        hybridMatrixSubtraction(n, n, T, U, resultMatrix + n * (pitch +
            1));

        free(U);  // deallocate temp matrices
        free(T);
        for (int i = 6; i >= 0; i--) {
            free(P[i]);
        }
    }
}

#endif // !HYBRID_MULTIPLICATION
