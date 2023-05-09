#ifndef STRASSEN_MULTIPLICATION_H
#define STRASSEN_MULTIPLICATION_H

// resultMatrix = leftMatrix * rightMatrix
void matrixMultiplication(const int matrixSize, const int leftMatrixPitch,
    const double leftMatrix[],
    const int rightMatrixPitch, const double rightMatrix[], 
    const int resultMatrixPitch, double resultMatrix[]) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            double sum = 0.0;
            for (int k = 0; k < matrixSize; k++) {
                sum += leftMatrix[i * leftMatrixPitch + k] * rightMatrix[k
                    * rightMatrixPitch + j];
                resultMatrix[i * resultMatrixPitch + j] = sum;
            }
        }
    }
}

// resultMatrix = leftMatrix + rightMatrix
void matrixAddition(const int matrixSize, const int leftMatrixPitch,
    const double leftMatrix[], const int rightMatrixPitch,
    const double rightMatrix[], const int resultMatrixPitch,
    double resultMatrix[]) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            resultMatrix[i * resultMatrixPitch + j] = leftMatrix[i *
                leftMatrixPitch + j] + rightMatrix[i * rightMatrixPitch +
                j];
        }
    }
}


//  resultMatrix = leftMatrix - rightMatrix
void matrixSubtraction(const int matrixSize,const int leftMatrixPitch,
    const double leftMatrix[], const int rightMatrixPitch,
    const double rightMatrix[], const int resultMatrixPitch,
    double resultMatrix[]) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            resultMatrix[i * resultMatrixPitch + j] = leftMatrix[i *
                leftMatrixPitch + j] - rightMatrix[i * rightMatrixPitch +
                j];
        }
    }
}

void printStrassen(const int rows,const int pitch, const double matrix[]) {
    std::cout << "Strassen Multiplication Result: " << std::endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < pitch; j++) {
            std::cout << std::fixed << std::setprecision(kResultPrecision)
                << matrix[i * pitch + j] << " ";;
        }
        std::cout << std::endl;
    }
}

/*
 * Multiplies two matrices using Strassen's algorithm, which recursively
 * splits the matrices into smaller submatrices and uses seven matrix
 * multiplications to compute the result. This algorithm has a theoretical
 * runtime of O(N^2.807).
 *
 * This implementation assumes that both matrices have dimensions NxN,
 * where N is a power of two. The algorithm splits each matrix into four
 * smaller (N/2)x(N/2) submatrices as follows:
 *
 *                   _    _                    _   _
 *     leftMatrix = | A  B |    rightMatrix = | E F |
 *                  | C  D |                  | G H |
 *                   -    -                    -   -
 *
 * The following seven matrices are then computed using seven (N/2)x(N/2)
 * matrix multiplications:
 *
 *     P0 = A*(F - H);
 *     P1 = (A + B)*H
 *     P2 = (C + D)*E
 *     P3 = D*(G - E);
 *     P4 = (A + D)*(E + H)
 *     P5 = (B - D)*(G + H)
 *     P6 = (A - C)*(E + F)
 *
 * Finally, the result is computed as follows:
 *
 *                   _                                            _
 *   resultMatrix = | (P3 + P4) + (P5 - P1)   P0 + P1              |
 *                  | P2 + P3                 (P0 + P4) - (P2 + P6)|
 *                   -                                            -
 *
 * @param matrixSize: The dimension of the matrices.
 * @param leftMatrixPitch: The pitch of the left matrix (the distance
 *                         between elements at (i,j) and (i+1,j), in
 *                         doubles).
 * @param leftMatrix: A pointer to the first matrix to multiply.
 * @param rightMatrixPitch: The pitch of the right matrix (the distance
 *                         between elements at (i,j) and (i+1,j), in
 *                         doubles).
 * @param rightMatrix: A pointer to the second matrix to multiply.
 * @param resultMatrixPitch: The pitch of the result matrix (the distance
 *                         between elements at (i,j) and (i+1,j), in
 *                         doubles).
 * @param resultMatrix: A pointer to the matrix to store the result of the
 *                      multiplication in. This matrix must be 
 *                      pre-allocated with enough space to hold the
 *                      result.
 */

void strassenMultiplication(const int matrixSize, const int leftMatrixPitch,
    const double leftMatrix[], const int rightMatrixPitch,
    const double rightMatrix[], const int resultMatrixPitch,
    double resultMatrix[]) {
   
    if (matrixSize == 1) {
        matrixMultiplication(matrixSize, leftMatrixPitch, leftMatrix,
            rightMatrixPitch, rightMatrix, resultMatrixPitch,
            resultMatrix);
        return;
    }

    const int n = matrixSize / 2;    // size of sub-matrices

    const double* A = leftMatrix;    // A-D matrices embedded in leftMatrix
    const double* B = leftMatrix + n;
    const double* C = leftMatrix + n * leftMatrixPitch;
    const double* D = C + n;

    const double* E = rightMatrix;   // E-H matrices embeded in rightMatrix
    const double* F = rightMatrix + n;
    const double* G = rightMatrix + n * rightMatrixPitch;
    const double* H = G + n;

    double* P[7];                    // allocate temp matrices off heap
    const int sz = n * n * sizeof(double);
    for (int i = 0; i < 7; i++) {
        P[i] = (double*)malloc(sz);
    }
    double* T = (double*)malloc(sz);
    double* U = (double*)malloc(sz);

    // P0 = A * (F - H);
    matrixSubtraction(n, rightMatrixPitch, F, rightMatrixPitch, H, n, T);
    strassenMultiplication(n, leftMatrixPitch, A, n, T, n, P[0]);

    // P1 = (A + B) * H
    matrixAddition(n, leftMatrixPitch, A, leftMatrixPitch, B, n, T);
    strassenMultiplication(n, n, T, rightMatrixPitch, H, n, P[1]);

    // P2 = (C + D) * E
    matrixAddition(n, leftMatrixPitch, C, leftMatrixPitch, D, n, T);
    strassenMultiplication(n, n, T, rightMatrixPitch, E, n, P[2]);

    // P3 = D * (G - E);
    matrixSubtraction(n, rightMatrixPitch, G, rightMatrixPitch, E, n, T);
    strassenMultiplication(n, leftMatrixPitch, D, n, T, n, P[3]);

    // P4 = (A + D) * (E + H)
    matrixAddition(n, leftMatrixPitch, A, leftMatrixPitch, D, n, T);
    matrixAddition(n, rightMatrixPitch, E, rightMatrixPitch, H, n, U);
    strassenMultiplication(n, n, T, n, U, n, P[4]);

    // P5 = (B - D) * (G + H)
    matrixSubtraction(n, leftMatrixPitch, B, leftMatrixPitch, D, n, T);
    matrixAddition(n, rightMatrixPitch, G, rightMatrixPitch, H, n, U);
    strassenMultiplication(n, n, T, n, U, n, P[5]);

    // P6 = (A - C) * (E + F)
    matrixSubtraction(n, leftMatrixPitch, A, leftMatrixPitch, C, n, T);
    matrixAddition(n, rightMatrixPitch, E, rightMatrixPitch, F, n, U);
    strassenMultiplication(n, n, T, n, U, n, P[6]);

    // resultMatrix upper left = (P3 + P4) + (P5 - P1)
    matrixAddition(n, n, P[4], n, P[3], n, T);
    matrixSubtraction(n, n, P[5], n, P[1], n, U);
    matrixAddition(n, n, T, n, U, resultMatrixPitch, resultMatrix);

    // resultMatrix lower left = P2 + P3
    matrixAddition(n, n, P[2], n, P[3], resultMatrixPitch, resultMatrix +
        n * resultMatrixPitch);

    // resultMatrix upper right = P0 + P1
    matrixAddition(n, n, P[0], n, P[1], resultMatrixPitch, resultMatrix +
        n);

    // resultMatrix lower right = (P0 + P4) - (P2 + P6)
    matrixAddition(n, n, P[0], n, P[4], n, T);
    matrixAddition(n, n, P[2], n, P[6], n, U);
    matrixSubtraction(n, n, T, n, U, resultMatrixPitch, resultMatrix + n *
        (resultMatrixPitch + 1));

    free(U);  // deallocate temp matrices
    free(T);
    for (int i = 6; i >= 0; i--) {
        free(P[i]);
    }
}

#endif // STRASSEN_MULTIPLICATION_H
