#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>
#include <windows.h>

# include "alman-williamMultiplication.h"
#include "config.h"
#include "fillMatrices.h"
#include "hybridMultiplication.h"
#include "naiveMultiplication.h"
#include "strassen.h"
#include "timer.h"
#include "winograd.h"

struct CalculationTime {
    bool displayFlag;
    double time;
    std::string name;
};

void displayPerCalculationTime(float startTime, float endTime,
    std::vector<float>& times) {
    if (kDisplayPerCalculationTime) {
        times.push_back((endTime - startTime) / 1000.0);
        for (int i = 0; i < times.size(); ++i) {
            std::cout << "per calculation time: " <<
                std::fixed << std::setprecision(kTimePrecision)
                << times[i] << "\n\n";
        }
        times.clear();
    }
}

void displayResult(const int rows,const int pitch, const double matrix[],
    std::string name) {
    if (name == "LeftMatrix" || name == "RightMatrix")
        std::cout << name << "(" << rows << "x" << pitch << "): " << std::endl;
    else { 
        std::cout << name << " multiplication result: " << std::endl;
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < pitch; j++) {
            std::cout << std::fixed << std::setprecision(kResultPrecision)
                << matrix[i * pitch + j] << " ";;
        }
        std::cout << std::endl;
    }
}

/*
 * Program that tests the performance of the
 * various algorithm for matrix multiplication.
 *
 * Robert Hallmark  HALLMA2K@Outlook.com
 *
 * Conventional algorithm time complexity of O(N ^ 3).
 * Winograd algorithm time complexity of  O(N^2.373)
 * Strassen algorithm time complexity of O(N^2.807).
 *
 * @param leftMatrix: A pointer to the first matrix to multiply.
 * @param rightMatrix: A pointer to the second matrix to multiply.
 * @param resultMatrix: A pointer to the matrix to store the result of the
 * multiplication in. This matrix must be pre-allocated with enough space
 * to hold the result.
 */
int main() {

    time_t seconds;
    time(&seconds);
    srand((unsigned int)seconds);
    Timer timer;
    std::vector<float> times;

    double naiveTotalTime = 0.0, winogradTotalTime = 0.0,
        strassenTotalTime = 0.0, hybridTotalTime = 0.0,
        almanWilliamsTotalTime = 0.0;
    
    double* leftMatrix = new double[kMatrixRows * kMatrixCols];
    double* rightMatrix = new double[kMatrixRows * kMatrixCols];
    double* resultMatrix = new double[kMatrixRows * kMatrixCols];
    
    std::cout << "Starting Calculations: This may take some time "
        "depending on matrices sizes." << std::endl;
    std::cout << "-------------------------------------------------------"
        "--------------------" << std::endl << std::endl;
    
    try {
        for (int i = 0; i < kCalculationCount; i++) {

            generateMatrices(kMatrixRows, leftMatrix, rightMatrix);
        
            if (kDisplayMatrix) {
                std::string name = "LeftMatrix";
                displayResult(kMatrixRows, kMatrixCols, leftMatrix, name);
                std::cout << std::endl;
                name = "RightMatrix";
                displayResult(kMatrixRows, kMatrixCols, rightMatrix, name);
                std::cout << std::endl;
            }
        
            if (kMatrixRows != kMatrixCols) {
                throw std::invalid_argument("The matrices to be "
                    "multiplied must be square, i.e. the number of rows and "
                    "columns in the left and right matrices must be equal. "
                    "Please update the dimensions of your matrices in the "
                    "config.cpp file and ensure that the 'kMatrixRows' and "
                    "'kMatrixCols' macros have the same value.");
            }
            else
            {
                if (kRunNaive) {
                    std::string name = "Naive";
                    double naiveStartTime = timer.milliseconds();
                    naiveMultiplication(kMatrixRows, kMatrixCols,
                        leftMatrix,
                        rightMatrix, resultMatrix);
                    double naiveEndTime = timer.milliseconds();
                    displayResult(kMatrixRows, kMatrixCols, resultMatrix,
                        name);
                    displayPerCalculationTime(naiveStartTime, naiveEndTime,
                        times);
                    naiveTotalTime += naiveEndTime - naiveStartTime;
                }
                if (kRunWinograd) {
                    std::string name = "Winograd";
                    double winogradStartTime = timer.milliseconds();
                    winogradMultiplication(kMatrixRows, kMatrixCols,
                        leftMatrix, rightMatrix, resultMatrix);
                    double winogradEndTime = timer.milliseconds();
                    displayResult(kMatrixRows, kMatrixCols, resultMatrix,
                        name);
                    displayPerCalculationTime(winogradStartTime,
                        winogradEndTime, times);
                    winogradTotalTime += (winogradEndTime -
                        winogradStartTime);
                }
                if (kRunStrassen) {
                    std::string name = "Strassen";
                    double strassenStartTime = timer.milliseconds();
                    strassenMultiplication(kMatrixRows, kMatrixCols, leftMatrix,
                        kMatrixCols, rightMatrix, kMatrixCols, resultMatrix);
                    double strassenEndTime = timer.milliseconds();
                    displayResult(kMatrixRows, kMatrixCols, resultMatrix, name);
                    displayPerCalculationTime(strassenStartTime, strassenEndTime, times);
                    strassenTotalTime += (strassenEndTime - strassenStartTime);
                }
                if (kRunHybrid) {
                    std::string name = "Hybrid";
                    double hybridStartTime = timer.milliseconds();
                    hybridMultiplication(kMatrixRows, kMatrixCols, leftMatrix,
                        rightMatrix, resultMatrix);
                    double hybridEndTime = timer.milliseconds();
                    displayResult(kMatrixRows, kMatrixCols, resultMatrix, name);
                    displayPerCalculationTime(hybridStartTime, hybridEndTime, times);
                    hybridTotalTime += (hybridEndTime - hybridStartTime);
                }
                if (kRunAlmanWilliams)
                {
                    std::string name = "Alman-Williams";
                    double almanWilliamsStartTime = timer.milliseconds();
                    almanWilliamsMultiplication(kMatrixRows, kMatrixCols, leftMatrix,
                        rightMatrix, resultMatrix);
                    double almanWilliamsEndTime = timer.milliseconds();
                    displayResult(kMatrixRows, kMatrixCols, resultMatrix, name);
                    displayPerCalculationTime(almanWilliamsStartTime, almanWilliamsEndTime, times);
                    almanWilliamsTotalTime += (almanWilliamsEndTime - almanWilliamsStartTime);
                }
            }
        }

        std::vector<CalculationTime> calculationTimes = {
            { kDisplayTotalNaiveTime, naiveTotalTime / 1000.0F, "Naive" },
            { kDisplayTotalWinogradTime, winogradTotalTime / 1000.0F, "Winograd" },
            { kDisplayTotalStrassenTime, strassenTotalTime / 1000.0F, "Strassen" },
            { kDisplayTotalHybridTime, hybridTotalTime / 1000.0F, "Hybrid" },
            { kDisplayTotalAlmanWilliamsTime, almanWilliamsTotalTime / 1000.0F, "Alman-Williams" }
        };

        auto displayFunc = [](const CalculationTime& ct) {
            return ct.displayFlag;
        };
        auto displayTimes = std::vector<CalculationTime>();
        std::copy_if(calculationTimes.begin(), calculationTimes.end(),
            std::back_inserter(displayTimes), displayFunc);

        for (const auto& calcTime : displayTimes) {
            std::cout << "total " << calcTime.name << " time: " << std::fixed
                << std::setprecision(kTimePrecision) << calcTime.time <<
                std::endl;
        }

        delete[] leftMatrix;
        delete[] rightMatrix;
        delete[] resultMatrix;

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred while generating matrices: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception occurred." << std::endl;
        return 1;
    }
}
