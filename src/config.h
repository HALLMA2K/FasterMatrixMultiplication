#ifndef CONFIG_H
#define CONFIG_H

// Strassen settings
constexpr bool kRunStrassen = true;
constexpr bool kDisplayTotalStrassenTime = true;

// Naive settings
constexpr bool kRunNaive = true;
constexpr bool kDisplayTotalNaiveTime = true;
// Winograd settings
constexpr bool kRunWinograd = true;
constexpr bool kDisplayTotalWinogradTime = true;

// Hybrid settings
constexpr bool kRunHybrid = true;
constexpr int kHybridNaiveThreshold = 16;
constexpr int kHybridWinogradThreshold = 2000;
constexpr bool kDisplayTotalHybridTime = true;

// Alma Williams settings
constexpr bool kRunAlmanWilliams = true;
constexpr bool kDisplayTotalAlmanWilliamsTime = true;

// Results and time display settings
constexpr int kTimePrecision = 7;
constexpr int kResultPrecision = 2;
constexpr bool kDisplayResults = true;
constexpr bool kDisplayPerCalculationTime = true;

// Display Matrices to multiply
constexpr bool kDisplayMatrix = true;

// Matrix dimensions
constexpr int kMatrixRows = 2;
constexpr int kMatrixCols = 2;

// Number of times to run the calculation
constexpr int kCalculationCount = 1;

// Random number generator bounds
constexpr int kRNGLowerBound = 1;
constexpr int kRNGUpperBound = 10;

#endif // !CONFIG_H
