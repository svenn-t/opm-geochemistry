#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <iostream>
#include <vector>

#include "RandomMatrices.hpp"
#include <opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp>

static constexpr double abs_tolerance_ = 1.0e-8;


TEST_CASE("Test identity matrix...")
{
    const int matrixSize = 4;

    double** matrix = allocateMemoryForMatrix(matrixSize);
    fillIdentityMatrix(matrix, matrixSize);

    double determinant = -1;
    double** matrixInverse = allocateMemoryForMatrix(matrixSize);
    invert_matrix(matrix, matrixSize, matrixInverse, determinant);

    CHECK_THAT(determinant, Catch::Matchers::WithinAbs(1.0, 1.0e-10));
    for(int i=0; i < matrixSize; ++i){
        for(int j=0; j < matrixSize; ++j){
            CHECK_THAT(matrix[i][j], Catch::Matchers::WithinAbs(matrixInverse[i][j], 1.0e-10));
        }
    }

    freeMatrixMemory(matrix, matrixSize);
    freeMatrixMemory(matrixInverse, matrixSize);
}

TEST_CASE("Test invert_matrix() and LU_invert_matrix() for an example 3x3 matrix...")
{
    // Used NumPy to randomly generate a 3x3 matrix with elements in the
    // range (0, 1) (rounded to 2 digits)

    double** matrix = allocateMemoryForMatrix(3);
    fillMatrix(matrix, 3, 0.0);

    matrix[0][0] = 0.94;
    matrix[0][1] = 0.95;
    matrix[0][2] = 0.47;
    matrix[1][0] = 0.91;
    matrix[1][1] = 0.55;
    matrix[1][2] = 0.22;
    matrix[2][0] = 0.11;
    matrix[2][1] = 0.55;
    matrix[2][2] = 0.29;

    double determinant = -1;
    double** matrixInverse = allocateMemoryForMatrix(3);

    // Original implementation [1-based index in ludcmp() and lubksb()]
    invert_matrix(matrix, 3, matrixInverse, determinant);
    CHECK_THAT(determinant, Catch::Matchers::WithinAbs(0.015275, abs_tolerance_));

    CHECK_THAT(matrixInverse[0][0], Catch::Matchers::WithinAbs(2.52045827, abs_tolerance_));
    CHECK_THAT(matrixInverse[0][1], Catch::Matchers::WithinAbs(-1.11292962, abs_tolerance_));
    CHECK_THAT(matrixInverse[0][2], Catch::Matchers::WithinAbs(-3.2405892, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][0], Catch::Matchers::WithinAbs(-15.69230769, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][1], Catch::Matchers::WithinAbs(14.46153846, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][2], Catch::Matchers::WithinAbs(14.46153846, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][0], Catch::Matchers::WithinAbs(28.80523732, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][1], Catch::Matchers::WithinAbs(-27.00490998, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][2], Catch::Matchers::WithinAbs(-22.74959083, abs_tolerance_));

    // New implementation [0-based index in LU_decomposition and LU_backsubstitution]
    LU_invert_matrix(matrix, 3, matrixInverse, determinant);
    CHECK_THAT(determinant, Catch::Matchers::WithinAbs(0.015275, abs_tolerance_));

    CHECK_THAT(matrixInverse[0][0], Catch::Matchers::WithinAbs(2.52045827, abs_tolerance_));
    CHECK_THAT(matrixInverse[0][1], Catch::Matchers::WithinAbs(-1.11292962, abs_tolerance_));
    CHECK_THAT(matrixInverse[0][2], Catch::Matchers::WithinAbs(-3.2405892, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][0], Catch::Matchers::WithinAbs(-15.69230769, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][1], Catch::Matchers::WithinAbs(14.46153846, abs_tolerance_));
    CHECK_THAT(matrixInverse[1][2], Catch::Matchers::WithinAbs(14.46153846, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][0], Catch::Matchers::WithinAbs(28.80523732, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][1], Catch::Matchers::WithinAbs(-27.00490998, abs_tolerance_));
    CHECK_THAT(matrixInverse[2][2], Catch::Matchers::WithinAbs(-22.74959083, abs_tolerance_));

    freeMatrixMemory(matrix, 3);
    freeMatrixMemory(matrixInverse, 3);
}

TEST_CASE("Check that independent implementations of LU gives the same matrix output")
{
    auto id = GENERATE(range(1, 12));

    // Note: Below, it doesn't matter that we overwrite the determinant and
    //       some of the other "out parameters", all we care about is the
    //       modified matrix.
    //
    int matrixSize;
    double expectedDeterminant;
    double** matrix_new = getMatrixAndDeterminant(id, matrixSize, expectedDeterminant);
    // Make a "copy" that has an extra row and an extra column:
    double** matrix_old = allocateMemoryForMatrix(matrixSize+1);
    for(int i=0; i < matrixSize; ++i)
    {
        for(int j=0; j < matrixSize; ++j)
        {
            matrix_old[i+1][j+1] = matrix_new[i][j];
        }
    }

    // CALCULATE...
    auto vec_wksp = std::vector<double>(matrixSize+1);
    auto vec_idx  = std::vector<int>(matrixSize+1);
    double even_or_odd = 0.0;

    ludcmp(matrix_old, vec_wksp.data(), matrixSize, vec_idx.data(), even_or_odd);
    LU_decomposition(matrix_new, vec_wksp.data(), matrixSize, vec_idx.data(), even_or_odd);

    // COMPARE...
    for(int i=0; i < matrixSize; ++i)
    {
        for(int j=0; j < matrixSize; ++j)
        {
            const double aij_old = matrix_old[i+1][j+1];
            const double aij_new = matrix_new[i][j];
            CHECK_THAT(aij_new, Catch::Matchers::WithinAbs(aij_old, abs_tolerance_));
        }
    }

    // CALCULATE...
    auto vec_sol = std::vector<double>(matrixSize+1);
    lubksb(matrix_old, matrixSize, vec_idx.data(), vec_sol.data());
    LU_backsubstitution(matrix_new, matrixSize, vec_idx.data(), vec_sol.data());

    // COMPARE...
    for(int i=0; i < matrixSize; ++i)
    {
        for(int j=0; j < matrixSize; ++j)
        {
            const double aij_old = matrix_old[i+1][j+1];
            const double aij_new = matrix_new[i][j];
            CHECK_THAT(aij_new, Catch::Matchers::WithinAbs(aij_old, abs_tolerance_));
        }
    }

    freeMatrixMemory(matrix_old, matrixSize+1);
    freeMatrixMemory(matrix_new, matrixSize);
}
