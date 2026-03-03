/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef LINEAR_SOLVER_ROUTINES_DEF_H
#define LINEAR_SOLVER_ROUTINES_DEF_H

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

// ****************************************************************************************
//					                    HELPER FUNCTIONS
// ****************************************************************************************
double** allocateMemoryForMatrix(int matrixSize);
void freeMatrixMemory(double** matrix, int matrixSize);
void fillMatrix(double** matrix, int matrixSize, double value);
void fillIdentityMatrix(double** matrix, int matrixSize);

void printSquareMatrix(double** matrix, int matrixSize, std::ostream& out=std::cout);

// ****************************************************************************************
//					LINEAR SOLVER ROUTINES (LU-DECOMPOSITION)
// ****************************************************************************************

/**
* Solves a set of n linear equations: A*x = b.
* The routine takes into account the possibility that b will begin with many zero elements,
* hence it is efficient for use in matrix inversion.
*
* IMPORTANT: The input matrix is NOT the full matrix!!
*
* @param[in] A LU decomposition of matrix (can be computed with the separate routine ludcmp).
* @param[in] n Dimension of matrix.
* @param[in] index Permutation of original matrix computed by LU-decomposition, e.g., by ludcmp().
* @param[in,out] b Originally the right-hand side vector, at the end of the function holds the solution.
*/
void lubksb(double** A, int n, const int* index, double* b);

/**
* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition
* of a rowwise permutation of itself.
*
* This routine is used in combination with lubksb to solve linear equations, or invert a matrix.
*
* @param[in,out] A Matrix which is modified [see equation (2.3.14)]
* @param[out] wksp A 1D workspace array of at least the same size as the matrix (can also be larger).
* @param[in] n Dimension of matrix.
* @param[in,out] index Is modified to store the row permutation resulting from the partial pivoting.
* @param[in,out] sign Depending on whether the number of row interchanges is even or odd, set to 1 (even) or -1 (odd).
*
*/
int ludcmp(double** A, double* wksp, int n, int* index, double& sign);

/**
* Inverts a matrix, and calculates the determinant.
* The inverse is calculated by means of a LU decomposition.
*
* NB: Should not be used to calculate determinant of a large matrix!
*
* @param[in] A The matrix to be inverted (NOT modified)
* @param[in] n The dimension of the matrix.
* @param[in,out] a_inv The inverse of the input matrix.
* @param[in,out] determinant The determinant of the input matrix.
 */
void invert_matrix(double** A, int n, double** A_inv, double& determinant);

// 10/12-2022:
//	Attempt to create alternate impl. of the above functions that use 0-based indexing.
//
void LU_backsubstitution(double** A, int n, const int* index, double* b);
int LU_decomposition(double** A, double* wksp, int n, int* index, double& sign);
void LU_invert_matrix(double** A, int n, double** A_inv, double& determinant);

#endif