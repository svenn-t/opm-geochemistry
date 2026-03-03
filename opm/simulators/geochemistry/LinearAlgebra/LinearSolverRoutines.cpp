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
#include <opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp>

// ****************************************************************************************
//					                    HELPER FUNCTIONS
// ****************************************************************************************
double** allocateMemoryForMatrix(int matrixSize)
{
    double** matrix = new double*[matrixSize];
    for(int i=0; i < matrixSize; ++i)
    {
        matrix[i] = new double[matrixSize];
    }
    return matrix;
}

void freeMatrixMemory(double** matrix, int matrixSize)
{
    for(int i=0; i < matrixSize; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void fillMatrix(double** matrix, int matrixSize, double value)
{
    for (int i = 0; i < matrixSize; ++i)
    {
        for (int j = 0; j < matrixSize; ++j)
        {
            matrix[i][j] = value;
        }
    }
}

void fillIdentityMatrix(double** matrix, int matrixSize)
{
    for (int i = 0; i < matrixSize; ++i)
    {
        for (int j = 0; j < matrixSize; ++j)
        {
            if (i == j) matrix[i][j] = 1.0;
            else        matrix[i][j] = 0.0;
        }
    }
}

void printSquareMatrix(double** matrix, int matrixSize, std::ostream& out)
{
    for(int i = 0; i < matrixSize; ++i)
    {
        out << matrix[i][0];
        for(int j = 1; j < matrixSize; ++j)
        {
            out << "\t" << matrix[i][j];
        }
        out << "\n";
    }
}

// ****************************************************************************************
//					LINEAR SOLVER ROUTINES (LU-DECOMPOSITION)
// ****************************************************************************************

void lubksb(double** A, int n, const int* index, double* b)
{
    // When ii is set to a positive value, it will become the index
    // of the first nonvanishing element of b.
    int ii = 0;

    // We now do the forward substitution.
    // The only new wrinkle is to unscramble the permutation as we go.
    double sum;
    for(int i=1; i <= n; ++i)
    {
        const int ip    = index[i];
        sum             = b[ip];
        b[ip]           = b[i];
        if(ii > 0)
        {
            for(int j=ii; j < i; ++j)
            {
                sum -= A[i][j]*b[j];
            }
        }
        else if(sum != 0)
        {
            // A nonzero element was encountered, so from now on we will
            // have to do the sums in the loop above...
            ii = i;
        }
        b[i] = sum;
    }

    // Finally, do the backsubstitution.
    for (int i=n; i >= 1; --i)
    {
        sum = b[i];
        for (int j=i+1; j <= n; ++j)
        {
            sum -= A[i][j]*b[j];
        }
        b[i] = sum / A[i][i];  // Store a component of the solution vector.
    }
}

int ludcmp(double** A, double* wksp, int n, int* index, double& sign)
{
    int imax = -1;

    sign = 1.0;  // No row interchanges yet.
    for(int i=1; i <= n; ++i)
    {   // Loop over rows to get the implicit scaling information.
        double big = 0.0;
        for(int j=1; j <= n; ++j)
        {
            const double temp = std::fabs(A[i][j]);
            if(temp > big)
            {
                big = temp;
            }
        }

        if(big == 0.0)
        {
            printf("Singular matrix in routine ludcmp\n");
            return 1;
        }
        wksp[i] = 1.0 / big;  // Store the scaling of each row.
    }

    for(int j=1; j <= n; ++j)
    {
        for(int i=1; i < j; ++i)
        {
            double sum = A[i][j];
            for(int k=1; k < i; ++k)
            {
                sum -= A[i][k]*A[k][j];
            }
            A[i][j] = sum;
        }

        double big = 0.0;  // Initialize for the search for largest pivot element.
        for(int i=j; i <= n; ++i)
        {   // This is i = j of equation (2.3.12) and i = j+1. . .N.
            // of equation (2.3.13).
            double sum = A[i][j];
            for (int k=1; k < j; ++k)
            {
                sum -= A[i][k]*A[k][j];
            }

            A[i][j] = sum;

            const double dum = wksp[i]*std::fabs(sum);
            if(dum >= big)
            {
                // Is the figure of merit for the pivot better than the best so far?
                big  = dum;
                imax = i;
            }
        }
        assert(imax > -1);

        if(j != imax)
        {   // Do we need to interchange rows?
            for(int k=1; k <= n; ++k)
            {   // Yes, do so...
                const double dum    = A[imax][k];
                A[imax][k]          = A[j][k];
                A[j][k]             = dum;
            }
            sign        = -sign;    // ...and change the parity...
            wksp[imax]  = wksp[j];  // ...and also interchange the scale factor.
        }
        index[j] = imax;

        // If the pivot element is zero, the matrix is singular (at least to the precision of the
        // algorithm). For some applications, it is desirable to substitute TINY for zero...
        if(A[j][j] == 0.0)
        {
            A[j][j] = 1.0e-20;
        }

        if(j != n)
        {   // Now, nally, divide by the pivot element.
            const double dum = 1.0 / A[j][j];
            for (int i=j+1; i <= n; ++i)
            {
                A[i][j] *= dum;
            }
        }
    }  // Go back for the next column in the reduction.

    return 0;
}

void invert_matrix(double** A, int n, double** A_inv, double& determinant)
{
    double** A_tmp = allocateMemoryForMatrix(n+1);
    double* wksp   = new double[n+1];

    // Make a copy of the original matrix
    for(int i=0; i < n; ++i)
    {
        for(int j=0; j < n; ++j)
        {
            A_tmp[i+1][j+1] = A[i][j];
        }
    }

    int* index = new int[n + 1];
    ludcmp(A_tmp, wksp, n, index, determinant); // So far, either 1 or -1.

    for(int j=1; j<= n; ++j)
    {
        determinant *= A_tmp[j][j];
    }

    if(std::fabs(determinant) < 1.0e-20)
    {
        printf("Singular matrix, no inverse was found.\n");
        std::exit(1);
    }

    double* col = new double[n + 1];
    for(int j=1; j <= n; ++j)
    {
        for(int i=1; i <= n; ++i)
        {
            col[i] = 0.0;
        }
        col[j] = 1.0;

        lubksb(A_tmp, n, index, col);

        for(int i=1; i <= n; ++i)
        {
            A_inv[i-1][j-1] = col[i];
        }
    }

    freeMatrixMemory(A_tmp, n+1);
    delete[] wksp;
    delete[] index;
    delete[] col;
}

void LU_backsubstitution(double** A, int n, const int* index, double* b)
{
    // When ii is set to a non-negative value, it will become the index
    // of the first nonvanishing element of b.
    int ii = -1;

    // We now do the forward substitution.
    // The only new wrinkle is to unscramble the permutation as we go.
    double sum;
    for(int i=0;  i < n;  ++i)
    {
        const int ip    = index[i];
        sum             = b[ip];
        b[ip]           = b[i];
        if(ii >= 0)
        {
            for(int j=ii; j < i; ++j)
            {
                sum -= A[i][j]*b[j];
            }
        }
        else if(sum != 0)
        {
            // A nonzero element was encountered, so from now on we will
            // have to do the sums in the loop above...
            ii = i;
        }
        b[i] = sum;
    }

    // Finally, do the backsubstitution.
    for (int i=n-1; i >= 0; --i)
    {
        sum = b[i];
        for (int j=i+1; j < n; ++j)
        {
            sum -= A[i][j]*b[j];
        }
        b[i] = sum / A[i][i];  // Store a component of the solution vector.
    }
}

int LU_decomposition(double** A, double* wksp, int n, int* index, double& sign)
{
    int imax = 0;

    sign = 1.0;  // No row interchanges yet.
    for(int i=0; i < n; ++i)
    {   // Loop over rows to get the implicit scaling information.
        double big = 0.0;
        for (int j=0; j < n; ++j)
        {
            const double temp = std::fabs(A[i][j]);
            if(temp > big)
            {
                big = temp;
            }
        }

        if(big == 0.0)
        {
            printf("Singular matrix in routine LU_decomposition\n");
            return 1;
        }
        wksp[i] = 1.0 / big;  // Store the scaling of each row.
    }

    for(int j=0; j < n; ++j)
    {
        for(int i=0; i < j; ++i)
        {
            double sum = A[i][j];
            for(int k=0; k < i; ++k)
            {
                sum -= A[i][k]*A[k][j];
            }
            A[i][j] = sum;
        }

        double big = 0.0;  // Initialize for the search for largest pivot element.
        for(int i=j; i < n; ++i)
        {   // This is i = j of equation (2.3.12) and i = j+1. . .N.
            // of equation (2.3.13).
            double sum = A[i][j];
            for (int k=0; k < j; ++k)
            {
                sum -= A[i][k]*A[k][j];
            }

            A[i][j] = sum;

            const double dum = wksp[i]*std::fabs(sum);
            if(dum >= big)
            {
                // Is the figure of merit for the pivot better than the best so far?
                big  = dum;
                imax = i;
            }
        }
        //assert(imax > -1);

        if(j != imax)
        {   // Do we need to interchange rows?
            for(int k=0; k < n; ++k)
            {   // Yes, do so...
                const double dum    = A[imax][k];
                A[imax][k]          = A[j][k];
                A[j][k]             = dum;
            }
            sign        = -sign;    // ...and change the parity...
            wksp[imax]  = wksp[j];  // ...and also interchange the scale factor.
        }
        index[j] = imax;

        // If the pivot element is zero, the matrix is singular (at least to the precision of the
        // algorithm). For some applications, it is desirable to substitute TINY for zero...
        if(A[j][j] == 0.0)
        {
            A[j][j] = 1.0e-20;
        }

        if(j != n)
        {   // Now, nally, divide by the pivot element.
            const double dum = 1.0 / A[j][j];
            for(int i=j+1; i < n; ++i)
            {
                A[i][j] *= dum;
            }
        }
    }  // Go back for the next column in the reduction.

    return 0;
}

void LU_invert_matrix(double** A, int n, double** A_inv, double& determinant)
{
    double** A_tmp = allocateMemoryForMatrix(n);
    double* wksp   = new double[n];

    // Make a copy of the original matrix
    for(int i=0; i < n; ++i)
    {
        for(int j=0; j < n; ++j)
        {
            A_tmp[i][j] = A[i][j];
        }
    }

    int* index = new int[n];
    LU_decomposition(A_tmp, wksp, n, index, determinant); // So far, either 1 or -1.

    for(int j=0; j < n; ++j)
    {
        determinant *= A_tmp[j][j];
    }

    if(std::fabs(determinant) < 1.0e-20)
    {
        printf("Singular matrix, no inverse was found.\n");
        std::exit(1);
    }

    double* col = new double[n];
    for(int j=0; j < n; ++j)
    {
        for(int i=0; i < n; ++i)
        {
            col[i] = 0.0;
        }
        col[j] = 1.0;

        LU_backsubstitution(A_tmp, n, index, col);

        for(int i=0; i < n; ++i)
        {
            A_inv[i][j] = col[i];
        }
    }

    freeMatrixMemory(A_tmp, n);
    delete[] wksp;
    delete[] index;
    delete[] col;
}