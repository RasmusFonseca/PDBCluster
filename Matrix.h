/*
    PDBCluster - k-medoids and k-centers clustering of molecular structures.
    Copyright (C) 2014 Rasmus Fonseca

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Rasmus Fonseca - fonseca.rasmus@gmail.com
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

using namespace std;

/**
 * @brief Convenient wrapper for matrix data used by SVD and for multiplication.
 * The matrix data is stored in an array of doubles using column-major format which
 * is convenient for the SVD routine.
 */
class Matrix
{
public:
    /**
     * @brief Construct an mxn matrix with all 0-entries.
     * @param m The number of rows
     * @param n The number of columns
     */
    Matrix(int m, int n);

    /**
     * @brief Construct a "deep" copy of m
     * All values are copied using memcpy.
     * @param m The matrix to copy
     */
	Matrix(Matrix& m);

	~Matrix();

	double operator () (int i, int j) { return get(i,j); }

    /**
     * @brief Get the value of the (i,j)th entry
     * @param i The row index
     * @param j The column index
     * @return A double value
     */
	double get(int i,int j) const;

    /**
     * @brief Set the value of the (i,j)th entry
     * @param i The row index
     * @param j The column index
     * @param v A double value
     */
	void set(int i,int j, double v);

    /** @brief The number of rows in this matrix. */
	unsigned int rows() const;

    /** @brief The number of columns in this matrix. */
	unsigned int cols() const;

    double getDifferenceSum(Matrix*);

private:
	Matrix(double*, int m, int n);
	double *data;
    int m,n;

	friend class SVD;
};

ostream& operator<<(ostream& os, const Matrix& m);
ostream& operator<<(ostream& os, const Matrix* m);


#endif // MATRIX_H
