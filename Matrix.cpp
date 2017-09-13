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

#include "Matrix.h"

#include <assert.h>
#include <new>
#include <math.h>
#include <string.h>

Matrix::Matrix(int _m, int _n): data(new double[_m*_n]), m(_m), n(_n)
{
    memset(data, 0, m*n*sizeof(double));
}

Matrix::Matrix(double* newData, int _m, int _n): data(newData), m(_m), n(_n)
{}
Matrix::Matrix(Matrix& M): data(new double[M.m*M.n]), m(M.m), n(M.n)
{
    memcpy( data, M.data, m*n*sizeof(double) );
}

Matrix::~Matrix(){
	delete[] data;
}

unsigned int Matrix::rows() const { return m; }
unsigned int Matrix::cols() const { return n; }
double Matrix::get(int i,int j) const{
    assert ( i>=0 );
    assert ( i<m );
    assert ( j>=0 );
    assert ( j<n );
    return data[j*m+i];
}

void Matrix::set(int i,int j, double v){
    assert ( i>=0 );
    assert ( i<m );
    assert ( j>=0 );
    assert ( j<n );
    data[j*m+i] = v;
}

double Matrix::getDifferenceSum(Matrix* A)
{
    assert(m==A->m);
    assert(n==A->n);
	double ret = 0;
    for(int i=0;i<m*n;i++){
        double diff = data[i]-A->data[i];
        ret+=diff*diff;
    }
    return sqrt(ret/(m*n));
}

ostream& operator<<(ostream& os, const Matrix& m){
	for (unsigned int i=0; i<m.rows(); ++i) {
		for (unsigned int j=0; j<m.cols(); ++j){
            os << m.get(i,j);
			if(j == (m.cols()-1))
                os << endl;
			else
                os << '\t';
		}
	}
	return os;
}

ostream& operator<<(ostream& os, const Matrix* m){
	return os<<(*m);
}

