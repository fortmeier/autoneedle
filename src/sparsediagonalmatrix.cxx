/**
 * Copyright (C) 2014 Dirk Fortmeier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE
 *
 */

#include <iomanip>
#include <exception> 

#include "sparsediagonalmatrix.h"

template<typename Real>
SparseDiagonalMatrix<Real>::SparseDiagonalMatrix( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{
  if( m < bandwidth_2 ) throw std::runtime_error("bandwitdh / 2 greater than matrix size");
  values = new Real[m*bandwidth];
  zero();

}

template<typename Real>
SparseDiagonalMatrix<Real>::~SparseDiagonalMatrix()
{
	delete[] values;
}

template<typename Real>
void SparseDiagonalMatrix<Real>::zero()
{
  for( int i = 0; i < bandwidth; i++ )
  {
    for( int j = 0; j < size; j++ )
    {
      _at(i,j) = 0;
    }
  }
}

template<typename Real>
void SparseDiagonalMatrix<Real>::zeroRow( int j)
{
  for( int i = 0; i < bandwidth; i++ )
  {
    _at(i,j) = 0;
  }
}

template<typename Real>
std::ostream& SparseDiagonalMatrix<Real>::print ( std::ostream &out ) const
{
  for( int j = 0; j < size; j++ )
  {
    int remaining = size;
    out << "[ ";
    for( int k = bandwidth_2; k < j && k < size - bandwidth_2-1; k++) 
    {
      out << "   . ";
      remaining--;
    }

    for( int i = 0; i < bandwidth; i++ )
    {
      out << std::setw( 4 ) << _at(i,j) << " ";
      remaining--;
    }

    for( int k = 0; k < remaining; k++ ) 
    {
      out << "   . ";
    }

    out << " ]\n";
  }
  return out;
}

template<typename Real>
Real& SparseDiagonalMatrix<Real>::_at( int i, int j ) const
{
  return values[ i + bandwidth*j ];
}

template<typename Real>
Real& SparseDiagonalMatrix<Real>::operator() ( int i, int j ) const
{
  if( j <= bandwidth_2) return _at(i,j);
  else if( j >= size - bandwidth_2) return _at(i-size+bandwidth,j);
  return _at(i-j+bandwidth_2,j);
}

template<typename Real>
cml::vectord SparseDiagonalMatrix<Real>::sumRows() const
{
  cml::vectord r(getSize());
  for( int j = 0; j < getSize(); j++ )
  {
    for( int i = 0; i < bandwidth; i++)
    {
      r[i] += _at(i, j);
    }
  }
  return r;
}

template<typename Real>
cml::vectord SparseDiagonalMatrix<Real>::operator* (const cml::vectord& x) const
{
  cml::vectord r(x.size());

  for(int j = 0; j < bandwidth_2+1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += (*this)(i,j) * x[i];
    }

  }

  for(int j = bandwidth_2+1; j < size - bandwidth_2-1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+j-bandwidth_2];
    }

  }

  for(int j = size - bandwidth_2-1; j < size; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+size-bandwidth];
    }

  }

  return r;
}

template<typename Real>
int SparseDiagonalMatrix<Real>::getSize() const
{
  return size;
}

template<typename Real>
int SparseDiagonalMatrix<Real>::getBandwidth() const
{
  return bandwidth;
}

template<typename Real>
void SparseDiagonalMatrix<Real>::multiplyWith( const cml::vectord& x, cml::vectord& r ) const
{
  r = *this * x;
}

//template std::ostream& operator<< ( std::ostream &out, const BandMatrixInterface<double> &matrix );

template class SparseDiagonalMatrix<double>;