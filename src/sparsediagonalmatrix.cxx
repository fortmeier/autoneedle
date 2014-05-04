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

SparseDiagonalMatrix::SparseDiagonalMatrix( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{
  if( m < bandwidth_2 ) throw std::runtime_error("bandwitdh / 2 greater than matrix size");
  values = new double[m*bandwidth];
  zero();

}

SparseDiagonalMatrix::~SparseDiagonalMatrix()
{
	delete[] values;
}

void SparseDiagonalMatrix::zero()
{
  for( int i = 0; i < bandwidth; i++ )
  {
    for( int j = 0; j < size; j++ )
    {
      _at(i,j) = 0;
    }
  }
}

std::ostream& operator<< ( std::ostream &out, const SparseDiagonalMatrix &matrix )
{
  for( int j = 0; j < matrix.size; j++ )
  {
    int remaining = matrix.size;
    out << "[ ";
    for( int k = matrix.bandwidth_2; k < j && k < matrix.size - matrix.bandwidth_2-1; k++) 
    {
      out << "   . ";
      remaining--;
    }

    for( int i = 0; i < matrix.bandwidth; i++ )
    {
      out << std::setw( 4 ) << matrix._at(i,j) << " ";
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

double& SparseDiagonalMatrix::_at( int i, int j ) const
{
  return values[ i + bandwidth*j ];
}

double& SparseDiagonalMatrix::operator() ( int i, int j ) const
{
  if( j <= bandwidth_2) return _at(i,j);
  else if( j >= size - bandwidth_2) return _at(i-size+bandwidth,j);
  return _at(i-j+bandwidth_2,j);
}

cml::vectord SparseDiagonalMatrix::operator* (const cml::vectord& x) const
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

int SparseDiagonalMatrix::getSize()
{
  return size;
}


int SparseDiagonalMatrix::getBandwidth()
{
  return bandwidth;
}

void SparseDiagonalMatrix::multiplyWith( const cml::vectord& x, cml::vectord& r )
{
  // first, check if size of matrix and vectors are right:
  if( size != x.size() )
    throw std::runtime_error("size of x must be the same as matrix size");

  if( size != r.size() )
    throw std::runtime_error("size of r must be the same as matrix size");


  if (this->bandwidth % 4 != 0 ) 
    throw std::runtime_error("fast multiplication only works with matrix with bandwidth is multiple of 4");

  // now, do the multiplication

  // do the start non-optimized
  for(int j = 0; j < bandwidth_2+1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += (*this)(i,j) * x[i];
    }

  }

  // do the middle part fully optimized
  for(int j = bandwidth_2+1; j < size - bandwidth_2-1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i+=4 )
    {
      r[j] += _at(i+0,j) * x[i+0+j-bandwidth_2];
      r[j] += _at(i+1,j) * x[i+1+j-bandwidth_2];
      r[j] += _at(i+2,j) * x[i+2+j-bandwidth_2];
      r[j] += _at(i+3,j) * x[i+3+j-bandwidth_2];
    }

  }

  // do the end part non-optimized as well
  for(int j = size - bandwidth_2-1; j < size; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+size-bandwidth];
    }

  }

}