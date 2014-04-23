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

void testMatrix() {
  SparseDiagonalMatrix A( 10, 5 );
  A(0,0) = 1;
  A(1,1) = 2;
  A(0,2) = 3;
  A(3,3) = 4;
  A(4,4) = 5;
  A(5,5) = 6;
  A(6,6) = 7;
  A(9,7) = 8;
  A(8,8) = 9;
  A(9,9) = 3;

  cml::vectord r(10);
  r[0] = 11;
  r[1] = 22;
  r[2] = 33;
  r[3] = 44;
  r[4] = 55;
  r[5] = 66;
  r[6] = 77;
  r[7] = 88;
  r[8] = 99;
  r[9] = 11;
  std::cout<<A<<std::endl;
  std::cout<<r<<std::endl;
  std::cout<<A*r<<std::endl;
  exit(0);

}

int SparseDiagonalMatrix::getSize()
{
  return size;
}


int SparseDiagonalMatrix::getBandwidth()
{
  return bandwidth;
}

